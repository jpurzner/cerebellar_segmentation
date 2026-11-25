"""
Unfold a labelled P7 cerebellum section into a 2D (arclength, depth) ribbon.

Input: per-pixel label TIFF with the scheme
  0 bg, 1 cereb, 2 iEGL, 3 oEGL, 4 IGL, 5 ML, 6 DWL, 7 PCL, (8 DCN, 9 anomaly)

For every pixel in the cortical ribbon (EGL+ML+PCL+IGL) we compute:
  d = Euclidean distance from the pial surface (outer boundary of tissue)
  s = arclength of the nearest pial point, measured along the pia
The cortex is then re-rasterized in (s, d) space — bands should appear as
horizontal stripes (oEGL → iEGL → ML → PCL → IGL from d=0 inward).

This is the empirical target the parametric ribbon-growth model needs to match.
"""

from pathlib import Path
import numpy as np
import cv2
from scipy import ndimage as ndi
from PIL import Image
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

LBL = dict(BG=0, CEREB=1, iEGL=2, oEGL=3, IGL=4, ML=5, DWL=6, PCL=7, DCN=8, ANOM=9)
CORTEX_LABELS = (LBL["iEGL"], LBL["oEGL"], LBL["IGL"], LBL["ML"], LBL["PCL"])
TISSUE_LABELS = tuple(range(1, 10))

# Colors aligned with HANDOFF.md
COLORS = {
    0: (0.00, 0.00, 0.00),  # bg
    1: (0.10, 0.10, 0.40),  # cereb dark blue
    2: (0.40, 0.70, 0.95),  # iEGL sky blue
    3: (0.00, 0.80, 0.90),  # oEGL cyan
    4: (1.00, 0.95, 0.20),  # IGL yellow
    5: (0.55, 0.90, 0.55),  # ML light green
    6: (1.00, 0.55, 0.10),  # DWL orange
    7: (0.95, 0.10, 0.10),  # PCL red
    8: (0.85, 0.20, 0.85),  # DCN magenta
    9: (0.40, 0.40, 0.40),  # anomaly gray
}


def load_labels(path: Path) -> np.ndarray:
    img = np.array(Image.open(path))
    return img.astype(np.int32)


def downsample(arr: np.ndarray, factor: int) -> np.ndarray:
    # Nearest-neighbor downsample for label maps.
    return arr[::factor, ::factor]


def tissue_mask(lbl: np.ndarray) -> np.ndarray:
    m = np.isin(lbl, TISSUE_LABELS)
    # Take largest connected component to drop stray junk.
    n, comp, stats, _ = cv2.connectedComponentsWithStats(m.astype(np.uint8), 8)
    if n <= 1:
        return m
    biggest = 1 + np.argmax(stats[1:, cv2.CC_STAT_AREA])
    return comp == biggest


def pial_contour(tissue: np.ndarray) -> np.ndarray:
    """Outermost closed contour of the tissue (longest one)."""
    contours, _ = cv2.findContours(
        tissue.astype(np.uint8), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE
    )
    if not contours:
        raise RuntimeError("no tissue contour")
    contour = max(contours, key=lambda c: c.shape[0])
    return contour[:, 0, :]  # (N, 2) as (x, y)


def smooth_contour(xy: np.ndarray, win: int = 15) -> np.ndarray:
    # Circular moving average — pia is closed.
    k = np.ones(win) / win
    pad = np.concatenate([xy[-win:], xy, xy[:win]], axis=0)
    sx = np.convolve(pad[:, 0], k, mode="same")[win:-win]
    sy = np.convolve(pad[:, 1], k, mode="same")[win:-win]
    return np.stack([sx, sy], axis=1)


def arclength_param(xy: np.ndarray) -> np.ndarray:
    d = np.linalg.norm(np.diff(xy, axis=0, append=xy[:1]), axis=1)
    s = np.concatenate([[0.0], np.cumsum(d)[:-1]])
    return s


def build_ribbon(
    lbl: np.ndarray, tissue: np.ndarray, pia_xy: np.ndarray, pia_s: np.ndarray
):
    """For every cortex pixel, return (s, d, layer) arrays."""
    cortex = np.isin(lbl, CORTEX_LABELS) & tissue

    # Pia rasterized as a binary mask the same size as the image.
    H, W = lbl.shape
    pia_mask = np.zeros((H, W), dtype=np.uint8)
    px = np.clip(pia_xy[:, 0].astype(int), 0, W - 1)
    py = np.clip(pia_xy[:, 1].astype(int), 0, H - 1)
    pia_mask[py, px] = 1

    # Distance to pia + index of nearest pia pixel.
    inv = (pia_mask == 0).astype(np.uint8)
    dist, (iy, ix) = ndi.distance_transform_edt(inv, return_indices=True)

    # Map each cortex pixel's nearest pia (y, x) back to its index in pia_xy.
    pia_idx_grid = -np.ones((H, W), dtype=np.int64)
    pia_idx_grid[py, px] = np.arange(pia_xy.shape[0])
    nn_pia_idx = pia_idx_grid[iy, ix]

    ys, xs = np.where(cortex)
    d = dist[ys, xs]
    s = pia_s[nn_pia_idx[ys, xs]]
    layer = lbl[ys, xs]
    return s, d, layer, cortex


def rasterize_unrolled(s, d, layer, s_bins=2000, d_bins=200):
    s_max = s.max()
    d_max = d.max()
    # Bin counts per (s, d, label).
    s_i = np.clip((s / s_max * (s_bins - 1)).astype(int), 0, s_bins - 1)
    d_i = np.clip((d / d_max * (d_bins - 1)).astype(int), 0, d_bins - 1)
    out = np.zeros((d_bins, s_bins), dtype=np.int32)
    # Take the modal (majority) label per bin — quick proxy.
    flat_idx = d_i * s_bins + s_i
    order = np.argsort(flat_idx, kind="stable")
    fi_sorted = flat_idx[order]
    lay_sorted = layer[order]
    # Per-cell majority via lexsort trick.
    boundaries = np.concatenate(
        [[0], np.where(np.diff(fi_sorted) != 0)[0] + 1, [len(fi_sorted)]]
    )
    for a, b in zip(boundaries[:-1], boundaries[1:]):
        vals, cnts = np.unique(lay_sorted[a:b], return_counts=True)
        winner = vals[cnts.argmax()]
        cell = fi_sorted[a]
        out.flat[cell] = winner
    return out, s_max, d_max


def label_rgb(label_img: np.ndarray) -> np.ndarray:
    rgb = np.zeros((*label_img.shape, 3), dtype=np.float32)
    for v, c in COLORS.items():
        rgb[label_img == v] = c
    return rgb


def main():
    src = Path(
        "/Users/jpurzner/Dropbox/imaging_analysis/cerebellar_segmentation/python/labelled/2018_05_22_s1_2_p27-0021_labelled.tiff"
    )
    out_dir = Path(
        "/Users/jpurzner/Dropbox/Q_research/py_projects/cerebellar_atlas_labels/outputs/ribbon_unfold"
    )
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"loading {src.name}")
    lbl_full = load_labels(src)
    print(f"  full shape: {lbl_full.shape}, dtype {lbl_full.dtype}")

    # Downsample for tractable distance transform (full res is 4788x4860).
    DS = 4
    lbl = downsample(lbl_full, DS)
    print(f"  downsampled by {DS}: {lbl.shape}")
    for v in np.unique(lbl):
        print(f"    label {v}: {(lbl==v).sum()} px")

    tissue = tissue_mask(lbl)
    print(f"tissue pixels: {tissue.sum()}")

    pia_raw = pial_contour(tissue)
    pia_xy = smooth_contour(pia_raw.astype(np.float32), win=21)
    pia_s = arclength_param(pia_xy)
    print(f"pia contour: {pia_xy.shape[0]} pts, total length {pia_s[-1]:.1f}px")

    s, d, layer, cortex_mask = build_ribbon(lbl, tissue, pia_xy, pia_s)
    print(f"cortex pixels: {cortex_mask.sum()}")
    print(f"depth range: {d.min():.1f} – {d.max():.1f} px")
    for v in np.unique(layer):
        sel = layer == v
        print(
            f"  label {v}: n={sel.sum()}, "
            f"d mean={d[sel].mean():.1f}, d p10={np.percentile(d[sel],10):.1f}, "
            f"d p90={np.percentile(d[sel],90):.1f}"
        )

    print("rasterising unrolled ribbon ...")
    unrolled, s_max, d_max = rasterize_unrolled(s, d, layer, s_bins=2400, d_bins=180)
    print(f"  unrolled shape: {unrolled.shape} (d x s)")
    print(f"  s_max={s_max:.1f}px, d_max={d_max:.1f}px")

    # --- Figure ---
    fig, axes = plt.subplots(2, 1, figsize=(14, 9), constrained_layout=True)

    ax = axes[0]
    ax.imshow(label_rgb(lbl))
    ax.plot(pia_xy[:, 0], pia_xy[:, 1], "w-", lw=0.6, alpha=0.8)
    # Tick marks at 8 evenly spaced arclengths so the reader can map s→space.
    n_ticks = 8
    tick_idx = np.linspace(0, pia_xy.shape[0] - 1, n_ticks).astype(int)
    for i, ti in enumerate(tick_idx):
        ax.plot(pia_xy[ti, 0], pia_xy[ti, 1], "o", color="white", ms=6)
        ax.annotate(
            f"{pia_s[ti]/pia_s[-1]:.2f}",
            (pia_xy[ti, 0], pia_xy[ti, 1]),
            color="white",
            fontsize=8,
            xytext=(4, 4),
            textcoords="offset points",
        )
    ax.set_title("native space (pia in white, ticks = fractional arclength)")
    ax.set_axis_off()

    ax = axes[1]
    ax.imshow(
        label_rgb(unrolled),
        aspect="auto",
        extent=[0, s_max, d_max, 0],
        interpolation="nearest",
    )
    ax.set_xlabel("arclength s along pia (px)")
    ax.set_ylabel("depth from pia d (px)")
    ax.set_title("unrolled (s, d) ribbon — layers should appear as horizontal bands")

    out_png = out_dir / "ribbon_unfold.png"
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    print(f"wrote {out_png}")

    # Also dump the raw (s, d, layer) arrays for downstream use.
    np.savez_compressed(
        out_dir / "ribbon_coords.npz",
        s=s.astype(np.float32),
        d=d.astype(np.float32),
        layer=layer.astype(np.uint8),
        pia_xy=pia_xy.astype(np.float32),
        pia_s=pia_s.astype(np.float32),
        s_max=s_max,
        d_max=d_max,
        downsample=DS,
    )
    print(f"wrote {out_dir/'ribbon_coords.npz'}")


if __name__ == "__main__":
    main()
