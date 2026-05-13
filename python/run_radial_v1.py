"""Run v1 smooth-angle radial detector on s1_2 using the port's masks."""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import tifffile
from radial_spokes import build_spokes, PIAL_SAMPLE_UM

INPUT = Path(
    "/Users/jpurzner/Dropbox/images/edu_repeat/p27/"
    "2018_05_22_s1_2_p27-0021/2018_05_22_s1_2_p27-0021_fused_crop.tif"
)
PORT_MASKS = Path(__file__).resolve().parent / "figs" / "v_port" / "port_masks.npz"
OUT = Path(__file__).resolve().parent / "figs" / "v_radial_v1"
OUT.mkdir(parents=True, exist_ok=True)
PX_UM = 0.5119049
DRAW_EVERY = 8


def main():
    print(f"loading {INPUT.name}")
    stack = tifffile.imread(INPUT)
    p27_raw, neun_raw, dapi_raw = stack[0], stack[1], stack[2]
    def to_unit(im):
        im = im.astype(np.float32); lo, hi = im.min(), im.max()
        return (im - lo) / max(hi - lo, 1e-9)
    rgb = np.stack([to_unit(p27_raw), to_unit(neun_raw), to_unit(dapi_raw)], axis=-1)
    rgb /= max(rgb.max(), 1e-9)
    H, W = p27_raw.shape

    print(f"loading port masks {PORT_MASKS}")
    z = np.load(PORT_MASKS)
    cereb = z['mask'].astype(bool)
    igl   = z['c_final'].astype(bool)
    print(f"  cereb area = {cereb.sum()*(PX_UM/1000)**2:.3f} mm^2")
    print(f"  IGL area   = {igl.sum()*(PX_UM/1000)**2:.3f} mm^2")

    print("building spokes")
    sl = build_spokes(cereb, igl, pixel_size_um=PX_UM)
    print(f"  total spokes: {len(sl.pia_points)}")
    print(f"  valid (reach IGL): {sl.valid.sum()}/{len(sl.valid)} "
          f"({100*sl.valid.mean():.1f}%)")
    if sl.valid.any():
        print(f"  ray length um: median={np.nanmedian(sl.lengths_um):.0f}  "
              f"p10={np.nanpercentile(sl.lengths_um,10):.0f}  "
              f"p90={np.nanpercentile(sl.lengths_um,90):.0f}")
    print(f"  angle deflection: median={np.median(np.abs(sl.angle_deflection_deg)):.1f} deg, "
          f"max={np.max(np.abs(sl.angle_deflection_deg)):.1f} deg")

    rgb_dim = rgb * 0.45

    # plot 1: full overview
    fig, ax = plt.subplots(1, 1, figsize=(16, 14))
    ax.imshow(rgb_dim)
    ax.contour(igl, levels=[0.5], colors=["cyan"], linewidths=0.6, alpha=0.85)
    ax.contour(cereb, levels=[0.5], colors=["yellow"], linewidths=0.4, alpha=0.6)
    n = len(sl.pia_points)
    for i in range(0, n, DRAW_EVERY):
        py, px = sl.pia_points[i]
        if sl.valid[i]:
            ye, xe = sl.igl_endpoints[i]
            ax.plot([px, xe], [py, ye], "-", color="magenta", lw=0.8, alpha=0.9)
            ax.plot(px, py, ".", color="orange", markersize=2.5)
            ax.plot(xe, ye, ".", color="lime", markersize=2.5)
        else:
            ax.plot(px, py, "x", color="red", markersize=4)
    ax.set_title(f"v1 radial spokes: pia (orange) -> IGL (lime)  "
                 f"every {PIAL_SAMPLE_UM*DRAW_EVERY:.0f} um  "
                 f"valid={sl.valid.sum()}/{n}")
    ax.axis("off")
    fig.tight_layout(); fig.savefig(OUT / "01_spokes_full.png", dpi=160)
    plt.close(fig)

    # plot 2: zoom on a deep fold area (centroid of cereb)
    cy, cx = np.array(np.where(cereb)).mean(axis=1).astype(int)
    half = 800
    y0, y1 = max(0, cy - half), min(H, cy + half)
    x0, x1 = max(0, cx - half), min(W, cx + half)
    fig, ax = plt.subplots(1, 1, figsize=(14, 14))
    ax.imshow(rgb_dim[y0:y1, x0:x1])
    ax.contour(igl[y0:y1, x0:x1], levels=[0.5], colors=["cyan"], linewidths=0.7,
                alpha=0.85)
    for i in range(0, n, max(1, DRAW_EVERY//2)):
        py, px = sl.pia_points[i]
        if not (y0 <= py < y1 and x0 <= px < x1):
            continue
        if sl.valid[i]:
            ye, xe = sl.igl_endpoints[i]
            ax.plot([px-x0, xe-x0], [py-y0, ye-y0], "-", color="magenta",
                    lw=1.0, alpha=0.95)
            ax.plot(px-x0, py-y0, ".", color="orange", markersize=4)
            ax.plot(xe-x0, ye-y0, ".", color="lime", markersize=4)
        else:
            ax.plot(px-x0, py-y0, "x", color="red", markersize=6)
    ax.set_title("zoom: pia -> IGL spokes (deep folds visible)")
    ax.axis("off")
    fig.tight_layout(); fig.savefig(OUT / "02_spokes_zoom.png", dpi=160)
    plt.close(fig)

    # plot 3: stats — length distribution + angle deflection
    fig, ax = plt.subplots(1, 2, figsize=(14, 5))
    L = sl.lengths_um[sl.valid]
    ax[0].hist(L, bins=60); ax[0].set_xlabel("ray length to IGL (um)")
    ax[0].set_title(f"valid spokes n={len(L)}, median={np.median(L):.0f}um")
    ax[1].hist(sl.angle_deflection_deg, bins=np.linspace(-60, 60, 31))
    ax[1].set_xlabel("angle deflection from local normal (deg)")
    ax[1].set_title(f"saturated at +/-60: {(np.abs(sl.angle_deflection_deg)>=58).sum()}")
    fig.tight_layout(); fig.savefig(OUT / "03_stats.png", dpi=120)
    plt.close(fig)

    print(f"figures saved to {OUT}")


if __name__ == "__main__":
    main()
