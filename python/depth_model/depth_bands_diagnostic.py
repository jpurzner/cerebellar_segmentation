"""
Step 1 of the depth-based re-classifier: visual diagnostic of how cleanly
the cortical layers separate in (s, ψ) space.

Reads a `<slide>_ribbon_coords.npz` produced by `ribbon_unfold_laplace.py`
and produces three plots:

  A. Marginal ψ distribution per cortical layer — should show 5 separated
     peaks if the depth signature is clean. Overlaid with global-median
     interface positions (dashed lines) for context.

  B. Per-layer (s, ψ) density heatmap (one row per cortical layer). Each
     layer should appear as a continuous horizontal band; outlier swarms
     above or below the band are misannotation candidates.

  C. Smooth layer-interface curves on a single axis. Thin = raw per-bin
     midpoints; thick = running-median smoothed. Smooth, well-separated
     curves → re-classifier will work. Wobbly / crossing curves → not yet.

Also writes a `<slide>_depth_bands.npz` with the per-bin quantile profiles
and smoothed interface curves, which the re-classifier (next step) will
consume.
"""

import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Re-use palette + label dict from the main pipeline.
sys.path.insert(0, str(Path(__file__).parent))
from ribbon_unfold_laplace import COLORS, LBL, CORTEX_LABELS, OUTPUT_ROOT  # noqa: E402

LAYER_ORDER_PIA_TO_CORE = [
    LBL["oEGL"],   # 3
    LBL["iEGL"],   # 2
    LBL["ML"],     # 5
    LBL["PCL"],    # 7
    LBL["IGL"],    # 4
]
LAYER_NAMES = {3: "oEGL", 2: "iEGL", 5: "ML", 7: "PCL", 4: "IGL"}


def running_median(x: np.ndarray, w: int) -> np.ndarray:
    """Window-w running median that ignores NaNs."""
    n = len(x)
    out = np.full(n, np.nan, dtype=np.float64)
    half = w // 2
    for i in range(n):
        a = max(0, i - half)
        b = min(n, i + half + 1)
        window = x[a:b]
        window = window[~np.isnan(window)]
        if len(window) >= 3:
            out[i] = np.median(window)
    return out


def per_bin_quantiles(s, psi, layer, lbl_value, edges,
                      qs=(0.1, 0.25, 0.5, 0.75, 0.9), min_count=8):
    """Per s-bin ψ quantiles for one label."""
    sel = layer == lbl_value
    s_sel = s[sel]
    psi_sel = psi[sel]
    n_bins = len(edges) - 1
    out = np.full((n_bins, len(qs)), np.nan, dtype=np.float64)
    counts = np.zeros(n_bins, dtype=np.int64)
    if len(s_sel) == 0:
        return out, counts
    idx = np.clip(np.searchsorted(edges, s_sel, side="right") - 1,
                  0, n_bins - 1)
    for b in range(n_bins):
        m = idx == b
        c = int(m.sum())
        if c >= min_count:
            counts[b] = c
            out[b] = np.quantile(psi_sel[m], qs)
    return out, counts


def main(npz_path: Path, n_s_bins: int = 200, smoothing_window: int = 25):
    print(f"loading {npz_path}")
    d = np.load(npz_path, allow_pickle=False)
    layer = d["layer"]
    psi = d["psi"]
    s = d["s"]
    reached = d["reached"]
    s_max = float(d["s_max"])
    comp_ranges = d["comp_ranges"]
    slide_id = npz_path.stem.replace("_ribbon_coords", "")
    out_dir = npz_path.parent
    print(f"slide: {slide_id}")

    valid = (reached
             & np.isfinite(s)
             & np.isfinite(psi)
             & np.isin(layer, list(CORTEX_LABELS)))
    s = s[valid]
    psi = psi[valid]
    layer = layer[valid]
    print(f"valid cortex pixels: {len(s)}")

    edges = np.linspace(0.0, s_max, n_s_bins + 1)
    centres = 0.5 * (edges[:-1] + edges[1:])
    qs = (0.10, 0.25, 0.50, 0.75, 0.90)

    profiles = {}
    counts = {}
    for L in LAYER_ORDER_PIA_TO_CORE:
        prof, cnt = per_bin_quantiles(s, psi, layer, L, edges, qs)
        profiles[L] = prof
        counts[L] = cnt
        n_pix = int((layer == L).sum())
        n_active_bins = int((cnt > 0).sum())
        med_global = float(np.nanmedian(psi[layer == L])) if n_pix else float("nan")
        print(f"  layer {L} ({LAYER_NAMES[L]:>4}): n={n_pix:>7d}, "
              f"active bins={n_active_bins:>3d}, "
              f"global ψ p50={med_global:.3f}")

    # Smooth interface curves between adjacent layers.
    interfaces = []
    for L_above, L_below in zip(LAYER_ORDER_PIA_TO_CORE[:-1],
                                LAYER_ORDER_PIA_TO_CORE[1:]):
        med_above = profiles[L_above][:, 2]
        med_below = profiles[L_below][:, 2]
        iface = np.where(np.isfinite(med_above) & np.isfinite(med_below),
                         0.5 * (med_above + med_below),
                         np.nan)
        iface_smooth = running_median(iface, smoothing_window)
        med_glob = np.nanmedian(iface_smooth)
        p10 = np.nanpercentile(iface_smooth, 10)
        p90 = np.nanpercentile(iface_smooth, 90)
        gap_above = med_glob - float(np.nanmedian(profiles[L_above][:, 2]))
        gap_below = float(np.nanmedian(profiles[L_below][:, 2])) - med_glob
        interfaces.append((L_above, L_below, iface, iface_smooth))
        print(f"interface {LAYER_NAMES[L_above]:>4} | {LAYER_NAMES[L_below]:<4} "
              f"med ψ={med_glob:.3f}  spread {p10:.3f}–{p90:.3f}  "
              f"margins ±({gap_above:.3f}, {gap_below:.3f})")

    # ----------------------------- figure ----------------------------- #
    fig = plt.figure(figsize=(14, 17), constrained_layout=True)
    gs = fig.add_gridspec(7, 1, height_ratios=[3, 1.2, 1.2, 1.2, 1.2, 1.2, 4])
    ax_hist = fig.add_subplot(gs[0])
    ax_layers = [fig.add_subplot(gs[i + 1]) for i in range(5)]
    ax_iface = fig.add_subplot(gs[6])

    # Pane A — marginal histograms per cortical layer.
    bins = np.linspace(0.0, 1.0, 101)
    bcenters = 0.5 * (bins[:-1] + bins[1:])
    for L in LAYER_ORDER_PIA_TO_CORE:
        psi_L = psi[layer == L]
        if len(psi_L) == 0:
            continue
        h, _ = np.histogram(psi_L, bins=bins, density=True)
        ax_hist.plot(bcenters, h, color=COLORS[L], lw=1.6,
                     label=f"{L} {LAYER_NAMES[L]} (n={len(psi_L)})")
    for L_above, L_below, _, iface_smooth in interfaces:
        med = np.nanmedian(iface_smooth)
        ax_hist.axvline(med, color="black", lw=0.5, alpha=0.4, ls="--")
    ax_hist.set_xlabel("Laplace depth ψ (0 = pia, 1 = WM)")
    ax_hist.set_ylabel("density")
    ax_hist.set_title("A — marginal ψ distribution per cortical layer  "
                      "(dashed lines = median interface ψ)")
    ax_hist.legend(loc="upper right", fontsize=8)
    ax_hist.grid(alpha=0.3)
    ax_hist.set_xlim(0, 1)

    # Panes B — per-layer (s, ψ) density heatmap with overlaid band.
    psi_edges = np.linspace(0.0, 1.0, 121)
    for ax, L in zip(ax_layers, LAYER_ORDER_PIA_TO_CORE):
        s_L = s[layer == L]
        psi_L = psi[layer == L]
        if len(s_L) == 0:
            ax.set_axis_off()
            continue
        H, _, _ = np.histogram2d(s_L, psi_L, bins=[edges, psi_edges])
        H_log = np.log1p(H.T)
        ax.imshow(H_log, origin="upper", aspect="auto",
                  extent=[0.0, s_max, 1.0, 0.0],
                  cmap="magma", interpolation="nearest")
        prof = profiles[L]
        ax.plot(centres, prof[:, 2], color=COLORS[L], lw=1.0, alpha=0.95)
        ax.fill_between(centres, prof[:, 1], prof[:, 3],
                        color=COLORS[L], alpha=0.25, linewidth=0)
        # Bounding interfaces for this layer (above and below).
        for L_above, L_below, _, iface_smooth in interfaces:
            if L == L_above or L == L_below:
                ax.plot(centres, iface_smooth, color="white",
                        lw=0.7, alpha=0.85)
        for cid, lo, hi, sz in comp_ranges:
            ax.axvline(hi, color="white", lw=0.3, alpha=0.5)
        ax.set_ylim(1.0, 0.0)
        ax.set_xlim(0, s_max)
        ax.set_ylabel(f"ψ — {LAYER_NAMES[L]} ({L})", fontsize=9)
        if L != LAYER_ORDER_PIA_TO_CORE[-1]:
            ax.set_xticklabels([])
        else:
            ax.set_xlabel("arclength s along oEGL components (px)")
    ax_layers[0].set_title("B — (s, ψ) density per cortical layer  "
                           "(log scale; overlay = p50 line, p25–p75 shaded; "
                           "white = layer interface)", loc="left", fontsize=9)

    # Pane C — all interface curves on a single axis.
    for L_above, L_below, iface_raw, iface_smooth in interfaces:
        c = COLORS[L_below]
        ax_iface.plot(centres, iface_raw, color=c, lw=0.4, alpha=0.35)
        ax_iface.plot(centres, iface_smooth, color=c, lw=1.6,
                      label=f"{LAYER_NAMES[L_above]} | {LAYER_NAMES[L_below]}")
    for cid, lo, hi, sz in comp_ranges:
        ax_iface.axvline(hi, color="gray", lw=0.3, alpha=0.5)
    ax_iface.set_ylim(1.0, 0.0)
    ax_iface.set_xlim(0, s_max)
    ax_iface.set_xlabel("arclength s along oEGL components (px)")
    ax_iface.set_ylabel("interface ψ")
    ax_iface.set_title("C — smooth layer-interface curves "
                       "(thick = running-median smoothed, thin = raw "
                       "per-bin midpoint)", loc="left", fontsize=9)
    ax_iface.legend(loc="upper right", fontsize=8, ncol=2)
    ax_iface.grid(alpha=0.3)

    fig.suptitle(f"depth-band diagnostic — {slide_id}",
                 fontsize=12, fontweight="bold")
    out_png = out_dir / f"{slide_id}_depth_bands.png"
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    print(f"wrote {out_png}")

    # Persist the band profiles + interfaces for the re-classifier.
    save_dict = {
        "bin_centers": centres.astype(np.float32),
        "smoothing_window": smoothing_window,
        "n_s_bins": n_s_bins,
        "quantiles": np.array(qs, dtype=np.float32),
    }
    for L in LAYER_ORDER_PIA_TO_CORE:
        save_dict[f"profile_{LAYER_NAMES[L]}"] = profiles[L].astype(np.float32)
        save_dict[f"count_{LAYER_NAMES[L]}"] = counts[L].astype(np.int64)
    for L_above, L_below, raw, smooth in interfaces:
        key = f"interface_{LAYER_NAMES[L_above]}_{LAYER_NAMES[L_below]}"
        save_dict[key + "_raw"] = raw.astype(np.float32)
        save_dict[key + "_smooth"] = smooth.astype(np.float32)
    np.savez_compressed(out_dir / f"{slide_id}_depth_bands.npz", **save_dict)
    print(f"wrote {out_dir / (slide_id + '_depth_bands.npz')}")


if __name__ == "__main__":
    if len(sys.argv) >= 2:
        npz = Path(sys.argv[1])
    else:
        sid = "2018_05_22_s1_2_p27-0021_iter1deep"
        npz = OUTPUT_ROOT / sid / f"{sid}_ribbon_coords.npz"
    main(npz)
