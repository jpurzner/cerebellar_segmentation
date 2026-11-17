"""
Pool (s_norm, ψ, layer) data from a set of reference slides into a single
depth-band model that describes the typical cerebellar layer contours.

Selection strategies:
  - "best_of_section": top-scoring slide per distinct physical section
  - "top_n":           top N globally by quality_score
  - "all_above":       all slides with quality_score >= threshold

The model stores per-(s_norm bin, layer) ψ quantiles plus smoothed
layer-interface curves. `apply_reference_model.py` then uses this model
to predict each pixel's cortical layer from its (s_norm, ψ) coordinate.
"""

import sys
import re
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from ribbon_unfold_laplace import COLORS, LBL, OUT_ROOT, OUTPUT_ROOT  # noqa: E402
from depth_bands_diagnostic import running_median  # noqa: E402

# --------------------------------------------------------------------------- #
# Layer set
# --------------------------------------------------------------------------- #

# Cortex layers in canonical depth order from pia (ψ=0) to WM (ψ=1).
CORTEX_ORDER = [LBL["oEGL"], LBL["iEGL"], LBL["ML"], LBL["PCL"], LBL["IGL"]]
# Include DWL so we can disambiguate "deepest cortex" vs "white matter" at ψ≈1.
MODEL_LAYERS = CORTEX_ORDER + [LBL["DWL"]]
LAYER_NAMES = {
    LBL["oEGL"]: "oEGL",
    LBL["iEGL"]: "iEGL",
    LBL["ML"]: "ML",
    LBL["PCL"]: "PCL",
    LBL["IGL"]: "IGL",
    LBL["DWL"]: "DWL",
}

# --------------------------------------------------------------------------- #
# Defaults / paths
# --------------------------------------------------------------------------- #

DEFAULT_OUTPUTS = OUTPUT_ROOT
OUTPUTS: Path = DEFAULT_OUTPUTS
SCORECARD: Path = DEFAULT_OUTPUTS / "batch_scorecard.csv"
MODEL_DIR = OUT_ROOT / "reference_models"
N_S_BINS_DEFAULT = 200
SMOOTH_WINDOW_DEFAULT = 21
MIN_PIXELS_PER_BIN = 10


# --------------------------------------------------------------------------- #
# Slide selection
# --------------------------------------------------------------------------- #

def section_id_from_slide(slide_id: str) -> str:
    m = re.match(r"(\d{4}_\d{2}_\d{2}_s\d+_\d+_p27(?:-\d+)?)", slide_id)
    return m.group(1) if m else slide_id


def select_slides(df: pd.DataFrame, strategy: str, top_n: int = 8,
                  thresh: float = 0.85) -> pd.DataFrame:
    df = df.copy()
    df["section"] = df["slide_id"].apply(section_id_from_slide)
    if strategy == "best_of_section":
        best = (
            df.sort_values("quality_score", ascending=False)
              .groupby("section")
              .head(1)
              .sort_values("quality_score", ascending=False)
              .reset_index(drop=True)
        )
        return best
    elif strategy == "top_n":
        return df.sort_values("quality_score", ascending=False).head(top_n).reset_index(drop=True)
    elif strategy == "all_above":
        return (df[df["quality_score"] >= thresh]
                .sort_values("quality_score", ascending=False)
                .reset_index(drop=True))
    elif strategy == "all":
        return df.sort_values("quality_score", ascending=False).reset_index(drop=True)
    else:
        raise ValueError(f"unknown strategy: {strategy}")


# --------------------------------------------------------------------------- #
# Pool + bin
# --------------------------------------------------------------------------- #

def pool_pixel_observations(slides_df: pd.DataFrame):
    s_norm_chunks, psi_chunks, layer_chunks, slide_chunks = [], [], [], []
    used_slides = []
    for i, row in slides_df.iterrows():
        slide_id = row["slide_id"]
        p = OUTPUTS / slide_id / f"{slide_id}_ribbon_coords.npz"
        if not p.exists():
            print(f"  WARN: missing {p}, skipping")
            continue
        d = np.load(p, allow_pickle=False)
        layer = d["layer"]
        psi = d["psi"]
        s = d["s"]
        reached = d["reached"]
        s_max = float(d["s_max"])
        valid = (reached
                 & np.isfinite(s)
                 & np.isfinite(psi)
                 & np.isin(layer, MODEL_LAYERS))
        s_norm = s[valid] / s_max
        s_norm_chunks.append(s_norm.astype(np.float32))
        psi_chunks.append(psi[valid].astype(np.float32))
        layer_chunks.append(layer[valid].astype(np.uint8))
        slide_chunks.append(np.full(int(valid.sum()), len(used_slides),
                                    dtype=np.int32))
        used_slides.append(slide_id)
        print(f"  pooled {slide_id}: {int(valid.sum())} px "
              f"(quality={row['quality_score']:.3f})")
    s_norm = np.concatenate(s_norm_chunks)
    psi = np.concatenate(psi_chunks)
    layer = np.concatenate(layer_chunks)
    slide_idx = np.concatenate(slide_chunks)
    return s_norm, psi, layer, slide_idx, used_slides


def per_bin_layer_quantiles(s_norm, psi, layer, n_bins, qs):
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    centres = 0.5 * (edges[:-1] + edges[1:])
    bins = np.clip(np.searchsorted(edges, s_norm, side="right") - 1,
                   0, n_bins - 1)
    profiles = {}
    counts = {}
    for L in MODEL_LAYERS:
        sel = layer == L
        prof = np.full((n_bins, len(qs)), np.nan, dtype=np.float64)
        cnt = np.zeros(n_bins, dtype=np.int64)
        if sel.sum() > 0:
            bins_L = bins[sel]
            psi_L = psi[sel].astype(np.float64)
            for b in range(n_bins):
                m = bins_L == b
                c = int(m.sum())
                if c >= MIN_PIXELS_PER_BIN:
                    cnt[b] = c
                    prof[b] = np.quantile(psi_L[m], qs)
        profiles[L] = prof
        counts[L] = cnt
        active = int((cnt > 0).sum())
        n_pix = int(sel.sum())
        print(f"  layer {L} ({LAYER_NAMES[L]:>4}): n={n_pix:>8d}  "
              f"active bins={active:>4d}/{n_bins}")
    return centres, profiles, counts


# --------------------------------------------------------------------------- #
# Interface curves
# --------------------------------------------------------------------------- #

def interface_curves(profiles, smoothing_window):
    """Smooth midpoint between adjacent cortex layers' p50 medians."""
    cortex = CORTEX_ORDER
    out = []
    for La, Lb in zip(cortex[:-1], cortex[1:]):
        med_a = profiles[La][:, 2]
        med_b = profiles[Lb][:, 2]
        raw = np.where(
            np.isfinite(med_a) & np.isfinite(med_b),
            0.5 * (med_a + med_b),
            np.nan,
        )
        smooth = running_median(raw, smoothing_window)
        out.append({
            "layer_above": La,
            "layer_below": Lb,
            "raw": raw,
            "smooth": smooth,
        })
        med = float(np.nanmedian(smooth))
        print(f"  interface {LAYER_NAMES[La]:>4} | {LAYER_NAMES[Lb]:<4}: "
              f"median ψ = {med:.3f}")
    return out


# --------------------------------------------------------------------------- #
# Plot the model
# --------------------------------------------------------------------------- #

def plot_model(centres, profiles, counts, interfaces, used_slides, out_path):
    fig = plt.figure(figsize=(15, 11), constrained_layout=True)
    gs = fig.add_gridspec(3, 1, height_ratios=[4, 3, 0.5])
    ax_bands = fig.add_subplot(gs[0])
    ax_counts = fig.add_subplot(gs[1], sharex=ax_bands)
    ax_legend = fig.add_subplot(gs[2])

    # Pane 1: median + IQR band for each layer; interfaces as black curves.
    for L in MODEL_LAYERS:
        col = COLORS[L]
        prof = profiles[L]
        ax_bands.fill_between(centres, prof[:, 1], prof[:, 3],
                              color=col, alpha=0.30, linewidth=0,
                              label=f"{LAYER_NAMES[L]} ({L}) p25–p75")
        ax_bands.plot(centres, prof[:, 2], color=col, lw=1.6,
                      label=f"{LAYER_NAMES[L]} median")
    for iface in interfaces:
        ax_bands.plot(centres, iface["smooth"], color="black", lw=1.2,
                      alpha=0.85)
        # label at right edge
        last_val = next((v for v in iface["smooth"][::-1] if np.isfinite(v)),
                        None)
        if last_val is not None:
            ax_bands.text(
                1.005, last_val,
                f"{LAYER_NAMES[iface['layer_above']]}|"
                f"{LAYER_NAMES[iface['layer_below']]}",
                fontsize=8, va="center", ha="left", color="black",
            )
    ax_bands.set_ylim(1.0, 0.0)
    ax_bands.set_xlim(0, 1)
    ax_bands.set_xlabel("normalised arclength s / s_max")
    ax_bands.set_ylabel("Laplace depth ψ")
    ax_bands.set_title(
        f"pooled depth-band reference — {len(used_slides)} slides, "
        f"median + IQR per layer; black curves = smooth interfaces"
    )
    ax_bands.grid(alpha=0.3)
    ax_bands.legend(loc="lower right", fontsize=7, ncol=2)

    # Pane 2: counts per bin per layer (how many pooled pixels at each s_norm)
    for L in MODEL_LAYERS:
        ax_counts.plot(centres, counts[L], color=COLORS[L], lw=1.0,
                       label=LAYER_NAMES[L])
    ax_counts.set_yscale("log")
    ax_counts.set_xlim(0, 1)
    ax_counts.set_ylabel("pooled pixels per bin")
    ax_counts.set_xlabel("normalised arclength s / s_max")
    ax_counts.set_title("per-bin sample sizes (log) — drives the trust in "
                        "the band edges")
    ax_counts.grid(alpha=0.3, which="both")
    ax_counts.legend(loc="lower right", fontsize=7, ncol=3)

    # Bottom strip: list which slides contributed.
    ax_legend.set_axis_off()
    ax_legend.text(0.0, 0.5,
                   "Slides pooled: " + ", ".join(used_slides),
                   fontsize=7, va="center", ha="left",
                   wrap=True)

    fig.savefig(out_path, dpi=140)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #

def main():
    global OUTPUTS, SCORECARD
    ap = argparse.ArgumentParser()
    ap.add_argument("--strategy", default="best_of_section",
                    choices=["best_of_section", "top_n", "all_above", "all"])
    ap.add_argument("--top_n", type=int, default=8)
    ap.add_argument("--thresh", type=float, default=0.85)
    ap.add_argument("--tag", default=None,
                    help="tag used in output filename; default = strategy")
    ap.add_argument("--n_s_bins", type=int, default=N_S_BINS_DEFAULT)
    ap.add_argument("--smoothing_window", type=int,
                    default=SMOOTH_WINDOW_DEFAULT)
    ap.add_argument("--ribbon_root", type=Path, default=DEFAULT_OUTPUTS,
                    help="directory containing per-slide ribbon_coords.npz "
                         "+ batch_scorecard.csv")
    args = ap.parse_args()
    tag = args.tag or args.strategy
    OUTPUTS = args.ribbon_root
    SCORECARD = OUTPUTS / "batch_scorecard.csv"

    MODEL_DIR.mkdir(parents=True, exist_ok=True)
    print(f"loading scorecard: {SCORECARD}")
    df = pd.read_csv(SCORECARD)
    slides = select_slides(df, args.strategy, args.top_n, args.thresh)
    print(f"strategy={args.strategy}: {len(slides)} slides selected")
    print(slides[["slide_id", "quality_score"]].to_string(index=False))

    print("\npooling per-pixel (s_norm, ψ, layer) ...")
    s_norm, psi, layer, slide_idx, used = pool_pixel_observations(slides)
    print(f"total pooled pixels: {len(psi)}")

    print(f"\nbinning into {args.n_s_bins} s_norm bins ...")
    qs = (0.10, 0.25, 0.50, 0.75, 0.90)
    centres, profiles, counts = per_bin_layer_quantiles(
        s_norm, psi, layer, args.n_s_bins, qs
    )

    print(f"\nsmoothing interfaces (window={args.smoothing_window}) ...")
    interfaces = interface_curves(profiles, args.smoothing_window)

    # ---- save ----
    save = {
        "centres": centres.astype(np.float32),
        "quantiles": np.array(qs, dtype=np.float32),
        "n_s_bins": int(args.n_s_bins),
        "smoothing_window": int(args.smoothing_window),
        "strategy": args.strategy,
        "slides_used": np.array(used, dtype="U200"),
        "model_layers": np.array(MODEL_LAYERS, dtype=np.int32),
    }
    for L in MODEL_LAYERS:
        nm = LAYER_NAMES[L]
        save[f"profile_{nm}"] = profiles[L].astype(np.float32)
        save[f"count_{nm}"] = counts[L].astype(np.int64)
    for iface in interfaces:
        a = LAYER_NAMES[iface["layer_above"]]
        b = LAYER_NAMES[iface["layer_below"]]
        save[f"interface_{a}_{b}_smooth"] = iface["smooth"].astype(np.float32)
        save[f"interface_{a}_{b}_raw"] = iface["raw"].astype(np.float32)

    out_npz = MODEL_DIR / f"depth_model_{tag}.npz"
    np.savez_compressed(out_npz, **save)
    print(f"\nwrote {out_npz}")

    out_png = MODEL_DIR / f"depth_model_{tag}.png"
    plot_model(centres, profiles, counts, interfaces, used, out_png)
    print(f"wrote {out_png}")


if __name__ == "__main__":
    main()
