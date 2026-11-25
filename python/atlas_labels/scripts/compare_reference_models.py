"""
Side-by-side comparison of two reference depth-band models. Shows where
their per-layer medians and IQR bands differ along the normalised
arclength axis.
"""

import sys
import argparse
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from ribbon_unfold_laplace import COLORS  # noqa: E402
from build_reference_model import LAYER_NAMES, CORTEX_ORDER, MODEL_DIR  # noqa: E402


def load_model(tag):
    p = MODEL_DIR / f"depth_model_{tag}.npz"
    return np.load(p, allow_pickle=False), p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", default="best_of_section")
    ap.add_argument("--b", default="10x_all")
    args = ap.parse_args()

    m_a, p_a = load_model(args.a)
    m_b, p_b = load_model(args.b)
    print(f"A: {p_a}")
    print(f"B: {p_b}")

    centres_a = m_a["centres"]
    centres_b = m_b["centres"]

    fig, axes = plt.subplots(2, 1, figsize=(13, 9), constrained_layout=True,
                             sharex=True)
    ax_med, ax_diff = axes

    print("\nlayer-by-layer median ψ comparison (global medians across s):")
    print(f"  {'layer':>5}  {args.a:>20}  {args.b:>20}  shift")
    for L in CORTEX_ORDER:
        nm = LAYER_NAMES[L]
        prof_a = m_a[f"profile_{nm}"]
        prof_b = m_b[f"profile_{nm}"]
        med_a = prof_a[:, 2]
        med_b = prof_b[:, 2]
        col = COLORS[L]

        ax_med.plot(centres_a, med_a, color=col, lw=1.5, ls="-",
                    label=f"{nm} ({L}) — {args.a}")
        ax_med.plot(centres_b, med_b, color=col, lw=1.5, ls="--",
                    label=f"{nm} ({L}) — {args.b}")
        ax_med.fill_between(centres_a, prof_a[:, 1], prof_a[:, 3],
                            color=col, alpha=0.10, linewidth=0)
        ax_med.fill_between(centres_b, prof_b[:, 1], prof_b[:, 3],
                            color=col, alpha=0.10, linewidth=0)

        diff = med_b - med_a
        med_a_g = float(np.nanmedian(med_a))
        med_b_g = float(np.nanmedian(med_b))
        shift = med_b_g - med_a_g
        print(f"  {nm:>5}  {med_a_g:>20.3f}  {med_b_g:>20.3f}  {shift:+.3f}")

        ax_diff.plot(centres_a, diff, color=col, lw=1.0,
                     label=f"{nm} ({L})")

    ax_med.set_ylim(1.0, 0.0)
    ax_med.set_xlim(0, 1)
    ax_med.set_xlabel("normalised arclength s/s_max")
    ax_med.set_ylabel("ψ — median (lines) + IQR (shaded)")
    ax_med.set_title(f"per-layer ψ bands: solid = {args.a}, "
                     f"dashed = {args.b}")
    ax_med.legend(loc="lower right", fontsize=7, ncol=2)
    ax_med.grid(alpha=0.3)

    ax_diff.axhline(0, color="black", lw=0.5, alpha=0.5)
    ax_diff.set_xlabel("normalised arclength s/s_max")
    ax_diff.set_ylabel(f"ψ shift  ({args.b} − {args.a})")
    ax_diff.set_title("per-bin median shift  "
                      "(negative = 10x model puts the layer shallower)")
    ax_diff.legend(loc="lower right", fontsize=7, ncol=3)
    ax_diff.grid(alpha=0.3)

    out = MODEL_DIR / f"compare_{args.a}_vs_{args.b}.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
