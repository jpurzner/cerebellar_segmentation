"""
Apply a reference depth-band model to every slide in a ribbon-unfold output
root, then aggregate the per-slide overwrite counts and transition matrices
into a single summary report.

Outputs (under <ribbon_root>/../reference_models/applied/):
  apply_summary_<model>.csv   — per-slide overwrite stats
  apply_summary_<model>.png   — per-slide bars + pooled transition heatmap
"""

import sys
import argparse
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from ribbon_unfold_laplace import COLORS  # noqa: E402

SCRIPT_DIR = Path(__file__).parent
APPLIED_DIR = Path(
    "/Users/jpurzner/Dropbox/Q_research/py_projects/cerebellar_atlas_labels/"
    "outputs/reference_models/applied"
)

# pia → core ordering for the transition matrix axes
CORTEX_ORDER = [3, 2, 5, 7, 4]
LAYER_NAMES = {3: "oEGL", 2: "iEGL", 5: "ML", 7: "PCL", 4: "IGL"}


def run_apply(slide_id: str, model: str, ribbon_root: Path, extra_args):
    cmd = [
        "/usr/bin/python3", str(SCRIPT_DIR / "apply_reference_model.py"),
        "--slide", slide_id,
        "--model", model,
        "--ribbon_root", str(ribbon_root),
    ] + list(extra_args)
    res = subprocess.run(cmd, capture_output=True, text=True)
    return res.returncode, res.stdout, res.stderr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="10x_all")
    ap.add_argument("--ribbon_root", type=Path, required=True,
                    help="ribbon-unfold output root holding per-slide dirs + "
                         "batch_scorecard.csv")
    ap.add_argument("--conf_thresh", type=float, default=None)
    ap.add_argument("--margin_thresh", type=float, default=None)
    ap.add_argument("--prior_boost", type=float, default=None)
    args = ap.parse_args()

    scorecard = args.ribbon_root / "batch_scorecard.csv"
    slide_ids = pd.read_csv(scorecard)["slide_id"].tolist()
    print(f"applying model '{args.model}' to {len(slide_ids)} slides "
          f"from {args.ribbon_root}")

    extra = []
    if args.conf_thresh is not None:
        extra += ["--conf_thresh", str(args.conf_thresh)]
    if args.margin_thresh is not None:
        extra += ["--margin_thresh", str(args.margin_thresh)]
    if args.prior_boost is not None:
        extra += ["--prior_boost", str(args.prior_boost)]

    rows = []
    pooled_trans = np.zeros((len(CORTEX_ORDER), len(CORTEX_ORDER)),
                            dtype=np.int64)
    order_idx = {L: i for i, L in enumerate(CORTEX_ORDER)}

    for i, sid in enumerate(slide_ids):
        print(f"[{i+1:>2}/{len(slide_ids)}] {sid} ...", flush=True)
        rc, out, err = run_apply(sid, args.model, args.ribbon_root, extra)
        if rc != 0:
            print(f"   !! apply failed (rc={rc})")
            print("   " + err.strip().splitlines()[-1] if err.strip() else "")
            continue

        pred_p = APPLIED_DIR / sid / f"{sid}_predicted.npz"
        if not pred_p.exists():
            print(f"   !! no predicted npz at {pred_p}")
            continue
        d = np.load(pred_p, allow_pickle=False)
        cur = d["layer_current"]
        prd = d["layer_predicted"]
        ow = d["overwrote"]
        n_valid = len(cur)
        n_disagree = int((cur != prd).sum())
        n_ow = int(ow.sum())
        rows.append({
            "slide_id": sid,
            "n_valid": n_valid,
            "n_disagree": n_disagree,
            "pct_disagree": 100.0 * n_disagree / max(1, n_valid),
            "n_overwrite": n_ow,
            "pct_overwrite": 100.0 * n_ow / max(1, n_valid),
        })
        # Pool overwrite transitions.
        for c, p in zip(cur[ow], prd[ow]):
            if int(c) in order_idx and int(p) in order_idx:
                pooled_trans[order_idx[int(c)], order_idx[int(p)]] += 1
        print(f"   disagree {n_disagree} ({100*n_disagree/n_valid:.1f}%)  "
              f"overwrite {n_ow} ({100*n_ow/n_valid:.1f}%)", flush=True)

    df = pd.DataFrame(rows).sort_values("pct_overwrite",
                                        ascending=False).reset_index(drop=True)
    csv_out = APPLIED_DIR / f"apply_summary_{args.model}.csv"
    df.to_csv(csv_out, index=False)
    print(f"\nwrote {csv_out}")
    print("\n=== per-slide overwrite summary ===")
    print(df.to_string(index=False))

    # ----------------------------- figure ----------------------------- #
    fig, (ax_bar, ax_heat) = plt.subplots(
        1, 2, figsize=(15, max(5, 0.5 * len(df))), constrained_layout=True,
        gridspec_kw={"width_ratios": [2, 2]},
    )

    n = len(df)
    y = np.arange(n)
    ax_bar.barh(y, df["pct_disagree"], color="lightgray",
                label="disagree (any conf)")
    ax_bar.barh(y, df["pct_overwrite"], color="firebrick",
                label="overwrite (qualifying)")
    ax_bar.set_yticks(y)
    ax_bar.set_yticklabels(df["slide_id"], fontsize=8)
    ax_bar.invert_yaxis()
    ax_bar.set_xlabel("% of valid cortex pixels")
    ax_bar.set_title(f"per-slide relabel rate — model '{args.model}'")
    ax_bar.legend(fontsize=8, loc="lower right")
    ax_bar.grid(axis="x", alpha=0.3)

    # Transition heatmap (row-normalised log counts), pooled across slides.
    K = len(CORTEX_ORDER)
    logT = np.log1p(pooled_trans)
    im = ax_heat.imshow(logT, cmap="magma", aspect="equal")
    ax_heat.set_xticks(range(K))
    ax_heat.set_yticks(range(K))
    ax_heat.set_xticklabels([LAYER_NAMES[L] for L in CORTEX_ORDER])
    ax_heat.set_yticklabels([LAYER_NAMES[L] for L in CORTEX_ORDER])
    ax_heat.set_xlabel("predicted layer")
    ax_heat.set_ylabel("current (annotated) layer")
    ax_heat.set_title("pooled overwrite transitions (log count)")
    for r in range(K):
        for c in range(K):
            v = pooled_trans[r, c]
            if v > 0:
                ax_heat.text(c, r, str(v), ha="center", va="center",
                             fontsize=7,
                             color="white" if logT[r, c] < logT.max() * 0.6
                             else "black")
    fig.colorbar(im, ax=ax_heat, fraction=0.046, label="log(1+count)")

    png_out = APPLIED_DIR / f"apply_summary_{args.model}.png"
    fig.savefig(png_out, dpi=140)
    plt.close(fig)
    print(f"wrote {png_out}")

    # Pooled stats
    tot_valid = df["n_valid"].sum()
    tot_ow = df["n_overwrite"].sum()
    print(f"\npooled: {tot_ow} / {tot_valid} pixels overwritten "
          f"({100*tot_ow/tot_valid:.2f}%) across {len(df)} slides")
    print("\npooled transition matrix (rows=current, cols=predicted):")
    hdr = "       " + "".join(f"{LAYER_NAMES[L]:>7}" for L in CORTEX_ORDER)
    print(hdr)
    for i, L in enumerate(CORTEX_ORDER):
        line = f"{LAYER_NAMES[L]:>5}  " + "".join(
            f"{pooled_trans[i, j]:>7d}" for j in range(K))
        print(line)


if __name__ == "__main__":
    main()
