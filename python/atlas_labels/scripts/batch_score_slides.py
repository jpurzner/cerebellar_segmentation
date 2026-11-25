"""
Batch driver: run the full ribbon-unfold + depth-bands diagnostic on every
sandbox slide, then aggregate per-slide metrics into a scorecard that ranks
slides by their suitability as a depth-band reference.

Outputs:
  outputs/ribbon_unfold/batch_scorecard.csv          — per-slide metrics, sorted
  outputs/ribbon_unfold/batch_scorecard_ranked.png   — visual scorecard
  outputs/ribbon_unfold/batch_run.log                — per-slide stdout
  outputs/ribbon_unfold/batch_errors.log             — tracebacks (if any)
  outputs/ribbon_unfold/<slide_id>/...               — per-slide artefacts
"""
import io
import sys
import time
import traceback
import contextlib
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Make sibling scripts importable.
sys.path.insert(0, str(Path(__file__).parent))
from ribbon_unfold_laplace import main as run_ribbon_unfold, LBL  # noqa: E402
from depth_bands_diagnostic import (                                # noqa: E402
    main as run_diagnostic,
    LAYER_ORDER_PIA_TO_CORE,
    LAYER_NAMES,
)

DEFAULT_SANDBOX = Path(
    "/Users/jpurzner/Dropbox/imaging_analysis/cerebellar_segmentation/"
    "python/matlab_anomaly/sandbox"
)
DEFAULT_OUTPUTS = Path(
    "/Users/jpurzner/Dropbox/Q_research/py_projects/cerebellar_atlas_labels/"
    "outputs/ribbon_unfold"
)
# These globals are set by main() so compute_metrics() can locate per-slide
# artefacts. Keep the legacy capitalised names so existing helpers still work.
SANDBOX: Path = DEFAULT_SANDBOX
OUTPUTS: Path = DEFAULT_OUTPUTS
CSV_OUT = DEFAULT_OUTPUTS / "batch_scorecard.csv"
PNG_OUT = DEFAULT_OUTPUTS / "batch_scorecard_ranked.png"
RUN_LOG = DEFAULT_OUTPUTS / "batch_run.log"
ERR_LOG = DEFAULT_OUTPUTS / "batch_errors.log"


def slide_id_from_path(p: Path) -> str:
    s = p.stem
    for sfx in ("_native_labels", "_labelled"):
        if s.endswith(sfx):
            return s[: -len(sfx)]
    return s


# --------------------------------------------------------------------------- #
# Per-slide metric extraction
# --------------------------------------------------------------------------- #

def compute_metrics(slide_id: str) -> dict:
    slide_dir = OUTPUTS / slide_id
    coords_p = slide_dir / f"{slide_id}_ribbon_coords.npz"
    bands_p = slide_dir / f"{slide_id}_depth_bands.npz"

    coords = np.load(coords_p, allow_pickle=False)
    bands = np.load(bands_p, allow_pickle=False)

    layer = coords["layer"]
    psi = coords["psi"]
    reached = coords["reached"]

    m: dict = {"slide_id": slide_id}
    m["streamline_reach"] = float(reached.mean())

    for L in LAYER_ORDER_PIA_TO_CORE:
        nm = LAYER_NAMES[L]
        sel = (layer == L) & reached
        m[f"n_pix_{nm}"] = int(sel.sum())
        m[f"psi_med_{nm}"] = (float(np.nanmedian(psi[sel]))
                              if sel.sum() > 0 else float("nan"))

    n_bins = int(bands["n_s_bins"])
    m["n_s_bins"] = n_bins
    for L in LAYER_ORDER_PIA_TO_CORE:
        nm = LAYER_NAMES[L]
        cnt = bands[f"count_{nm}"]
        m[f"coverage_{nm}"] = float((cnt > 0).mean())

    # Depth ordering: count of correctly-ordered adjacent layer medians
    # (oEGL < iEGL < ML < PCL < IGL).  Score 0..4.
    meds = [m[f"psi_med_{LAYER_NAMES[L]}"] for L in LAYER_ORDER_PIA_TO_CORE]
    ok = 0
    for a, b in zip(meds[:-1], meds[1:]):
        if np.isfinite(a) and np.isfinite(b) and a < b:
            ok += 1
    m["depth_order_correct"] = ok

    # Interface separations (signed |Δmedian|).
    seps = []
    for a, b in zip(meds[:-1], meds[1:]):
        seps.append(abs(b - a) if (np.isfinite(a) and np.isfinite(b))
                    else float("nan"))
    m["sep_oEGL_iEGL"] = seps[0]
    m["sep_iEGL_ML"] = seps[1]
    m["sep_ML_PCL"] = seps[2]
    m["sep_PCL_IGL"] = seps[3]
    m["sep_total"] = float(np.nansum(seps))

    # Layer-collapse saturation: fraction of PCL / IGL pixels at ψ ≈ 1.
    psi_pcl = psi[(layer == LBL["PCL"]) & reached]
    psi_igl = psi[(layer == LBL["IGL"]) & reached]
    m["pcl_saturation"] = (float((psi_pcl >= 0.99).mean())
                           if len(psi_pcl) else float("nan"))
    m["igl_saturation"] = (float((psi_igl >= 0.99).mean())
                           if len(psi_igl) else float("nan"))

    # oEGL thickness: fraction of arclength bins with thickness > 0.
    thk = coords["oegl_thickness_p95_ds_px"]
    m["oegl_thickness_coverage"] = float((thk > 0).mean())

    return m


def quality_score(row: pd.Series) -> float:
    """Composite reference-quality score, higher = better.

    Weights are tuned to penalise the two main failure modes we saw:
      - sparse cortex labelling (low coverage) — weight 0.35
      - depth-collapsed PCL / IGL (high saturation) — weight 0.20
    """
    cov = np.nanmean(
        [row[f"coverage_{nm}"] for nm in ["oEGL", "iEGL", "ML", "PCL", "IGL"]]
    )
    order_ok = row["depth_order_correct"] / 4.0
    sep = min(1.0, row["sep_total"])
    sat = 1.0 - 0.5 * (
        np.nan_to_num(row["pcl_saturation"], nan=1.0)
        + np.nan_to_num(row["igl_saturation"], nan=1.0)
    )
    reach = row["streamline_reach"]
    return float(
        0.35 * cov + 0.20 * order_ok + 0.15 * sep + 0.20 * sat + 0.10 * reach
    )


# --------------------------------------------------------------------------- #
# Summary figure
# --------------------------------------------------------------------------- #

def render_scorecard(df: pd.DataFrame, out_path: Path):
    n = len(df)
    fig_h = max(12, 0.30 * n + 2)
    fig, axes = plt.subplots(
        1, 2, figsize=(15, fig_h), constrained_layout=True,
        gridspec_kw={"width_ratios": [1.2, 4]},
    )
    ax_score, ax_heat = axes

    # Composite score on the left.
    colors = plt.get_cmap("viridis")(
        (df["quality_score"] - df["quality_score"].min())
        / max(1e-9, df["quality_score"].max() - df["quality_score"].min())
    )
    ax_score.barh(range(n), df["quality_score"], color=colors)
    ax_score.set_yticks(range(n))
    ax_score.set_yticklabels(df["slide_id"], fontsize=7)
    ax_score.invert_yaxis()
    ax_score.set_xlabel("composite quality score (higher = better reference)")
    ax_score.set_title("ranking", loc="left")
    ax_score.grid(axis="x", alpha=0.3)
    ax_score.set_xlim(left=0)

    # Metric heatmap on the right.
    metric_cols = [
        "coverage_oEGL", "coverage_iEGL", "coverage_ML",
        "coverage_PCL", "coverage_IGL",
        "depth_order_correct", "sep_total",
        "streamline_reach",
        "pcl_saturation", "igl_saturation",
        "oegl_thickness_coverage",
    ]
    data = df[metric_cols].to_numpy(dtype=float)
    # Normalise each column to 0–1, then invert "lower-is-better" columns.
    lo = np.nanmin(data, axis=0)
    hi = np.nanmax(data, axis=0)
    rng = np.where(hi - lo < 1e-9, 1.0, hi - lo)
    norm = (data - lo) / rng
    invert = {metric_cols.index(c) for c in ("pcl_saturation", "igl_saturation")}
    for j in invert:
        norm[:, j] = 1.0 - norm[:, j]

    im = ax_heat.imshow(norm, aspect="auto", cmap="RdYlGn",
                        interpolation="nearest", vmin=0, vmax=1)
    ax_heat.set_xticks(range(len(metric_cols)))
    ax_heat.set_xticklabels(metric_cols, rotation=45, ha="right", fontsize=8)
    ax_heat.set_yticks([])
    ax_heat.set_title("per-metric score (green = better; saturation cols "
                      "inverted)", loc="left")
    # annotate with raw values
    for i in range(n):
        for j in range(len(metric_cols)):
            v = data[i, j]
            if np.isfinite(v):
                ax_heat.text(j, i, f"{v:.2f}", ha="center", va="center",
                             fontsize=5.5, color="black")
    fig.colorbar(im, ax=ax_heat, fraction=0.02, label="normalised score")

    fig.suptitle(f"sandbox depth-band scorecard — {n} slides sorted by "
                 f"composite quality",
                 fontsize=12, fontweight="bold")
    fig.savefig(out_path, dpi=130)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #

def main():
    import argparse
    global SANDBOX, OUTPUTS, CSV_OUT, PNG_OUT, RUN_LOG, ERR_LOG
    ap = argparse.ArgumentParser()
    ap.add_argument("--sandbox", type=Path, default=DEFAULT_SANDBOX,
                    help="directory containing *_native_labels.tif files")
    ap.add_argument("--output_root", type=Path, default=DEFAULT_OUTPUTS,
                    help="root directory for batch + per-slide outputs")
    ap.add_argument("--ds", type=int, default=4,
                    help="downsample factor for ribbon_unfold")
    args = ap.parse_args()
    SANDBOX = args.sandbox
    OUTPUTS = args.output_root
    CSV_OUT = OUTPUTS / "batch_scorecard.csv"
    PNG_OUT = OUTPUTS / "batch_scorecard_ranked.png"
    RUN_LOG = OUTPUTS / "batch_run.log"
    ERR_LOG = OUTPUTS / "batch_errors.log"

    # The ribbon-unfold script also writes into a per-slide directory under
    # its own OUTPUT_ROOT constant. Point that at our chosen output root.
    import ribbon_unfold_laplace as ru
    ru.OUTPUT_ROOT = OUTPUTS

    OUTPUTS.mkdir(parents=True, exist_ok=True)
    files = sorted(SANDBOX.glob("*_native_labels.tif"))
    print(f"found {len(files)} slides under {SANDBOX}")
    print(f"output root: {OUTPUTS}")
    print(f"downsample factor: {args.ds}")
    print(f"writing per-slide stdout to {RUN_LOG}")

    errors = []
    t_start = time.time()
    with open(RUN_LOG, "w") as run_log:
        run_log.write(f"# batch started at {time.ctime()}\n\n")
        for i, f in enumerate(files):
            slide_id = slide_id_from_path(f)
            print(f"[{i+1:>2}/{len(files)}] {slide_id} ...", flush=True)
            t0 = time.time()
            buf = io.StringIO()
            try:
                with contextlib.redirect_stdout(buf):
                    run_ribbon_unfold(f, ds=args.ds)
                    coords_p = OUTPUTS / slide_id / f"{slide_id}_ribbon_coords.npz"
                    if coords_p.exists():
                        run_diagnostic(coords_p)
            except Exception as e:
                tb = traceback.format_exc()
                errors.append((slide_id, tb))
                buf.write(f"\nERROR: {e}\n{tb}\n")
                print(f"   !! ERROR: {e}", flush=True)
            elapsed = time.time() - t0
            run_log.write(f"\n=========================================\n"
                          f"=== {slide_id} (elapsed {elapsed:.1f}s) ===\n"
                          f"=========================================\n")
            run_log.write(buf.getvalue())
            run_log.flush()
            print(f"   done in {elapsed:.1f}s  "
                  f"(total {(time.time() - t_start)/60:.1f} min)", flush=True)

    if errors:
        with open(ERR_LOG, "w") as fh:
            for sid, tb in errors:
                fh.write(f"\n=== {sid} ===\n{tb}\n")
        print(f"\n{len(errors)} slides failed; tracebacks in {ERR_LOG}")
    else:
        print("\nall slides processed without errors")

    # --- Aggregate metrics ---
    print("\nAggregating metrics ...")
    rows = []
    for f in files:
        sid = slide_id_from_path(f)
        try:
            rows.append(compute_metrics(sid))
        except Exception as e:
            print(f"  skip metrics for {sid}: {e}")
    df = pd.DataFrame(rows)
    df["quality_score"] = df.apply(quality_score, axis=1)
    df = df.sort_values("quality_score", ascending=False).reset_index(drop=True)
    df.to_csv(CSV_OUT, index=False)
    print(f"\nwrote {CSV_OUT}  ({len(df)} slides scored)")

    # --- Top 5 summary ---
    print("\nTop 5 reference candidates:")
    cols_to_show = ["slide_id", "quality_score",
                    "coverage_oEGL", "coverage_iEGL", "coverage_ML",
                    "coverage_PCL", "coverage_IGL",
                    "depth_order_correct", "sep_total",
                    "streamline_reach",
                    "pcl_saturation", "igl_saturation"]
    print(df[cols_to_show].head(5).to_string(index=False))

    # --- Visual scorecard ---
    if len(df) > 0:
        render_scorecard(df, PNG_OUT)
        print(f"\nwrote {PNG_OUT}")

    print(f"\nbatch done in {(time.time() - t_start)/60:.1f} min")


if __name__ == "__main__":
    main()
