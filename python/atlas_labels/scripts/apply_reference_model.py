"""
Apply a trained depth-band reference model to a target slide.

For every cortex pixel in the target slide, predict the cortical layer
from its (s_norm, ψ) coordinate using a per-layer Gaussian likelihood
(median ± robust std derived from the reference model's per-bin IQR).
Disagreements between the prediction and the current label are flagged;
high-confidence disagreements are optionally overwritten.

Usage:
    python apply_reference_model.py \
        --slide <slide_id>  \
        --model <model_tag> \
        [--conf_thresh 0.5] [--overwrite]

Outputs:
    outputs/reference_models/applied/<slide_id>/
        <slide_id>_predicted.npz       per-pixel pred + confidence
        <slide_id>_filled_labels.tif   downsampled label image w/ overwrites
        <slide_id>_filled_labels_native.tif  nearest-upsampled to native
        <slide_id>_diagnostic.png      4-pane diagnostic figure
"""

import sys
import argparse
from pathlib import Path

import numpy as np
from PIL import Image
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from ribbon_unfold_laplace import (                                # noqa: E402
    COLORS, LBL, load_labels, downsample, label_rgb
)
from build_reference_model import (                                # noqa: E402
    LAYER_NAMES, MODEL_LAYERS, CORTEX_ORDER, MODEL_DIR, OUTPUTS,
)

APPLIED_DIR = Path(
    "/Users/jpurzner/Dropbox/Q_research/py_projects/cerebellar_atlas_labels/"
    "outputs/reference_models/applied"
)
SANDBOX_40X = Path(
    "/Users/jpurzner/Dropbox/imaging_analysis/cerebellar_segmentation/"
    "python/matlab_anomaly/sandbox"
)
SANDBOX_10X = Path(
    "/Users/jpurzner/Dropbox/imaging_analysis/cerebellar_segmentation/"
    "python/matlab_anomaly_10x/sandbox"
)


# --------------------------------------------------------------------------- #
# Predictor
# --------------------------------------------------------------------------- #

_EGL_PAIR = {LBL["iEGL"], LBL["oEGL"]}


def gaussian_log_lik_per_layer(s_norm, psi, model):
    """For each pixel, return log p(ψ | layer, s_norm) for every layer."""
    centres = model["centres"]
    n_bins = int(model["n_s_bins"])
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    bin_idx = np.clip(np.searchsorted(edges, s_norm, side="right") - 1,
                      0, n_bins - 1)

    n_pix = len(psi)
    n_lay = len(CORTEX_ORDER)
    log_lik = np.full((n_pix, n_lay), -np.inf, dtype=np.float64)

    for k, L in enumerate(CORTEX_ORDER):
        prof = model[f"profile_{LAYER_NAMES[L]}"]   # (n_bins, 5)
        med = prof[:, 2]
        iqr = prof[:, 3] - prof[:, 1]
        std = np.where(np.isfinite(iqr), iqr / 1.349, np.nan)
        std = np.where((std > 0.01) & np.isfinite(std), std, 0.01)

        med_pix = med[bin_idx]
        std_pix = std[bin_idx]
        valid = np.isfinite(med_pix)
        z = (psi - med_pix) / std_pix
        ll = -0.5 * z * z - np.log(std_pix * np.sqrt(2 * np.pi))
        log_lik[:, k] = np.where(valid, ll, -np.inf)

    return log_lik


def posterior_and_pred(log_lik):
    """Softmax over layers → posterior, argmax → predicted layer index."""
    # Subtract row-wise max for stability.
    row_max = np.max(log_lik, axis=1, keepdims=True)
    lik = np.exp(log_lik - row_max)
    # Replace -inf rows (all layers invalid) with uniform.
    row_sum = lik.sum(axis=1, keepdims=True)
    bad = (row_sum.flatten() == 0) | ~np.isfinite(row_sum.flatten())
    if bad.any():
        lik[bad] = 1.0
        row_sum[bad] = float(lik.shape[1])
    post = lik / row_sum

    pred_idx = post.argmax(axis=1)
    pred_layer = np.array(CORTEX_ORDER)[pred_idx]
    # Confidence = top posterior; margin = top minus second-best.
    top = post.max(axis=1)
    second = np.partition(post, -2, axis=1)[:, -2]
    margin = top - second
    return pred_layer, top, margin, post


# --------------------------------------------------------------------------- #
# Build the corrected label image
# --------------------------------------------------------------------------- #

def apply_overwrite(orig_labels, cy, cx, pred_layer, qualify, do_overwrite):
    """Return a new label image, overwriting exactly the pixels in `qualify`.

    `qualify` already encodes every decision rule (disagreement, confidence,
    margin, EGL-swap block), so this function must NOT re-derive any
    criterion — it just stamps the predicted label at the qualifying pixels.
    """
    out = orig_labels.copy()
    changed = np.zeros(len(pred_layer), dtype=bool)
    if not do_overwrite:
        return out, changed
    out[cy[qualify], cx[qualify]] = pred_layer[qualify]
    changed = qualify.copy()
    return out, changed


# --------------------------------------------------------------------------- #
# Native-resolution upsampling (nearest-neighbour for labels)
# --------------------------------------------------------------------------- #

def nearest_upsample(label_ds: np.ndarray, ds_factor: int,
                     native_shape: tuple) -> np.ndarray:
    H, W = native_shape
    out = np.zeros((H, W), dtype=label_ds.dtype)
    H_ds, W_ds = label_ds.shape
    # Map each native pixel to its source ds pixel index by floor-division.
    ys = np.minimum(np.arange(H) // ds_factor, H_ds - 1)
    xs = np.minimum(np.arange(W) // ds_factor, W_ds - 1)
    out = label_ds[np.ix_(ys, xs)]
    return out


# --------------------------------------------------------------------------- #
# Diagnostic figure
# --------------------------------------------------------------------------- #

def render_diagnostic(lbl_ds, pred_lbl_ds, changed_mask_ds, conf_ds,
                      slide_id, model_tag, conf_thresh, n_changed,
                      out_path):
    fig, axes = plt.subplots(2, 2, figsize=(15, 13), constrained_layout=True)

    # 1. Original
    ax = axes[0, 0]
    ax.imshow(label_rgb(lbl_ds))
    ax.set_title("original labels (downsampled)")
    ax.set_axis_off()

    # 2. After overwrite
    ax = axes[0, 1]
    ax.imshow(label_rgb(pred_lbl_ds))
    ax.set_title(f"after overwrite ({n_changed} pixels changed; "
                 f"conf ≥ {conf_thresh})")
    ax.set_axis_off()

    # 3. Disagreement / changed mask
    ax = axes[1, 0]
    ax.imshow(label_rgb(lbl_ds), alpha=0.4)
    if changed_mask_ds.any():
        ax.imshow(np.ma.masked_where(~changed_mask_ds, changed_mask_ds),
                  cmap="autumn", alpha=0.85, interpolation="nearest")
    ax.set_title("changed pixels (yellow/red)")
    ax.set_axis_off()

    # 4. Confidence map
    ax = axes[1, 1]
    im = ax.imshow(conf_ds, cmap="magma", vmin=0, vmax=1)
    ax.set_title("max posterior (1 = certain)")
    ax.set_axis_off()
    fig.colorbar(im, ax=ax, fraction=0.025)

    fig.suptitle(f"depth-model apply — slide={slide_id}  model={model_tag}",
                 fontsize=12, fontweight="bold")
    fig.savefig(out_path, dpi=140)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--slide", required=True,
                    help="slide_id (e.g. 2018_05_22_s1_2_p27-0021_iter1deep)")
    ap.add_argument("--model", default="best_of_section",
                    help="model tag in outputs/reference_models/")
    ap.add_argument("--conf_thresh", type=float, default=0.7,
                    help="overwrite if max posterior ≥ this")
    ap.add_argument("--margin_thresh", type=float, default=0.3,
                    help="ALSO require posterior margin ≥ this for overwrite")
    ap.add_argument("--prior_boost", type=float, default=1.5,
                    help="log-prior boost added to current label's likelihood "
                         "(higher = more conservative; 1.5 ≈ 4.5x linear bias)")
    ap.add_argument("--allow_egl_swap", action="store_true", default=False,
                    help="permit iEGL↔oEGL overwrites (off by default — those "
                         "two layers are not depth-distinguishable)")
    ap.add_argument("--overwrite", action="store_true", default=True,
                    help="actually overwrite (default true; --no-overwrite "
                         "to suppress)")
    ap.add_argument("--no-overwrite", dest="overwrite", action="store_false")
    ap.add_argument("--ribbon_root", type=Path, default=None,
                    help="directory containing per-slide ribbon_coords.npz "
                         "(auto-detect 40x vs 10x if omitted)")
    ap.add_argument("--sandbox", type=Path, default=None,
                    help="directory containing *_native_labels.tif "
                         "(auto-detect 40x vs 10x if omitted)")
    args = ap.parse_args()

    slide_id = args.slide
    model_path = MODEL_DIR / f"depth_model_{args.model}.npz"
    assert model_path.exists(), f"missing model: {model_path}"

    # Auto-detect ribbon_coords + native label location across 40x/10x.
    candidate_roots = ([args.ribbon_root] if args.ribbon_root else
                       [OUTPUTS,
                        Path("/Users/jpurzner/Dropbox/Q_research/py_projects/"
                             "cerebellar_atlas_labels/outputs/ribbon_unfold_10x")])
    coords_path = None
    for r in candidate_roots:
        p = r / slide_id / f"{slide_id}_ribbon_coords.npz"
        if p.exists():
            coords_path = p
            break
    assert coords_path is not None, (
        f"missing ribbon coords for {slide_id}; tried: "
        + ", ".join(str(r / slide_id) for r in candidate_roots)
    )

    candidate_sandboxes = ([args.sandbox] if args.sandbox
                           else [SANDBOX_40X, SANDBOX_10X])
    sandbox_tif = None
    for sb in candidate_sandboxes:
        p = sb / f"{slide_id}_native_labels.tif"
        if p.exists():
            sandbox_tif = p
            break
    assert sandbox_tif is not None, (
        f"no native_labels.tif for {slide_id} in any known sandbox"
    )

    print(f"slide:  {slide_id}")
    print(f"model:  {model_path}")
    print(f"coords: {coords_path}")

    model = np.load(model_path, allow_pickle=False)
    coords = np.load(coords_path, allow_pickle=False)

    cy = coords["cy"]
    cx = coords["cx"]
    layer_cur = coords["layer"]
    psi = coords["psi"]
    s = coords["s"]
    reached = coords["reached"]
    s_max = float(coords["s_max"])
    ds_factor = int(coords["downsample"])

    valid = (reached
             & np.isfinite(s)
             & np.isfinite(psi)
             & np.isin(layer_cur, CORTEX_ORDER))
    print(f"valid cortex pixels: {int(valid.sum())} / {len(psi)}")

    cy_v = cy[valid]; cx_v = cx[valid]
    s_norm_v = s[valid] / s_max
    psi_v = psi[valid]
    layer_cur_v = layer_cur[valid]

    log_lik = gaussian_log_lik_per_layer(s_norm_v, psi_v, model)

    # Prior boost: add log-boost to the CURRENT label's log-likelihood so the
    # alternative has to be substantially more likely to dethrone it.
    if args.prior_boost > 0:
        for k, L in enumerate(CORTEX_ORDER):
            log_lik[layer_cur_v == L, k] += args.prior_boost

    pred, post_top, margin, post = posterior_and_pred(log_lik)

    # Decision logic
    qualify = ((pred != layer_cur_v)
               & (post_top >= args.conf_thresh)
               & (margin >= args.margin_thresh))
    # Block iEGL↔oEGL swaps unless explicitly allowed.
    if not args.allow_egl_swap:
        egl_swap = (np.isin(layer_cur_v, list(_EGL_PAIR))
                    & np.isin(pred, list(_EGL_PAIR)))
        qualify = qualify & ~egl_swap
    n_disagree = int((pred != layer_cur_v).sum())
    n_overwrite = int(qualify.sum())
    print(f"disagreements (any conf):       {n_disagree}  "
          f"({100*n_disagree/len(pred):.1f}%)")
    print(f"qualifying for overwrite:       {n_overwrite}  "
          f"({100*n_overwrite/len(pred):.1f}%)")

    # Breakdown of overwrite transitions
    print("\noverwrite transitions (current → predicted):")
    transitions = {}
    for cur, prd in zip(layer_cur_v[qualify], pred[qualify]):
        key = (int(cur), int(prd))
        transitions[key] = transitions.get(key, 0) + 1
    for (a, b), n in sorted(transitions.items(), key=lambda kv: -kv[1]):
        print(f"  {LAYER_NAMES[a]:>4} ({a}) → {LAYER_NAMES[b]:>4} ({b}):  {n}")

    # Reload original label image (downsampled, just for visualisation).
    # We also need the native shape to upsample.
    native = load_labels(sandbox_tif)
    native_shape = native.shape
    lbl_ds = downsample(native, ds_factor)
    print(f"\nlabels: native shape={native_shape}, ds shape={lbl_ds.shape}")

    # Apply overwrite on the downsampled image, using the full qualify mask
    # (disagreement + confidence + margin + EGL-swap block) so the written
    # labels match the reported "qualifying for overwrite" count exactly.
    pred_lbl_ds, changed = apply_overwrite(
        lbl_ds, cy_v, cx_v, pred, qualify, args.overwrite
    )

    # Build a changed-mask raster (ds resolution).
    changed_mask_ds = np.zeros(lbl_ds.shape, dtype=bool)
    changed_mask_ds[cy_v[changed], cx_v[changed]] = True

    # Build a confidence raster (ds resolution) for the diagnostic.
    conf_ds = np.zeros(lbl_ds.shape, dtype=np.float32)
    conf_ds[cy_v, cx_v] = post_top.astype(np.float32)

    # Output dir
    out_dir = APPLIED_DIR / slide_id
    out_dir.mkdir(parents=True, exist_ok=True)

    # Save per-pixel predictions npz
    np.savez_compressed(
        out_dir / f"{slide_id}_predicted.npz",
        cy=cy_v.astype(np.int32),
        cx=cx_v.astype(np.int32),
        s_norm=s_norm_v.astype(np.float32),
        psi=psi_v.astype(np.float32),
        layer_current=layer_cur_v.astype(np.uint8),
        layer_predicted=pred.astype(np.uint8),
        posterior_top=post_top.astype(np.float32),
        posterior_margin=margin.astype(np.float32),
        posterior_full=post.astype(np.float32),
        overwrote=changed,
        conf_thresh=float(args.conf_thresh),
        margin_thresh=float(args.margin_thresh),
        model_tag=args.model,
    )

    # Save downsampled filled label image
    Image.fromarray(pred_lbl_ds.astype(np.uint8)).save(
        out_dir / f"{slide_id}_filled_labels.tif"
    )

    # Upsample (nearest) to native resolution and save
    pred_lbl_native = nearest_upsample(pred_lbl_ds, ds_factor, native_shape)
    # Re-stamp anomaly + DCN from the ORIGINAL native (they were untouched
    # by the model and should not be lost during the round-trip).
    pred_lbl_native = np.where(
        (native == LBL["ANOMALY"]) | (native == LBL["DCN"]),
        native, pred_lbl_native
    )
    Image.fromarray(pred_lbl_native.astype(np.uint8)).save(
        out_dir / f"{slide_id}_filled_labels_native.tif"
    )

    render_diagnostic(
        lbl_ds, pred_lbl_ds, changed_mask_ds, conf_ds,
        slide_id, args.model, args.conf_thresh,
        n_overwrite, out_dir / f"{slide_id}_diagnostic.png"
    )
    print(f"\nwrote outputs to {out_dir}")


if __name__ == "__main__":
    main()
