"""
Use the MATLAB segments output for s1_2 as the EGL ground-truth, then run the
v4-style multivariate matched-filter kernel on the s1_2 input image and compare.

Step 1: load segments.tif, identify EGL labels by cluster-color analysis.
Step 2: run kernel (deep-pair + outer-band) on s1_2.
Step 3: compute IoU, recall, precision, and per-pixel overlay.

Baseline comparison: which fraction of MATLAB EGL does our kernel hit?
                     which fraction of kernel pixels are inside MATLAB EGL?
"""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage as ndi
import tifffile
from scipy.signal import fftconvolve
from skimage import exposure, filters, morphology
from skimage.morphology import disk
from sklearn.cluster import KMeans

INPUT = Path(
    "/Users/jpurzner/Dropbox/images/edu_repeat/p27/"
    "2018_05_22_s1_2_p27-0021/2018_05_22_s1_2_p27-0021_fused_crop.tif"
)
SEG_TIF = Path(
    "/Users/jpurzner/Dropbox/imaging_analysis/cerebellar_segmentation/_test_run/"
    "20x/2018_05_22_s1_2_p27-0021/"
    "2018_05_22_s1_2_p27-0021_fused_crop_segments.tif"
)
A_IEGL_TIF = Path(
    "/Users/jpurzner/Dropbox/imaging_analysis/cerebellar_segmentation/_test_run/"
    "20x/2018_05_22_s1_2_p27-0021/"
    "2018_05_22_s1_2_p27-0021_fused_crop_a_iegl_only.tif"
)

OUT = Path(__file__).resolve().parent / "figs" / "v_baseline_s1_2"
OUT.mkdir(parents=True, exist_ok=True)
PX_UM = 0.5119049
CH_P27, CH_NEUN, CH_DAPI = 0, 1, 2

# Pair-template params (from v4)
KERNEL_TICK_UM    = 1.0
PAIR_HALFLEN_UM   = 70.0
GAP_HALFWIDTH_UM  = 12.0
LOBE_HALFWIDTH_UM = 18.0
ML_FLANK_WIDTH    = 20.0
KERNEL_THICK_UM   = 5.0
N_THETA           = 12
W_DAPI = 1.0
W_P27  = 1.5

PAIR_THRESH_PCT       = 95
EGL_BAND_UM           = 50.0     # outer-band geometric prior
EGL_INTENSITY_PCT     = 88       # within outer band, top X% by intensity


def to_unit(img):
    img = img.astype(np.float32)
    lo, hi = img.min(), img.max()
    return (img - lo) / max(hi - lo, 1e-9)


def clahe(img, clip=0.01):
    return exposure.equalize_adapthist(img, clip_limit=clip, nbins=256)


def cerebellum_mask(p27, dapi):
    s = (p27 + dapi) / 2; s = s / s.max()
    raw = s > 0.10
    closed = morphology.closing(raw, footprint=disk(int(round(5/PX_UM))))
    filled = ndi.binary_fill_holes(closed)
    lbl, n = ndi.label(filled); sizes = np.bincount(lbl.ravel()); sizes[0] = 0
    if n == 0: return filled
    return lbl == sizes.argmax()


def extract_egl_from_segments(seg_rgb):
    """Identify iEGL and oEGL pixels in the MATLAB segments TIF.

    The MATLAB code writes label2rgb(set_bin, 'jet', 'k') where set_bin has values
    0-7 (background, all_cerebellum, iEGL, oEGL, IGL, ML, DWL, PCL). Empirically
    the colors are jet(8) palette:
       label 2 (iEGL) -> (0, 0, 255)   blue
       label 3 (oEGL) -> (0, 128, 255) azure
       label 4 (IGL)  -> (0, 255, 255) cyan
       label 5 (ML)   -> (128,255,128) lt green
       label 6 (DWL)  -> (255,128, 0)  orange
       label 7 (PCL)  -> (255,  0, 0)  red

    We identify each layer by RGB color match (tolerance 16 per channel).
    """
    h, w, _ = seg_rgb.shape
    s = seg_rgb.astype(int)
    def col(target):
        return np.all(np.abs(s - np.array(target)) < 16, axis=-1)
    # Empirical mapping verified via overlap with a_iegl_only.tif:
    #   azure (0,128,255) is iEGL (79% overlap with a_iegl_only)
    #   cyan  (0,255,255) is oEGL (outermost thin band)
    #   dark blue (0,0,255) is all_cerebellum (unlabelled fallback inside mask)
    #   red is PCL (small scattered)
    masks = {
        "iEGL": col([0, 128, 255]),
        "oEGL": col([0, 255, 255]),
        "IGL":  col([255,255,  0]),
        "ML":   col([128,255,128]),
        "DWL":  col([255,128,  0]),
        "PCL":  col([255,  0,  0]),
        "all_cerebellum": col([0,   0, 255]),
    }
    for name, m in masks.items():
        print(f"  {name:14s}  {m.sum():>9d} px  {m.sum()*(PX_UM/1000)**2:.3f} mm^2")
    egl = masks["iEGL"] | masks["oEGL"]
    return masks, egl


def make_pair_template():
    half_n = int(PAIR_HALFLEN_UM / KERNEL_TICK_UM)
    d = np.arange(-half_n, half_n + 1) * KERNEL_TICK_UM
    def build(lobe_w):
        t = np.zeros_like(d, dtype=np.float32)
        t[np.abs(d) <= GAP_HALFWIDTH_UM] = -1.0
        for sign in (-1, +1):
            ctr = sign * (GAP_HALFWIDTH_UM + LOBE_HALFWIDTH_UM)
            t[(d >= ctr - lobe_w) & (d <= ctr + lobe_w)] = +1.0
        for sign in (-1, +1):
            base = sign * (GAP_HALFWIDTH_UM + LOBE_HALFWIDTH_UM + lobe_w)
            if sign > 0: t[(d > base) & (d <= base + ML_FLANK_WIDTH)] = -0.7
            else:        t[(d < base) & (d >= base - ML_FLANK_WIDTH)] = -0.7
        t -= t.mean(); t /= max(np.sqrt((t**2).sum()), 1e-9)
        return t
    dapi_t = build(LOBE_HALFWIDTH_UM)
    p27_t  = build(LOBE_HALFWIDTH_UM - 4)
    thick = max(1, int(round(KERNEL_THICK_UM / KERNEL_TICK_UM)))
    def to_pix(t):
        K = np.tile(t.reshape(-1, 1), (1, thick))
        K = ndi.zoom(K, (1.0/PX_UM, 1.0/PX_UM), order=1).astype(np.float32)
        if K.shape[0] % 2 == 0: K = K[:-1, :]
        if K.shape[1] % 2 == 0: K = K[:, :-1]
        return K
    return to_pix(dapi_t), to_pix(p27_t)


def matched_pair_response(dapi_sm, p27_sm, K_d, K_p):
    H, W = dapi_sm.shape
    best = np.full((H, W), -np.inf, dtype=np.float32)
    d32 = dapi_sm.astype(np.float32); p32 = p27_sm.astype(np.float32)
    for i in range(N_THETA):
        ang = i * 180.0 / N_THETA
        Kd = ndi.rotate(K_d, ang, order=1, reshape=True, mode="constant",
                        cval=0.0).astype(np.float32)
        Kp = ndi.rotate(K_p, ang, order=1, reshape=True, mode="constant",
                        cval=0.0).astype(np.float32)
        Kd -= Kd.mean(); Kp -= Kp.mean()
        rd = fftconvolve(d32, Kd[::-1, ::-1], mode="same")
        rp = fftconvolve(p32, Kp[::-1, ::-1], mode="same")
        np.maximum(best, W_DAPI*rd + W_P27*rp, out=best)
    return best


def main():
    print(f"loading input  {INPUT.name}")
    stack = tifffile.imread(INPUT)
    if stack.shape[0] not in (3, 4) and stack.shape[-1] in (3, 4):
        stack = np.moveaxis(stack, -1, 0)
    p27  = clahe(to_unit(stack[CH_P27]))
    neun = clahe(to_unit(stack[CH_NEUN]))
    dapi = clahe(to_unit(stack[CH_DAPI]))
    rgb = np.stack([p27, neun, dapi], axis=-1); rgb = rgb / max(rgb.max(), 1e-9)
    H, W = p27.shape

    print(f"loading segments {SEG_TIF.name}")
    seg = tifffile.imread(SEG_TIF)
    print(f"  segments shape {seg.shape}, dtype {seg.dtype}")
    if seg.shape[:2] != (H, W):
        print(f"  WARN: shape mismatch input ({H},{W}) vs seg {seg.shape[:2]}")

    print("extracting EGL from segments via color clustering")
    masks_by_color, egl_truth = extract_egl_from_segments(seg)
    egl_area_mm2 = egl_truth.sum() * (PX_UM/1000)**2
    print(f"  EGL ground truth area = {egl_area_mm2:.3f} mm^2")

    # Try the dedicated a_iegl_only.tif as a tighter iEGL-only reference
    if A_IEGL_TIF.exists():
        a_iegl = tifffile.imread(A_IEGL_TIF)
        print(f"  a_iegl tif shape {a_iegl.shape}, dtype {a_iegl.dtype}")
        if a_iegl.ndim == 3:
            a_iegl = a_iegl.sum(axis=-1)
        iegl_truth = a_iegl > 0
        print(f"  iEGL-only area = {iegl_truth.sum()*(PX_UM/1000)**2:.3f} mm^2")
    else:
        iegl_truth = None

    print("cerebellum mask (python)")
    mask = cerebellum_mask(p27, dapi)
    print(f"  cerebellum area = {mask.sum()*(PX_UM/1000)**2:.3f} mm^2")

    print("kernel: deep-pair response (FFT)")
    pre_um = 4.0
    dapi_sm = filters.gaussian(dapi, sigma=pre_um/PX_UM, preserve_range=True) * mask
    p27_sm  = filters.gaussian(p27,  sigma=pre_um/PX_UM, preserve_range=True) * mask
    K_d, K_p = make_pair_template()
    pair_resp = matched_pair_response(dapi_sm, p27_sm, K_d, K_p) * mask
    inside = pair_resp[mask]
    pair_thr = np.percentile(inside, PAIR_THRESH_PCT)
    pair_centerlines = (pair_resp > pair_thr) & mask

    # convert centerlines to actual EGL pixels: dilate by lobe distance + intersect
    # with high-DAPI tissue
    lobe_dist_um = GAP_HALFWIDTH_UM + LOBE_HALFWIDTH_UM
    dilate_r = int(round(lobe_dist_um / PX_UM))
    pair_band = morphology.binary_dilation(pair_centerlines,
                                            footprint=disk(dilate_r)) & mask
    # within the band, restrict to high-intensity (DAPI+p27 above local median)
    dapi_p27 = (dapi_sm + p27_sm) / 2
    in_band_intensity = dapi_p27[pair_band]
    if len(in_band_intensity) > 0:
        thr_band = np.percentile(in_band_intensity, 50)
        pair_egl = pair_band & (dapi_p27 > thr_band)
    else:
        pair_egl = pair_band

    print("kernel: outer geometric band + intensity refinement")
    dist_to_bg_um = ndi.distance_transform_edt(mask) * PX_UM
    band_outer = mask & (dist_to_bg_um <= EGL_BAND_UM)
    inside = dapi_p27[band_outer]
    thr_outer = np.percentile(inside, EGL_INTENSITY_PCT)
    outer_egl = band_outer & (dapi_p27 > thr_outer)

    egl_kernel = pair_egl | outer_egl
    egl_kernel = morphology.remove_small_objects(egl_kernel, min_size=int(80/PX_UM**2))
    egl_kernel = morphology.closing(egl_kernel, disk(int(round(3/PX_UM))))

    # === metrics ===
    inter = (egl_kernel & egl_truth).sum()
    union = (egl_kernel | egl_truth).sum()
    precision = inter / max(egl_kernel.sum(), 1)
    recall    = inter / max(egl_truth.sum(),  1)
    iou       = inter / max(union, 1)
    print(f"\n=== kernel vs MATLAB EGL ground truth ===")
    print(f"  GT area       = {egl_truth.sum()*(PX_UM/1000)**2:.3f} mm^2")
    print(f"  kernel area   = {egl_kernel.sum()*(PX_UM/1000)**2:.3f} mm^2")
    print(f"  intersection  = {inter*(PX_UM/1000)**2:.3f} mm^2")
    print(f"  precision     = {precision:.3f}  (kernel pixels in GT)")
    print(f"  recall        = {recall:.3f}     (GT pixels caught by kernel)")
    print(f"  IoU           = {iou:.3f}")

    # === plots ===
    print("plot 01: side-by-side overview")
    rgb_dim = rgb * 0.45
    fig, ax = plt.subplots(2, 3, figsize=(20, 12))
    ax[0, 0].imshow(rgb); ax[0, 0].set_title("RGB"); ax[0, 0].axis("off")
    ax[0, 1].imshow(rgb_dim)
    ax[0, 1].contour(egl_truth, levels=[0.5], colors=["lime"], linewidths=0.5)
    ax[0, 1].set_title(f"MATLAB EGL ground truth ({egl_area_mm2:.2f} mm^2)")
    ax[0, 1].axis("off")
    ax[0, 2].imshow(rgb_dim)
    ax[0, 2].contour(egl_kernel, levels=[0.5], colors=["yellow"], linewidths=0.5)
    ax[0, 2].set_title(f"kernel EGL ({egl_kernel.sum()*(PX_UM/1000)**2:.2f} mm^2)")
    ax[0, 2].axis("off")
    # diff overlay: green = TP, red = FN (missed), yellow = FP (extra)
    tp = egl_truth & egl_kernel
    fn = egl_truth & ~egl_kernel
    fp = egl_kernel & ~egl_truth
    diff = np.zeros((H, W, 3), dtype=np.float32)
    diff[..., 0] = (rgb_dim[..., 0] + fp.astype(np.float32) * 0.9 +
                    fn.astype(np.float32) * 0.9)
    diff[..., 1] = (rgb_dim[..., 1] + tp.astype(np.float32) * 0.9 +
                    fp.astype(np.float32) * 0.7)
    diff[..., 2] = rgb_dim[..., 2]
    diff = np.clip(diff, 0, 1)
    ax[1, 0].imshow(diff)
    ax[1, 0].set_title(f"agreement: TP green, FN red, FP yellow  IoU={iou:.2f}")
    ax[1, 0].axis("off")
    ax[1, 1].imshow(np.where(mask, pair_resp, np.nan), cmap="hot")
    ax[1, 1].set_title("pair-kernel response")
    ax[1, 1].axis("off")
    ax[1, 2].imshow(rgb_dim)
    ax[1, 2].contour(pair_egl,  levels=[0.5], colors=["red"],   linewidths=0.5,
                      alpha=0.85)
    ax[1, 2].contour(outer_egl, levels=[0.5], colors=["yellow"], linewidths=0.4,
                      alpha=0.7)
    ax[1, 2].set_title("kernel decomposition: red=pair-EGL, yellow=outer-band")
    ax[1, 2].axis("off")
    fig.tight_layout(); fig.savefig(OUT / "01_overview.png", dpi=130); plt.close(fig)

    # zoom into deep folds
    cy, cx = np.array(np.where(mask)).mean(axis=1).astype(int)
    half = 800
    y0, y1 = max(0, cy - half), min(H, cy + half)
    x0, x1 = max(0, cx - half), min(W, cx + half)
    fig, ax = plt.subplots(1, 3, figsize=(22, 8))
    ax[0].imshow(rgb[y0:y1, x0:x1]); ax[0].set_title("RGB (zoom)"); ax[0].axis("off")
    ax[1].imshow(rgb_dim[y0:y1, x0:x1])
    ax[1].contour(egl_truth[y0:y1, x0:x1], levels=[0.5], colors=["lime"], linewidths=0.7)
    ax[1].set_title("MATLAB EGL GT"); ax[1].axis("off")
    ax[2].imshow(diff[y0:y1, x0:x1])
    ax[2].set_title(f"TP green, FN red, FP yellow (IoU={iou:.2f})"); ax[2].axis("off")
    fig.tight_layout(); fig.savefig(OUT / "02_zoom.png", dpi=140); plt.close(fig)

    print(f"done. figures in {OUT}")


if __name__ == "__main__":
    main()
