"""Parameter sweep against the s1_2 baseline.

Cache the expensive matched-filter response on disk, then sweep cheap thresholds
and morphology to find the best (precision, recall, IoU) tradeoff.
"""
from pathlib import Path
import json
import numpy as np
import scipy.ndimage as ndi
import tifffile
import matplotlib.pyplot as plt
from scipy.signal import fftconvolve
from skimage import exposure, filters, morphology
from skimage.morphology import disk

INPUT = Path(
    "/Users/jpurzner/Dropbox/images/edu_repeat/p27/"
    "2018_05_22_s1_2_p27-0021/2018_05_22_s1_2_p27-0021_fused_crop.tif"
)
SEG_TIF = Path(
    "/Users/jpurzner/Dropbox/imaging_analysis/cerebellar_segmentation/_test_run/"
    "20x/2018_05_22_s1_2_p27-0021/"
    "2018_05_22_s1_2_p27-0021_fused_crop_segments.tif"
)
OUT = Path(__file__).resolve().parent / "figs" / "v_param_sweep"
OUT.mkdir(parents=True, exist_ok=True)
CACHE = OUT / "cache.npz"
PX_UM = 0.5119049
CH_P27, CH_NEUN, CH_DAPI = 0, 1, 2

KERNEL_TICK_UM = 1.0
PAIR_HALFLEN_UM = 70.0
GAP_HALFWIDTH_UM = 12.0
LOBE_HALFWIDTH_UM = 18.0
ML_FLANK_WIDTH = 20.0
KERNEL_THICK_UM = 5.0
N_THETA = 12
W_DAPI = 1.0
W_P27 = 1.5


def to_unit(im):
    im = im.astype(np.float32); lo, hi = im.min(), im.max()
    return (im - lo) / max(hi - lo, 1e-9)


def clahe(im, clip=0.01):
    return exposure.equalize_adapthist(im, clip_limit=clip, nbins=256)


def cerebellum_mask(p27, dapi):
    s = (p27 + dapi) / 2; s = s / s.max()
    raw = s > 0.10
    closed = morphology.closing(raw, footprint=disk(int(round(5/PX_UM))))
    filled = ndi.binary_fill_holes(closed)
    lbl, n = ndi.label(filled); sizes = np.bincount(lbl.ravel()); sizes[0] = 0
    return (lbl == sizes.argmax()) if n else filled


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
    dt = build(LOBE_HALFWIDTH_UM)
    pt = build(LOBE_HALFWIDTH_UM - 4)
    thick = max(1, int(round(KERNEL_THICK_UM / KERNEL_TICK_UM)))
    def to_pix(t):
        K = np.tile(t.reshape(-1, 1), (1, thick))
        K = ndi.zoom(K, (1.0/PX_UM, 1.0/PX_UM), order=1).astype(np.float32)
        if K.shape[0] % 2 == 0: K = K[:-1, :]
        if K.shape[1] % 2 == 0: K = K[:, :-1]
        return K
    return to_pix(dt), to_pix(pt)


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


def extract_egl_gt(seg_rgb):
    s = seg_rgb.astype(int)
    def col(target): return np.all(np.abs(s - np.array(target)) < 16, axis=-1)
    iEGL = col([0, 128, 255])
    oEGL = col([0, 255, 255])
    return iEGL | oEGL, iEGL, oEGL


def expensive_compute():
    print("loading input")
    stack = tifffile.imread(INPUT)
    if stack.shape[0] not in (3,4) and stack.shape[-1] in (3,4):
        stack = np.moveaxis(stack, -1, 0)
    p27 = clahe(to_unit(stack[CH_P27]))
    dapi = clahe(to_unit(stack[CH_DAPI]))
    mask = cerebellum_mask(p27, dapi)
    pre_um = 4.0
    dapi_sm = filters.gaussian(dapi, sigma=pre_um/PX_UM, preserve_range=True) * mask
    p27_sm  = filters.gaussian(p27,  sigma=pre_um/PX_UM, preserve_range=True) * mask
    print("computing matched-filter response (FFT)")
    K_d, K_p = make_pair_template()
    pair_resp = matched_pair_response(dapi_sm, p27_sm, K_d, K_p) * mask
    dist_to_bg_um = ndi.distance_transform_edt(mask).astype(np.float32) * PX_UM
    dapi_p27 = (dapi_sm + p27_sm) / 2

    seg = tifffile.imread(SEG_TIF)
    egl_gt, iegl_gt, oegl_gt = extract_egl_gt(seg)

    print("caching")
    np.savez_compressed(CACHE,
        p27=p27.astype(np.float32),
        dapi=dapi.astype(np.float32),
        mask=mask,
        pair_resp=pair_resp,
        dist_to_bg_um=dist_to_bg_um,
        dapi_p27=dapi_p27.astype(np.float32),
        egl_gt=egl_gt,
        iegl_gt=iegl_gt,
        oegl_gt=oegl_gt,
    )
    return dict(p27=p27, dapi=dapi, mask=mask, pair_resp=pair_resp,
                dist_to_bg_um=dist_to_bg_um, dapi_p27=dapi_p27,
                egl_gt=egl_gt, iegl_gt=iegl_gt, oegl_gt=oegl_gt)


def kernel_egl(d, pair_thr_pct, band_um, band_int_pct, dilate_um=12):
    """Build EGL mask from cached arrays + cheap params."""
    mask = d["mask"]
    pair_resp = d["pair_resp"]
    inside = pair_resp[mask]
    pair_thr = np.percentile(inside, pair_thr_pct)
    pair_centerlines = (pair_resp > pair_thr) & mask
    dilate_r = max(1, int(round(dilate_um / PX_UM)))
    pair_band = morphology.dilation(pair_centerlines, footprint=disk(dilate_r)) & mask
    in_band_int = d["dapi_p27"][pair_band]
    if len(in_band_int) > 0:
        thr_band = np.percentile(in_band_int, 50)
        pair_egl = pair_band & (d["dapi_p27"] > thr_band)
    else:
        pair_egl = pair_band

    band_outer = mask & (d["dist_to_bg_um"] <= band_um)
    inside_int = d["dapi_p27"][band_outer]
    thr_int = np.percentile(inside_int, band_int_pct)
    outer_egl = band_outer & (d["dapi_p27"] > thr_int)

    egl = pair_egl | outer_egl
    egl = morphology.remove_small_objects(egl, min_size=int(80/PX_UM**2))
    egl = morphology.closing(egl, disk(int(round(3/PX_UM))))
    return egl, pair_egl, outer_egl


def metrics(pred, gt):
    inter = (pred & gt).sum()
    pred_n = pred.sum(); gt_n = gt.sum()
    union = pred_n + gt_n - inter
    p = inter / max(pred_n, 1); r = inter / max(gt_n, 1)
    iou = inter / max(union, 1)
    return dict(precision=float(p), recall=float(r), iou=float(iou),
                pred_mm2=float(pred_n*(PX_UM/1000)**2),
                gt_mm2=float(gt_n*(PX_UM/1000)**2),
                inter_mm2=float(inter*(PX_UM/1000)**2))


def main():
    if CACHE.exists():
        print(f"loading cache {CACHE}")
        z = np.load(CACHE)
        d = {k: z[k] for k in z.files}
    else:
        d = expensive_compute()
    print(f"GT EGL area = {d['egl_gt'].sum()*(PX_UM/1000)**2:.3f} mm^2")
    print(f"cerebellum  = {d['mask'].sum()*(PX_UM/1000)**2:.3f} mm^2")

    # === sweep ===
    pair_pcts  = [90, 95, 97]
    band_ums   = [30, 40, 50]
    band_pcts  = [85, 90, 95]
    print("\n  pair_pct  band_um  band_pct  pred_mm2  precision  recall   IoU")
    results = []
    for pp in pair_pcts:
        for bu in band_ums:
            for bp in band_pcts:
                egl, _, _ = kernel_egl(d, pp, bu, bp)
                m = metrics(egl, d["egl_gt"])
                results.append({"pair_pct": pp, "band_um": bu, "band_pct": bp, **m})
                print(f"  {pp:3d}     {bu:3.0f}     {bp:3d}      "
                      f"{m['pred_mm2']:5.2f}     {m['precision']:.3f}     "
                      f"{m['recall']:.3f}    {m['iou']:.3f}")

    # best by IoU
    best = max(results, key=lambda r: r["iou"])
    print(f"\nBest IoU: {best['iou']:.3f} at "
          f"pair_pct={best['pair_pct']} band_um={best['band_um']} "
          f"band_pct={best['band_pct']}  (P={best['precision']:.3f} R={best['recall']:.3f})")
    with open(OUT / "sweep_results.json", "w") as f:
        json.dump(results, f, indent=2)

    # plot the best vs current parameters
    egl_best, _, _ = kernel_egl(d, best["pair_pct"], best["band_um"], best["band_pct"])
    egl_curr, _, _ = kernel_egl(d, 95, 50, 88)  # current default (poor IoU=0.38)
    rgb = np.stack([d["p27"], np.zeros_like(d["p27"]), d["dapi"]], axis=-1)
    rgb /= max(rgb.max(), 1e-9); rgb_dim = rgb * 0.45

    def diff_overlay(pred, gt, base):
        H, W = pred.shape
        d_im = np.copy(base)
        tp = pred & gt; fn = gt & ~pred; fp = pred & ~gt
        d_im[..., 0] = np.clip(d_im[..., 0] + fn*0.9 + fp*0.9, 0, 1)
        d_im[..., 1] = np.clip(d_im[..., 1] + tp*0.9 + fp*0.7, 0, 1)
        return d_im

    fig, ax = plt.subplots(1, 3, figsize=(22, 8))
    ax[0].imshow(rgb_dim)
    ax[0].contour(d["egl_gt"], levels=[0.5], colors=["lime"], linewidths=0.4)
    ax[0].set_title(f"GT ({d['egl_gt'].sum()*(PX_UM/1000)**2:.2f} mm^2)")
    ax[0].axis("off")
    ax[1].imshow(diff_overlay(egl_curr, d["egl_gt"], rgb_dim.copy()))
    m_c = metrics(egl_curr, d["egl_gt"])
    ax[1].set_title(f"current params P={m_c['precision']:.2f} R={m_c['recall']:.2f} IoU={m_c['iou']:.2f}")
    ax[1].axis("off")
    ax[2].imshow(diff_overlay(egl_best, d["egl_gt"], rgb_dim.copy()))
    ax[2].set_title(f"best params P={best['precision']:.2f} R={best['recall']:.2f} IoU={best['iou']:.2f}\n"
                    f"pair={best['pair_pct']} band={best['band_um']} bp={best['band_pct']}")
    ax[2].axis("off")
    fig.tight_layout(); fig.savefig(OUT / "best_vs_current.png", dpi=140); plt.close(fig)
    print(f"saved best_vs_current.png in {OUT}")


if __name__ == "__main__":
    main()
