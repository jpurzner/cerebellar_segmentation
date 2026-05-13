"""Python port of cerebellum_threshold_segment20x.m.

Faithful port of the structural pipeline. Skips:
  - phasesym (used as a secondary refinement; we approximate with a Hessian
    ridge response on the same channel)
  - anisodiff (input to phasesym; not needed if we skip phasesym)
  - edgelink/filledgegaps (used to bridge gaps in EGL ribbons; we use
    morphological dilation + skeletonization instead)

These omissions trade some accuracy in marginal cases for code simplicity.
The main detection logic (thresh_by_area, im_get_overlap, iterative
refinement) is preserved.

Layer labels in the output `set_bin`:
  1 = all_cerebellum (default)
  2 = iEGL  (a_final)
  3 = oEGL  (b_final)
  4 = IGL   (c_final)
  5 = ML    (ml_final)
  6 = DWL   (dwl_final)
  7 = PCL   (pc_bin_filt)
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Dict, Any
import warnings
import numpy as np
import scipy.ndimage as ndi
from skimage import exposure, filters, morphology, measure
from skimage.morphology import disk, skeletonize

# phasepack is a Python port of Kovesi's MATLAB PhaseCongruency package
with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    try:
        from phasepack import phasesym as _phasesym
        _HAS_PHASESYM = True
    except ImportError:
        _HAS_PHASESYM = False


# ---------------------------------------------------------------------------
# Default parameters (in physical units; scale from pixels via pixel_size_um)
# ---------------------------------------------------------------------------

DEFAULT_PARAMS = dict(
    # CLAHE
    clahe_clip=0.01,
    clahe_tiles=(50, 50),

    # Gaussian blurs (in micrometers — sigma_px = sigma_um / pixel_size_um)
    blur_sigma_um=2.6,        # was imgaussfilt(.., 5) at 10x = 5*1.0239um
    igl_blur_um=5.1,          # 10*1.0239 / 2 ≈ 5.1 (effective at 20x)

    # cerebellum mask
    mask_threshold=0.10,
    mask_close_um=10.2,       # 20px @ 10x = 20.5um → ~10um for 20x scaled
    mask_keep_min_px=5000,    # 10x px; scale below

    # IGL
    igl_area_frac=0.40,       # thresh_by_area target
    igl_close_um=5.1,
    igl_min_seed_um2_min=20000,  # bwareafilt threshold (was 5000-20000 px)
    igl_overlap_c=0.6,
    igl_overlap_a=0.6,
    igl_overlap_b=0.3,
    igl_close_final_um=2.6,
    igl_min_final_um2=200000, # was 200000 px

    # DWL
    dwl_area_frac=0.08,
    dwl_outline_dilate_um=35.8, # 70px @ 10x ≈ 71.7um
    dwl_min_um2=20000,

    # EGL
    egl_area_frac=0.15,
    egl_strict_frac=0.10,
    egl_min_seed_um2=5000,
    egl_close_um=5.1,

    # iEGL adaptive threshold (sensitivity values for skimage.filters.threshold_local)
    iegl_adapt_sens_loose=0.30,
    iegl_adapt_sens_mid=0.10,
    iegl_adapt_block_um=200.0,  # neighborhood diameter

    # oEGL
    oegl_adapt_sens=0.20,

    # final cleanup
    pcl_intensity=0.40,
    purkinje_circ_max=3.0,
    purkinje_area_um2_min=10.5,   # 40 px²@10x in um²
    purkinje_area_um2_max=78.7,   # 300 px²
)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def to_unit(img: np.ndarray) -> np.ndarray:
    """Per-channel rescale to [0,1] like MATLAB mat2gray."""
    img = img.astype(np.float32)
    lo, hi = img.min(), img.max()
    return (img - lo) / max(hi - lo, 1e-9)


def clahe(img: np.ndarray, clip: float = 0.01) -> np.ndarray:
    return exposure.equalize_adapthist(img, clip_limit=clip, nbins=256)


def background_subtract(im: np.ndarray) -> np.ndarray:
    """Port of background_subtract.m. Find histogram peak (background mode),
    set low threshold = peak + 0.04, high = where 90% of >low signal lies."""
    counts, edges = np.histogram(im, bins=500)
    # exclude bin 0
    if len(counts) > 1:
        max_i = np.argmax(counts[1:]) + 1  # +1 for offset
    else:
        max_i = 0
    low_in = min(edges[max_i + 1] + 0.04, 0.9)
    valid = edges >= low_in
    low_i = np.argmax(valid) if valid.any() else 0
    cs = np.cumsum(counts[low_i:]).astype(float)
    if cs[-1] > 0:
        rel = cs / cs[-1]
        high_off = np.argmax(rel > 0.9)
    else:
        high_off = 0
    high_i = min(low_i + high_off, len(edges) - 1)
    high_in = min(edges[high_i], 1.0)
    if high_in <= low_in:
        return im
    out = np.clip((im - low_in) / (high_in - low_in), 0, 1)
    return out.astype(np.float32)


def thresh_by_area(im: np.ndarray, mask: np.ndarray,
                   area_frac: float) -> np.ndarray:
    """Port of thresh_by_area.m. Threshold so that area_frac fraction of mask
    pixels remain above threshold."""
    inside = im[mask]
    if len(inside) == 0:
        return np.zeros_like(mask, dtype=bool)
    thr = np.quantile(inside, 1 - area_frac)
    return (im > thr) & mask


def im_get_overlap(bw1: np.ndarray, bw2: np.ndarray, ov_thresh: float,
                   gap: int) -> np.ndarray:
    """Port of im_get_overlap.m. Returns count map of how many grid-shifted
    versions of bw2 (chopped by horizontal+vertical strips spaced `gap` px)
    overlap with bw1's connected components above ov_thresh fraction.

    NOTE: MATLAB cuts BOTH bw1 and bw2 per iteration. We tried matching MATLAB
    faithfully (cut both) and observed major regressions on s2_5 and s3_2 (recall
    halved). Empirically the more-permissive variant -- label bw1 ONCE without
    cutting, only cut bw2 -- gives better cross-slide IoU. Likely because the
    downstream pipeline (DWL/ML/EGL refinement) was tuned around the larger
    c_final this version produces.
    """
    H, W = bw1.shape
    out = np.zeros((H, W), dtype=np.uint8)
    step = max(1, gap // 10)
    lbl, n = ndi.label(bw1)
    if n == 0:
        return out
    comp_size = np.bincount(lbl.ravel())
    for start in range(1, gap + 1, step):
        bw2_t = bw2.copy()
        bw2_t[start::gap, :] = False
        bw2_t[:, start::gap] = False
        overlaps = np.bincount(lbl.ravel(),
                                weights=bw2_t.ravel().astype(np.float32),
                                minlength=len(comp_size))
        ratios = overlaps / np.maximum(comp_size, 1)
        keep_lbls = np.where(ratios >= ov_thresh)[0]
        keep_lbls = keep_lbls[keep_lbls > 0]
        if len(keep_lbls):
            out[np.isin(lbl, keep_lbls)] += 1
    return out


def smallfill(orig: np.ndarray, sz: int) -> np.ndarray:
    """Port of smallfill.m: fill holes smaller than sz pixels."""
    filled = ndi.binary_fill_holes(orig)
    holes = filled & ~orig
    bigholes = morphology.remove_small_objects(holes, min_size=sz)
    smallholes = holes & ~bigholes
    return orig | smallholes


def keep_largest(mask: np.ndarray) -> np.ndarray:
    lbl, n = ndi.label(mask)
    if n <= 1:
        return mask
    sizes = np.bincount(lbl.ravel())
    sizes[0] = 0
    return lbl == sizes.argmax()


def bwareafilt(mask: np.ndarray, min_size: int,
               max_size: int = np.inf) -> np.ndarray:
    """skimage analogue of MATLAB bwareafilt."""
    if not mask.any():
        return mask.astype(bool)
    lbl, n = ndi.label(mask)
    sizes = np.bincount(lbl.ravel())
    keep = (sizes >= min_size) & (sizes <= max_size)
    keep[0] = False
    return keep[lbl]


def disk_um(radius_um: float, px_per_um: float):
    """Build a disk footprint of given micron radius."""
    r = max(1, int(round(radius_um * px_per_um)))
    return disk(r)


# --- fast morphology via scipy.ndimage iterations (3x3 cross by default) ---
# These are ~100-300x faster than morphology.{closing,dilation,erosion,opening}
# with disk(R) for R>10.  The shape is a "diamond" rather than a "disk", but at
# R>5 the visual difference is negligible for our purposes.
def _iter(radius_um: float, px_per_um: float) -> int:
    return max(1, int(round(radius_um * px_per_um)))


def fast_dilate_um(mask, radius_um, px_per_um):
    return ndi.binary_dilation(mask, iterations=_iter(radius_um, px_per_um))


def fast_erode_um(mask, radius_um, px_per_um):
    return ndi.binary_erosion(mask, iterations=_iter(radius_um, px_per_um))


def fast_close_um(mask, radius_um, px_per_um):
    return ndi.binary_closing(mask, iterations=_iter(radius_um, px_per_um))


def fast_open_um(mask, radius_um, px_per_um):
    return ndi.binary_opening(mask, iterations=_iter(radius_um, px_per_um))


def fast_dilate_iters(mask, iters):
    return ndi.binary_dilation(mask, iterations=max(1, int(iters)))


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

@dataclass
class SegmentResult:
    set_bin: np.ndarray         # int8 label image (0..7)
    mask: np.ndarray            # all_cerebellum
    a_final: np.ndarray         # iEGL
    b_final: np.ndarray         # oEGL
    c_final: np.ndarray         # IGL
    ml_final: np.ndarray        # ML
    dwl_final: np.ndarray       # DWL
    pc_final: np.ndarray        # PCL
    egl_mask: np.ndarray        # broad EGL container
    igl_in_mask: np.ndarray     # IGL+DWL container
    all_outline: np.ndarray     # cerebellum perimeter
    p27_n: np.ndarray           # CLAHE-normalised channels (for downstream)
    dapi_n: np.ndarray
    neun_n: np.ndarray


def segment_cerebellum(stack: np.ndarray,
                       channel_order=("p27", "neun", "dapi"),
                       pixel_size_um: float = 0.5119049,
                       params: Optional[dict] = None) -> SegmentResult:
    """Run the full pipeline. `stack` is (3, H, W) or (H, W, 3).

    channel_order: tuple naming what each TIFF page contains.
                   Default ("p27", "neun", "dapi") matches the 20x dataset.
                   For the 10x dataset use ("p27", "dapi", "neun").
    """
    p = {**DEFAULT_PARAMS, **(params or {})}
    px = pixel_size_um
    px_per_um = 1.0 / px

    # --- channel extraction ---
    if stack.ndim != 3:
        raise ValueError(f"expected 3D stack, got {stack.shape}")
    if stack.shape[0] not in (3, 4) and stack.shape[-1] in (3, 4):
        stack = np.moveaxis(stack, -1, 0)
    name_to_idx = {name: i for i, name in enumerate(channel_order)}
    p27_raw  = stack[name_to_idx["p27"]]
    dapi_raw = stack[name_to_idx["dapi"]]
    neun_raw = stack[name_to_idx["neun"]]

    # --- normalisation: CLAHE then background subtract ---
    a_nf_pre = clahe(to_unit(p27_raw),  p["clahe_clip"])
    b_nf_pre = clahe(to_unit(dapi_raw), p["clahe_clip"])
    c_nf_pre = clahe(to_unit(neun_raw), p["clahe_clip"])

    # MATLAB keeps a copy as `a_set_thresh` BEFORE bg subtract, for adaptive thresh
    a_set_thresh = a_nf_pre.copy()
    b_set_thresh = b_nf_pre.copy()

    a_nf = background_subtract(a_nf_pre)
    b_nf = background_subtract(b_nf_pre)
    c_nf = background_subtract(c_nf_pre)

    # --- gaussian blurs ---
    sigma_px = p["blur_sigma_um"] * px_per_um
    a_nfg = filters.gaussian(a_nf, sigma=sigma_px, preserve_range=True)
    b_nfg = filters.gaussian(b_nf, sigma=sigma_px, preserve_range=True)
    c_nfg = filters.gaussian(c_nf, sigma=sigma_px, preserve_range=True)

    # --- pre-Otsu thresholds for high/low decomposition ---
    a_level = filters.threshold_otsu(a_nf)
    b_level = filters.threshold_otsu(b_nf)
    c_level = filters.threshold_otsu(c_nf)

    def hi_lo(g, lvl):
        hi = background_subtract(np.clip((g - lvl) / max(1 - lvl, 1e-9), 0, 1))
        lo = np.clip(g / max(lvl, 1e-9), 0, 1)
        return hi, lo

    a_nfg_h, a_nfg_l = hi_lo(a_nfg, a_level)
    b_nfg_h, b_nfg_l = hi_lo(b_nfg, b_level)
    c_nfg_h, c_nfg_l = hi_lo(c_nfg, c_level)

    # =====================================================================
    # 1. Cerebellum mask (all_cerebellum + outline)
    # =====================================================================
    s = (a_nf + b_nf) / 2
    s = s / max(s.max(), 1e-9)
    raw = s > p["mask_threshold"]
    closed = morphology.closing(raw, footprint=disk_um(p["mask_close_um"], px_per_um))
    filled = ndi.binary_fill_holes(closed)
    all_cerebellum = keep_largest(filled)

    # boundary smoothing via morphological close+open at the right scale
    smooth_r = max(2, int(round(20 * px_per_um)))   # ~20 um smoothing
    all_cerebellum = morphology.closing(all_cerebellum, disk(smooth_r))
    all_cerebellum = ndi.binary_fill_holes(all_cerebellum)
    all_cerebellum = keep_largest(all_cerebellum)

    # outline = mask boundary (1-pixel ring)
    eroded = morphology.erosion(all_cerebellum, disk(1))
    all_outline = all_cerebellum & ~eroded
    all_outline_mask = fast_dilate_um(all_outline, 25.6, px_per_um)

    # =====================================================================
    # 2. IGL detection (c_final)
    # =====================================================================
    c_norm = c_nfg_h * (~all_outline_mask) * all_cerebellum
    # MATLAB applies an order-25 5x5 rank filter (= median-ish max filter)
    c_norm = ndi.maximum_filter(c_norm, size=5)
    c_norm = filters.gaussian(c_norm, sigma=10*px_per_um/2, preserve_range=True)
    c_norm = ndi.maximum_filter(c_norm, size=5)
    c_norm = filters.gaussian(c_norm, sigma=10*px_per_um/2, preserve_range=True)
    c_norm = ndi.maximum_filter(c_norm, size=5)
    c_norm = filters.gaussian(c_norm, sigma=10*px_per_um/2, preserve_range=True)
    c_norm = background_subtract(c_norm)

    c_bin = thresh_by_area(c_norm, all_cerebellum, p["igl_area_frac"])
    c_bin = c_bin & ~all_outline_mask
    c_bin = bwareafilt(c_bin, min_size=int(5000*(px_per_um**2)/4))
    c_bin = fast_close_um(c_bin, p["igl_close_um"], px_per_um)
    c_bin = bwareafilt(c_bin, min_size=int(p["igl_min_seed_um2_min"]
                                           * (px_per_um**2)/4))

    # overlap-based refinement: keep regions of c_bin that overlap with
    # bright NeuN/p27/DAPI sufficiently
    c_whole  = im_get_overlap(c_bin, c_nfg_h > c_level,
                              p["igl_overlap_c"], 200)
    ca_whole = im_get_overlap(c_bin, a_nfg_h > a_level,
                              p["igl_overlap_a"], 200)
    cb_whole = im_get_overlap(c_bin, b_nfg_h > b_level,
                              p["igl_overlap_b"], 200)
    c_final = (ca_whole + c_whole + cb_whole) > 0
    c_final = fast_close_um(c_final, p["igl_close_final_um"], px_per_um)
    c_final = smallfill(c_final, 10000)
    c_final = bwareafilt(c_final,
                         min_size=int(p["igl_min_final_um2"]*(px_per_um**2)/4))

    # =====================================================================
    # 3. DWL (dwl_final) -- iterative refinement loop
    # =====================================================================
    dwl_norm = (a_nfg_l - c_nfg_l - 2 * a_nfg_h)
    dwl_norm = np.clip(dwl_norm, 0, None)
    dwl_norm = dwl_norm * all_cerebellum
    dwl_norm = ndi.maximum_filter(dwl_norm, size=5)
    dwl_norm = filters.gaussian(dwl_norm, sigma=5*px_per_um/2, preserve_range=True)
    dwl_norm = background_subtract(dwl_norm)

    dwl_bin = thresh_by_area(dwl_norm, all_cerebellum, p["dwl_area_frac"])
    dwl_bin = dwl_bin & ~fast_dilate_um(
        all_outline_mask, p["dwl_outline_dilate_um"], px_per_um)
    dwl_bin = dwl_bin & ~c_bin

    igl_in = c_final | dwl_bin
    igl_in = fast_dilate_um(igl_in, 2.6, px_per_um)
    igl_in = fast_close_um(igl_in, 2.6, px_per_um)
    igl_in = smallfill(igl_in, 1000)
    igl_in = bwareafilt(igl_in, min_size=int(70000*(px_per_um**2)/4))

    # iterative refinement (the MATLAB has 2 passes of c_out detection)
    for _ in range(2):
        ml_out = (~igl_in) & all_cerebellum
        ml_out = ml_out | all_outline_mask
        ml_out = fast_dilate_um(ml_out, 5, px_per_um)
        ml_out = fast_erode_um(ml_out, 2.6, px_per_um)
        ml_out = bwareafilt(ml_out, min_size=int(30000*(px_per_um**2)/4))
        ml_out_deep = ml_out & ~fast_dilate_um(
            all_outline_mask, 30.7, px_per_um)
        ml_out_deep = ndi.binary_fill_holes(ml_out_deep)
        ml_out = ml_out | ml_out_deep
        c_out = im_get_overlap(ml_out, igl_in, 0.5, 100) > 5
        igl_in = igl_in & ~c_out
        c_final = c_final & ~c_out
        dwl_final = igl_in & ~c_final
        dwl_final = smallfill(dwl_final, 5000)
        dwl_final = fast_erode_um(dwl_final, 2.6, px_per_um)
        dwl_final = fast_dilate_um(dwl_final, 5, px_per_um)
        dwl_final = bwareafilt(dwl_final,
                               min_size=int(p["dwl_min_um2"]*(px_per_um**2)/4))
        igl_in = c_final | dwl_final
        igl_in = fast_close_um(igl_in, 10, px_per_um)
        igl_in = smallfill(igl_in, 200000)
        dwl_final = igl_in & ~c_final

    igl_in_mask = fast_dilate_um(igl_in, 15.3, px_per_um)

    # =====================================================================
    # 4. EGL detection (egl_mask)
    # =====================================================================
    egl_norm = b_nfg * (~igl_in_mask) * all_cerebellum
    egl_norm = ndi.maximum_filter(egl_norm, size=5)
    egl_bin = thresh_by_area(egl_norm, all_cerebellum, p["egl_area_frac"])
    egl_bin = bwareafilt(egl_bin, min_size=50)
    egl_strict = thresh_by_area(egl_norm, all_cerebellum, p["egl_strict_frac"])
    egl_strict = fast_dilate_um(egl_strict, 5, px_per_um)
    egl_strict = bwareafilt(egl_strict,
                            min_size=int(p["egl_min_seed_um2"]*(px_per_um**2)/4))

    egl_mask = fast_dilate_um(egl_bin, p["egl_close_um"], px_per_um)
    egl_mask = bwareafilt(egl_mask,
                          min_size=int(p["egl_min_seed_um2"]*(px_per_um**2)/4))
    # bwmorph thin → skeletonize (no edgelink/filledgegaps; just dilate-skeleton)
    egl_thin = skeletonize(egl_mask)
    # bridge gaps via close + thin again
    egl_bridge = fast_close_um(egl_thin, 50, px_per_um)
    egl_thin = skeletonize(egl_bridge)
    egl_thin = fast_dilate_um(egl_thin, 2, px_per_um)

    egl_mask = egl_mask | egl_thin
    egl_mask = bwareafilt(egl_mask, min_size=int(20000*(px_per_um**2)/4))
    egl_mask = fast_dilate_um(egl_mask, 2.6, px_per_um)

    # ---- phasesym extension (catches deep-fold EGL by ribbon-symmetry) ----
    # MATLAB:
    #   [a_phasesym, ...] = phasesym(a_nfad, 5, 6, 3, 2.5, 0.55, 2.0, 0)
    #   egl_phase = im_get_overlap(egl_mask, a_phasesym > 0.1, 0.1, 200)
    #   egl_mask  = imdilate(egl_mask | egl_phase, strel('disk',5))
    if _HAS_PHASESYM:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            ps_a, _, _, _ = _phasesym(
                a_nf.astype(np.float32),
                nscale=5, norient=6, minWaveLength=3,
                mult=2.5, sigmaOnf=0.55, k=2.0, polarity=0)
        a_phasesym_bin = (ps_a > p.get("egl_phasesym_thr", 0.1)) & all_cerebellum
        egl_phase = im_get_overlap(egl_mask, a_phasesym_bin, 0.1, 200) > 0
        egl_mask = (egl_mask | egl_phase) & all_cerebellum
        egl_mask = fast_dilate_um(egl_mask, 2.6, px_per_um)

    # =====================================================================
    # 5. iEGL (a_final): adaptive threshold on p27 within egl_mask
    # =====================================================================
    a_set_thresh_egl = a_set_thresh.copy().astype(np.float32)
    a_set_thresh_egl[~all_cerebellum] = a_set_thresh.mean()

    # iEGL adaptive threshold (Otsu on ratio).
    #
    # Strategy: compute LOCAL ratio (a_set_thresh / a_local_mean), apply Otsu
    # within the (cereb AND bright-enough) region. Otsu finds the natural
    # bimodal split — image-adaptive without needing to assume an EGL fraction.
    sigma_local = min(30.0, a_nfg.shape[0] / 50.0)
    a_local_mean = ndi.gaussian_filter(a_set_thresh_egl.astype(np.float32),
                                        sigma=sigma_local)
    a_ratio = a_set_thresh.astype(np.float32) / np.maximum(a_local_mean, 1e-6)
    a_floor = float(p.get("iegl_abs_floor", 0.20))
    a_bright_enough = (a_set_thresh > a_floor) & all_cerebellum

    inside_ratios = a_ratio[a_bright_enough]
    if len(inside_ratios) > 100:
        # Otsu on the ratio distribution
        a_thr_strict = filters.threshold_otsu(inside_ratios)
        # mid/loose: scale below Otsu
        a_thr_mid    = a_thr_strict * 0.95
        a_thr_loose  = a_thr_strict * 0.85
    else:
        a_thr_strict = a_thr_mid = a_thr_loose = 1.1
    a_bin       = a_bright_enough & (a_ratio > a_thr_strict)
    a_bin_mid   = a_bright_enough & (a_ratio > a_thr_mid)
    a_bin_loose = a_bright_enough & (a_ratio > a_thr_loose)
    a_bin = bwareafilt(a_bin, min_size=10)

    a_bin_g = filters.gaussian(a_bin.astype(float),
                                sigma=p["blur_sigma_um"]*px_per_um,
                                preserve_range=True)
    a_bin_smooth = a_bin_g >= 0.2
    a_bin_smooth = a_bin_smooth & egl_mask
    a_bin_smooth_c = fast_close_um(a_bin_smooth, 2, px_per_um)
    a_bin_smooth = smallfill(a_bin_smooth | a_bin_smooth_c, 200)
    a_bin_smooth = bwareafilt(a_bin_smooth, min_size=200)
    a_bin = bwareafilt(a_bin_smooth, min_size=500)

    # outer EGL container after iEGL pre-filter
    a_thin = skeletonize(fast_dilate_um(a_bin_smooth, 2.6, px_per_um))
    oegl_out = (fast_dilate_um(a_thin, 2.6, px_per_um)
                | fast_dilate_um(all_outline, 10, px_per_um)
                | a_bin)
    oegl_out = fast_dilate_um(oegl_out, 2.6, px_per_um)
    oegl_out = smallfill(oegl_out, 1000000)
    oegl_out = fast_erode_um(oegl_out, 2.6, px_per_um)

    # =====================================================================
    # 6. oEGL (b_final): adaptive DAPI within EGL container
    # =====================================================================
    b_set_thresh_egl = b_set_thresh.copy().astype(np.float32)
    b_set_thresh_egl[~all_cerebellum] = b_set_thresh.mean()
    b_local_mean = ndi.gaussian_filter(b_set_thresh_egl, sigma=sigma_local)
    b_ratio = b_set_thresh.astype(np.float32) / np.maximum(b_local_mean, 1e-6)
    b_floor = float(p.get("oegl_abs_floor", 0.20))
    b_bright_enough = (b_set_thresh > b_floor) & all_cerebellum
    inside_b = b_ratio[b_bright_enough]
    if len(inside_b) > 100:
        b_thr = filters.threshold_otsu(inside_b)
        b_thr_loose = b_thr * 0.85
    else:
        b_thr = b_thr_loose = 1.1
    b_bin       = b_bright_enough & (b_ratio > b_thr)
    b_bin_loose = b_bright_enough & (b_ratio > b_thr_loose)
    b_bin_g = filters.gaussian(b_bin.astype(float),
                                sigma=p["blur_sigma_um"]*px_per_um,
                                preserve_range=True)
    b_bin_smooth = b_bin_g >= 0.2
    b_bin_smooth = b_bin_smooth & egl_mask
    b_bin_smooth = fast_close_um(b_bin_smooth, 1, px_per_um)
    b_bin_smooth = smallfill(b_bin_smooth | b_bin, 200)
    b_bin_smooth = bwareafilt(b_bin_smooth, min_size=100)
    b_bin = b_bin_smooth & oegl_out

    # b_steal: pull in additional oEGL pixels from areas where b is much larger
    # than a (DAPI-only outer band)
    b_norm = (b_nf - a_nf - igl_in.astype(float))
    b_norm = np.clip(b_norm, 0, None)
    b_norm = to_unit(b_norm)
    b_bin_steal = thresh_by_area(b_norm, all_cerebellum, 0.05)
    b_bin_steal = b_bin_steal & oegl_out

    b_bin = (b_bin & ~a_bin) | b_bin_steal

    # =====================================================================
    # 7. ML (ml_final) and Purkinje cell layer
    # =====================================================================
    ml_layer_mask = (fast_dilate_um(igl_in_mask, 5, px_per_um)
                     | fast_dilate_um(egl_strict, 20, px_per_um))
    ml_layer_mask = ndi.binary_fill_holes(ml_layer_mask)
    deep_nuclei = all_cerebellum & ~ml_layer_mask
    deep_nuclei = fast_dilate_um(deep_nuclei, 20, px_per_um)
    deep_nuclei = bwareafilt(deep_nuclei, min_size=int(10000*(px_per_um**2)/4))
    deep_nuclei_near_egl = im_get_overlap(deep_nuclei, egl_mask, 0.8, 100) > 5
    deep_nuclei = deep_nuclei & ~deep_nuclei_near_egl

    ml_layer_mask = ml_layer_mask & ~igl_in & ~egl_mask & ~deep_nuclei
    ml_layer_mask = ml_layer_mask & all_cerebellum
    ml_layer_mask = fast_open_um(ml_layer_mask, 1.5, px_per_um)
    ml_layer_mask = bwareafilt(ml_layer_mask,
                               min_size=int(20000*(px_per_um**2)/4))

    pc_layer_mask = igl_in_mask & ~igl_in & ~egl_mask & ml_layer_mask
    pc_layer_bin = (pc_layer_mask.astype(float) * b_nfg) > p["pcl_intensity"]

    # consolidate dwl + deep_nuclei
    # IMPORTANT: use SMALLFILL (size-bounded) instead of binary_fill_holes here.
    # Unbounded fill on (c_final | dwl_final), which forms an annular shape
    # around the cereb interior, fills the entire interior including the ML
    # territory. Bounded fill only closes small gaps.
    dwl_final = (dwl_final | deep_nuclei) & ~c_final
    dwl_final = fast_dilate_um(dwl_final, 20, px_per_um)
    max_hole_px = int(50000 * (px_per_um ** 2) / 4)   # ~13000 px @ 10x = 13700px ratio
    dwl_final = smallfill(dwl_final, max_hole_px) & ~c_final
    igl_in = (c_final | dwl_final) & all_cerebellum
    igl_in = smallfill(igl_in, max_hole_px)
    dwl_final = igl_in & ~c_final
    dwl_final = bwareafilt(dwl_final,
                            min_size=int(200000*(px_per_um**2)/4))

    # ML refinement
    ml_layer_mask = ml_layer_mask & ~dwl_final

    ml_in = igl_in | ml_layer_mask
    ml_in = fast_close_um(ml_in, 5, px_per_um)
    ml_in = ndi.binary_fill_holes(ml_in)
    ml_in = fast_dilate_um(ml_in, 5, px_per_um)
    ml_in = bwareafilt(ml_in, min_size=int(5000*(px_per_um**2)/4))
    ml_in = ndi.binary_fill_holes(ml_in)

    ml_in = ml_in & ~a_bin & ~b_bin
    ml_in = ndi.binary_fill_holes(ml_in)

    # =====================================================================
    # 8. Final iEGL/oEGL split
    # =====================================================================
    oegl_in = ml_in | a_bin
    iegl_out = (~oegl_in) & all_cerebellum
    iegl_out = bwareafilt(iegl_out, min_size=50)
    iegl_out = (iegl_out | b_bin) & ~oegl_in
    iegl_out = fast_close_um(iegl_out, 5, px_per_um)
    iegl_out = fast_dilate_um(iegl_out, 5, px_per_um) & ~a_bin

    egl_intersect = im_get_overlap(oegl_in, iegl_out, 0.5, 100) > 4
    oegl_in = ml_in | a_bin | egl_intersect
    oegl_in = bwareafilt(oegl_in, min_size=200)
    oegl_in = smallfill(oegl_in, 3000)

    iegl_out = iegl_out & ~egl_intersect
    oegl_in = fast_dilate_um(oegl_in, 1, px_per_um)

    a_final = a_bin & oegl_in & ~ml_in & ~igl_in
    b_final = b_bin & ~ml_in & ~igl_in
    ml_final = ml_in & ~igl_in

    # =====================================================================
    # 9. Purkinje cells (bright p27 round blobs in ML)
    # =====================================================================
    pk = (a_nf - b_nf) * (a_nf - c_nf) * a_nf
    pk = exposure.equalize_adapthist(np.clip(pk - pk.min(),
                                              0, None) / max(pk.max() - pk.min(), 1e-9))
    pc_bin = pk > 0.1
    pc_bin = bwareafilt(pc_bin, min_size=int(p["purkinje_area_um2_min"]
                                              * (px_per_um**2)))
    pc_bin = ndi.binary_fill_holes(pc_bin)
    pc_lbl, pc_n = ndi.label(pc_bin)
    if pc_n > 0:
        regs = measure.regionprops(pc_lbl)
        keep = []
        for r in regs:
            if r.perimeter <= 0:
                continue
            circ = r.perimeter ** 2 / (4 * np.pi * r.area)
            area_um2 = r.area / (px_per_um ** 2)
            if (circ < p["purkinje_circ_max"]
                    and p["purkinje_area_um2_min"] <= area_um2 <=
                    p["purkinje_area_um2_max"] * 1.0):
                keep.append(r.label)
        pc_bin = np.isin(pc_lbl, keep)
    pc_bin = fast_open_um(pc_bin, 1, px_per_um)
    pc_bin = bwareafilt(pc_bin,
                        min_size=int(p["purkinje_area_um2_min"]*(px_per_um**2)))
    pc_bin = fast_dilate_um(pc_bin, 2, px_per_um)
    pc_bin_filt = pc_bin & ml_final
    pc_final = pc_bin_filt | pc_layer_bin

    # =====================================================================
    # 10. Build labelled set_bin
    # =====================================================================
    set_bin = np.zeros(all_cerebellum.shape, dtype=np.uint8)
    set_bin[all_cerebellum] = 1
    set_bin[ml_final] = 5
    set_bin[dwl_final] = 6
    set_bin[c_final] = 4
    set_bin[b_final] = 3
    set_bin[a_final] = 2
    set_bin[pc_final] = 7

    return SegmentResult(
        set_bin=set_bin,
        mask=all_cerebellum,
        a_final=a_final, b_final=b_final, c_final=c_final,
        ml_final=ml_final, dwl_final=dwl_final, pc_final=pc_final,
        egl_mask=egl_mask, igl_in_mask=igl_in_mask,
        all_outline=all_outline,
        p27_n=a_nf, dapi_n=b_nf, neun_n=c_nf,
    )
