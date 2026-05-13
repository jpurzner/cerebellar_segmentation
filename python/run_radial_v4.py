"""v4: pia spokes UNION IGL spokes for full deep-fold coverage.

  pia spokes: pia -> IGL (good for gyral crowns)
  IGL spokes: IGL -> pia (good for deep folds where pia is hard to find)

For each set, run 1D layer detection along the spoke. Back-project EGL labels
from each, then union them. Compare against the port's EGL and the MATLAB GT.
"""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage as ndi
import tifffile
from scipy.spatial import cKDTree
from skimage import exposure
from radial_spokes import build_spokes
from igl_spokes import build_igl_spokes
from spoke_profiles import sample_along_spoke, detect_layers_1d

INPUT = Path(
    "/Users/jpurzner/Dropbox/images/edu_repeat/p27/"
    "2018_05_22_s1_2_p27-0021/2018_05_22_s1_2_p27-0021_fused_crop.tif"
)
PORT_MASKS = Path(__file__).resolve().parent / "figs" / "v_port" / "port_masks.npz"
SEG_GT_TIF = Path(
    "/Users/jpurzner/Dropbox/imaging_analysis/cerebellar_segmentation/_test_run/"
    "20x/2018_05_22_s1_2_p27-0021/"
    "2018_05_22_s1_2_p27-0021_fused_crop_segments.tif"
)
OUT = Path(__file__).resolve().parent / "figs" / "v_radial_v4"
OUT.mkdir(parents=True, exist_ok=True)
PX_UM = 0.5119049
SAMPLE_STEP_UM = 1.0


def to_unit(im):
    im = im.astype(np.float32); lo, hi = im.min(), im.max()
    return (im - lo) / max(hi - lo, 1e-9)


def clahe(im, clip=0.01):
    return exposure.equalize_adapthist(im, clip_limit=clip, nbins=256)


def extract_egl_gt(seg_rgb):
    s = seg_rgb.astype(int)
    def col(t): return np.all(np.abs(s - np.array(t)) < 16, axis=-1)
    return col([0,128,255]) | col([0,255,255])


def back_project_voronoi(spokes_pts, spoke_dirs, lengths,
                          boundaries_egl, boundaries_iegl_s, boundaries_iegl_e,
                          boundaries_other, cereb_mask, sense="pia",
                          PX_UM=PX_UM):
    """Generic voronoi back-projection.

    sense='pia':  spokes start at pia (s=0) and extend to IGL.
                  EGL is for s in [0, egl_end].
                  iEGL is the high-p27 region within EGL.
    sense='igl':  spokes start at IGL (s=0) and extend OUTWARD to pia.
                  EGL is detected at the OUTER end (s in [egl_start, L]).
                  Boundaries here mean: egl_start_um, iegl_s_um, iegl_e_um (relative
                  to the IGL anchor going outward).
    """
    H, W = cereb_mask.shape
    yy, xx = np.where(cereb_mask)
    pixel_pts = np.column_stack([yy, xx]).astype(np.float32)
    tree = cKDTree(spokes_pts)
    _, idx = tree.query(pixel_pts, k=1)

    pia = spokes_pts[idx]
    dr  = spoke_dirs[idx]
    v = pixel_pts - pia
    s_px = (v * dr).sum(axis=1) * PX_UM
    L = lengths[idx]
    ee = boundaries_egl[idx]
    ie_s = boundaries_iegl_s[idx]
    ie_e = boundaries_iegl_e[idx]

    egl_pred = np.zeros(len(pixel_pts), dtype=bool)
    iegl_pred = np.zeros(len(pixel_pts), dtype=bool)
    if sense == "pia":
        # EGL is s in [0, ee]
        is_egl = (s_px >= 0) & (s_px < ee)
        is_iegl = is_egl & (s_px >= ie_s) & (s_px <= ie_e)
    else:
        # IGL spokes: EGL is at the OUTER end -- s in [ee, L]
        # Here ee is "egl_start_um" measured from the IGL anchor outward.
        # iegl_s, iegl_e are also measured from IGL anchor outward.
        in_range = (s_px >= 0) & (s_px <= L + 5)
        is_egl  = in_range & (s_px >= ee) & (s_px <= L + 5)
        is_iegl = is_egl & (s_px >= ie_s) & (s_px <= ie_e)

    egl_pred[is_egl] = True
    iegl_pred[is_iegl] = True

    egl_2d = np.zeros((H, W), dtype=bool)
    iegl_2d = np.zeros((H, W), dtype=bool)
    egl_2d[yy, xx] = egl_pred
    iegl_2d[yy, xx] = iegl_pred
    return egl_2d, iegl_2d


def detect_pia_spoke_boundaries(chans, sl, PX_UM=PX_UM):
    """1D layer detection per pia spoke.
    Returns arrays: ee, ie_s, ie_e (all um from pia)."""
    n = len(sl.pia_points)
    ee = np.full(n, np.nan); ie_s = np.full(n, np.nan); ie_e = np.full(n, np.nan)
    for i in range(n):
        if not sl.valid[i]: continue
        prof = sample_along_spoke(chans, sl.pia_points[i], sl.igl_endpoints[i],
                                   SAMPLE_STEP_UM, PX_UM)
        layers = detect_layers_1d(prof, SAMPLE_STEP_UM,
                                   channels={"p27": 0, "neun": 1, "dapi": 2})
        ee[i] = layers["egl_end_um"]
        ie_s[i] = layers["iegl_start_um"]
        ie_e[i] = layers["iegl_end_um"]
    # smooth along arc length to remove sudden jumps
    sigma_pts = 5
    for arr in (ee, ie_s, ie_e):
        arr_valid = ~np.isnan(arr)
        if arr_valid.sum() > 5:
            arr[arr_valid] = ndi.gaussian_filter1d(
                arr[arr_valid], sigma=sigma_pts, mode='nearest')
    return ee, ie_s, ie_e


def detect_igl_spoke_boundaries(chans, isl, PX_UM=PX_UM):
    """1D layer detection per IGL spoke (going OUTWARD from IGL).

    Layer order along the spoke: IGL -> ML -> iEGL -> oEGL -> pia
    Detection: find DAPI peak in OUTER region (the far end), define EGL from
    egl_start (some um before peak) to L (the pia).
    """
    n = len(isl.igl_anchor)
    ee = np.full(n, np.nan); ie_s = np.full(n, np.nan); ie_e = np.full(n, np.nan)
    for i in range(n):
        if not isl.valid[i]: continue
        prof = sample_along_spoke(chans, isl.igl_anchor[i], isl.pia_endpoint[i],
                                   SAMPLE_STEP_UM, PX_UM)
        L_um = prof.shape[1] * SAMPLE_STEP_UM
        if L_um < 25: continue   # spoke too short
        # We want the EGL: high DAPI region near the OUTER end (last ~60um)
        sigma = max(0.5, 4.0 / SAMPLE_STEP_UM)
        dapi = ndi.gaussian_filter1d(prof[2], sigma=sigma)
        p27  = ndi.gaussian_filter1d(prof[0], sigma=sigma)

        # search the outer 80um for the DAPI peak (=EGL center)
        outer_um = 80.0
        outer_steps = min(prof.shape[1] - 1, int(outer_um / SAMPLE_STEP_UM))
        outer_start = max(0, prof.shape[1] - outer_steps)
        outer_seg = dapi[outer_start:]
        if len(outer_seg) < 5: continue

        peak_local = int(np.argmax(outer_seg))
        peak_step = outer_start + peak_local
        peak_val = outer_seg[peak_local]
        # walk INWARD from peak (toward IGL) -> EGL_start at first DROP below 50%
        thr = max(0.15, peak_val * 0.5)
        below = dapi < thr
        kernel_len = max(3, int(5.0 / SAMPLE_STEP_UM))
        sustained = ndi.binary_opening(below, structure=np.ones(kernel_len))
        cand = np.where(sustained[:peak_step])[0]
        if len(cand) > 0:
            egl_start_step = cand[-1]
        else:
            egl_start_step = max(0, peak_step - int(40 / SAMPLE_STEP_UM))

        # iEGL within EGL: high-p27 region
        p27_seg = p27[egl_start_step:]
        if len(p27_seg) > 0:
            pthr = max(0.20, p27_seg.max() * 0.55)
            in_iegl = p27_seg >= pthr
            idx = np.where(in_iegl)[0]
            if len(idx) > 0:
                iegl_s_step = egl_start_step + idx[0]
                iegl_e_step = egl_start_step + idx[-1]
            else:
                iegl_s_step = egl_start_step
                iegl_e_step = prof.shape[1] - 1
        else:
            iegl_s_step = egl_start_step
            iegl_e_step = prof.shape[1] - 1

        ee[i] = egl_start_step * SAMPLE_STEP_UM
        ie_s[i] = iegl_s_step * SAMPLE_STEP_UM
        ie_e[i] = iegl_e_step * SAMPLE_STEP_UM

    sigma_pts = 5
    for arr in (ee, ie_s, ie_e):
        arr_valid = ~np.isnan(arr)
        if arr_valid.sum() > 5:
            arr[arr_valid] = ndi.gaussian_filter1d(
                arr[arr_valid], sigma=sigma_pts, mode='nearest')
    return ee, ie_s, ie_e


def main():
    print(f"loading {INPUT.name}")
    stack = tifffile.imread(INPUT)
    p27_n  = clahe(to_unit(stack[0]))
    neun_n = clahe(to_unit(stack[1]))
    dapi_n = clahe(to_unit(stack[2]))
    chans = np.stack([p27_n, neun_n, dapi_n], axis=0)
    H, W = p27_n.shape
    rgb = np.stack([p27_n, neun_n, dapi_n], axis=-1); rgb /= max(rgb.max(), 1e-9)

    z = np.load(PORT_MASKS)
    cereb = z['mask'].astype(bool); igl = z['c_final'].astype(bool)
    port_egl = z['a_final'].astype(bool) | z['b_final'].astype(bool)

    print("building pia spokes")
    sl = build_spokes(cereb, igl, pixel_size_um=PX_UM, boundary_smooth_um=3.0)
    keep_p = sl.valid & (sl.lengths_um <= 250) & (sl.lengths_um >= 30)
    pia_pts = sl.pia_points[keep_p]; pia_dirs = sl.spoke_dirs[keep_p]
    pia_igl = sl.igl_endpoints[keep_p]; pia_L = sl.lengths_um[keep_p]
    print(f"  pia spokes: {keep_p.sum()}/{len(keep_p)}")

    # build a temp SpokeLayout-like dataclass for the boundary detection helper
    class TmpSL:
        pass
    tmp = TmpSL()
    tmp.pia_points = pia_pts
    tmp.igl_endpoints = pia_igl
    tmp.valid = np.ones(len(pia_pts), dtype=bool)
    print("detecting pia-spoke boundaries")
    ee_p, ies_p, iee_p = detect_pia_spoke_boundaries(chans, tmp)

    print("building IGL spokes")
    isl = build_igl_spokes(cereb, igl, pixel_size_um=PX_UM)
    keep_i = isl.valid & (isl.lengths_um <= 250) & (isl.lengths_um >= 30)
    igl_pts = isl.igl_anchor[keep_i]
    igl_dirs = isl.spoke_dirs[keep_i]
    igl_pia = isl.pia_endpoint[keep_i]
    igl_L = isl.lengths_um[keep_i]
    print(f"  IGL spokes: {keep_i.sum()}/{len(keep_i)}")

    class TmpISL:
        pass
    tmp_i = TmpISL()
    tmp_i.igl_anchor = igl_pts
    tmp_i.pia_endpoint = igl_pia
    tmp_i.valid = np.ones(len(igl_pts), dtype=bool)
    print("detecting IGL-spoke boundaries")
    ee_i, ies_i, iee_i = detect_igl_spoke_boundaries(chans, tmp_i)

    print("back-projecting via voronoi")
    egl_pia, iegl_pia = back_project_voronoi(
        pia_pts, pia_dirs, pia_L, ee_p, ies_p, iee_p, None, cereb, sense="pia")
    egl_igl, iegl_igl = back_project_voronoi(
        igl_pts, igl_dirs, igl_L, ee_i, ies_i, iee_i, None, cereb, sense="igl")

    egl_radial_union = (egl_pia | egl_igl) & cereb
    egl_full_union   = (egl_pia | egl_igl | port_egl) & cereb

    # === GT comparison ===
    seg = tifffile.imread(SEG_GT_TIF)
    egl_gt = extract_egl_gt(seg)

    def metrics(pred, gt):
        inter = (pred & gt).sum()
        union = (pred | gt).sum()
        return (inter / max(pred.sum(), 1),
                inter / max(gt.sum(), 1),
                inter / max(union, 1))

    p_pia, r_pia, iou_pia = metrics(egl_pia, egl_gt)
    p_igl, r_igl, iou_igl = metrics(egl_igl, egl_gt)
    p_ru,  r_ru,  iou_ru  = metrics(egl_radial_union, egl_gt)
    p_port,r_port,iou_port = metrics(port_egl, egl_gt)
    p_fu,  r_fu,  iou_fu  = metrics(egl_full_union, egl_gt)

    print(f"\n  source              area_mm2  P     R     IoU")
    print(f"  pia spokes only     {egl_pia.sum()*(PX_UM/1000)**2:.3f}     "
          f"{p_pia:.2f}  {r_pia:.2f}  {iou_pia:.3f}")
    print(f"  IGL spokes only     {egl_igl.sum()*(PX_UM/1000)**2:.3f}     "
          f"{p_igl:.2f}  {r_igl:.2f}  {iou_igl:.3f}")
    print(f"  radial union        {egl_radial_union.sum()*(PX_UM/1000)**2:.3f}     "
          f"{p_ru:.2f}  {r_ru:.2f}  {iou_ru:.3f}")
    print(f"  port EGL            {port_egl.sum()*(PX_UM/1000)**2:.3f}     "
          f"{p_port:.2f}  {r_port:.2f}  {iou_port:.3f}")
    print(f"  port + radial union {egl_full_union.sum()*(PX_UM/1000)**2:.3f}     "
          f"{p_fu:.2f}  {r_fu:.2f}  {iou_fu:.3f}")
    print(f"  GT total            {egl_gt.sum()*(PX_UM/1000)**2:.3f}")

    # === plots ===
    rgb_dim = rgb * 0.45
    fig, ax = plt.subplots(2, 3, figsize=(22, 14))
    ax = ax.ravel()
    panels = [
        ("RGB",                rgb,       None,           None),
        ("MATLAB GT EGL",      rgb_dim,   egl_gt,         "lime"),
        ("port EGL (a|b_final)", rgb_dim, port_egl,       "yellow"),
        ("pia spokes only",    rgb_dim,   egl_pia,        "magenta"),
        ("IGL spokes only",    rgb_dim,   egl_igl,        "orange"),
        ("port + radial UNION",rgb_dim,   egl_full_union, "red"),
    ]
    for k, (title, im, m, c) in enumerate(panels):
        ax[k].imshow(im)
        if m is not None:
            ax[k].contour(m, levels=[0.5], colors=[c], linewidths=0.5)
        ax[k].set_title(title); ax[k].axis("off")
    fig.tight_layout(); fig.savefig(OUT / "01_egl_sources.png", dpi=140)
    plt.close(fig)
    print(f"figures saved to {OUT}")


if __name__ == "__main__":
    main()
