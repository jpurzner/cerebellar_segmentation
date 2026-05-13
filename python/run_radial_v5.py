"""v5: ribbon-stamp back-projection (instead of voronoi).

For each valid spoke, paint the EGL ribbon as a segment of length egl_end_um
along the spoke direction, dilated by half spoke spacing. Union pia spokes
+ IGL spokes. Compare against MATLAB GT.

Should give a thinner, more accurate ribbon than voronoi-cell back-projection.
"""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage as ndi
import tifffile
from skimage import exposure
from skimage.draw import line as draw_line
from radial_spokes import build_spokes, PIAL_SAMPLE_UM
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
OUT = Path(__file__).resolve().parent / "figs" / "v_radial_v5"
OUT.mkdir(parents=True, exist_ok=True)
PX_UM = 0.5119049
SAMPLE_STEP_UM = 1.0
RIBBON_THICKEN_UM = 6.0   # dilate the spoke-segment ribbon by this much


def to_unit(im):
    im = im.astype(np.float32); lo, hi = im.min(), im.max()
    return (im - lo) / max(hi - lo, 1e-9)


def clahe(im, clip=0.01):
    return exposure.equalize_adapthist(im, clip_limit=clip, nbins=256)


def extract_egl_gt(seg_rgb):
    s = seg_rgb.astype(int)
    def col(t): return np.all(np.abs(s - np.array(t)) < 16, axis=-1)
    iEGL = col([0,128,255]); oEGL = col([0,255,255])
    return iEGL | oEGL, iEGL, oEGL


def stamp_segment(mask_2d, start_yx, dir_yx, length_um, px_per_um):
    """Stamp a line segment from start in direction dir, of length length_um, into mask_2d."""
    H, W = mask_2d.shape
    end_y = start_yx[0] + dir_yx[0] * length_um * px_per_um
    end_x = start_yx[1] + dir_yx[1] * length_um * px_per_um
    y0 = int(round(start_yx[0])); x0 = int(round(start_yx[1]))
    y1 = int(round(end_y));      x1 = int(round(end_x))
    y0 = np.clip(y0, 0, H-1); x0 = np.clip(x0, 0, W-1)
    y1 = np.clip(y1, 0, H-1); x1 = np.clip(x1, 0, W-1)
    rr, cc = draw_line(y0, x0, y1, x1)
    mask_2d[rr, cc] = True


def detect_pia_spoke_boundaries(chans, pia_pts, igl_pts, valid):
    n = len(pia_pts)
    ee = np.full(n, np.nan); ie_s = np.full(n, np.nan); ie_e = np.full(n, np.nan)
    for i in range(n):
        if not valid[i]: continue
        prof = sample_along_spoke(chans, pia_pts[i], igl_pts[i],
                                   SAMPLE_STEP_UM, PX_UM)
        layers = detect_layers_1d(prof, SAMPLE_STEP_UM,
                                   channels={"p27": 0, "neun": 1, "dapi": 2})
        ee[i] = layers["egl_end_um"]
        ie_s[i] = layers["iegl_start_um"]
        ie_e[i] = layers["iegl_end_um"]
    sigma_pts = 5
    for arr in (ee, ie_s, ie_e):
        arr_valid = ~np.isnan(arr)
        if arr_valid.sum() > 5:
            arr[arr_valid] = ndi.gaussian_filter1d(
                arr[arr_valid], sigma=sigma_pts, mode='nearest')
    return ee, ie_s, ie_e


def detect_igl_spoke_boundaries(chans, igl_pts, pia_pts, valid):
    """Detect EGL_start (um from IGL anchor) per IGL spoke."""
    n = len(igl_pts)
    egl_start_um = np.full(n, np.nan)   # how far from IGL the EGL begins
    iegl_s_um = np.full(n, np.nan)
    iegl_e_um = np.full(n, np.nan)
    for i in range(n):
        if not valid[i]: continue
        prof = sample_along_spoke(chans, igl_pts[i], pia_pts[i],
                                   SAMPLE_STEP_UM, PX_UM)
        L_steps = prof.shape[1]
        if L_steps < 25: continue
        sigma = max(0.5, 4.0 / SAMPLE_STEP_UM)
        dapi = ndi.gaussian_filter1d(prof[2], sigma=sigma)
        p27  = ndi.gaussian_filter1d(prof[0], sigma=sigma)
        # EGL = high DAPI region near outer end (last 80um)
        outer_um = 80
        outer_steps = min(L_steps - 1, int(outer_um / SAMPLE_STEP_UM))
        outer_start = max(0, L_steps - outer_steps)
        outer_seg = dapi[outer_start:]
        if len(outer_seg) < 5: continue
        peak_local = int(np.argmax(outer_seg))
        peak_step = outer_start + peak_local
        peak_val = outer_seg[peak_local]
        thr = max(0.15, peak_val * 0.5)
        below = dapi < thr
        kernel_len = max(3, int(5.0 / SAMPLE_STEP_UM))
        sustained = ndi.binary_opening(below, structure=np.ones(kernel_len))
        cand = np.where(sustained[:peak_step])[0]
        if len(cand) > 0:
            egl_start_step = cand[-1]
        else:
            egl_start_step = max(0, peak_step - int(40 / SAMPLE_STEP_UM))
        egl_start_um[i] = egl_start_step * SAMPLE_STEP_UM
        # iEGL within EGL
        p27_seg = p27[egl_start_step:]
        if len(p27_seg) > 0:
            pthr = max(0.20, p27_seg.max() * 0.55)
            in_iegl = p27_seg >= pthr
            idx = np.where(in_iegl)[0]
            if len(idx) > 0:
                iegl_s_um[i] = (egl_start_step + idx[0]) * SAMPLE_STEP_UM
                iegl_e_um[i] = (egl_start_step + idx[-1]) * SAMPLE_STEP_UM
    sigma_pts = 5
    for arr in (egl_start_um, iegl_s_um, iegl_e_um):
        arr_valid = ~np.isnan(arr)
        if arr_valid.sum() > 5:
            arr[arr_valid] = ndi.gaussian_filter1d(
                arr[arr_valid], sigma=sigma_pts, mode='nearest')
    return egl_start_um, iegl_s_um, iegl_e_um


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

    px_per_um = 1.0 / PX_UM

    print("building pia spokes")
    sl = build_spokes(cereb, igl, pixel_size_um=PX_UM, boundary_smooth_um=3.0)
    keep_p = sl.valid & (sl.lengths_um <= 250) & (sl.lengths_um >= 30)
    pia_pts = sl.pia_points[keep_p]; pia_dirs = sl.spoke_dirs[keep_p]
    pia_igl = sl.igl_endpoints[keep_p]; pia_L = sl.lengths_um[keep_p]
    print(f"  {keep_p.sum()} pia spokes")

    print("detecting pia-spoke boundaries")
    ee_p, ies_p, iee_p = detect_pia_spoke_boundaries(chans, pia_pts, pia_igl,
                                                      np.ones(len(pia_pts), bool))

    print("building IGL spokes")
    isl = build_igl_spokes(cereb, igl, pixel_size_um=PX_UM)
    keep_i = isl.valid & (isl.lengths_um <= 250) & (isl.lengths_um >= 30)
    igl_pts = isl.igl_anchor[keep_i]; igl_dirs = isl.spoke_dirs[keep_i]
    igl_pia = isl.pia_endpoint[keep_i]; igl_L = isl.lengths_um[keep_i]
    print(f"  {keep_i.sum()} IGL spokes")

    print("detecting IGL-spoke boundaries")
    egl_s_i, ies_i, iee_i = detect_igl_spoke_boundaries(chans, igl_pts, igl_pia,
                                                         np.ones(len(igl_pts), bool))

    print("ribbon-stamp back-projection")
    egl_pia_thin = np.zeros((H, W), dtype=bool)
    for i in range(len(pia_pts)):
        if np.isnan(ee_p[i]): continue
        # paint segment from pia along spoke for length egl_end_um
        stamp_segment(egl_pia_thin, pia_pts[i], pia_dirs[i], ee_p[i], px_per_um)

    egl_igl_thin = np.zeros((H, W), dtype=bool)
    for i in range(len(igl_pts)):
        if np.isnan(egl_s_i[i]): continue
        # IGL spoke EGL: from egl_s_i along the IGL->pia direction for (L - egl_s_i) um
        L_um = igl_L[i]
        # the EGL covers s in [egl_s_i, L]; start point = igl_anchor + dir * egl_s_i
        start_y = igl_pts[i, 0] + igl_dirs[i, 0] * egl_s_i[i] * px_per_um
        start_x = igl_pts[i, 1] + igl_dirs[i, 1] * egl_s_i[i] * px_per_um
        seg_len_um = max(1.0, L_um - egl_s_i[i])
        stamp_segment(egl_igl_thin, (start_y, start_x), igl_dirs[i], seg_len_um,
                       px_per_um)

    # dilate ribbons by half spoke spacing
    dil_iters = max(1, int(round(RIBBON_THICKEN_UM * px_per_um)))
    egl_pia_thick = ndi.binary_dilation(egl_pia_thin, iterations=dil_iters) & cereb
    egl_igl_thick = ndi.binary_dilation(egl_igl_thin, iterations=dil_iters) & cereb
    egl_radial   = egl_pia_thick | egl_igl_thick

    # Biological prior: EGL is within ~70um of the cerebellum boundary.
    # Filter the RADIAL contributions by this constraint to drop spurious FPs in
    # the deep interior (the port itself doesn't need this filter as it already
    # uses DAPI thresholding).
    dist_to_bg_um = ndi.distance_transform_edt(cereb).astype(np.float32) * PX_UM
    near_pia = dist_to_bg_um <= 70.0
    egl_radial_filt = egl_radial & near_pia
    egl_full     = egl_radial_filt | port_egl

    # === GT comparison ===
    seg = tifffile.imread(SEG_GT_TIF)
    egl_gt, _, _ = extract_egl_gt(seg)

    def m(p, g):
        i = (p & g).sum(); u = (p | g).sum()
        return (i / max(p.sum(), 1), i / max(g.sum(), 1), i / max(u, 1))

    P, R, IoU = {}, {}, {}
    for name, pred in [("pia ribbon", egl_pia_thick), ("IGL ribbon", egl_igl_thick),
                        ("radial union", egl_radial),
                        ("radial UNION (filt 70um)", egl_radial_filt),
                        ("port", port_egl),
                        ("port + radial(filt)", egl_full)]:
        p, r, iou = m(pred, egl_gt)
        P[name] = p; R[name] = r; IoU[name] = iou
        print(f"  {name:18s}  area={pred.sum()*(PX_UM/1000)**2:5.3f} mm^2  "
              f"P={p:.3f}  R={r:.3f}  IoU={iou:.3f}")
    print(f"  GT total            area={egl_gt.sum()*(PX_UM/1000)**2:5.3f} mm^2")

    # === plots ===
    rgb_dim = rgb * 0.45
    fig, ax = plt.subplots(2, 3, figsize=(22, 14)); ax = ax.ravel()
    panels = [
        ("RGB",                rgb,       None,          None),
        ("MATLAB GT",          rgb_dim,   egl_gt,        "lime"),
        ("port EGL",           rgb_dim,   port_egl,      "yellow"),
        ("pia ribbon",         rgb_dim,   egl_pia_thick, "magenta"),
        ("IGL ribbon",         rgb_dim,   egl_igl_thick, "orange"),
        ("port + radial union", rgb_dim,  egl_full,      "red"),
    ]
    for k, (t, im, mm, c) in enumerate(panels):
        ax[k].imshow(im)
        if mm is not None:
            ax[k].contour(mm, levels=[0.5], colors=[c], linewidths=0.5)
        ax[k].set_title(t); ax[k].axis("off")
    fig.tight_layout(); fig.savefig(OUT / "01_egl_sources.png", dpi=140)
    plt.close(fig)
    print(f"figures saved to {OUT}")


if __name__ == "__main__":
    main()
