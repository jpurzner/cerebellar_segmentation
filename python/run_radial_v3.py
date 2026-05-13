"""v3 radial: cap ray length + voronoi-style back-projection.

Improvements over v2:
  - filter spokes longer than RAY_MAX_REASONABLE_UM (default 250 um)
  - for each cerebellum pixel: find nearest spoke, project pixel onto spoke
    direction to get depth, look up layer from spoke's boundary table.
  - this replaces the brittle "draw spoke line + dilate" approach.
"""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage as ndi
import tifffile
from scipy.spatial import cKDTree
from skimage import exposure
from radial_spokes import build_spokes
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
OUT = Path(__file__).resolve().parent / "figs" / "v_radial_v3"
OUT.mkdir(parents=True, exist_ok=True)

PX_UM = 0.5119049
RAY_MAX_REASONABLE_UM = 250.0   # drop spokes longer than this
SAMPLE_STEP_UM = 1.0


def to_unit(im):
    im = im.astype(np.float32); lo, hi = im.min(), im.max()
    return (im - lo) / max(hi - lo, 1e-9)


def clahe(im, clip=0.01):
    return exposure.equalize_adapthist(im, clip_limit=clip, nbins=256)


def extract_egl_gt(seg_rgb):
    s = seg_rgb.astype(int)
    def col(target): return np.all(np.abs(s - np.array(target)) < 16, axis=-1)
    iEGL = col([0, 128, 255]); oEGL = col([0, 255, 255])
    return iEGL | oEGL, iEGL, oEGL


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
    cereb = z['mask'].astype(bool)
    igl   = z['c_final'].astype(bool)
    dwl   = z['dwl_final'].astype(bool)
    port_egl_a = z['a_final'].astype(bool)
    port_egl_b = z['b_final'].astype(bool)

    sl = build_spokes(cereb, igl, pixel_size_um=PX_UM, boundary_smooth_um=3.0)
    n_total = len(sl.pia_points)

    # filter to "reasonable" spokes
    keep = sl.valid & (sl.lengths_um <= RAY_MAX_REASONABLE_UM) & \
           (sl.lengths_um >= 30.0)  # also drop very short spokes
    print(f"spokes total {n_total}; valid {sl.valid.sum()}; "
          f"reasonable (30..{RAY_MAX_REASONABLE_UM:.0f}um): {keep.sum()}")

    pia_pts  = sl.pia_points[keep]
    igl_pts  = sl.igl_endpoints[keep]
    dirs     = sl.spoke_dirs[keep]
    lengths  = sl.lengths_um[keep]

    # --- 1D layer detection per spoke ---
    print("sampling 1D profiles + detecting layer boundaries")
    boundaries = np.full((len(pia_pts), 4), np.nan)  # egl_end, iegl_start, iegl_end, igl_start (um)
    for i in range(len(pia_pts)):
        prof = sample_along_spoke(chans, pia_pts[i], igl_pts[i],
                                   SAMPLE_STEP_UM, PX_UM)
        layers = detect_layers_1d(prof, SAMPLE_STEP_UM,
                                   channels={"p27": 0, "neun": 1, "dapi": 2})
        boundaries[i] = [layers['egl_end_um'], layers['iegl_start_um'],
                          layers['iegl_end_um'], layers['igl_start_um']]

    # --- arc-length smoothing of boundaries (no sudden jumps) ---
    # The spokes are ordered along arc length, so we can smooth each column
    # with a 1D gaussian over indices.
    sigma_pts = 5  # neighbours = 5 spokes ~ 40um arc length
    boundaries_sm = np.zeros_like(boundaries)
    for c in range(4):
        boundaries_sm[:, c] = ndi.gaussian_filter1d(boundaries[:, c],
                                                     sigma=sigma_pts, mode='wrap')

    # --- back-projection via nearest-spoke + projected depth ---
    print("back-projecting via voronoi (nearest spoke + projected depth)")
    set_bin = np.zeros((H, W), dtype=np.uint8)

    # cereb pixel indices
    yy, xx = np.where(cereb)
    pixel_pts = np.column_stack([yy, xx]).astype(np.float32)
    tree = cKDTree(pia_pts)
    dist, idx = tree.query(pixel_pts, k=1)

    # project pixel onto its nearest spoke direction
    pia_for_px = pia_pts[idx]
    dir_for_px = dirs[idx]
    v = pixel_pts - pia_for_px      # (M, 2)
    s_px = (v * dir_for_px).sum(axis=1) * PX_UM  # depth (um), can be negative
    # also compute perpendicular distance to the spoke line
    perp = (v[:, 0] * dir_for_px[:, 1] - v[:, 1] * dir_for_px[:, 0]) * PX_UM

    # spoke length per pixel
    L_for_px = lengths[idx]
    # bounds per pixel
    egl_end, iegl_s, iegl_e, igl_s = boundaries_sm[idx].T

    # Strategy:
    #   - Port handles everything inside the cerebellum (its IGL/DWL/PCL/ML output).
    #     We default each cereb pixel to the port's label.
    #   - Radial detector ONLY refines the EGL (oEGL + iEGL) at the periphery,
    #     where the spoke's projected depth s ∈ [0, egl_end_um]. This is
    #     where the radial 1D profile is reliable.

    s = s_px; ee = egl_end; ie_s = iegl_s; ie_e = iegl_e

    # default label per pixel from the port
    px_labels = np.full(len(pixel_pts), 5, dtype=np.uint8)  # default = ML (broad layer)
    igl_pix = igl[yy, xx]
    dwl_pix = dwl[yy, xx]
    px_labels[dwl_pix] = 6
    px_labels[igl_pix] = 4

    # radial overrides ONLY for EGL (s ∈ [0, egl_end])
    is_egl  = (s >= 0) & (s < ee)
    is_iegl = is_egl & (s >= ie_s) & (s <= ie_e)
    is_oegl = is_egl & ~is_iegl
    px_labels[is_oegl] = 3
    px_labels[is_iegl] = 2

    set_bin[yy, xx] = px_labels

    # === GT comparison ===
    seg_rgb = tifffile.imread(SEG_GT_TIF)
    egl_gt, iegl_gt, oegl_gt = extract_egl_gt(seg_rgb)
    egl_pred = (set_bin == 2) | (set_bin == 3)
    inter = (egl_pred & egl_gt).sum()
    union = (egl_pred | egl_gt).sum()
    p = inter / max(egl_pred.sum(), 1); r = inter / max(egl_gt.sum(), 1)
    iou = inter / max(union, 1)
    print(f"\n=== EGL match vs MATLAB GT ===")
    print(f"  GT     = {egl_gt.sum()*(PX_UM/1000)**2:.3f} mm^2")
    print(f"  pred   = {egl_pred.sum()*(PX_UM/1000)**2:.3f} mm^2")
    print(f"  P      = {p:.3f}  R = {r:.3f}  IoU = {iou:.3f}")

    LABEL_COLORS = np.array([
        [0,   0,   0],  [0,   0, 128], [0, 128, 255], [0, 255, 255],
        [255, 255, 0], [128, 255, 128],[255, 128, 0], [255,   0, 0],
    ], dtype=np.uint8)
    seg_color = LABEL_COLORS[set_bin]

    fig, ax = plt.subplots(1, 3, figsize=(22, 9))
    ax[0].imshow(rgb); ax[0].set_title("RGB"); ax[0].axis("off")
    ax[1].imshow(seg_color); ax[1].set_title(f"v3 radial (n_spokes={keep.sum()})")
    ax[1].axis("off")
    ax[2].imshow(rgb * 0.5)
    ax[2].contour(egl_gt,  levels=[0.5], colors=["lime"],   linewidths=0.4, alpha=0.85)
    ax[2].contour(egl_pred, levels=[0.5], colors=["yellow"], linewidths=0.4, alpha=0.85)
    ax[2].set_title(f"GT (lime) vs pred (yellow)  IoU={iou:.2f}")
    ax[2].axis("off")
    fig.tight_layout(); fig.savefig(OUT / "01_seg_vs_gt.png", dpi=140); plt.close(fig)

    print(f"figures saved to {OUT}")


if __name__ == "__main__":
    main()
