"""
Re-draw the IGL of an IF section from its NeuN channel, anchored on the white matter.

    from igl_refit import refit_igl
    new_labels, info = refit_igl(labels, neun, px_um)

Labels: 9-class MATLAB scheme (0 bg, 1 cereb, 2 iEGL, 3 oEGL, 4 IGL, 5 ML, 6 DWL,
7 PCL, 8 DCN, 9 anomaly). Only cortex pixels between the EGL and the white matter
(current labels 1, 4, 5, 7) can change; EGL, white matter, DCN, tears stay.

  1. NeuN smoothed (sigma 3 um), divided by the median NeuN of the section's current ML.
  2. Threshold = geometric mean of the current IGL and ML medians (the section's own
     contrast, so acquisition batches with different brightness need no tuning).
  3. Bright pixels: opening 2 um, closing 4 um, holes < 300 um2 filled.
  4. Bright pieces are kept if they come within 30 um of the white matter or are
     larger than 5000 um2 (granule cells sit in a shell around each white-matter branch).
  5. Kept bright -> IGL. Pixels that were IGL and are no longer bright -> ML;
     current PCL pixels that are not bright stay PCL (Purkinje cells are NeuN-dim).
"""
import numpy as np
from scipy import ndimage as ndi
from skimage import morphology

DEFAULTS = dict(sigma_um=3.0, open_um=2.0, close_um=4.0, hole_um2=300.0, near_wm_um=30.0,
                big_um2=5000.0)
ST = np.ones((3, 3), bool)


def refit_igl(lab, neun, px, prm=None):
    prm = {**DEFAULTS, **(prm or {})}
    lab = lab.copy()
    wm = np.isin(lab, (6, 8))
    zone = np.isin(lab, (1, 4, 5, 7))
    ml, igl = lab == 5, lab == 4
    n = ndi.gaussian_filter(neun.astype(np.float32), prm["sigma_um"] / px)
    ref = np.median(n[ml]) if ml.any() else np.median(n[zone])
    n = n / max(ref, 1e-6)
    t = float(np.sqrt(np.median(n[igl]) * np.median(n[ml])))
    bright = zone & (n > t)
    bright = ndi.binary_opening(bright, structure=morphology.disk(max(1, int(round(prm["open_um"] / px)))))
    rc = max(1, int(round(prm["close_um"] / px)))
    bright = ndi.binary_closing(np.pad(bright, rc + 1), structure=morphology.disk(rc))[rc + 1:-rc - 1, rc + 1:-rc - 1]
    holes = ndi.binary_fill_holes(bright) & ~bright
    hl, nh = ndi.label(holes, structure=ST)
    if nh:
        ha = np.bincount(hl.ravel()) * px * px
        bright |= np.isin(hl, np.nonzero(ha < prm["hole_um2"])[0][1:]) & zone
    bright &= zone
    cc, nc = ndi.label(bright, structure=ST)
    keep = np.zeros_like(bright)
    if nc:
        area = np.bincount(cc.ravel()) * px * px
        dW = ndi.distance_transform_edt(~wm) * px
        near = np.unique(cc[(dW <= prm["near_wm_um"]) & bright])
        ok = set(near[near > 0].tolist()) | set(np.nonzero(area >= prm["big_um2"])[0][1:].tolist())
        keep = np.isin(cc, list(ok))
    out = lab.copy()
    out[keep] = 4
    lost = igl & ~keep
    out[lost] = 5
    info = dict(threshold=round(t, 3), ref_ml_neun=float(ref),
                igl_before_mm2=round(float(igl.sum() * px * px / 1e6), 3),
                igl_after_mm2=round(float((out == 4).sum() * px * px / 1e6), 3),
                gained_mm2=round(float((keep & ~igl).sum() * px * px / 1e6), 3),
                lost_mm2=round(float(lost.sum() * px * px / 1e6), 3))
    return out, info


def igl_pieces(lab, px, min_um2=1000.0):
    cc, n = ndi.label(lab == 4, structure=ST)
    a = np.bincount(cc.ravel())[1:] * px * px
    a = a[a >= min_um2]
    return int(a.size), round(float(100 * a.max() / a.sum()), 1) if a.size else 0.0
