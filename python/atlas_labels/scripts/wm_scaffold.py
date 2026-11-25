"""
White-matter scaffold for P7 cerebellum label maps (9-class MATLAB scheme).

Every fissure runs between two neighbouring folia, i.e. between two white-matter
branches. So inside the tissue, a fissure lies on the midline of (tissue - WM)
where white matter is on BOTH sides: the distance to white matter (dW) is smaller
than the distance to the outside surface (dO). The midline of a single folium's
cortex band has the surface on one side and fails that test.

    fissure_lines(lab, px) -> dict(lines, dW, dO, filled, wm)

The fused EGL of a closed fissure sits on that line, from the depth where the EGL
starts (dW = distance from white matter to the EGL inner edge, read locally off
the existing labels) outwards.
"""
import numpy as np
from scipy import ndimage as ndi
from skimage import morphology

ST = np.ones((3, 3), bool)
DEFAULTS = dict(wm_min_um2=2000.0, wm_smooth_um=0.0, min_do_um=40.0, min_dw_um=100.0,
                side_margin_um=None, min_line_um=0.0, tears_outside=True)


def wm_mask(lab, px, prm):
    wm = np.isin(lab, (6, 8))
    r = int(round(prm["wm_smooth_um"] / px))
    if r > 0:
        wm = ndi.binary_closing(np.pad(wm, r + 1), structure=morphology.disk(r))[r + 1:-r - 1, r + 1:-r - 1]
        wm = ndi.binary_opening(wm, structure=morphology.disk(r))
    return morphology.remove_small_objects(wm, max_size=int(prm["wm_min_um2"] / px / px))


def fissure_lines(lab, px, prm=None):
    prm = {**DEFAULTS, **(prm or {})}
    tissue = (lab > 0) & ((lab != 9) if prm["tears_outside"] else True)
    filled = ndi.binary_fill_holes(tissue)
    wm = wm_mask(lab, px, prm) & filled
    dO = ndi.distance_transform_edt(filled) * px
    dW = ndi.distance_transform_edt(~wm) * px
    mid = morphology.skeletonize(filled & ~wm)
    lines = mid & (dO > prm["min_do_um"]) & (dW >= prm["min_dw_um"])
    if prm["side_margin_um"] is not None:          # white matter on both sides (stricter)
        lines &= dO > dW + prm["side_margin_um"]
    if prm["min_line_um"]:
        lines = morphology.remove_small_objects(lines, max_size=int(prm["min_line_um"] / px), connectivity=2)
    return dict(lines=lines, dW=dW, dO=dO, filled=filled, wm=wm, mid=mid)


def score_lines(pred, ref, px, tol_um=(20.0, 30.0)):
    dR = ndi.distance_transform_edt(~ref) * px
    dP = ndi.distance_transform_edt(~pred) * px
    out = {}
    for t in tol_um:
        out[f"recall@{t:.0f}"] = round(float((dP[ref] <= t).mean()), 3) if ref.any() else None
        out[f"precision@{t:.0f}"] = round(float((dR[pred] <= t).mean()), 3) if pred.any() else None
    out["pred_mm"] = round(float(pred.sum() * px / 1000), 1)
    out["ref_mm"] = round(float(ref.sum() * px / 1000), 1)
    return out
