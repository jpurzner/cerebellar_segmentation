"""
EGL repair (2026-11-21, final = version 5, signed off by the user after five review rounds). The EGL
label is not marked all the way along the surface, and the flattening leans on it. The IGL (current
labels) and the stains repair it. Everything runs on the 4 um label grid (flatten.labels codes:
0 bg, 1 unassigned, 2 EGL, 3 ML, 4 IGL, 5 WM, 6 anomaly); columns are sampled every 2 um.

Input: flatten.labels -> `best` (reviewed controls, the user's corrected cKO sections, MATLAB for
s1_5) and the stain features of igl_detect (nuclei, NeuN, p27; local and section-median versions).

1. stain_maps. Per section, Youden thresholds on its own labels: stain EGL = bright nuclei + low NeuN,
   stain IGL = high NeuN + bright nuclei, stain ML = dim nuclei + low NeuN (EGL / ML also need absolute
   NeuN below the EGL-IGL midpoint). `egl_loose` = nuclei above the IGL / ML split with low NeuN.
   `nonegl` (meninges, pia) = p27 < 0.6 x and NeuN < 0.7 x the section's ML level.
2. cuts -> the reference line. Fissure midlines down the fused stain EGL (flatten.midline), each
   continued from its deep end along its own direction (extend_ends: <= 30 deg per step, 45 total,
   stops 30 um from the IGL or at another cut; may cross WM labels only where the stain shows EGL).
   Fallback: cuts halfway between facing IGL bands (igl_midline), kept only for fissures the stain did
   not find, with >= 20 % stain EGL and < 10 % WM beside them. The tissue cut along these lines gives
   the reference contour (flatten.reference) and the columns (flatten.streamlines).
   v6, pia_snap: in a closed fissure two EGLs face each other across a thin pia / arachnoid. As the R
   layer_quant pipeline splits its fissure profiles at the p27 minimum, each cut point moves to the
   p27 + NeuN minimum across the fissure nearest it (EGL on both sides, below the ML level), smoothed
   along the cut. The cut before the move is kept as mid_pre.
   v8 (the user's flatten v5 notes, deep furrows missed): flatten.single_line ate 24 um off every cut end at each
   of its 3-4 passes, so the cuts stopped short of the furrow base. Now an extension that ran into the ML is cut
   back (trim_ends), each deep end is walked on down the middle of the fused EGL while it looks like two EGLs back
   to back (walk_base), and the last passes (pia_snap, flatten.prune_cuts) keep the ends (flatten.final_prm). Tiny
   holes inside the tissue are filled for cutting (flatten.cut_tissue); depth for dropping shallow cuts is measured
   past narrow open furrows (flatten.cut_depth).
   v9 (the user's flatten v6 notes: the line followed arachnoid, nonsense perimeters): the tissue mask (MATLAB: p27 +
   DAPI thresholded, holes filled) holds arachnoid flaps and the meninges packing open fissures. pia_trim removes them
   before cutting: low p27 + NeuN (the antibody background stops at the pia, as in the R pia_prunner), not DAPI-bright
   (the EGL is), not labelled ML / IGL / WM, touching the outside, mostly against EGL (not a white-matter cut surface).
   Cuts and the reference line run on tissue & ~arach (export: final/<short>_pia.tif, pia_tissue(short)); EGL labels
   in the trimmed tissue become unassigned (not on the gold).
   v10 (the user's flatten v8 notes, cKO: the line ran round the meninges filling an open fissure, s2_2 at 8.6, s3_2 at
   9.8): where the cKO EGL is missing the meninges lie on ML or IGL, so a piece also counts as arachnoid when EGL / ML /
   IGL (not white matter, the cut-surface sign) make up half its labelled border and unassigned / anomaly cells at most
   half of the whole border (on the gold s1_2 and in s4_A / s5_A those are traced meninges), unless more than 30 % of its cells
   reach the background only across EGL or trimming it cuts EGL off the section (s5_C at 3.3: a dim strip between a
   folium's ML and its EGL).
3. relabel (not on the gold s1_2). Within 60 um of the line: IGL -> EGL where `egl_loose`, in pieces
   touching EGL or the surface. Deeper than 24 um: EGL -> IGL where the stain shows IGL. IGL -> ML only on
   IGL the earlier repairs added to controls.
4. Column order EGL -> ML -> IGL -> WM (not on gold). IGL starting within 60 um of the line, with
   >= 12 um of mostly ML (no WM) and IGL again below it, is inner EGL -> EGL (igl_above_ml); leftover
   IGL islands < 2000 um2 lying within 60 um take the nearest other label.
5. Column fit: the log nuclei profile (averaged +-8 um along the line) fitted with three flat segments
   EGL / ML / IGL; z1 = EGL bottom, z2 = IGL top held -24 / +50 um around the labelled IGL top.
   Smoothed along the line (40 um). EGL present = EGL - ML >= 0.8, EGL - IGL >= -0.4, z1 >= 12 um.
   Labelled IGL above a fitted IGL top >= 16 um deeper becomes ML unless the stain shows IGL (igl_in_ml).
6. EGL add (add-only). Gap = EGL label missing or < 12 um near the line, unless labelled EGL lies
   further down the column above its IGL (outer surface only: something covers it). Fill = gap with
   EGL present and IGL below; other gaps <= 300 um are bridged from the thickness either side.
   Thickness capped at min(60, 1.5 x median labelled EGL), never past the first IGL, >= 16 um from the
   IGL. Which cells may turn EGL depends on the column's start: outer surface -> unassigned, ML, WM,
   anomaly, but not `nonegl`; cut -> unassigned, ML; cut extension -> those plus WM the stain calls EGL.
   Added pieces < 200 um2 not touching EGL are dropped.
7. Deep EGL (not on gold): EGL > 90 um below the line, or specks < 2000 um2 more than 24 um down,
   take the nearest non-EGL label; with the `nonegl` signature they become WM.

    python scripts/egl_repair.py s5_B s2_4      # some sections -> outputs/egl_repair/<short>.npz
    python scripts/egl_repair.py                # all, + summary.json
    python scripts/egl_repair.py export         # final label grids -> outputs/egl_repair/final/

Use the result through repaired(short) -> (labels, tissue), the same grid and codes as flatten.labels.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from skimage import morphology
from skimage.graph import MCP_Geometric

sys.path.insert(0, str(Path(__file__).parent))
import flatten as F  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/egl_repair"
G, STEP = F.G, F.STEP
GOLD = "s1_2"
PRM = dict(F.PRM, igl_min_um2=2000.0, seed_min_um=30.0, seed_max_um=200.0, seed_bg_um=30.0, seed_cos=-0.5,
           seed_jump_um=40.0, seed_layer_um=24.0, fallback_um=24.0, ext_near_um=80.0, ext_end_bg_um=60.0, ext_back_um=60.0, ext_max_um=400.0, turn_deg=30.0, max_turn_deg=45.0,
           egl_cap_factor=1.5, egl_cap_max_um=60.0, ml_room_um=16.0, deep_egl_um=90.0, speck_egl_um2=2000.0, egl_feature=3, neun_feature=10, neun_abs_feature=11,
           egl_min_um2=1000.0, igl2egl_um=60.0, egl2igl_um=24.0, igl2ml_um=120.0, swap_min_um2=800.0, egl2igl_min_um2=2000.0,
           bridge_egl_um=300.0,
           prof_avg_um=8.0, z1_min_um=6.0, ml_min_um=8.0, z2_lo_um=24.0, z2_hi_um=50.0, z2_free_max_um=200.0,
           igl_seg_um=80.0, fit_depth_um=240.0, z1_max_um=150.0, contrast=0.8, egl_vs_igl=-0.4, z1_present_um=12.0,
           thin_um=12.0, z1_med_um=40.0,
           # column order EGL -> ML -> IGL -> WM (v4)
           order_start_um=60.0, order_run_um=60.0, order_ml_um=12.0, order_deep_um=16.0, order_ml_frac=0.7,
           ml_gap_um=16.0, ml_keep_um=8.0, igl_island_um2=2000.0,
           # v5: meninges, fallback cuts, furrow extensions
           p27_abs_feature=17, nonegl_p27=0.6, nonegl_neun=0.7, nonegl_min_um2=200.0,
           fallback_egl_frac=0.2, fallback_wm_frac=0.1, use_egl_under=True,
           # v6: the cut moves onto the pia between the two fused EGLs (R layer_quant split)
           pia_snap=True, pia_win_um=32.0, pia_side_um=24.0, pia_along_um=8.0, pia_depth=0.2, pia_abs_margin=0.1,
           pia_bg_um=20.0, pia_smooth_um=16.0, pia_egl_frac=0.3,
           # v7 (the user's flatten v4 notes): the pia dip is measured against the peaks on either side (the two
           # inner EGLs, high p27), no ML-level test; the IGL the user traced in the controls is kept; anomaly cells
           # the stain calls EGL may become EGL from a cut; open furrows continue as a cut along EGL hidden under labels
           pia_flank="peak", pia_abs=False, protect_traced=True, cut_takes_anom_egl=True,
           tip_ext=False, tip_narrow_um=60.0, tip_back_um=60.0, tip_gap_um=160.0, tip_egl_frac=0.3, tip_cut_near_um=16.0,
           tip_bg_steps=10,
           # v8 (the user's flatten v5 notes, deep furrows missed): cut ends kept (flatten keep_ends) and ended at the
           # furrow base
           end_trim=True, end_frac=0.5, end_max_um=12.0, end_win_um=80.0,
           base_bg_um=20.0, base_back_um=40.0, base_turn_deg=30.0, base_max_turn_deg=60.0, base_frac=0.5, base_max_um=300.0, base_side=2.5, base_hi=1.6,
           # v9 (the user's flatten v6 notes: the line follows arachnoid): the pial surface is the outer face of the EGL;
           # tissue outside it is trimmed before cutting (pia_trim)
           pia_trim=True, nuc_abs_feature=4, arach_thr=-0.6, arach_bright=0.5, arach_close=1, arach_contact=0.5,
           arach_island_um2=20000.0, arach_island_cortex_um2=2000.0,
           # v10 (the user's flatten v8 notes on the cKO sections): meninges on ML / IGL where the EGL is missing
           arach_relax=True, arach_relax_egl_um2=200.0, arach_behind=0.3, arach_relax_unl=0.5)


def _connect(seeds, cost, bg):
    """join every piece of seeds that does not touch the background to it by the cheapest path."""
    paths = np.zeros_like(seeds)
    starts = np.argwhere(np.isfinite(cost) & ndi.binary_dilation(bg))
    if not len(starts) or not seeds.any():
        return paths
    m = MCP_Geometric(cost, fully_connected=True)
    cum, _ = m.find_costs(starts)
    comp, _ = ndi.label(seeds, structure=np.ones((3, 3), bool))
    touch = ndi.binary_dilation(bg)
    for k, sl in enumerate(ndi.find_objects(comp), start=1):
        cells = np.argwhere(comp[sl] == k) + [sl[0].start, sl[1].start]
        if touch[cells[:, 0], cells[:, 1]].any():
            continue
        cvals = cum[cells[:, 0], cells[:, 1]]
        if not np.isfinite(cvals).any():
            continue
        for r_, c_ in m.traceback(tuple(cells[int(np.argmin(cvals))])):
            paths[r_, c_] = True
    return paths


def _prune_shallow(mid, tissue, prm=PRM):
    comp, nc = ndi.label(mid, structure=np.ones((3, 3), bool))
    if not nc:
        return mid
    deep = np.asarray(ndi.maximum(F.cut_depth(tissue, prm), comp, index=np.arange(1, nc + 1)))
    return np.r_[False, deep >= prm["min_cut_um"]][comp]


def igl_midline(igl, zone, tissue, prm=PRM):
    """closed-fissure midlines from facing IGL bands: cells of the outer zone whose nearest IGL cell lies in
    the opposite direction from their neighbour's, at least seed_min_um from the IGL (room for the ML),
    with ML or EGL within seed_layer_um, joined to the background along the middle of the zone."""
    bg = ~tissue
    dist, ind = ndi.distance_transform_edt(~igl, return_indices=True)
    dbg = ndi.distance_transform_edt(tissue)
    rr, cc = np.indices(igl.shape)
    v = np.stack([ind[0] - rr, ind[1] - cc]).astype(np.float32)
    cand = zone & (dist * G >= prm["seed_min_um"]) & (dist * G <= prm["seed_max_um"]) & (dbg * G >= prm["seed_bg_um"])
    seeds = np.zeros_like(cand)
    for dy, dx in ((0, 1), (1, 0), (1, 1), (1, -1)):
        ca, cb, sa, sb = F._pair(cand, dy, dx)
        va, vb, _, _ = F._pair(v, dy, dx)
        ia, ib, _, _ = F._pair(ind.astype(np.float32), dy, dx)
        dot = (va * vb).sum(0)
        na, nb = np.hypot(*va), np.hypot(*vb)
        hit = ca & cb & (dot < prm["seed_cos"] * na * nb) & (np.hypot(*(ia - ib)) * G > prm["seed_jump_um"])
        seeds[sa] |= hit
        seeds[sb] |= hit
    seeds = morphology.remove_small_objects(seeds, max_size=2, connectivity=2)
    cost = np.where(zone | seeds, 1.0 / (1.0 + dist) ** 2, np.inf)
    cost[seeds] = np.minimum(cost[seeds], 0.02)
    mid = F.single_line(seeds | _connect(seeds, cost, bg), bg, prm) & tissue
    return _prune_shallow(mid, tissue, prm), dist


def _youden(a, b):
    """threshold best separating a (higher) from b."""
    qs = np.quantile(np.r_[a, b], np.linspace(0.02, 0.98, 97))
    return max((np.mean(a > q) - np.mean(b > q), q) for q in qs)[1]


def _clean(m, min_um2):
    return morphology.remove_small_objects(ndi.binary_opening(m), max_size=int(min_um2 / G ** 2), connectivity=2)


def stain_maps(short, lab, tissue, prm=PRM):
    """EGL / IGL / ML from the stains: nuclei against the surrounding 150 um (smoothed 12 um) separates the
    bright EGL and IGL from the dim ML; NeuN against the surrounding 150 um separates the IGL (NeuN+) from
    the EGL. Thresholds are each section's own best split of its labels (most labels are right)."""
    Fx = F.D.load(short)[0]
    nuc, neun = Fx[prm["egl_feature"]], Fx[prm["neun_feature"]]
    E, M, I = lab == 2, lab == 3, lab == 4
    t = dict(nuc_egl_ml=float(_youden(nuc[E], nuc[M])), nuc_igl_ml=float(_youden(nuc[I], nuc[M])),
             neun_igl_egl=float(_youden(neun[I], neun[E])))
    # absolute NeuN (against the section median): a cell only counts as EGL or ML if its NeuN is below the
    # midpoint of the section's EGL and IGL levels, and as IGL if above it, so a relative dip in good IGL is
    # not taken for EGL (bleached IGL is below that midpoint, and there the nuclei stain alone decides)
    absn = Fx[prm["neun_abs_feature"]]
    t["neun_abs_mid"] = float((np.median(absn[E]) + np.median(absn[I])) / 2)
    lowN, highN = absn < t["neun_abs_mid"], absn >= t["neun_abs_mid"]
    egl = tissue & (nuc > t["nuc_egl_ml"]) & (neun < t["neun_igl_egl"]) & lowN
    igl = tissue & (neun >= t["neun_igl_egl"]) & (nuc > t["nuc_igl_ml"]) & highN
    ml = tissue & (nuc < min(t["nuc_egl_ml"], t["nuc_igl_ml"])) & (neun < t["neun_igl_egl"]) & lowN
    # granule-layer nuclei (brighter than the IGL / ML split) without NeuN: EGL the relative nuclei feature
    # misses next to very bright EGL (v5, s3_2)
    egl_loose = tissue & (nuc > t["nuc_igl_ml"]) & (neun < t["neun_igl_egl"]) & lowN
    # not cortex: p27 and NeuN both well below the section's ML (meninges, pia; v5). The EGL has more p27 than the ML.
    p27a = Fx[prm["p27_abs_feature"]]
    t["p27_ml"], t["neun_ml"] = float(np.median(p27a[M])), float(np.median(absn[M]))
    nonegl = tissue & (p27a < t["p27_ml"] + np.log(prm["nonegl_p27"])) & (absn < t["neun_ml"] + np.log(prm["nonegl_neun"]))
    nonegl = morphology.remove_small_objects(ndi.binary_opening(nonegl), max_size=int(prm["nonegl_min_um2"] / G ** 2), connectivity=2)
    return dict(egl=_clean(egl, prm["egl_min_um2"]), igl=_clean(igl, prm["egl_min_um2"]), ml=_clean(ml, prm["egl_min_um2"]),
                egl_loose=_clean(egl_loose, prm["swap_min_um2"]), nonegl=nonegl, nuc=nuc, neun=neun, p27a=p27a, neun_a=absn,
                nuc_a=Fx[prm["nuc_abs_feature"]], t=t)


def pia_trim(lab, tissue, S, prm=PRM, gold=False):
    """v9: tissue outside the pial surface (arachnoid, meninges, the fill of open fissures the tissue mask closed
    over), which the reference line would otherwise run around. The pia is the outer face of the EGL. p27 and NeuN
    give the tissue a background that stops at the pia (the drop the R layer_quant pia_prunner looks for), while the
    EGL is bright on DAPI. Arachnoid cells: p27 + NeuN (log, each against the section's ML) below arach_thr, nuclei
    (DAPI against the section median) below arach_bright of the way from the ML to the EGL level, not labelled ML,
    IGL or white matter (on the gold, not labelled at all). Closed by arach_close cells (scattered p27+ pial cells),
    in pieces that touch the background and whose border with the rest of the tissue is at least arach_contact EGL
    (label or stain): arachnoid lies on the EGL, while low-background tissue against white matter is a cut surface.
    Tissue this cuts off from the main body goes too unless it is >= arach_island_um2 or holds
    >= arach_island_cortex_um2 of ML / IGL / white matter. Returns the trimmed cells."""
    k8 = np.ones((3, 3), bool)
    M, E = lab == 3, lab == 2
    nm, pm = float(np.median(S["neun_a"][M])), float(np.median(S["p27a"][M]))
    um, ue = float(np.median(S["nuc_a"][M])), float(np.median(S["nuc_a"][E]))
    nb = (S["neun_a"] - nm) + (S["p27a"] - pm) < prm["arach_thr"]
    bright = S["nuc_a"] >= um + prm["arach_bright"] * (ue - um)
    cand = tissue & nb & ~bright & ~np.isin(lab, (2, 3, 4, 5) if gold else (3, 4, 5))
    if prm["arach_close"]:
        cand = ndi.binary_closing(cand, k8, iterations=int(prm["arach_close"])) & tissue & ~np.isin(lab, (2, 3, 4, 5) if gold else (3, 4, 5))
    comp, nc = ndi.label(cand, structure=k8)
    if not nc:
        return np.zeros_like(tissue)
    bgd = ndi.binary_dilation(~tissue, k8)
    touch = np.unique(comp[bgd & cand])
    eglish = E | S["egl"]
    keep, relax = [], []
    for k, sl in enumerate(ndi.find_objects(comp), start=1):
        if k not in touch:
            continue
        sl = tuple(slice(max(a.start - 2, 0), a.stop + 2) for a in sl)
        m = comp[sl] == k
        ring = ndi.binary_dilation(m, k8) & ~m & tissue[sl] & ~cand[sl]
        if not ring.any() or eglish[sl][ring].mean() >= prm["arach_contact"]:
            keep.append(k)
        elif prm.get("arach_relax"):
            # v10 (the user's flatten v8 notes, cKO: s2_2 at 8.6, s3_2 at 9.8): where the cKO EGL is missing, the meninges
            # filling an open fissure lie on ML or even IGL, so the EGL share of the border falls below arach_contact.
            # Only white matter marks a cut surface: kept as arachnoid when EGL / ML / IGL make up arach_contact or more
            # of the labelled border (unassigned / anomaly not counted), unless trimming it cuts EGL off the section.
            pial = (eglish[sl] | np.isin(lab[sl], (3, 4)))[ring].sum()
            cutf = ((lab[sl] == 5) & ~eglish[sl])[ring].sum()
            unl = (np.isin(lab[sl], (0, 1, 6)) & ~eglish[sl])[ring].mean()   # unassigned / anomaly border: not this case
            if pial >= prm["arach_contact"] * max(pial + cutf, 1) and unl <= prm["arach_relax_unl"]:
                relax.append(k)
    arach = np.isin(comp, keep)
    cx_lab = np.isin(lab, (3, 4, 5))

    def islands(a):
        rest = tissue & ~a
        comp2, n2 = ndi.label(rest, structure=k8)
        if n2 <= 1:
            return np.zeros_like(a)
        idx = np.arange(1, n2 + 1)
        size = np.asarray(ndi.sum(rest, comp2, index=idx)) * G * G
        cx = np.asarray(ndi.sum(cx_lab, comp2, index=idx)) * G * G
        ok = (size >= prm["arach_island_um2"]) | (cx >= prm["arach_island_cortex_um2"])
        return rest & ~np.r_[False, ok][comp2]
    base_lost = islands(arach)
    if relax:
        from skimage.draw import line as draw_line
        _, (br, bc) = ndi.distance_transform_edt(tissue, return_indices=True)
    for k in relax:
        # outside the EGL: the straight way from each cell to the nearest background crosses no EGL for at least
        # 1 - arach_behind of the cells (a dim strip behind the EGL band is not arachnoid)
        pts = np.argwhere(comp == k)
        behind = 0
        for r_, c_ in pts:
            rr, cc = draw_line(r_, c_, br[r_, c_], bc[r_, c_])
            behind += bool((eglish[rr, cc] & ~cand[rr, cc]).any())
        if behind > prm["arach_behind"] * len(pts):
            continue
        trial = arach | (comp == k)
        lost = islands(trial) & ~base_lost
        if (lost & E).sum() * G * G <= prm["arach_relax_egl_um2"]:
            arach = trial
            base_lost = base_lost | lost
    return arach | base_lost


def pia_snap(mid, lab, tissue, S, prm=PRM):
    """move a fissure cut onto the pia between the two fused EGLs, as the R layer_quant pipeline splits its
    fissure profiles: p27 and NeuN have no tissue background in the pia / arachnoid, so the pia is the
    minimum of p27 + NeuN across the fissure nearest the centre. At each skeleton point of the cut >= pia_bg_um
    from the background (closer, the fissure is open), the profile across the cut (+-pia_win_um, averaged over
    +-pia_along_um along it) is searched for local minima on EGL / unassigned / anomaly cells with EGL on
    both sides (>= pia_egl_frac within 4-pia_side_um), at least pia_depth (log) below both flanks (8 um and
    further out) and with p27 + NeuN below the section's ML level (+ pia_abs_margin); the one nearest the cut
    wins. The shifts are averaged along the cut over pia_smooth_um and taper where under half the points
    found the pia; points without one stay put. Returns the new cut and the per-point record."""
    from scipy.spatial import cKDTree
    from skimage.draw import line as draw_line
    p27, neun = S["p27a"], S["neun_a"]
    ml = lab == 3
    ml_level = float(np.median(p27[ml]) + np.median(neun[ml]))
    q = p27 + neun
    eglish = (lab == 2) | S["egl"]
    bandok = np.isin(lab, (1, 2, 6)) | S["egl"]
    bg = ~tissue
    dbg = ndi.distance_transform_edt(tissue) * G
    sk = morphology.skeletonize(mid)
    pts = np.argwhere(sk)
    info = dict(found=np.zeros(len(pts), bool), deep=np.zeros(len(pts), bool), moved_um=np.zeros(len(pts)))
    if len(pts) < 4:
        return mid, info
    tree = cKDTree(pts)
    w, sd = prm["pia_win_um"], prm["pia_side_um"]
    T = np.arange(-w - sd, w + sd + 1, 1.0)
    alongs = np.arange(-prm["pia_along_um"], prm["pia_along_um"] + 0.1, G) / G
    H, W = lab.shape
    shift = np.zeros((len(pts), 2))
    for k, p in enumerate(pts):
        if dbg[tuple(p)] < prm["pia_bg_um"]:
            continue
        nb = pts[tree.query_ball_point(p, 3.5)]
        if len(nb) < 4:
            continue
        _, _, vt = np.linalg.svd(nb - nb.mean(0), full_matrices=False)
        tan = vt[0]
        nrm = np.array([-tan[1], tan[0]])
        yy = p[0] + alongs[:, None] * tan[0] + nrm[0] * T[None] / G
        xx = p[1] + alongs[:, None] * tan[1] + nrm[1] * T[None] / G
        prof = ndi.map_coordinates(q, [yy, xx], order=1, mode="nearest").mean(0)
        yi = np.clip(np.round(p[0] + nrm[0] * T / G).astype(int), 0, H - 1)
        xi = np.clip(np.round(p[1] + nrm[1] * T / G).astype(int), 0, W - 1)
        if bg[yi, xi][np.abs(T) <= w].any():
            continue                                         # the fissure is open here: the cut reaches the gap
        info["deep"][k] = True
        E, B = eglish[yi, xi], bandok[yi, xi]
        qs = ndi.gaussian_filter1d(prof, 2)
        best = None
        for i in range(1, len(T) - 1):
            t = T[i]
            if abs(t) > w or not B[i] or not (qs[i] <= qs[i - 1] and qs[i] <= qs[i + 1]):
                continue
            l, r = (T >= t - sd) & (T <= t - 4), (T >= t + 4) & (T <= t + sd)
            if E[l].mean() < prm["pia_egl_frac"] or E[r].mean() < prm["pia_egl_frac"]:
                continue
            if prm.get("pia_flank", "mean") == "peak":              # v7: the two inner EGLs either side of the pia
                depth = min(qs[l].max(), qs[r].max()) - qs[i]
            else:
                depth = min(qs[l & (T <= t - 8)].mean(), qs[r & (T >= t + 8)].mean()) - qs[i]
            if depth < prm["pia_depth"] or (prm.get("pia_abs", True) and prof[i] >= ml_level + prm["pia_abs_margin"]):
                continue
            if best is None or abs(t) < abs(best):
                best = t
        if best is not None:
            info["found"][k] = True
            shift[k] = nrm * best / G
    sm = np.zeros_like(shift)
    for k, p in enumerate(pts):
        nb = tree.query_ball_point(p, prm["pia_smooth_um"] / G)
        f = info["found"][nb]
        if f.any():
            sm[k] = shift[nb][f].mean(0) * min(1.0, f.mean() / 0.5)
    info["moved_um"] = np.hypot(*sm.T) * G
    ni = np.clip(np.round(pts + sm).astype(int), 0, [H - 1, W - 1])
    out = np.zeros_like(mid)
    out[ni[:, 0], ni[:, 1]] = True
    for a, b in tree.query_pairs(1.5):                       # keep the line connected after the move
        rr, cc = draw_line(ni[a, 0], ni[a, 1], ni[b, 0], ni[b, 1])
        out[rr, cc] = True
    new = F.single_line((out & tissue) | (mid & ndi.binary_dilation(bg, iterations=2)), bg, F.final_prm(prm)) & tissue
    return new, info


def extend_ends(ms, zone_ok, dist_igl, tissue, prm=PRM):
    """continue each stain-based cut from its deep end along the fissure's own direction (the last
    back_um of the cut): steps of one cell, turning at most turn_deg per step and max_turn_deg overall,
    choosing among straight / left / right the step furthest from the IGL (the middle between the two
    folia); stops within seed_min_um of the IGL, outside the outer zone, after ext_max_um, or next to
    another cut (it never joins two fissures)."""
    sk = morphology.skeletonize(ms)
    dbg = ndi.distance_transform_edt(tissue) * G
    deg = ndi.convolve(sk.astype(int), np.ones((3, 3), int), mode="constant") - 1
    ends = np.argwhere(sk & (deg == 1) & (dbg > prm["ext_end_bg_um"]))
    comp, _ = ndi.label(sk, structure=np.ones((3, 3), bool))
    ext = np.zeros_like(ms)
    H, W = ms.shape
    back = int(prm["ext_back_um"] / G)
    rot = lambda v, a: np.array([v[0] * np.cos(a) - v[1] * np.sin(a), v[0] * np.sin(a) + v[1] * np.cos(a)])
    other_cut = ndi.distance_transform_edt(~sk)
    for e in ends:
        # walk back along the skeleton to get the direction
        seen, front, far = {tuple(e)}, [tuple(e)], tuple(e)
        for _ in range(back):
            nxt = []
            for (r, c) in front:
                for dr in (-1, 0, 1):
                    for dc in (-1, 0, 1):
                        q = (r + dr, c + dc)
                        if 0 <= q[0] < H and 0 <= q[1] < W and sk[q] and q not in seen:
                            seen.add(q)
                            nxt.append(q)
            if not nxt:
                break
            front, far = nxt, nxt[0]
        d0 = np.array(e, float) - np.array(far, float)
        if np.hypot(*d0) < back / 3:
            continue
        d0 /= np.hypot(*d0)
        pos, d = np.array(e, float), d0.copy()
        own = comp[tuple(e)]
        for k in range(int(prm["ext_max_um"] / G)):
            best = None
            for a in (0.0, np.radians(prm["turn_deg"]), -np.radians(prm["turn_deg"])):
                dd = rot(d, a)
                if np.degrees(np.arccos(np.clip((dd * d0).sum(), -1, 1))) > prm["max_turn_deg"]:
                    continue
                q = np.round(pos + dd).astype(int)
                if not (0 <= q[0] < H and 0 <= q[1] < W) or not zone_ok[tuple(q)] or dist_igl[tuple(q)] * G < prm["seed_min_um"]:
                    continue
                if k > 3 and other_cut[tuple(q)] <= 2 and comp[tuple(np.clip(q, 0, [H - 1, W - 1]))] != own:
                    continue
                score = dist_igl[tuple(q)] + (0.5 if a == 0 else 0.0)
                if best is None or score > best[0]:
                    best = (score, dd, q)
            if best is None:
                break
            _, d, q = best
            pos = pos + d
            ext[tuple(q)] = True
    return ext


def furrow_tips(tissue, mid, prm=PRM):
    """deepest points of open furrows (narrow channels of background between two folia) that no cut leaves
    from: end points of the background's medial axis within tip_narrow_um of the tissue, on a stretch of axis at
    least 0.6 x tip_back_um long that stays that narrow, more than tip_cut_near_um from any cut. Returns
    [(point, direction into the tissue)]."""
    bg = ~tissue
    skel, dist = morphology.medial_axis(bg, return_distance=True)
    dist_um = dist * G
    sk = skel & (dist_um <= prm["tip_narrow_um"])
    deg = ndi.convolve(sk.astype(int), np.ones((3, 3), int), mode="constant") - 1
    dcut = ndi.distance_transform_edt(~mid) * G if mid.any() else np.full(tissue.shape, np.inf)
    back = int(prm["tip_back_um"] / G)
    H, W = tissue.shape
    out = []
    for e in np.argwhere(sk & (deg == 1)):
        if dcut[tuple(e)] <= prm["tip_cut_near_um"]:
            continue
        seen, cur, path = {tuple(e)}, tuple(e), [tuple(e)]
        for _ in range(back):
            nxt = [(cur[0] + dr, cur[1] + dc) for dr in (-1, 0, 1) for dc in (-1, 0, 1)]
            nxt = [q for q in nxt if 0 <= q[0] < H and 0 <= q[1] < W and sk[q] and q not in seen]
            if not nxt:
                break
            cur = nxt[0]
            seen.add(cur)
            path.append(cur)
        if len(path) < 0.6 * back:
            continue                                         # a notch in the surface, not a furrow
        # the furrow's axis: principal direction of the narrow channel within 100 um of the tip (the medial
        # axis forks at a tip, so its last steps can point sideways), oriented from the channel into the tissue
        r0 = int(100.0 / G)
        y0, y1, x0, x1 = max(e[0] - r0, 0), min(e[0] + r0 + 1, H), max(e[1] - r0, 0), min(e[1] + r0 + 1, W)
        yy, xx = np.nonzero(bg[y0:y1, x0:x1] & (dist_um[y0:y1, x0:x1] <= prm["tip_narrow_um"]))
        pts = np.c_[yy + y0, xx + x0].astype(float)
        pts = pts[np.hypot(*(pts - e).T) <= r0]
        if len(pts) < 10:
            continue
        c = pts.mean(0)
        _, _, vt = np.linalg.svd(pts - c, full_matrices=False)
        d0 = vt[0] if ((e - c) * vt[0]).sum() >= 0 else -vt[0]
        out.append((np.array(e), d0))
    return out


def walk_tip(e, d0, ok_mask, dist_igl, egl, other_cut, prm=PRM):
    """from a furrow tip, continue the furrow into the tissue: steps of one cell, turning at most turn_deg per
    step and max_turn_deg overall, preferring stain EGL and the middle between the folia (furthest from the IGL);
    any tissue that is not IGL may be crossed, but the walk ends tip_gap_um after the last stain-EGL cell and is
    cut back to it; kept if at least tip_egl_frac of it is stain EGL. Never within seed_min_um of the IGL or
    next to another cut."""
    H, W = ok_mask.shape
    rot = lambda v, a: np.array([v[0] * np.cos(a) - v[1] * np.sin(a), v[0] * np.sin(a) + v[1] * np.cos(a)])
    pos, d = np.array(e, float), d0.copy()
    path, last, n_egl, entered = [], 0, 0, False
    gap = int(prm["tip_gap_um"] / G)
    for k in range(int(prm["ext_max_um"] / G)):
        best = None
        for a in (0.0, np.radians(prm["turn_deg"]), -np.radians(prm["turn_deg"])):
            dd = rot(d, a)
            if np.degrees(np.arccos(np.clip((dd * d0).sum(), -1, 1))) > prm["max_turn_deg"]:
                continue
            q = tuple(np.round(pos + dd).astype(int))
            if not (0 <= q[0] < H and 0 <= q[1] < W):
                continue
            if not ok_mask[q] and (entered or k >= prm["tip_bg_steps"]):
                continue                                     # crossing the channel to the tissue first
            if dist_igl[q] * G < prm["seed_min_um"] or (k > 3 and other_cut[q] <= 2):
                continue
            score = dist_igl[q] + (0.5 if a == 0 else 0.0) + (3.0 if egl[q] else 0.0)
            if best is None or score > best[0]:
                best = (score, dd, q)
        if best is None:
            break
        _, d, q = best
        pos = pos + d
        if not ok_mask[q] and not entered:
            continue
        entered = True
        path.append(q)
        if egl[q]:
            last, n_egl = len(path), n_egl + 1
        if len(path) - last > gap:
            break
    path = path[:last]
    if not path or n_egl < prm["tip_egl_frac"] * len(path):
        return []
    return path


def trim_ends(mid, cortex, tissue, prm=PRM):
    """v8: end each cut at the furrow base. With the ends no longer eaten by single_line (flatten keep_ends), an
    extension that ran on past the fused EGL into the ML stays; so each deep end is cut back, one skeleton cell at
    a time, while it lies closer to the cortex (ML / IGL / white matter the stain does not call EGL) than end_frac
    x the half-width of the fused EGL behind it (the median of that distance along the cut within end_win_um, over the
    cut cells inside the EGL, so a long run into the ML does not pull it to zero), and
    never more than end_max_um: the run into the ML and the last few um before it go. At the base the two EGLs wrap
    round, so the base is still about a half-width from the ML."""
    from scipy.spatial import cKDTree
    bg = ~tissue
    sk = morphology.skeletonize(mid)
    if not sk.any():
        return mid
    dcx = ndi.distance_transform_edt(~cortex) * G
    pts = np.argwhere(sk)
    tree = cKDTree(pts)
    hw = np.zeros(sk.shape)
    vals = dcx[pts[:, 0], pts[:, 1]]
    for k, nb in enumerate(tree.query_ball_point(pts, prm["end_win_um"] / G)):
        v = vals[nb][vals[nb] > G]
        hw[tuple(pts[k])] = min(prm["end_frac"] * np.median(v), prm["end_max_um"]) if len(v) else prm["end_max_um"]
    keep = ndi.binary_dilation(bg, iterations=3)
    for _ in range(int(prm["ext_max_um"] / G)):
        bad = F.skel_ends(sk) & ~keep & (dcx < hw)
        if not bad.any():
            break
        sk &= ~bad
    return (ndi.binary_dilation(sk, structure=np.ones((2, 2), bool)) | (mid & keep)) & tissue


def walk_base(mid, eglish, cortex, tissue, prm=PRM):
    """v8: continue each deep end of a cut down the fused EGL to the furrow base (the user's flatten v5 notes: the
    cuts stopped short of the base, so the columns at the end of a cut ran down the rest of the fused EGL, the sharp
    EGL peaks in the flattened strips). From every end more than base_bg_um from the background, in the direction
    of the last base_back_um of the cut: steps of one cell, turning at most base_turn_deg per step and
    base_max_turn_deg overall, through EGL-like cells only (EGL label or stain, unassigned, anomaly, white matter the
    stain calls EGL), choosing the step furthest from the cortex (the middle of the fused EGL, where the pia is); it
    stops where that distance falls below base_frac x the fused EGL's half-width behind the end (the median along the
    cut within end_win_um): there the two EGLs wrap round the base. It never comes within 2 cells of another cut.
    Each step must also look like two EGLs back to back (the double-EGL test): across the walk, cortex on both sides
    within base_side x the section's EGL half-width t (the median distance to cortex along all its cuts), with no
    background before it, and the step itself no more than base_hi x t from the cortex; a thick EGL mass (cKO), an
    open surface or the image edge stops it."""
    from scipy.spatial import cKDTree
    sk = morphology.skeletonize(mid)
    out = np.zeros_like(mid)
    if not sk.any():
        return out
    H, W = mid.shape
    dbg = ndi.distance_transform_edt(tissue) * G
    dcx = ndi.distance_transform_edt(~cortex) * G
    comp, _ = ndi.label(sk, structure=np.ones((3, 3), bool))
    dsk, (ir, ic) = ndi.distance_transform_edt(~sk, return_indices=True)
    pts = np.argwhere(sk)
    tree = cKDTree(pts)
    vals = dcx[pts[:, 0], pts[:, 1]]
    back = int(prm["base_back_um"] / G)
    sv = vals[vals > G]
    if not len(sv):
        return out
    t_sec = float(np.median(sv))
    side = np.arange(1, int(prm["base_side"] * t_sec / G) + 2)
    bg = ~tissue

    def two_sided(q, d):
        n = np.array([-d[1], d[0]])
        for sg in (1, -1):
            pp = np.round(np.array(q, float)[None, :] + sg * n[None, :] * side[:, None]).astype(int)
            inb = (pp[:, 0] >= 0) & (pp[:, 0] < H) & (pp[:, 1] >= 0) & (pp[:, 1] < W)
            if not inb.all():
                pp = pp[inb]
                if not len(pp):
                    return False
            hit_c, hit_b = cortex[pp[:, 0], pp[:, 1]], bg[pp[:, 0], pp[:, 1]]
            if not hit_c.any() or (hit_b.any() and np.argmax(hit_b) < np.argmax(hit_c)):
                return False
        return True
    rot = lambda v, a: np.array([v[0] * np.cos(a) - v[1] * np.sin(a), v[0] * np.sin(a) + v[1] * np.cos(a)])
    turns = [0.0] + [sg * np.radians(t) for t in (prm["base_turn_deg"] / 2, prm["base_turn_deg"]) for sg in (1, -1)]
    for e in np.argwhere(F.skel_ends(sk) & (dbg > prm["base_bg_um"])):
        seen, front, far = {tuple(e)}, [tuple(e)], tuple(e)
        for _ in range(back):
            nxt = [(r + dr, c + dc) for (r, c) in front for dr in (-1, 0, 1) for dc in (-1, 0, 1)
                   if 0 <= r + dr < H and 0 <= c + dc < W and sk[r + dr, c + dc] and (r + dr, c + dc) not in seen]
            if not nxt:
                break
            seen.update(nxt)
            front, far = nxt, nxt[0]
        d0 = np.array(e, float) - np.array(far, float)
        if np.hypot(*d0) < back / 3:
            continue
        d0 /= np.hypot(*d0)
        v = vals[tree.query_ball_point(e, prm["end_win_um"] / G)]
        v = v[v > G]
        if not len(v):
            continue
        stop = prm["base_frac"] * float(np.median(v))
        own = comp[tuple(e)]
        pos, d, path = np.array(e, float), d0.copy(), []
        for k in range(int(prm["base_max_um"] / G)):
            best = None
            for a in turns:
                dd = rot(d, a)
                if np.degrees(np.arccos(np.clip((dd * d0).sum(), -1, 1))) > prm["base_max_turn_deg"]:
                    continue
                q = tuple(np.round(pos + dd).astype(int))
                if not (0 <= q[0] < H and 0 <= q[1] < W) or not (tissue[q] and eglish[q]):
                    continue
                if dsk[q] <= 2 and comp[ir[q], ic[q]] != own:
                    continue
                score = dcx[q] + (0.5 if a == 0 else 0.0)
                if best is None or score > best[0]:
                    best = (score, dd, q)
            if best is None or dcx[best[2]] < stop or dcx[best[2]] > prm["base_hi"] * t_sec or not two_sided(best[2], best[1]):
                break
            _, d, q = best
            pos = pos + d
            path.append(q)
        for q in path:
            out[q] = True
    return ndi.binary_dilation(out, structure=np.ones((2, 2), bool)) & tissue


def cuts(lab, tissue, S, prm=PRM, protect=None):
    """fissure cuts: down the middle of the fused EGL the stains show (flatten.midline on labels whose
    EGL / ML in the outer zone come from the stains), continued from each deep end along the fissure's
    own direction towards the furrow base (extend_ends), plus fissures only the IGL shows (cut halfway
    between the facing IGL bands, more than ext_near_um from any stain-based cut)."""
    bg = ~tissue
    igl_ok = (lab == 4) & ~(S["egl"] & ~S["igl"])                # IGL the stain does not call EGL
    if protect is not None:
        igl_ok |= protect                                        # v7: IGL the user traced (bleached NeuN looks like EGL)
    igl_ok = morphology.remove_small_objects(igl_ok, max_size=int(prm["igl_min_um2"] / G ** 2), connectivity=2)
    zone = tissue & ~igl_ok & (lab != 5) & (lab != 6)
    lab_s = lab.copy()
    lab_s[zone & np.isin(lab, (2, 3, 4))] = 3
    lab_s[zone & S["egl"]] = 2
    ms = F.midline(lab_s, tissue, prm)
    layer = ndi.binary_dilation(np.isin(lab_s, (2, 3)), iterations=int(round(prm["seed_layer_um"] / G)))
    mi, dist_igl = igl_midline(igl_ok, zone & layer, tissue, prm)
    # v5: the extension may run through white-matter labels where the stain shows EGL (the fused EGL at a
    # furrow base labelled white matter, s4_A at 5); it still stops seed_min_um from the IGL
    wm_egl = (lab == 5) & ndi.binary_dilation(S["egl"], iterations=1)   # white-matter labels the stain calls EGL
    ext_ok, dist_ext = zone | np.isin(lab, (1, 6)) | wm_egl, dist_igl
    if protect is not None:
        # v7: a fissure cut may continue through traced IGL the stain calls EGL (the fused EGL a bleach bridge was
        # traced across, s4_C below the fissure at 7), but no new cut starts in it (zone above excludes it)
        through = protect & S["egl"]
        ext_ok = ext_ok | through
        dist_ext = ndi.distance_transform_edt(~(igl_ok & ~protect))      # traced IGL does not steer it
    ext = extend_ends(ms, ext_ok, dist_ext, tissue, prm) if ms.any() else np.zeros_like(ms)
    fb = np.zeros_like(mi)
    if mi.any():
        far = ndi.distance_transform_edt(~ms) * G > prm["ext_near_um"] if ms.any() else np.ones_like(mi)
        comp, nc = ndi.label(mi, structure=np.ones((3, 3), bool))
        idx = np.arange(1, nc + 1)
        frac = np.asarray(ndi.mean(far, comp, index=idx))
        # a closed fissure has EGL in it and no white matter beside it; the white-matter core between the two
        # IGL bands of one folium does not (s1_5 at 12-13)
        near = ndi.grey_dilation(comp, size=(7, 7)) * tissue
        egl_f = np.asarray(ndi.mean(S["egl"], near, index=idx))
        wm_f = np.asarray(ndi.mean(lab == 5, near, index=idx))
        ok = (frac > 0.8) & (egl_f >= prm["fallback_egl_frac"]) & (wm_f < prm["fallback_wm_frac"])
        fb = np.r_[False, ok][comp]                              # whole fissures the stain did not find
    if prm.get("tip_ext"):                                   # v7: open furrows continue as a cut (s5_F)
        oc = ndi.distance_transform_edt(~(ms | ext | fb)) if (ms | ext | fb).any() else np.full(tissue.shape, np.inf)
        ok_mask = tissue & ~igl_ok
        tips = np.zeros_like(ms)
        for e, d0 in furrow_tips(tissue, ms | ext | fb, prm):
            for q in walk_tip(e, d0, ok_mask, dist_igl, S["egl"], oc, prm):
                tips[q] = True
        ext = ext | tips
    mid = F.single_line(ms | ext | fb, bg, prm) & tissue
    if prm.get("end_trim"):                                  # v8: the deep ends run on to the furrow base
        cortex = tissue & np.isin(lab_s, (3, 4, 5)) & ~S["egl"]
        mid = trim_ends(mid, cortex, tissue, prm)              # an extension that ran into the ML is cut back,
        eglish = (lab_s == 2) | S["egl"] | np.isin(lab, (1, 6)) | wm_egl
        mid = mid | walk_base(mid, eglish, cortex, tissue, prm)  # then each deep end goes on down the fused EGL
    pre, pia = _prune_shallow(mid, tissue, prm), None
    if prm.get("pia_snap"):                                  # v6: onto the pia between the two EGLs
        mid, pia = pia_snap(mid, lab, tissue, S, prm)
    return _prune_shallow(mid, tissue, prm), ms, fb, ext, pre, pia


def relabel(lab, tissue, cut, S, prm=PRM, ml_allow=None, protect=None):
    """EGL / IGL swaps the stains support, near the surface (depth from the reference line):
    IGL -> EGL where the stain shows EGL (bright nuclei, low NeuN) within igl2egl_um of the surface, in
    pieces touching EGL or the surface; EGL -> IGL where the stain shows IGL (high NeuN) deeper than
    egl2igl_um; IGL -> ML where the stain shows ML (dim nuclei, low NeuN) within igl2ml_um, only inside ml_allow (repair() passes
    the IGL the earlier repairs added to the controls). EGL / ML need NeuN below the
    midpoint of the section's EGL and IGL levels, IGL above it (stain_maps)."""
    dref = ndi.distance_transform_edt(cut) * G
    piece = lambda m, um2: morphology.remove_small_objects(m, max_size=int(um2 / G ** 2), connectivity=2)
    i2e = piece((lab == 4) & S.get("egl_loose", S["egl"]) & (dref <= prm["igl2egl_um"]) & (True if protect is None else ~protect),
                prm["swap_min_um2"])
    if i2e.any():
        comp, nc = ndi.label(i2e, structure=np.ones((3, 3), bool))
        touch = ndi.binary_dilation((lab == 2) | (dref <= 2 * G))
        keep = np.asarray(ndi.maximum(touch, comp, index=np.arange(1, nc + 1))) > 0
        i2e = np.r_[False, keep][comp]
    e2i = piece((lab == 2) & S["igl"] & (dref > prm["egl2igl_um"]), prm["egl2igl_min_um2"])
    allow = np.ones_like(i2e) if ml_allow is None else ml_allow
    i2m = piece((lab == 4) & S["ml"] & ~S["igl"] & (dref <= prm["igl2ml_um"]) & ~i2e & allow, prm["swap_min_um2"]) if prm["igl2ml_um"] > 0 else np.zeros_like(i2e)
    new = lab.copy()
    new[i2e] = 2
    new[e2i] = 4
    new[i2m] = 3
    return new, dict(igl_to_egl=i2e, egl_to_igl=e2i, igl_to_ml=i2m)


def _profiles(S, valid, piece, avg):
    """log nuclei profile per column, averaged over +-avg columns within a piece."""
    P = np.log(np.maximum(S.astype(np.float32), 1.0))
    P = np.where(valid, P, np.nan)
    out = np.full_like(P, np.nan)
    for p in np.unique(piece):
        m = piece == p
        block = P[m]
        num = ndi.uniform_filter1d(np.nan_to_num(block), 2 * avg + 1, axis=0, mode="nearest")
        den = ndi.uniform_filter1d(np.isfinite(block).astype(np.float32), 2 * avg + 1, axis=0, mode="nearest")
        out[m] = np.where(den > 0.5, num / np.maximum(den, 1e-6), np.nan)
    return out


def fit_columns(P, valid, prior, prm=PRM, chunk=200):
    """three-segment fit (EGL / ML / IGL) of every column; returns z1, z2 (um), contrasts."""
    N = P.shape[0]
    D = min(P.shape[1], int(prm["fit_depth_um"] / STEP))
    z1o = np.arange(int(prm["z1_min_um"] / STEP), min(D - 8, int(prm["z1_max_um"] / STEP)))
    z2o = np.arange(z1o[0] + int(prm["ml_min_um"] / STEP), D - 4)
    seg = int(prm["igl_seg_um"] / STEP)
    res = {k: np.full(N, np.nan) for k in ("z1", "z2", "c_egl", "c_igl")}
    for a in range(0, N, chunk):
        b = min(N, a + chunk)
        X = P[a:b, :D]
        L = np.minimum(np.isfinite(X).cumprod(1).sum(1), D)          # usable length per column
        X = np.nan_to_num(X)
        S1 = np.concatenate([np.zeros((b - a, 1)), np.cumsum(X, 1)], 1)
        S2 = np.concatenate([np.zeros((b - a, 1)), np.cumsum(X * X, 1)], 1)

        # segment 1: [0, z1); segment 2: [z1, z2); segment 3: [z2, min(L, z2 + seg))
        Z1 = z1o[None, :, None]
        Z2 = z2o[None, None, :]
        e3 = np.minimum(L[:, None, None], Z2 + seg)
        n1, n2, n3 = Z1, Z2 - Z1, np.maximum(e3 - Z2, 1)
        I1 = np.broadcast_to(Z1, (b - a, len(z1o), len(z2o)))
        I2 = np.broadcast_to(Z2, (b - a, len(z1o), len(z2o)))
        I3 = np.broadcast_to(e3, (b - a, len(z1o), len(z2o)))
        take = lambda S, I: np.take_along_axis(S, I.reshape(b - a, -1), 1).reshape(I.shape)
        s1a, s1b, s1c = take(S1, I1), take(S1, I2), take(S1, I3)
        s2a, s2b, s2c = take(S2, I1), take(S2, I2), take(S2, I3)
        m1, m2, m3 = s1a / n1, (s1b - s1a) / n2, (s1c - s1b) / n3
        cost = (s2a - s1a * m1) + ((s2b - s2a) - (s1b - s1a) * m2) + ((s2c - s2b) - (s1c - s1b) * m3)
        ok = (Z2 + 4 <= L[:, None, None]) & (m1 > m2) & (m3 > m2) & (Z2 > Z1 + int(prm["ml_min_um"] / STEP) - 1)
        pr = prior[a:b]
        lo = np.where(np.isfinite(pr), (pr - prm["z2_lo_um"]) / STEP, z2o[0])
        hi = np.where(np.isfinite(pr), (pr + prm["z2_hi_um"]) / STEP, prm["z2_free_max_um"] / STEP)
        ok &= (Z2 >= lo[:, None, None]) & (Z2 <= hi[:, None, None])
        cost = np.where(ok, cost, np.inf)
        flat = cost.reshape(b - a, -1)
        best = np.argmin(flat, 1)
        good = np.isfinite(flat[np.arange(b - a), best])
        i1, i2 = np.unravel_index(best, cost.shape[1:])
        r = np.arange(b - a)
        res["z1"][a:b] = np.where(good, z1o[i1] * STEP, np.nan)
        res["z2"][a:b] = np.where(good, z2o[i2] * STEP, np.nan)
        res["c_egl"][a:b] = np.where(good, (m1 - m2)[r, i1, i2], np.nan)
        res["c_igl"][a:b] = np.where(good, (m3 - m2)[r, i1, i2], np.nan)
    return res


def smooth_runs(v, keep, win):
    out = np.full_like(v, np.nan)
    for a, b in F._runs(keep):
        out[a:b] = np.asarray(pd.Series(v[a:b]).rolling(win, center=True, min_periods=1).median())
    return out


def _sample(lab, X, valid):
    N, Dn = valid.shape
    Lb = ndi.map_coordinates(lab, X.reshape(-1, 2).T, order=0, mode="constant").reshape(N, Dn).astype(np.uint8)
    Lb[~valid] = 0
    return Lb


def _paint(mask_cols, X, shape):
    """grid cells the flagged column samples fall in: votes for (flagged) and total samples per cell."""
    gh, gw = shape
    r = np.clip(np.round(X[..., 0]).astype(int), 0, gh - 1)
    c = np.clip(np.round(X[..., 1]).astype(int), 0, gw - 1)
    yes = np.zeros(shape, np.int32)
    np.add.at(yes, (r[mask_cols], c[mask_cols]), 1)
    return yes


def igl_above_ml(Lb, valid, prm=PRM):
    """column samples of IGL that sit above an ML band with IGL below it (inner EGL called IGL): the
    first IGL run starts within order_start_um of the line and is at most order_run_um long, the next
    stretch is mostly ML (>= order_ml_um, no white matter), then IGL again (>= order_deep_um)."""
    N, Dn = Lb.shape
    out = np.zeros_like(valid)
    s0, rmax = int(prm["order_start_um"] / STEP), int(prm["order_run_um"] / STEP)
    mlmin, deep = int(prm["order_ml_um"] / STEP), int(prm["order_deep_um"] / STEP)
    I = ndi.binary_closing((Lb == 4) & valid, np.ones((1, 3), bool)) & valid
    cand = np.where(I[:, :s0 + 1].any(1))[0]
    for i in cand:
        runs = F._runs(I[i])
        a = runs[0][0]
        for j in range(len(runs) - 1):
            b, c, d = runs[j][1], runs[j + 1][0], runs[j + 1][1]
            if b - a > rmax:
                break
            mid = Lb[i, b:c]
            if (mid == 5).any():                                 # white matter: past the cortex
                break
            if (mid == 3).sum() >= mlmin and np.isin(mid, (3, 6)).mean() >= prm["order_ml_frac"] and d - c >= deep:
                out[i, a:b] = (Lb[i, a:b] == 4)
                break
    return out


def igl_in_ml(Lb, valid, fit, present, prm=PRM):
    """column samples labelled IGL between the EGL and the IGL top the column fit found, where the fit is
    ml_gap_um or more below the labelled IGL top (the labelled IGL covers the ML)."""
    N, Dn = Lb.shape
    d_um = (np.arange(Dn) + 0.5) * STEP
    I = (Lb == 4) & valid
    top = np.where(I.any(1), I.argmax(1) * STEP, np.nan)
    z1, z2 = fit["z1"], fit["z2"]
    ok = present & np.isfinite(z2) & np.isfinite(top) & (z2 - top >= prm["ml_gap_um"])
    lo = np.where(ok, np.maximum(z1, 0), np.inf)[:, None]
    hi = np.where(ok, z2 - prm["ml_keep_um"], -np.inf)[:, None]
    return I & (d_um[None, :] >= lo) & (d_um[None, :] < hi)


def repair(short, prm=PRM, save=True):
    t0 = time.time()
    sid = F.sid_of(short)
    auto, best, tissue = F.labels(short)
    S_ = stain_maps(short, best, tissue, prm)
    # v7: in the reviewed controls the IGL the user traced by hand (bleached NeuN, s4_C D4/E4) is IGL: no cut runs
    # through it and it is not turned into EGL (it looks like EGL to the stains: bright nuclei, low NeuN)
    traced = ((best == 4) & (auto != 4)) if (prm.get("protect_traced") and short in F.D.controls() and short != GOLD) else None
    # v9: the cuts and the reference line run on the tissue inside the pia (arachnoid trimmed)
    arach = pia_trim(best, tissue, S_, prm, gold=short == GOLD) if prm.get("pia_trim") else np.zeros_like(tissue)
    tis_p = tissue & ~arach
    mid, ms, mi, ext, mid_pre, pia = cuts(best, F.cut_tissue(tis_p, prm), S_, prm, protect=traced)
    me = F.midline(best, tissue, prm)                       # the old cuts, from the EGL labels, for comparison
    refs, cut, _ = F.reference(best, tis_p, prm, mid=mid)
    # IGL -> ML only on IGL the earlier repairs added to the controls (traces across bleached stretches);
    # on MATLAB's own IGL and in the cKO sections it removed good IGL in testing (eval_swaps)
    repaired = (best == 4) & (auto != 4) if short in F.D.controls() else np.zeros_like(tissue)
    if traced is not None:                                   # traced IGL right beside a cut extension is the fused EGL
        traced = traced & ~ndi.binary_dilation(ext & mid, iterations=3)
    lab1, swaps = relabel(best, tissue, cut, S_, prm, ml_allow=repaired, protect=traced)
    if short == GOLD:                                       # the gold standard's EGL / IGL are not swapped
        lab1, swaps = best.copy(), {k_: np.zeros_like(v) for k_, v in swaps.items()}
    Q = np.concatenate([r["Q"] for r in refs])
    n = np.concatenate([r["n"] for r in refs])
    piece = np.concatenate([np.full(len(r["Q"]), i) for i, r in enumerate(refs)])
    X, valid = F.streamlines(Q, n, cut, prm)
    N, Dn = valid.shape
    flat = X.reshape(-1, 2).T
    Lb = _sample(lab1, X, valid)
    piece_ok = lambda m: morphology.remove_small_objects(m, max_size=int(prm["swap_min_um2"] / G ** 2), connectivity=2)
    vote = lambda fl, allm: (lambda y, t: (y > 0) & (2 * y >= t))(_paint(fl, X, tissue.shape), _paint(allm, X, tissue.shape))
    # the layers come in the order EGL, ML, IGL, white matter along every column: IGL above an ML band
    # that has IGL below it is the inner EGL (not on the gold standard)
    order = np.zeros_like(tissue)
    if short != GOLD:
        order = piece_ok((lab1 == 4) & vote(igl_above_ml(Lb, valid, prm), (Lb == 4) & valid))
        lab1[order] = 2
        # what is left of such a band: small IGL islands lying wholly within order_start_um of the line take
        # the label of the nearest tissue that is not IGL
        dref0 = ndi.distance_transform_edt(cut) * G
        comp, nc = ndi.label(lab1 == 4, structure=np.ones((3, 3), bool))
        if nc:
            size = np.bincount(comp.ravel())[1:] * G * G
            maxd = np.asarray(ndi.maximum(dref0, comp, index=np.arange(1, nc + 1)))
            isl = np.r_[False, (size < prm["igl_island_um2"]) & (maxd <= prm["order_start_um"])][comp]
            if isl.any():
                other = tissue & ~np.isin(lab1, (0, 4))
                _, (ir, ic) = ndi.distance_transform_edt(~other, return_indices=True)
                lab1[isl] = lab1[ir[isl], ic[isl]]
                order |= isl
        Lb = _sample(lab1, X, valid)
    swaps["igl_above_ml"] = order
    R, k, (Hn, Wn) = F.stains_reduced(sid)
    gh, gw = tissue.shape
    red = (np.stack([flat[0] * (Hn - 1) / (gh - 1), flat[1] * (Wn - 1) / (gw - 1)]) + 0.5) / k - 0.5
    Sn = ndi.map_coordinates(R["nuclei"], red, order=1, mode="nearest").reshape(N, Dn)
    eb, top, bot = F.bounds(Lb, valid, prm)
    ftop, _ = F.fit_line(top, valid[:, 3], prm)
    prior = np.where(np.isfinite(top), top, ftop)
    fit = fit_columns(_profiles(Sn, valid, piece, int(prm["prof_avg_um"] / STEP)), valid, prior, prm)
    win = int(prm["z1_med_um"] / STEP) | 1
    sm = lambda v: np.concatenate([smooth_runs(np.nan_to_num(v[piece == p_], nan=-9.0), np.ones((piece == p_).sum(), bool), win)
                                   for p_ in np.unique(piece)])
    cs, ds, z1s = sm(fit["c_egl"]), sm(fit["c_egl"] - fit["c_igl"]), sm(fit["z1"])
    present = (cs >= prm["contrast"]) & (ds >= prm["egl_vs_igl"]) & (z1s >= prm["z1_present_um"])
    # labelled IGL over the ML: where the column fit puts the IGL top ml_gap_um or more below the labelled
    # IGL top, the labelled IGL between the EGL and the fitted top becomes ML (the stain must not show IGL)
    over_ml = np.zeros_like(tissue)
    if short != GOLD:
        fz = dict(fit, z1=z1s, z2=sm(fit["z2"]))
        over_ml = piece_ok((lab1 == 4) & ~S_["igl"] & vote(igl_in_ml(Lb, valid, fz, present, prm), (Lb == 4) & valid))
        lab1[over_ml] = 3
        Lb = _sample(lab1, X, valid)
    swaps["igl_over_ml"] = over_ml
    from igl_reference_page import egl_raw  # noqa: E402
    lab_has, lab_thick = egl_raw(Lb, valid)
    lab_ok = lab_has & (lab_thick >= prm["thin_um"])
    cortex = ((Lb == 4) & valid).any(1) | np.isfinite(prior)      # IGL under this point (or the line through it)
    d_um = (np.arange(Dn) + 0.5) * STEP
    Ib = (Lb == 4) & valid
    first_igl = np.where(Ib.any(1), Ib.argmax(1) * STEP, np.inf)
    # labelled EGL further down the column, above its IGL: the EGL is there, under something else (arachnoid,
    # a flap of tissue; s1_2 at 11 and 15), so this is not a gap (v5)
    egl_under = (((Lb == 2) & valid & (d_um[None, :] < first_igl[:, None])).sum(1) * STEP >= prm["thin_um"]) & ~lab_has
    qi = np.clip(np.round(Q).astype(int), 0, [gh - 1, gw - 1])
    dmid = ndi.distance_transform_edt(~mid)[qi[:, 0], qi[:, 1]]
    on_cut = dmid <= 1.5
    egl_under &= prm["use_egl_under"] & ~on_cut                    # arachnoid only covers the outer surface
    gap = ~lab_ok & ~egl_under                                    # EGL label missing or thin
    fill = present & gap & cortex
    # bridge the remaining gaps over IGL, up to bridge_egl_um, from the EGL thickness on either side
    cap = min(prm["egl_cap_max_um"], prm["egl_cap_factor"] * float(np.median(lab_thick[lab_ok]))) if lab_ok.any() else prm["egl_cap_max_um"]
    z1s = np.minimum(z1s, cap)
    th = np.where(fill, z1s, np.where(lab_ok, np.minimum(lab_thick, cap), np.nan))
    zb = np.full(N, np.nan)
    open_ = gap & cortex & ~fill
    for a, b in F._runs(open_):
        if a == 0 or b >= N or piece[a - 1] != piece[b] or (b - a) * STEP > prm["bridge_egl_um"]:
            continue
        if np.isfinite(th[a - 1]) and np.isfinite(th[b]):
            zb[a:b] = np.interp(np.arange(a, b), [a - 1, b], [th[a - 1], th[b]])
    bridge = np.isfinite(zb)
    depth_egl = np.where(fill, z1s, zb)
    # the added EGL never reaches past the first IGL in its column. Columns from a cut take unassigned and ML
    # cells; columns from a cut extension (a furrow base) also white-matter labels the stain calls EGL: the fused
    # EGL there is labelled white matter (s4_A at 5), and the extension stays 30 um from the IGL. Columns from the outer
    # surface also take white matter and anomaly above that IGL, but skip cells whose p27 and NeuN are both
    # well below the ML (arachnoid and pia, v5).
    on_ext = on_cut & (ndi.distance_transform_edt(~(ext & mid))[qi[:, 0], qi[:, 1]] <= 1.5) if (ext & mid).any() else np.zeros(N, bool)
    inside = valid & (fill | bridge)[:, None] & (d_um[None, :] < np.minimum(depth_egl, first_igl)[:, None])

    def _hit(cols):
        m = inside & cols[:, None]
        h = np.zeros(tissue.shape, bool)
        h[np.clip(np.round(X[..., 0][m]).astype(int), 0, gh - 1), np.clip(np.round(X[..., 1][m]).astype(int), 0, gw - 1)] = True
        return ndi.binary_closing(h, iterations=1) & cut
    room = ndi.distance_transform_edt(lab1 != 4) * G >= prm["ml_room_um"]       # leave the ML between EGL and IGL
    takes = np.isin(lab1, (1, 3))
    if prm.get("cut_takes_anom_egl"):                        # v7: anomaly the stain calls EGL lines the furrow (s5_E)
        takes |= (lab1 == 6) & S_["egl"]
    add = ((_hit(~on_cut) & np.isin(lab1, (1, 3, 5, 6)) & ~S_["nonegl"]) | (_hit(on_cut & ~on_ext) & takes)
           | (_hit(on_ext) & (takes | ((lab1 == 5) & ndi.binary_dilation(S_["egl"], iterations=1))))) & room
    # fragments the meninges exclusion leaves behind, unless they close a notch in the labelled EGL
    comp, nc = ndi.label(add, structure=np.ones((3, 3), bool))
    if nc:
        size = np.bincount(comp.ravel())[1:] * G * G
        touch = np.asarray(ndi.maximum(ndi.binary_dilation(lab1 == 2), comp, index=np.arange(1, nc + 1))) > 0
        add &= np.r_[False, (size >= prm["nonegl_min_um2"]) | touch][comp]
    new_lab = lab1.copy()
    new_lab[add] = 2
    # EGL that cannot be EGL: deeper than deep_egl_um below the line, or a small speck away from the surface;
    # it takes the label of the nearest tissue that is not EGL (white matter, IGL or ML). Not on the gold standard.
    deep = np.zeros_like(add)
    if short != GOLD:
        dref = ndi.distance_transform_edt(cut) * G
        E = new_lab == 2
        deep = E & (dref > prm["deep_egl_um"])
        comp, nc = ndi.label(E, structure=np.ones((3, 3), bool))
        if nc:
            size = np.bincount(comp.ravel()) * G * G
            mind = np.asarray(ndi.minimum(dref, comp, index=np.arange(1, nc + 1)))
            speck = np.r_[False, (size[1:] < prm["speck_egl_um2"]) & (mind > 24.0)]
            deep |= speck[comp]
        if deep.any():
            other = tissue & ~np.isin(new_lab, (0, 2)) & ~deep
            _, (ir, ic) = ndi.distance_transform_edt(~other, return_indices=True)
            new_lab[deep] = new_lab[ir[deep], ic[deep]]
            new_lab[deep & S_["nonegl"]] = 5             # no p27 or NeuN that deep: white matter (s1_5 at 12-13, v5)
    # v9: EGL labels outside the pia are arachnoid (not on the gold, whose labels the trim leaves alone)
    arach_egl = arach & (new_lab == 2)
    new_lab[arach_egl] = 1
    um = lambda m: int(m.sum() * STEP)
    a2 = lambda m: int(m.sum() * G * G)
    both = present & lab_ok
    diff = (z1s - lab_thick)[both]
    res = dict(section=short, group="control" if short in F.D.controls() else "cKO", thresholds={k_: round(v, 3) for k_, v in S_["t"].items()},
               cuts=dict(stain_um2=a2(ms), igl_um2=a2(mi), extension_um2=a2(ext), final_um2=a2(mid), old_egl_label_um2=a2(me),
                         old_near_new=round(float(((ndi.distance_transform_edt(~mid) * G)[me] <= prm["fallback_um"]).mean()), 3) if me.any() and mid.any() else None,
                         pia=None if pia is None else dict(
                             closed_um=int(pia["deep"].sum() * G), found_frac=round(float(pia["found"][pia["deep"]].mean()), 3) if pia["deep"].any() else None,
                             moved_over_4um_um=int((pia["moved_um"] > 4).sum() * G),
                             moved_med_um=round(float(np.median(pia["moved_um"][pia["moved_um"] > 4])), 1) if (pia["moved_um"] > 4).any() else 0.0)),
               swaps_um2={k_: a2(v) for k_, v in swaps.items()},
               line_um=int(N * STEP), pieces=len(refs),
               gap_over_igl_um=um(gap & cortex), gap_filled_um=um(fill), gap_bridged_um=um(bridge), gap_left_um=um(open_ & ~bridge), egl_under_um=um(egl_under),
               z1_vs_label=dict(median_diff_um=round(float(np.median(diff)), 1) if diff.size else None,
                                median_abs_diff_um=round(float(np.median(np.abs(diff))), 1) if diff.size else None,
                                within_8um=round(float(np.mean(np.abs(diff) <= 8)), 3) if diff.size else None),
               egl_added_um2=a2(add), egl_label_um2=a2(best == 2), egl_cap_um=round(cap, 1), deep_egl_relabelled_um2=a2(deep),
               arachnoid_um2=a2(arach), arachnoid_egl_um2=a2(arach_egl))
    if save:
        OUT.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(OUT / f"{short}.npz", arach=arach, mid=mid, mid_pre=mid_pre, mid_stain=ms, mid_igl=mi, mid_ext=ext, mid_egl=me, cut=cut,
                            egl_stain=S_["egl"], igl_stain=S_["igl"], ml_stain=S_["ml"], nonegl=S_["nonegl"], lab1=lab1, new_lab=new_lab, add=add,
                            **{f"swap_{k_}": v for k_, v in swaps.items()},
                            Q=Q.astype(np.float32), n=n.astype(np.float32), piece=piece, X=X.astype(np.float32), valid=valid, Lb=Lb,
                            S=np.round(Sn).astype(np.uint16), present=present, z1s=z1s.astype(np.float32), zb=zb.astype(np.float32),
                            fill=fill, bridge=bridge, gap=gap, cortex=cortex, deep_egl=deep, lab_has=lab_has, lab_ok=lab_ok,
                            lab_thick=lab_thick.astype(np.float32), prior=prior.astype(np.float32),
                            **{f"fit_{k2}": v.astype(np.float32) for k2, v in fit.items()})
    res["seconds"] = round(time.time() - t0, 1)
    print(json.dumps(res), flush=True)
    return res


def eval_swaps(prm=PRM):
    """the label swaps run on MATLAB's labels, scored against what the review (controls) and the user's
    round-1 edits (cKO) changed: IGL MATLAB had that is no longer IGL, and IGL added."""
    rows = []
    for short in F.D.training_sections():
        auto, best, tissue = F.labels(short)
        S_ = stain_maps(short, auto, tissue, prm)
        mid, *_ = cuts(auto, tissue, S_, prm)
        cut = (tissue & ~mid)
        new, sw = relabel(auto, tissue, cut, S_, prm)
        removed = (auto == 4) & (best != 4) & (best > 0)
        added = (auto != 4) & (best == 4)
        rem_rule = sw["igl_to_egl"] | sw["igl_to_ml"]
        typ = {k_: dict(um2=int(sw[k_].sum() * G * G), on_removed=round(float(removed[sw[k_]].mean()), 2) if sw[k_].any() else None)
               for k_ in ("igl_to_egl", "igl_to_ml")}
        r = dict(section=short, group="control" if short in F.D.controls() else "cKO",
                 truth_removed_um2=int(removed.sum() * G * G), rule_removed_um2=int(rem_rule.sum() * G * G),
                 recall=round(float((rem_rule & removed).sum() / max(removed.sum(), 1)), 2),
                 precision=round(float((rem_rule & removed).sum() / max(rem_rule.sum(), 1)), 2),
                 truth_added_um2=int(added.sum() * G * G), rule_egl_to_igl_um2=int(sw["egl_to_igl"].sum() * G * G),
                 egl_to_igl_on_truth_igl=round(float((sw["egl_to_igl"] & (best == 4)).sum() / max(sw["egl_to_igl"].sum(), 1)), 2), types=typ)
        rows.append(r)
        print(json.dumps(r), flush=True)
    (OUT / "swap_eval.json").write_text(json.dumps(rows, indent=1))
    return rows


def gold_pia_check(mid):
    """test of the fissure cuts on the gold s1_2: the user left the pia between two fused EGLs unassigned. Strips =
    unassigned cells with EGL on two opposite sides within 24 um (one of four directions), away from the
    background; returns (share of their centreline within 4 um of the cut, share of their cells under the cut)."""
    _, best, tissue = F.labels(GOLD)
    E, U = best == 2, (best == 1) & tissue

    def shift(a, dy, dx):
        out = np.zeros_like(a)
        H, W = a.shape
        out[max(-dy, 0):H + min(-dy, 0), max(-dx, 0):W + min(-dx, 0)] = a[max(dy, 0):H + min(dy, 0), max(dx, 0):W + min(dx, 0)]
        return out
    between = np.zeros_like(U)
    for dy, dx in ((0, 1), (1, 0), (1, 1), (1, -1)):
        plus, minus = np.zeros_like(U), np.zeros_like(U)
        for k in range(1, 7):
            plus |= shift(E, dy * k, dx * k)
            minus |= shift(E, -dy * k, -dx * k)
        between |= plus & minus
    gold = morphology.remove_small_objects(U & between & ~ndi.binary_dilation(~tissue, iterations=2), max_size=8, connectivity=2)
    sk = morphology.skeletonize(ndi.binary_closing(gold)) & gold
    d = ndi.distance_transform_edt(~mid) * G
    return round(float((d[sk] <= 4).mean()), 3), round(float(mid[gold].mean()), 3)


FINAL = OUT / "final"
CODES = {0: "background", 1: "tissue, unassigned", 2: "EGL", 3: "ML / PCL", 4: "IGL", 5: "white matter", 6: "anomaly"}


def export(shorts=None):
    """the repaired label grids (new_lab of each <short>.npz) -> final/<short>_labels.tif, the tissue mask
    -> final/<short>_tissue.tif, and manifest.json."""
    import tifffile
    FINAL.mkdir(parents=True, exist_ok=True)
    rows = {r["section"]: r for r in json.loads((OUT / "summary.json").read_text())["sections"]}
    man = dict(version="v10 (v9 + flatten v8 notes: meninges on ML / IGL where the cKO EGL is missing trimmed too; 2026-11-24)", grid_um=G, codes=CODES, script="scripts/egl_repair.py",
               note="same grid and codes as flatten.labels; the gold s1_2 gets no swaps or deep-EGL rules; <short>_pia.tif = the tissue "
                    "inside the pial surface (tissue less arachnoid / meninges), which the reference line follows", sections={})
    for short in shorts or F.all_shorts():
        lab = np.load(OUT / f"{short}.npz")["new_lab"].astype(np.uint8)
        _, best, tissue = F.labels(short)
        # outside the tissue mask nothing changes (the gold s1_2 has 72 labelled cells just outside it)
        assert lab.shape == best.shape and (lab[~tissue] == best[~tissue]).all() and (lab[tissue] > 0).all() and lab.max() <= 6
        tifffile.imwrite(FINAL / f"{short}_labels.tif", lab, compression="zlib")
        tifffile.imwrite(FINAL / f"{short}_tissue.tif", tissue.astype(np.uint8), compression="zlib")
        arach = np.load(OUT / f"{short}.npz")["arach"]
        tifffile.imwrite(FINAL / f"{short}_pia.tif", (tissue & ~arach).astype(np.uint8), compression="zlib")
        ch = lambda a, b: int(((best == a) & (lab == b)).sum() * G * G)
        man["sections"][short] = dict(group=rows[short]["group"], shape=list(lab.shape),
                                      egl_um2=dict(before=int((best == 2).sum() * G * G), after=int((lab == 2).sum() * G * G)),
                                      igl_um2=dict(before=int((best == 4).sum() * G * G), after=int((lab == 4).sum() * G * G)),
                                      changed_um2=int((lab != best).sum() * G * G), arachnoid_um2=int(arach.sum() * G * G),
                                      to_egl_um2={CODES[a]: ch(a, 2) for a in (1, 3, 4, 5, 6) if ch(a, 2)},
                                      from_egl_um2={CODES[b]: ch(2, b) for b in (1, 3, 4, 5, 6) if ch(2, b)})
        print(short, man["sections"][short], flush=True)
    (FINAL / "manifest.json").write_text(json.dumps(man, indent=1))
    return man


def repaired(short):
    """the final repaired labels of a section and its tissue mask (4 um grid, flatten.labels codes)."""
    import tifffile
    return tifffile.imread(FINAL / f"{short}_labels.tif"), tifffile.imread(FINAL / f"{short}_tissue.tif") > 0


def pia_tissue(short):
    """the tissue inside the pial surface (v9: the tissue mask less the arachnoid pia_trim found)."""
    import tifffile
    return tifffile.imread(FINAL / f"{short}_pia.tif") > 0


if __name__ == "__main__":
    if sys.argv[1:2] == ["eval"]:
        eval_swaps()
        sys.exit()
    if sys.argv[1:2] == ["export"]:
        export(sys.argv[2:] or None)
        sys.exit()
    shorts = sys.argv[1:] or F.all_shorts()
    rows = [repair(s) for s in shorts]
    if not sys.argv[1:]:
        (OUT / "summary.json").write_text(json.dumps(dict(prm=PRM, sections=rows), indent=1))
