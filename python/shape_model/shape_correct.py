"""Shape-based layer correction (post-processing of the MATLAB label map).

Idea: every slide calibrates its own shape model. The robust layers (the EGL ribbon
attached to the pial surface and the IGL) are used as anchors; every ML / PCL / DWL
pixel is described by its distance to those anchors; the slide's own labels give
P(class | distances); pixels whose current class is implausible for where they sit
are relabelled. See README.md in this folder for the evidence behind each rule.

Label index (canonical): 0 bg, 1 cereb, 2 iEGL, 3 oEGL, 4 IGL, 5 ML, 6 DWL, 7 PCL, 8 DCN, 9 anomaly.

Pipeline (correct()):
  1. egl_rescue      DCN that is EGL-dense (DAPI high, NeuN low), lies within 100 um of IGL and
                     is connected to the EGL ribbon goes back to EGL (undoes the iter-6
                     EGL->DCN over-conversion; 94% of these pixels are gold EGL on s1_2).
  2. relabel         x3 iterations of the anchor-distance model on {ML, PCL, DWL/DCN}.
  3. revert_thin_dwl new DWL that only forms a thin sliver (<15 um half-width) is undone.
  4. revert_dwl_islands new DWL not connected to existing DWL/DCN is undone (DWL is one tree).
  5. protect_dwl_tree DWL->ML changes that sever an IGL-wrapped white-matter core are undone.
"""
from __future__ import annotations
import numpy as np
from scipy import ndimage as ndi
from skimage import measure, morphology

EGL, ML, PCL, IGL, DWL, DCN = (2, 3), 5, 7, 4, 6, 8
OUT_CLS = np.array([ML, PCL, DWL], np.uint8)          # model classes 0,1,2


# ----------------------------------------------------------------------------- inputs
def smooth_channels(img, ds, order=("p27", "neun", "dapi"), win=9):
    """img: (3,H,W) or (H,W,3). order: channel roles of pages 0,1,2
    (20x = p27,neun,dapi; 10x = p27,dapi,neun). Returns dict of ~9 um box-smoothed maps at /ds."""
    if img.shape[0] == 3 and img.shape[-1] != 3:
        img = np.moveaxis(img, 0, -1)
    img = img[::ds, ::ds].astype(np.float32)
    return {n: ndi.uniform_filter(img[..., c], win) for c, n in enumerate(order)}


def _med(a, m, default=0.):
    return float(np.median(a[m])) if m.any() else default


def _dapi_rel(lab, ch):
    """DAPI scaled so slide-median DWL = 0 and slide-median EGL = 1."""
    egl = np.isin(lab, EGL)
    dE, dW = _med(ch["dapi"], egl), _med(ch["dapi"], lab == DWL)
    return (ch["dapi"] - dW) / max(dE - dW, 1e-6)


# ----------------------------------------------------------------------------- step 1
def egl_rescue(lab, ch, um, dapi_thr=0.6, neun_thr=0.5, max_igl_um=100., min_area=50,
               sources=(DCN,), dwl_dapi_thr=1.0, link_gap_um=20.):
    """Return DCN (and optionally DWL) pixels with EGL-like signal to the EGL ribbon."""
    lab = lab.copy()
    egl = np.isin(lab, EGL)
    if not egl.any():
        return lab, 0
    dapi_rel = _dapi_rel(lab, ch)
    nE, nI = _med(ch["neun"], egl), _med(ch["neun"], lab == IGL)
    neun_rel = (ch["neun"] - nE) / max(nI - nE, 1e-6)
    src = np.zeros(lab.shape, bool)
    if DCN in sources:
        src |= (lab == DCN) & (dapi_rel > dapi_thr)
    if DWL in sources:   # off by default: P7 white matter is cell-dense and fooled this (s1_4)
        src |= (lab == DWL) & (dapi_rel > dwl_dapi_thr)
    cand = src & (neun_rel < neun_thr)
    if max_igl_um:       # pial EGL has IGL beneath it; a cut white-matter edge does not
        cand &= ndi.distance_transform_edt(lab != IGL) * um <= max_igl_um
    cand = morphology.remove_small_objects(cand, max_size=min_area)
    if link_gap_um:   # same ribbon-continuity tolerance as the anchors: short gaps allowed
        rad = max(int(round(link_gap_um / 2 / um)), 1)
        ccE = measure.label(ndi.binary_dilation(egl | cand, iterations=rad), connectivity=2)
        ok = np.unique(ccE[egl])
        grown = cand & np.isin(ccE, ok[ok > 0])
    else:
        grown = ndi.binary_propagation(egl, mask=egl | cand) & cand
    pi, po = _med(ch["p27"], lab == 2), _med(ch["p27"], lab == 3)
    is_i = ch["p27"] > 0.5 * (pi + po)                # iEGL is the p27-high half
    lab[grown & is_i] = 2
    lab[grown & ~is_i] = 3
    return lab, int(grown.sum())


# ----------------------------------------------------------------------------- anchors
def egl_anchor(lab, d_bg, um, ch=None, surf_um=15., keep_area_um2=8000., keep_dapi=1.3,
               link_gap_um=20., neun_max=0.75):
    """EGL components trusted as real EGL:
       touching the outer surface | big ribbons (fissure EGL) | clearly EGL-dense by DAPI |
       chained to any of those by gaps < link_gap_um (ribbon continuity, short jumps);
       minus pixels whose local NeuN is IGL-like (neun_max on the EGL=0 / IGL=1 scale)."""
    egl = np.isin(lab, EGL)
    cc = measure.label(egl, connectivity=2)
    keep = set(np.unique(cc[egl & (d_bg < surf_um)]).tolist()) - {0}
    rel = _dapi_rel(lab, ch) if ch is not None else None
    for r in measure.regionprops(cc):
        if r.label in keep:
            continue
        a_um2 = r.area * um * um
        if a_um2 >= keep_area_um2:
            keep.add(r.label)
        elif rel is not None and a_um2 >= 1000. and np.median(rel[r.slice][r.image]) >= keep_dapi:
            keep.add(r.label)
    anc = np.isin(cc, list(keep)) if keep else np.zeros_like(egl)
    if link_gap_um and anc.any():
        rad = max(int(round(link_gap_um / 2 / um)), 1)
        ccE = measure.label(ndi.binary_dilation(egl, iterations=rad), connectivity=2)
        ok = np.unique(ccE[anc])
        anc = np.isin(ccE, ok[ok > 0]) & egl
    if ch is not None and neun_max is not None and anc.any():
        # "EGL" with IGL-level NeuN is mislabelled IGL (s4_C: median NeuN 0.95 on the
        # EGL=0 / IGL=1 scale); it must not anchor ML next to a white-matter trunk
        nE, nI = _med(ch["neun"], egl), _med(ch["neun"], lab == IGL)
        nrel = ndi.uniform_filter((ch["neun"] - nE) / max(nI - nE, 1e-6), 21)
        anc &= nrel < neun_max
    return anc


def betweenness(egl_a, igl):
    """cos of the angle between pixel->nearest EGL anchor and pixel->nearest IGL.
    ~-1: pixel lies between EGL and IGL (ML/PCL geometry); ~+1: both on one side (DWL).
    On s1_2 gold: 86% of ML but only 15% of DWL has cos < -0.3."""
    _, (ey, ex) = ndi.distance_transform_edt(~egl_a, return_indices=True)
    _, (iy, ix) = ndi.distance_transform_edt(~igl, return_indices=True)
    yy, xx = np.indices(igl.shape)
    a0, a1 = (ey - yy).astype(np.float32), (ex - xx).astype(np.float32)
    b0, b1 = (iy - yy).astype(np.float32), (ix - xx).astype(np.float32)
    return (a0 * b0 + a1 * b1) / np.maximum(np.hypot(a0, a1) * np.hypot(b0, b1), 1e-6)


def _hist_posterior(fa, fb, cur, ncls, bin_um, sig, dmax):
    """Smoothed 2-D histogram P(class | fa, fb) built from the slide's own labels."""
    NB = int(dmax / bin_um)
    ba = np.minimum((fa / bin_um).astype(np.int32), NB - 1)
    bb = np.minimum((fb / bin_um).astype(np.int32), NB - 1)
    m = cur >= 0
    H = np.bincount((ba[m] * NB + bb[m]) * ncls + cur[m], minlength=NB * NB * ncls)
    H = ndi.gaussian_filter(H.astype(np.float32).reshape(NB, NB, ncls), (sig, sig, 0))
    P = H / np.maximum(H.sum(-1, keepdims=True), 1e-6)
    return P[ba, bb]


# ----------------------------------------------------------------------------- step 2
def relabel(lab, um, ch=None, n_iter=3, bin_um=5., sig=1.5, dmax=300.,
            thr_new=0.7, thr_old=0.15, thr_mldwl=(0.6, 0.2),
            surf_um=15., cut_band_um=40., between_gate=-0.3, relabel_orphan_egl=False,
            anchor_kw=None):
    lab = lab.copy()
    log = []
    for it in range(n_iter):
        d_bg = ndi.distance_transform_edt(ndi.binary_fill_holes(lab > 0)) * um
        egl = np.isin(lab, EGL)
        egl_a = egl_anchor(lab, d_bg, um, ch, surf_um, **(anchor_kw or {}))
        d_egl = ndi.distance_transform_edt(~egl_a) * um
        d_surf = np.minimum(d_egl, d_bg)          # pia proxy where EGL is missing
        d_igl = ndi.distance_transform_edt(lab != IGL) * um
        cur = np.full(lab.shape, -1, np.int8)
        cur[lab == ML] = 0
        cur[lab == PCL] = 1
        cur[np.isin(lab, (DWL, DCN))] = 2
        post = _hist_posterior(d_surf, d_igl, cur, 3, bin_um, sig, dmax)
        tgt = post.argmax(-1)
        pt = post.max(-1)
        pc = np.take_along_axis(post, np.clip(cur, 0, 2)[..., None].astype(np.int64), -1)[..., 0]
        valid = cur >= 0
        flip = valid & (tgt != cur) & (pt >= thr_new) & (pc <= thr_old)
        if thr_mldwl is not None:                 # ML<->DWL is the main error; looser bar
            mldwl = ((cur == 0) & (tgt == 2)) | ((cur == 2) & (tgt == 0))
            flip |= valid & mldwl & (pt >= thr_mldwl[0]) & (pc <= thr_mldwl[1])
        flip &= ~((cur == 1) & (tgt == 2))        # never PCL -> DWL
        flip &= ~((cur == 2) & (d_bg < cut_band_um) & (d_egl > cut_band_um))   # cut surface
        if between_gate is not None:              # ML/PCL only created between EGL and IGL
            cosb = betweenness(egl_a, lab == IGL)
            flip &= ~((cur == 2) & (tgt != 2) & (cosb > between_gate))
        lab[flip] = OUT_CLS[tgt[flip]]
        n_orph = 0
        if relabel_orphan_egl:                    # off by default: can delete fissure EGL
            orphan = egl & ~egl_a & (pt >= thr_new)
            lab[orphan] = OUT_CLS[tgt[orphan]]
            n_orph = int(orphan.sum())
        log.append(dict(iter=it, flipped=int(flip.sum()), orphan_egl=n_orph))
    return lab, log


# ----------------------------------------------------------------------------- steps 3-4
def revert_thin_dwl(orig, lab, um, min_half_um=15., win_um=40.):
    dwl = np.isin(lab, (DWL, DCN))
    new = dwl & ~np.isin(orig, (DWL, DCN))
    if not new.any():
        return lab, 0
    half = ndi.distance_transform_edt(dwl) * um
    w = max(int(round(win_um / um)), 1)
    bad = new & (ndi.maximum_filter(half, size=2 * w + 1) < min_half_um)
    lab = lab.copy()
    lab[bad] = orig[bad]
    return lab, int(bad.sum())


def revert_dwl_islands(orig, lab):
    dwl = np.isin(lab, (DWL, DCN))
    old = np.isin(orig, (DWL, DCN))
    cc = measure.label(dwl, connectivity=2)
    ok = np.unique(cc[old & dwl])
    bad = dwl & ~old & ~np.isin(cc, ok[ok > 0])
    lab = lab.copy()
    lab[bad] = orig[bad]
    return lab, int(bad.sum())


def protect_dwl_tree(orig, lab, um, min_piece_um2=5000., min_cut_um2=10000., min_igl_ring=0.2,
                     n_rounds=5):
    """Undo DWL->cortex changes that cut the white-matter tree.

    For every significant pre-existing DWL/DCN component, the part that is still DWL should
    stay in one piece. Pieces cut off from it (>= 10000 um2) whose border is >= 20% IGL
    mark a cut; every connected blob of DWL->non-DWL changes touching such a piece is
    reverted. Repeat until stable.
    (s4_C 10x: a trunk running beside a fissure was otherwise turned into ML; its cut-off
    pieces had 48-64% IGL borders and were 77k-202k um2. The correct s1_2 cut had 0% IGL;
    the 1-2k um2 crumbs left on s2_2 after restoring ML must not trigger a revert.)"""
    old = np.isin(orig, (DWL, DCN))
    px = um * um
    occ = measure.label(old, connectivity=2)
    sizes = np.bincount(occ.ravel())
    big = np.flatnonzero(sizes * px >= min_piece_um2)
    big = big[big > 0]
    lab = lab.copy()
    total = 0
    for _ in range(n_rounds):
        now = np.isin(lab, (DWL, DCN))
        removed = old & ~now
        if not removed.any():
            break
        ncc = measure.label(now, connectivity=2)
        cut = np.zeros(lab.shape, bool)
        for o in big:
            keep = (occ == o) & now
            ids, cnt = np.unique(ncc[keep], return_counts=True)
            if len(ids) <= 1:
                continue
            main = ids[np.argmax(cnt)]
            for i, c in zip(ids, cnt):
                if i == main or c * px < min_cut_um2:
                    continue
                piece = keep & (ncc == i)
                ring = ndi.binary_dilation(piece, iterations=3) & ~piece
                # real white-matter core is wrapped in IGL; a DWL strip bordered only by
                # EGL/ML/background (s1_2 base) is not white matter and may be cut off
                if (lab[ring] == IGL).mean() >= min_igl_ring:
                    cut |= piece
        if not cut.any():
            break
        rcc = measure.label(removed, connectivity=2)
        hit = np.unique(rcc[ndi.binary_dilation(cut, iterations=2) & removed])
        undo = np.isin(rcc, hit[hit > 0])
        lab[undo] = orig[undo]
        total += int(undo.sum())
    return lab, total


# ----------------------------------------------------------------------------- driver
def correct(lab, um, ch=None, rescue=True, rescue_kw=None, thin_dwl_um=15.,
            dwl_islands=True, protect_tree=True, **relabel_kw):
    """lab: uint8 label map at working resolution `um` (um/px, ~1 um recommended).
    ch: smooth_channels() output at the same resolution (needed for rescue/anchors)."""
    orig = lab
    log = {}
    if rescue and ch is not None:
        lab, log["rescue"] = egl_rescue(lab, ch, um, **(rescue_kw or {}))
    lab, log["relabel"] = relabel(lab, um, ch=ch, **relabel_kw)
    if thin_dwl_um:
        lab, log["thin_dwl_revert"] = revert_thin_dwl(orig, lab, um, thin_dwl_um)
    if dwl_islands:
        lab, log["dwl_island_revert"] = revert_dwl_islands(orig, lab)
    if protect_tree:
        lab, log["dwl_tree_revert"] = protect_dwl_tree(orig, lab, um)
    return lab, log


def apply_at_full_res(full_lab, small_orig, small_new, ds):
    """Paste the low-res corrections back onto the native-resolution label map
    (only pixels whose low-res label changed are touched)."""
    if ds == 1:
        return small_new.copy()
    out = full_lab.copy()
    chg = small_orig != small_new
    up_chg = np.repeat(np.repeat(chg, ds, 0), ds, 1)[:out.shape[0], :out.shape[1]]
    up_new = np.repeat(np.repeat(small_new, ds, 0), ds, 1)[:out.shape[0], :out.shape[1]]
    H, W = up_chg.shape
    out[:H, :W][up_chg] = up_new[up_chg]
    return out
