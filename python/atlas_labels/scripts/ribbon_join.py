"""
Rejoin broken ribbon layers in P7 cerebellum label maps by restoring the layer order.

Labels: MATLAB 9-class scheme, 0 bg, 1 cereb, 2 iEGL, 3 oEGL, 4 IGL, 5 ML, 6 DWL,
7 PCL, 8 DCN, 9 anomaly.

From the pia inward the order is always  bg | EGL | ML | (PCL) | IGL | WM.
A broken ribbon shows up as a skipped layer at a contact. On the gold slide these
contacts are rare; on the fragmented 20x slides they are 3-10x more common:

  rule P  ML / IGL / PCL touching the outer surface  -> EGL is missing there.
          Pixels of those labels shallower than the local EGL depth (EGL within
          `win_um` along the surface) become EGL; oEGL vs iEGL by the local oEGL
          depth.
  rule M  IGL touching EGL                             -> ML is missing there.
          IGL pixels closer to the EGL than the local ML thickness, with no ML on
          the straight line to their nearest EGL pixel, become ML.

Both rules only act on pixels connected to the contact itself, and both take the
thickness to restore from the same layer nearby, so the bridge matches its flanks.
The template fit (fit_p7_template.py) is not used: where it disagrees with these
labels the gold sides with the labels 84 % of the time.
"""
import numpy as np
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
from skimage import morphology

ST = np.ones((3, 3), bool)
RIBBONS = {"EGL": dict(codes=(2, 3)), "IGL": dict(codes=(4,)), "ML": dict(codes=(5,))}
ORDER = ("EGL", "IGL", "ML")
DEFAULTS = dict(win_um=150.0, egl_cap_um=60.0, ml_cap_um=60.0, line_samples=12,
                egl_cover_um=10.0, max_run_um=300.0, max_wm_frac=0.2,
                max_run_nowm_um=800.0, max_wm_frac_long=0.02,
                w_margin_um=10.0, m_max_run_um=200.0,
                f_win_um=200.0, f_margin_um=20.0, f_margin_lo_um=3.0, f_grow_um=25.0, f_min_depth_um=30.0, f_wm_um=30.0,
                f_probe_extra_um=30.0, f_seed_um=15.0, f_ref_um=100.0, f_depth_weight=3.0,
                s_close_um=6.0, s_min_depth_um=20.0,
                fw_cover_um=25.0, fw_max_run_um=400.0, fw_min_flanks=1, fw_max_hw_um=30.0,
                fw_offset_um=25.0, fw_ml_margin_um=3.0,
                fw_win_um=300.0, fw_depth_tol_um=10.0,
                rules=("P", "M", "W", "F", "S"))


def _local_max(val, mask, px, win_um):
    """max of val over mask pixels within win_um/2 (0 where none)."""
    r = max(1, int(round(win_um / 2 / px)))
    return ndi.maximum_filter(np.where(mask, val, 0.0), size=2 * r + 1)


def _local_median_map(val, mask, px, win_um, step_um=40.0):
    """median of val over mask pixels in win_um boxes on a step_um lattice, nearest-filled."""
    H, W = mask.shape
    s = max(1, int(round(step_um / px)))
    h = max(1, int(round(win_um / 2 / px)))
    gy, gx = np.arange(0, H, s), np.arange(0, W, s)
    grid = np.full((len(gy), len(gx)), np.nan)
    for i, y in enumerate(gy):
        for j, x in enumerate(gx):
            sl = (slice(max(y - h, 0), y + h + 1), slice(max(x - h, 0), x + h + 1))
            v = val[sl][mask[sl]]
            if v.size >= 10:
                grid[i, j] = np.median(v)
    ok = np.isfinite(grid)
    if not ok.any():
        return np.zeros((H, W))
    idx = ndi.distance_transform_edt(~ok, return_distances=False, return_indices=True)
    grid = grid[tuple(idx)]
    return np.kron(grid, np.ones((s, s)))[:H, :W]


def surface_runs(lab, px, prm):
    """pia ring (tissue pixels touching the exterior) split into runs with / without EGL.
    Returns ring mask, ring-pixel run labels, and the set of run ids to fill:
    no EGL within egl_cover_um, and either ring length <= max_run_um with white
    matter on <= max_wm_frac of the run, or <= max_run_nowm_um with (almost) no
    white matter. Runs touching white matter longer than that are the base cut."""
    tissue = lab > 0
    ext = ~ndi.binary_fill_holes(tissue)
    ring = tissue & ndi.binary_dilation(ext, structure=ST)
    egl_cov = ndi.distance_transform_edt(~np.isin(lab, (2, 3))) * px <= prm["egl_cover_um"]
    wm = ndi.distance_transform_edt(~np.isin(lab, (6, 8))) * px <= prm["egl_cover_um"] / 2
    gap = ring & ~egl_cov
    runs, n = ndi.label(gap, structure=ST)
    ok = set()
    if n:
        length = ndi.sum(np.ones_like(runs), runs, range(1, n + 1)) * px
        wmf = ndi.mean(wm, runs, range(1, n + 1))
        ends = ndi.binary_dilation(gap, structure=ST) & ring & egl_cov
        for i in range(1, n + 1):
            if length[i - 1] <= prm["max_run_um"] and wmf[i - 1] <= prm["max_wm_frac"]:
                ok.add(i)
            elif length[i - 1] <= prm["max_run_nowm_um"] and wmf[i - 1] <= prm["max_wm_frac_long"]:
                ok.add(i)          # long run with no white matter: not the base cut
    return ext, ring, runs, ok


def rule_p(lab, px, prm):
    ext, ring, runs, ok = surface_runs(lab, px, prm)
    depth, (iy, ix) = ndi.distance_transform_edt(~ring, return_indices=True)
    depth = depth * px
    egl = np.isin(lab, (2, 3))
    near = depth <= prm["egl_cap_um"]
    t_egl = _local_max(depth, egl & near, px, prm["win_um"])
    t_o = _local_max(depth, (lab == 3) & near, px, prm["win_um"])
    # accepted runs, extended along the surface up to the flanking EGL
    acc_ring = np.isin(runs, list(ok) or [-1])
    acc_ring = ndi.binary_dilation(acc_ring, structure=ST,
                                   iterations=int(round(prm["egl_cover_um"] / px)) + 1, mask=ring)
    cand = np.isin(lab, (4, 5, 7)) & (depth <= t_egl) & (t_egl > 0) & acc_ring[iy, ix]
    yy, xx = np.nonzero(cand)        # no EGL between the pixel and the surface
    if yy.size:
        good = np.ones(yy.size, bool)
        for f in np.linspace(0.0, 1.0, prm["line_samples"]):
            sy = np.round(yy + f * (iy[yy, xx] - yy)).astype(int)
            sx = np.round(xx + f * (ix[yy, xx] - xx)).astype(int)
            good &= ~egl[sy, sx]
        cand[:] = False
        cand[yy[good], xx[good]] = True
    touch = cand & ring
    keep = morphology.reconstruction(touch.astype(np.uint8), cand.astype(np.uint8)).astype(bool) \
        if touch.any() else touch
    out = lab.copy()
    out[keep] = np.where(depth[keep] <= t_o[keep], 3, 2)
    return out, keep


def rule_m(lab, px, prm):
    egl, igl, ml = np.isin(lab, (2, 3)), lab == 4, lab == 5
    dE, (iy, ix) = ndi.distance_transform_edt(~egl, return_indices=True)
    dE = dE * px
    # ML thickness: distance to EGL + distance to IGL, at ML pixels; local median
    dI = ndi.distance_transform_edt(~igl) * px
    t_ml = np.minimum(_local_median_map(dE + dI, ml & (dE + dI <= 2 * prm["ml_cap_um"]),
                                        px, prm["win_um"]), prm["ml_cap_um"])
    src = igl | (lab == 7) if prm.get("m_pcl", True) else igl
    cand = src & (dE <= t_ml)
    yy, xx = np.nonzero(cand)
    if yy.size:
        ok = np.ones(yy.size, bool)
        for f in np.linspace(0.0, 1.0, prm["line_samples"]):
            sy = np.round(yy + f * (iy[yy, xx] - yy)).astype(int)
            sx = np.round(xx + f * (ix[yy, xx] - xx)).astype(int)
            ok &= ~ml[sy, sx]
        cand[:] = False
        cand[yy[ok], xx[ok]] = True
    touch = cand & ndi.binary_dilation(egl, structure=ST)
    if prm.get("m_max_run_um"):
        # contact runs along the EGL edge: EGL pixels touching IGL/PCL but not ML
        edge = egl & ndi.binary_dilation(src, structure=ST) & ~ndi.binary_dilation(ml, structure=ST)
        runs, n = ndi.label(ndi.binary_dilation(edge, structure=ST), structure=ST)
        if n:
            ln = ndi.sum(edge, runs, range(1, n + 1)) * px
            badids = set((1 + np.nonzero(ln > prm["m_max_run_um"])[0]).tolist())
            if prm.get("m_flanks", True):
                # a broken ML ribbon has normal EGL|ML contact at both ends of the run
                okedge = egl & ndi.binary_dilation(ml, structure=ST)
                for i, sl in enumerate(ndi.find_objects(runs), 1):
                    if i in badids:
                        continue
                    sl2 = tuple(slice(max(a.start - 3, 0), a.stop + 3) for a in sl)
                    fl = ndi.binary_dilation(runs[sl2] == i, structure=ST, iterations=2) & okedge[sl2]
                    if ndi.label(fl, structure=ST)[1] < 2:
                        badids.add(i)
            bad = np.isin(runs, list(badids) or [-1])
            touch &= ~ndi.binary_dilation(bad, structure=ST)
    keep = morphology.reconstruction(touch.astype(np.uint8), cand.astype(np.uint8)).astype(bool) \
        if touch.any() else touch
    out = lab.copy()
    out[keep] = 5
    return out, keep


def rule_w(lab, px, prm):
    """IGL missing between ML and WM: ML deeper than the local ML thickness that
    touches white matter, with no IGL on the line to the nearest WM -> IGL."""
    egl, igl, ml, wm = np.isin(lab, (2, 3)), lab == 4, lab == 5, np.isin(lab, (6, 8))
    dE = ndi.distance_transform_edt(~egl) * px
    dI = ndi.distance_transform_edt(~igl) * px
    dW, (iy, ix) = ndi.distance_transform_edt(~wm, return_indices=True)
    t_ml = np.minimum(_local_median_map(dE + dI, ml & (dE + dI <= 2 * prm["ml_cap_um"]),
                                        px, prm["win_um"]), prm["ml_cap_um"])
    cand = ml & (dE > t_ml + prm["w_margin_um"])
    yy, xx = np.nonzero(cand)
    if yy.size:
        ok = np.ones(yy.size, bool)
        for f in np.linspace(0.0, 1.0, prm["line_samples"]):
            sy = np.round(yy + f * (iy[yy, xx] - yy)).astype(int)
            sx = np.round(xx + f * (ix[yy, xx] - xx)).astype(int)
            ok &= ~igl[sy, sx]
        cand[:] = False
        cand[yy[ok], xx[ok]] = True
    touch = cand & ndi.binary_dilation(wm, structure=ST)
    keep = morphology.reconstruction(touch.astype(np.uint8), cand.astype(np.uint8)).astype(bool) \
        if touch.any() else touch
    out = lab.copy()
    out[keep] = 4
    return out, keep


def rule_f(lab, px, prm):
    """EGL missing inside a closed fissure. In a folium the EGL starts one ML
    thickness from the IGL; the local value of that is read off the existing
    ML|EGL boundary (distance to IGL there). ML farther than that from the IGL
    (+ f_margin_um), deeper than f_min_depth_um below the outer surface, not near
    white matter, in a region touching EGL -> EGL, but only where walking away
    from the pixel's own IGL reaches another IGL before any EGL (two folia meeting
    without EGL between them). Seeds need f_margin_um past the ML edge; they grow
    f_grow_um through the f_margin_lo_um band. iEGL / oEGL copied from the nearest
    flanking fissure EGL at the same distance from the IGL."""
    tissue = lab > 0
    ext = ~ndi.binary_fill_holes(tissue)
    dO = ndi.distance_transform_edt(~ext) * px
    egl, igl, ml = np.isin(lab, (2, 3)), lab == 4, lab == 5
    wm = np.isin(lab, (6, 8))
    dI = ndi.distance_transform_edt(~igl) * px
    dW = ndi.distance_transform_edt(~wm) * px
    ml_edge = ml & ndi.binary_dilation(egl, structure=ST)
    t_me = _local_median_map(dI, ml_edge & (dI <= prm["ml_cap_um"] * 2), px, prm["f_win_um"])
    io_edge = (lab == 2) & ndi.binary_dilation(lab == 3, structure=ST)
    t_io = _local_median_map(dI, io_edge & (dI <= prm["ml_cap_um"] * 3), px, prm["f_win_um"])
    band = ml & (t_me > 0) & (dI > t_me + prm["f_margin_lo_um"]) & (dO > prm["f_min_depth_um"]) \
        & (dW > prm["f_wm_um"])
    # two folia meet without EGL: walking away from the pixel's own IGL reaches
    # another IGL before any EGL or background
    _, (iy, ix) = ndi.distance_transform_edt(~igl, return_indices=True)
    yy, xx = np.nonzero(band)
    seed = np.zeros(lab.shape, bool)
    if yy.size:
        v = np.stack([yy - iy[yy, xx], xx - ix[yy, xx]], 1).astype(float)
        v /= np.maximum(np.hypot(v[:, 0], v[:, 1]), 1e-9)[:, None]
        lmax = (2 * t_me[yy, xx] + prm["f_probe_extra_um"]) / px
        H, W = lab.shape
        state = np.zeros(yy.size, np.int8)          # 0 running, 1 hit IGL, -1 blocked
        for k in range(1, int(lmax.max()) + 1):
            run = (state == 0) & (k <= lmax)
            if not run.any():
                break
            sy = np.clip(np.round(yy[run] + k * v[run, 0]).astype(int), 0, H - 1)
            sx = np.clip(np.round(xx[run] + k * v[run, 1]).astype(int), 0, W - 1)
            L = lab[sy, sx]
            st = np.where(L == 4, 1, np.where(np.isin(L, (0, 2, 3, 6, 8, 9)), -1, 0)).astype(np.int8)
            idx = np.nonzero(run)[0]
            state[idx] = st
        seed[yy[state == 1], xx[state == 1]] = True
    probe_ok = seed.copy()
    # seeds: well past the local ML edge (precise); then grow a limited distance
    # through probe-confirmed pixels at the loose margin to fill the EGL width
    seed &= (dI > t_me + prm["f_margin_um"]) & (ndi.distance_transform_edt(~egl) * px > prm["f_seed_um"])
    keep = seed.copy()
    if seed.any():
        keep = ndi.binary_dilation(seed, structure=ST, iterations=max(1, int(round(prm["f_grow_um"] / px))),
                                   mask=band)
    if keep.any():
        cc, n = ndi.label(keep, structure=ST)
        touch = np.unique(cc[ndi.binary_dilation(egl, structure=ST) & keep])
        keep = np.isin(cc, touch[touch > 0])
    out = lab.copy()
    if keep.any():
        # iEGL / oEGL: copy the nearest flanking fissure-EGL pixel at the same
        # distance from the IGL (distance weighted f_depth_weight x position)
        ref = egl & (dO > prm["f_min_depth_um"])
        ref &= ndi.binary_dilation(keep, structure=ST, iterations=int(round(prm["f_ref_um"] / px)))
        if not ref.any():
            ref = egl & ndi.binary_dilation(keep, structure=ST, iterations=int(round(prm["f_ref_um"] / px)))
        ky, kx = np.nonzero(keep)
        ry, rx = np.nonzero(ref)
        if ry.size:
            w = prm["f_depth_weight"]
            tree = cKDTree(np.stack([ry * px, rx * px, w * dI[ry, rx]], 1))
            _, nn = tree.query(np.stack([ky * px, kx * px, w * dI[ky, kx]], 1))
            out[ky, kx] = lab[ry[nn], rx[nn]]
        else:
            out[keep] = np.where(dI[keep] > np.maximum(t_io[keep], t_me[keep]), 3, 2)
    return out, keep


def rule_s(lab, px, prm):
    """thin seams in the EGL: the two folia's EGL layers in a fissure separated by a
    thin strip of ML / PCL / cereb (where they should be fused at the oEGL contact),
    or any EGL split by a thin strip of IGL. Pixels added by closing the EGL with a
    disk of s_close_um take the label of the nearest EGL pixel; ML / PCL / cereb
    only deeper than s_min_depth_um (the surface rim is left to rule P)."""
    tissue = lab > 0
    ext = ~ndi.binary_fill_holes(tissue)
    dO = ndi.distance_transform_edt(~ext) * px
    egl = np.isin(lab, (2, 3))
    r = max(1, int(round(prm["s_close_um"] / px)))
    pad = r + 2
    closed = ndi.binary_closing(np.pad(egl, pad), structure=morphology.disk(r))[pad:-pad, pad:-pad]
    keep = closed & ~egl & ((np.isin(lab, (1, 5, 7)) & (dO > prm["s_min_depth_um"])) | (lab == 4))
    out = lab.copy()
    if keep.any():
        _, (iy, ix) = ndi.distance_transform_edt(~egl, return_indices=True)
        out[keep] = lab[iy[keep], ix[keep]]
    return out, keep




def rule_fw(lab, px, prm):
    """[off by default: on s1_2 it adds nothing to P, M, W, F, S and costs 0.0005-0.0007
    on the clean slide; on the real slides it does not raise EGL continuity. The
    white-matter line is as uncertain (20-30 um) as the fused EGL is wide.]
    EGL missing along a fissure predicted from the white matter (wm_scaffold).
    Line pixels with no EGL within fw_cover_um form gap runs; runs <= fw_max_run_um
    with covered line (fissure EGL) at >= fw_min_flanks ends are extended to the
    flanks and painted EGL (any of ML / IGL / PCL) in a band as wide as the
    flanking fused EGL. iEGL / oEGL copied from the flanking EGL at the same
    distance from the line."""
    import wm_scaffold as WS
    sc = WS.fissure_lines(lab, px)
    L = sc["lines"]
    egl = np.isin(lab, (2, 3))
    dEg = ndi.distance_transform_edt(~egl) * px
    # EGL starts one cortex thickness out from the white matter: local distance from
    # WM to the EGL inner edge; line points closer to the WM are below the fissure bottom
    inner = egl & ndi.binary_dilation(np.isin(lab, (4, 5, 7)), structure=ST)
    t_we = _local_median_map(sc["dW"], inner & (sc["dW"] < 400), px, prm["fw_win_um"])
    deep_ok = (t_we > 0) & (sc["dW"] >= t_we - prm["fw_depth_tol_um"])
    cov = L & (dEg <= prm["fw_cover_um"])
    gap = L & ~cov & deep_ok
    igl, ml = lab == 4, lab == 5
    dI = ndi.distance_transform_edt(~igl) * px
    ml_edge = ml & ndi.binary_dilation(egl, structure=ST)
    t_me = _local_median_map(dI, ml_edge & (dI <= prm["ml_cap_um"] * 2), px, prm["f_win_um"])
    ml_beyond = (t_me > 0) & (dI > t_me + prm["fw_ml_margin_um"])
    runs, n = ndi.label(gap, structure=ST)
    out = lab.copy()
    keep = np.zeros(lab.shape, bool)
    if not n:
        return out, keep
    hw_map = ndi.distance_transform_edt(egl) * px     # at the line centre: fused-EGL half width
    ext = int(round(prm["fw_cover_um"] / px)) + 1
    paint_line = np.zeros(lab.shape, bool)
    for i, sl in enumerate(ndi.find_objects(runs), 1):
        m = runs[sl] == i
        if m.sum() * px > prm["fw_max_run_um"]:
            continue
        pad = ext + int(round(80 / px))
        sl2 = tuple(slice(max(a.start - pad, 0), a.stop + pad) for a in sl)
        mm = runs[sl2] == i
        Lm = L[sl2]
        grown = ndi.binary_dilation(mm, structure=ST, iterations=ext, mask=Lm)
        flank = ndi.binary_dilation(grown, structure=ST) & cov[sl2] & ~grown
        if ndi.label(flank, structure=ST)[1] < prm["fw_min_flanks"]:
            continue
        near = ndi.binary_dilation(grown, structure=ST, iterations=int(round(60 / px))) & cov[sl2]
        hw = np.median(hw_map[sl2][near]) if near.any() else 0.0
        if hw <= 0:
            continue
        paint_line[sl2] |= grown
        hw = min(hw, prm["fw_max_hw_um"])
        dl = ndi.distance_transform_edt(~grown) * px
        zone = (dl <= hw + prm["fw_offset_um"]) & deep_ok[sl2]
        sub = lab[sl2]
        # IGL / PCL across a fissure is wrong anywhere in the zone; ML only where it
        # lies beyond the local ML thickness from the IGL (the EGL's place)
        keep[sl2] |= zone & (np.isin(sub, (4, 7)) | ((sub == 5) & ml_beyond[sl2]))
    if keep.any():
        dline = ndi.distance_transform_edt(~paint_line) * px
        ref = egl & ndi.binary_dilation(keep, structure=ST, iterations=int(round(100 / px)))
        ky, kx = np.nonzero(keep)
        ry, rx = np.nonzero(ref)
        if ry.size:
            dl_ref = ndi.distance_transform_edt(~L) * px
            w = prm["f_depth_weight"]
            tree = cKDTree(np.stack([ry * px, rx * px, w * dl_ref[ry, rx]], 1))
            _, nn = tree.query(np.stack([ky * px, kx * px, w * dline[ky, kx]], 1))
            out[ky, kx] = lab[ry[nn], rx[nn]]
        else:
            out[keep] = 3
    return out, keep


RULES = {"P": rule_p, "M": rule_m, "W": rule_w, "F": rule_f, "Fw": rule_fw, "S": rule_s}


def join(lab, px, prm=None, log=None):
    prm = {**DEFAULTS, **(prm or {})}
    out = lab.copy()
    for r in prm["rules"]:
        out, k = RULES[r](out, px, prm)
        if log is not None:
            log.append(dict(rule=r, changed_um2=float(k.sum() * px * px)))
    return out, []
