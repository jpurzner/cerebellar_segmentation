"""
Flattened cerebellum (2026-11-21): unroll the cortical sheet along the outer face of the EGL, so
the layers become horizontal bands and errors in the IGL call show up as breaks in a band.

Reference line (4 um grid; labels from igl_detect: reviewed controls, the user's corrected cKO
sections, MATLAB for s1_5):
  - the outer face of the EGL: the tissue edge on open folia;
  - in closed fissures, where two EGL sheets are fused with no background between them, the
    fissure midline. Seeds: EGL cells whose nearest cortex cell (ML / IGL / WM) lies in the
    opposite direction from their neighbour's (the two folia). Each piece of seeds is joined to
    the open surface by the cheapest path through the EGL, cheap along the middle of the fused
    EGL, so the cut follows the fissure out to its mouth. Cuts shallower than 60 um are dropped
    (a short fused fissure does not bend the layers beneath it).
  - The tissue is cut along the midlines; the outer contour of the cut tissue runs down one side
    of every fissure and back up the other. Positions smoothed 8 um.
Columns: streamlines of the smoothed depth field (signed distance to the cut tissue edge), the
standard construction for cortical depth: they leave the reference at right angles, bend smoothly
around fissure ends and never cross. The direction comes from the field smoothed 8 um near the
surface, blending to 24 um by 40 um deep. A column is followed in 2 um steps to 300 um and
stops when it leaves the tissue, turns back, or passes the core of the folium (depth to the
nearest edge falls by 8 um, or stops rising: it is running along the core).

Strips (position along the reference x depth, 2 um x 2 um): nuclei (DAPI), NeuN, p27, and the
MATLAB and current label calls. Per column: EGL bottom, IGL top and bottom. A line is fitted
through the IGL top and bottom (rolling median 150 um, Gaussian 40 um) inside each stretch where
the IGL is present; gaps up to 300 um are bridged, longer ones are left as real discontinuities.
Columns flagged: IGL missing under a bridged line, IGL boundary off the line, IGL touching the
EGL with no ML between.

    python scripts/flatten.py s5_B s2_4      # some sections
    python scripts/flatten.py                # all 19 -> outputs/flatten/<short>.npz + summary.json
    python scripts/flatten.py --repaired     # on the final EGL repair (egl_repair.py): its fissure cuts (less
                                             # stretches with no folium on one side) and repaired labels
                                             # -> outputs/flatten_v11/ (v11, 2026-11-24, the user's v11 notes: the walk keeps to the
                                             # perimeter: corners cut by arcs, no jumps or hops over cuts, walk.smooth v2; no closed
                                             # slits, join_stranded;
                                             # v10, 2026-11-24, the user's v10 notes: the line is one walk from the
                                             # left end of the cut surface, clockwise, that may not turn sharply except at a furrow base or
                                             # a fissure mouth and turns round only at a furrow base; cuts ending inside a fused EGL are
                                             # continued; thin ML caps the stain calls EGL relabelled: walk.py;
                                             # v9, 2026-11-24, the user's walkthrough note: the EGL almost never
                                             # makes sharp turns; a stranded cut piece leaves along its own course, _smooth_ways, and the
                                             # line is shape-checked, shape_check.py; the grid it used is saved as lab_grid / pia;
                                             # v8, 2026-11-24: on EGL repair v10, meninges on ML / IGL trimmed in
                                             # the cKO sections; flatten_v7, 2026-11-23: cuts walked on down to the furrow base at every
                                             # sharp EGL peak of the strip, peak_cuts; stain EGL counts in prune_cuts' dead_noegl test;
                                             # flatten_v6: EGL repair v9, the line follows the pia, not the
                                             # arachnoid outside it: ER.pia_tissue; flatten_v5 = repair v8, cut ends kept; flatten_v4 = repair v7; flatten_v3 = repair v6,
                                             # cuts on the pia; flatten_v2 = repair v5). The page for v4 and later
                                             # is flat_fix_page.py (corrections on top of this run).
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from skimage import measure, morphology
from skimage.graph import MCP_Geometric
from skimage.transform import downscale_local_mean

sys.path.insert(0, str(Path(__file__).parent))
import igl_detect as D  # noqa: E402
import igl_gapfill_apply as GA  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/flatten"
G = D.G                       # 4 um label grid
STEP = 2.0                    # um per strip pixel, along the reference and in depth
PRM = dict(join_stranded=True, depth_um=300.0, pos_sigma_um=8.0, dir_sigma_um=16.0, near_sigma_um=8.0, far_sigma_um=24.0, blend_um=40.0,
           ridge_drop_um=8.0, stall_um=20.0, stall_frac=0.3, mid_jump_um=12.0, mid_cos=-0.3, min_cut_um=60.0, merge_um=8.0, spur_um=24.0, speck_um2=1500.0,
           cover_um=20.0, skip_um=40.0, cover_run_um=60.0, min_contour_um=300.0,
           min_igl_um=8.0, med_win_um=150.0, med_frac=0.3, fit_sigma_um=40.0, bridge_um=300.0,
           off_top_um=16.0, off_bot_um=30.0, no_ml_um=4.0, min_flag_um=8.0,
           dead_short_um=60.0, dead_frac=0.4, dead_win_um=40.0,          # prune_cuts (repaired mode)
           keep_ends=False, branch_um=24.0, final_keep_ends=True, furrow_close_um=40.0, fill_holes=True, dead_noegl=True,
           drop_sliver=True, sliver_um2=20000.0, sliver_cortex_um2=400.0, sliver_open=0.25,
           # v7 (the user's flatten v7 notes, furrows still missed): stain EGL at a column's top counts as EGL for the
           # dead_noegl test (s4_A at 7.5 mm: the fused EGL at the furrow base is labelled white matter) ...
           dead_stain=True,
           # ... and every sharp EGL peak left in the strip is followed down the fused EGL to the furrow base (peak_cuts)
           peak_cuts=True, pk_thr_um=24.0, pk_win_um=200.0, pk_sep_um=60.0, pk_start_k=5, pk_start_um=0.0,
           pk_turn_deg=30.0, pk_max_turn_deg=60.0, pk_look=3, pk_t_frac=0.5, pk_stop_frac=0.5, pk_stop_min_um=12.0,
           pk_side=2.5, pk_rays=24, pk_side_deg=45.0, pk_balance=2.5, pk_into_frac=0.5, pk_prog_k=4, pk_prog_frac=0.35,
           pk_max_fail=2, pk_pass_frac=0.7, pk_min_len_um=24.0, pk_egl_frac=0.6, pk_straight=0.75, pk_speck_cells=8, pk_check=True, pk_check_um=150.0, pk_short_um=40.0,
           # v9 (the user's walkthrough note, 2026-11-24: no sharp turns): a stranded cut piece leaves along its own course
           smooth_way=True, way_um=24.0, way_turn_deg=60.0, way_reach_um=600.0, way_alt_max_um=100.0, way_drop=True,
           # ... and the line is shape-checked (shape_check.py: meninges the line runs round trimmed, narrow detours filled
           # or cut off) before it is built
           shape_check=True, reprune=True,
           # v10 (the user's v10 notes: a mandatory movement, clockwise from the left end of the cut surface, no jumps, no
           # sharp turns, fix the small jogs): the line is walked (walk.py)
           walk=True)
           # v5 notes (2026-11-22): the last thinning passes keep the cuts' ends; depth measured past open furrows; tiny
           # holes inside the tissue (vessels) are not a surface a cut can end at


def sid_of(short):
    return next(s for s in GA.section_list() if GA.short(s) == short)


def labels(short):
    """grid codes (0 bg, 1 tissue, 2 EGL, 3 ML/PCL, 4 IGL, 5 WM, 6 anomaly): MATLAB and current best."""
    _, y, auto, tissue = D.load(short)
    if (y > 0).any():
        best = y.copy()
        fill = tissue & (best == 0)
        best[fill] = np.where(np.isin(auto[fill], (1, 6)), auto[fill], 1)
    else:
        best = auto.copy()
    return auto, best, tissue


def _pair(a, dy, dx):
    """views (a[r, c], a[r + dy, c + dx]) over the cells where both exist."""
    H, W = a.shape[-2:]
    y0, y1 = max(0, -dy), H - max(0, dy)
    x0, x1 = max(0, -dx), W - max(0, dx)
    return a[..., y0:y1, x0:x1], a[..., y0 + dy:y1 + dy, x0 + dx:x1 + dx], (slice(y0, y1), slice(x0, x1)), \
        (slice(y0 + dy, y1 + dy), slice(x0 + dx, x1 + dx))


def midline(lab, tissue, prm=PRM):
    """fissure midlines inside fused EGL, each joined to the open surface (see the module notes)."""
    bg = ~tissue
    cortex = np.isin(lab, (3, 4, 5))
    big = morphology.remove_small_objects(cortex, max_size=int(prm["speck_um2"] / G ** 2), connectivity=2)
    E = (lab == 2) | (cortex & ~big)                       # ML / IGL specks inside a ragged EGL count as EGL here
    cortex = big
    dist, ind = ndi.distance_transform_edt(~cortex, return_indices=True)
    rr, cc = np.indices(E.shape)
    v = np.stack([ind[0] - rr, ind[1] - cc]).astype(np.float32)
    cand = E & ~ndi.binary_dilation(bg) & (dist >= 1.5)
    seeds = np.zeros_like(E)
    jump = prm["mid_jump_um"] / G
    for dy, dx in ((0, 1), (1, 0), (1, 1), (1, -1)):
        ca, cb, sa, sb = _pair(cand, dy, dx)
        va, vb, _, _ = _pair(v, dy, dx)
        ia, ib, _, _ = _pair(ind.astype(np.float32), dy, dx)
        dot = (va * vb).sum(0)
        na, nb = np.hypot(*va), np.hypot(*vb)
        far = np.hypot(*(ia - ib)) > jump
        hit = ca & cb & (dot < prm["mid_cos"] * na * nb) & far
        seeds[sa] |= hit
        seeds[sb] |= hit
    seeds = morphology.remove_small_objects(seeds, max_size=2, connectivity=2)
    soft = np.isin(lab, (1, 6)) & ~seeds                     # unassigned tissue, tears / bleach lining a fissure
    passable = E | seeds | soft
    cost = np.where(passable, 1.0 / (1.0 + dist) ** 2, np.inf)
    cost[seeds] = np.minimum(cost[seeds], 0.02)
    cost[soft] = np.maximum(cost[soft], 0.3 / (1.0 + dist[soft]))
    starts = np.argwhere(passable & ndi.binary_dilation(bg))
    paths = np.zeros_like(seeds)
    if len(starts) and seeds.any():
        m = MCP_Geometric(cost, fully_connected=True)
        cum, _ = m.find_costs(starts)
        comp, _ = ndi.label(seeds, structure=np.ones((3, 3), bool))
        touch = ndi.binary_dilation(bg)
        for k, sl in enumerate(ndi.find_objects(comp), start=1):
            cells = np.argwhere(comp[sl] == k) + [sl[0].start, sl[1].start]
            if touch[cells[:, 0], cells[:, 1]].any():
                continue
            cv = cum[cells[:, 0], cells[:, 1]]
            if not np.isfinite(cv).any():
                continue
            for r_, c_ in m.traceback(tuple(cells[int(np.argmin(cv))])):
                paths[r_, c_] = True
    mid = single_line(seeds | paths, bg, prm) & tissue
    comp, nc = ndi.label(mid, structure=np.ones((3, 3), bool))
    if nc:
        deep = np.asarray(ndi.maximum(cut_depth(tissue, prm), comp, index=np.arange(1, nc + 1)))
        mid = np.r_[False, deep >= prm["min_cut_um"]][comp]
    return mid


def cut_depth(tissue, prm=PRM):
    """how deep a cut reaches (um), for dropping shallow ones: the distance to the background, with open furrows
    narrower than 2 x furrow_close_um filled in (v5 notes): a short fused stretch at the bottom of a deep open furrow
    (s5_F at 10.2 mm) is as deep as the furrow, while a fused notch in the open surface stays shallow."""
    r = prm.get("furrow_close_um", 0.0)
    if r <= 0:
        return ndi.distance_transform_edt(tissue) * G
    k = int(round(r / G))
    tc = ndi.binary_closing(np.pad(tissue, k + 1), structure=morphology.disk(k).astype(bool))[k + 1:-(k + 1), k + 1:-(k + 1)]
    return ndi.distance_transform_edt(tc | tissue) * G


def skel_ends(sk):
    """end cells of a skeleton by the crossing number (one run of line cells around the 8-neighbourhood): unlike a
    neighbour count, a corner cell of the little triangles skeletonize leaves at junctions counts as an end once
    its branch is gone, so trimming does not stall there."""
    p = np.pad(sk, 1)
    ring = [p[:-2, 1:-1], p[:-2, 2:], p[1:-1, 2:], p[2:, 2:], p[2:, 1:-1], p[2:, :-2], p[1:-1, :-2], p[:-2, :-2]]
    cn = sum((~ring[i]) & ring[(i + 1) % 8] for i in range(8)).astype(int)
    return sk & (cn <= 1)


def final_prm(prm):
    """the thinning for the last passes over a finished cut (pia snap, pruning): side spurs up to spur_um still go, the
    line's own ends are kept (keep_ends), so the cut is not shortened by spur_um at every pass."""
    return dict(prm, keep_ends=True, branch_um=prm["spur_um"]) if prm.get("final_keep_ends") else prm


def cut_tissue(tissue, prm=PRM):
    """the tissue as the cuts see it: tiny holes inside it (vessels, <= 400 um2 in every section) filled, so a cut is
    never joined to one instead of the open surface (s5_E: two cuts ending in such holes became closed slits, each a
    separate loop of the reference line placed at the end of the strip)."""
    return ndi.binary_fill_holes(tissue) if prm.get("fill_holes") else tissue


def single_line(m, bg, prm=PRM):
    """one thin cut line: merge strands closer than 2 x merge_um, thin to a skeleton, trim spurs
    shorter than spur_um (not where the line meets the open surface), widen to 2 cells.
    keep_ends (2026-11-22, the user's flatten v5 notes): the trimming eats spur_um off every free end, so each
    call shortened every cut's deep end by ~24 um, and a cut passes through this 3-4 times (s5_F at 2.4 mm: 70 um
    short of the furrow base, the sharp EGL peaks in the flattened strips). With keep_ends every free end is trimmed
    by branch_um, which removes side branches up to that long, and each end of what is left gets back the longest
    path through its trimmed cells: the line's own ends keep their full length, a short fork at an end is still
    dropped. Used for the last passes over a finished cut (final_prm); the first passes keep the old trimming, which
    also clears short stray branches (tried everywhere, it kept them: cuts 1-3 mm longer, some through the ML)."""
    r = max(1, int(round(prm["merge_um"] / G)))
    sk = morphology.skeletonize(ndi.binary_dilation(m, iterations=r) & ~bg)
    keep = ndi.binary_dilation(bg, iterations=r + 1)
    nb = np.ones((3, 3), int)
    sk0 = sk.copy()
    ke = prm.get("keep_ends", False)
    for _ in range(int(round((prm["branch_um"] if ke else prm["spur_um"]) / G))):
        if ke:
            ends = skel_ends(sk) & ~keep
        else:
            deg = ndi.convolve(sk.astype(int), nb, mode="constant") - 1
            ends = sk & (deg <= 1) & ~keep
        if not ends.any():
            break
        sk &= ~ends
    if ke:
        removed = sk0 & ~sk
        if removed.any() and sk.any():
            H, W = sk.shape
            nbrs = [(dr, dc) for dr in (-1, 0, 1) for dc in (-1, 0, 1) if dr or dc]
            back = np.zeros_like(sk)
            for t in np.argwhere(skel_ends(sk) & ~keep):
                # breadth-first through the trimmed cells from this end; restore the path to the farthest one
                par, front, last = {}, [], None
                for dr, dc in nbrs:
                    q = (t[0] + dr, t[1] + dc)
                    if 0 <= q[0] < H and 0 <= q[1] < W and removed[q] and q not in par:
                        par[q] = None
                        front.append(q)
                while front:
                    last = front[-1]
                    nxt = []
                    for (r_, c_) in front:
                        for dr, dc in nbrs:
                            q = (r_ + dr, c_ + dc)
                            if 0 <= q[0] < H and 0 <= q[1] < W and removed[q] and q not in par:
                                par[q] = (r_, c_)
                                nxt.append(q)
                    if nxt:
                        last = nxt[0]
                    front = nxt
                while last is not None:
                    back[last] = True
                    last = par[last]
            sk |= back
    sk |= m & keep                                     # the original cells at the mouth reach the background
    return ndi.binary_dilation(sk, structure=np.ones((2, 2), bool))


def _resample(P, closed, step):
    if closed:
        if np.allclose(P[0], P[-1]):
            P = P[:-1]
        P = np.vstack([P, P[:1]])
    L = np.r_[0, np.cumsum(np.hypot(*np.diff(P, axis=0).T))]
    n = max(int(L[-1] / step), 2)
    t = np.linspace(0, L[-1], n, endpoint=not closed)
    return np.c_[np.interp(t, L, P[:, 0]), np.interp(t, L, P[:, 1])]


def reference(lab, tissue, prm=PRM, mid=None):
    """contours of the cut tissue (grid units, 2 um apart), inward normals, closed flags. mid: fissure
    midlines to cut along (default: found from the EGL labels by midline())."""
    if mid is None:
        mid = midline(lab, tissue, prm)
    cut = tissue & ~mid
    cut = morphology.remove_small_objects(cut, max_size=int(20000 / G ** 2))
    phi = ndi.distance_transform_edt(cut) - ndi.distance_transform_edt(~cut)
    phi = ndi.gaussian_filter(phi.astype(np.float32), prm["near_sigma_um"] / G)
    gy, gx = np.gradient(phi)
    out = []
    for C in measure.find_contours(np.pad(cut, 1).astype(np.float32), 0.5):
        C = C - 1.0
        closed = bool(np.allclose(C[0], C[-1]))
        if np.hypot(*np.diff(C, axis=0).T).sum() * G < prm["min_contour_um"]:
            continue
        Q = _resample(C, closed, STEP / G)
        Q = ndi.gaussian_filter1d(Q, prm["pos_sigma_um"] / STEP, axis=0, mode="wrap" if closed else "nearest")
        t = ndi.gaussian_filter1d(Q, prm["dir_sigma_um"] / STEP, axis=0, order=1, mode="wrap" if closed else "nearest")
        n = np.c_[-t[:, 1], t[:, 0]]
        n /= np.maximum(np.hypot(*n.T), 1e-6)[:, None]
        g = np.stack([ndi.map_coordinates(gy, Q.T, order=1), ndi.map_coordinates(gx, Q.T, order=1)], 1)
        if np.mean((n * g).sum(1)) < 0:                      # point into the tissue
            n = -n
        out.append(dict(Q=Q, n=n, closed=closed))
    return out, cut, mid


def streamlines(Q, n0, cut, prm=PRM):
    """column sample positions (N, Dn, 2) in grid units along the depth field, and validity."""
    phi = (ndi.distance_transform_edt(cut) - ndi.distance_transform_edt(~cut)).astype(np.float32)
    fields = []
    for sig in (prm["near_sigma_um"], prm["far_sigma_um"]):
        f = ndi.gaussian_filter(phi, sig / G)
        fields.append(np.gradient(f))

    def direction(x, prev, depth_um):
        a = min(depth_um / prm["blend_um"], 1.0)
        d = np.zeros_like(x)
        for w, (gy, gx) in ((1 - a, fields[0]), (a, fields[1])):
            if w == 0:
                continue
            g = np.stack([ndi.map_coordinates(gy, x.T, order=1, mode="nearest"),
                          ndi.map_coordinates(gx, x.T, order=1, mode="nearest")], 1)
            d += w * g / np.maximum(np.hypot(*g.T), 1e-9)[:, None]
        m = np.hypot(*d.T)
        d = d / np.maximum(m, 1e-9)[:, None]
        d[m < 1e-6] = prev[m < 1e-6]
        return d

    N, Dn = len(Q), int(prm["depth_um"] / STEP)
    h = STEP / G
    X = np.zeros((N, Dn, 2), np.float32)
    valid = np.zeros((N, Dn), bool)
    x = Q.astype(np.float64).copy()
    prev = n0.astype(np.float64)
    alive = np.ones(N, bool)
    top = np.zeros(N)
    hist = []
    lag = int(prm["stall_um"] / STEP)
    for k in range(Dn):
        step = h / 2 if k == 0 else h
        depth = (k + 0.5) * STEP
        d1 = direction(x, prev, depth)
        d2 = direction(x + d1 * step / 2, d1, depth)
        xn = x + d2 * step
        ph = ndi.map_coordinates(phi, xn.T, order=1, mode="constant", cval=-1.0)
        turned = (d2 * prev).sum(1) < 0
        stalled = (ph - hist[-lag]) < prm["stall_frac"] * lag * h if len(hist) >= lag else np.zeros(N, bool)
        alive &= ((ph > -0.5) | (k < 3)) & (ph >= top - prm["ridge_drop_um"] / G) & ~turned & ~(stalled & (k * STEP > 2 * prm["stall_um"]))
        top = np.maximum(top, ph)
        hist.append(ph)
        X[:, k] = xn
        valid[:, k] = alive
        x = np.where(alive[:, None], xn, x)
        prev = np.where(alive[:, None], d2, prev)
    return X, valid


def coverage(L, valid, prm=PRM):
    """columns with EGL at the top: within cover_um of the first sample that is not unassigned
    tissue or anomaly (a tear or bleach patch lining a fissure, meninges)."""
    k = int(prm["cover_um"] / STEP)
    skip = np.isin(L, (1, 6)) & valid
    lead = np.argmin(skip, axis=1)                         # first sample past the leading 1 / 6 cells
    lead = np.minimum(lead, int(prm["skip_um"] / STEP))
    idx = lead[:, None] + np.arange(k)[None, :]
    idx = np.minimum(idx, L.shape[1] - 1)
    rows = np.arange(L.shape[0])[:, None]
    c = ((L[rows, idx] == 2) & valid[rows, idx]).any(1)
    r = max(1, int(prm["cover_run_um"] / STEP / 2))
    c = ndi.binary_closing(np.pad(c, r, mode="edge"), np.ones(2 * r + 1, bool))[r:-r]
    c = ndi.binary_opening(np.pad(c, r, mode="edge"), np.ones(2 * r + 1, bool))[r:-r]
    return c


def bounds(L, valid, prm=PRM):
    """per column (um): EGL bottom, IGL top, IGL bottom (NaN when absent)."""
    N, Dn = L.shape
    isI = (L == 4) & valid
    isI = ndi.binary_closing(isI, np.ones((1, 4), bool)) & valid
    minr = int(prm["min_igl_um"] / STEP)
    eb, top, bot = (np.full(N, np.nan) for _ in range(3))
    for i in range(N):
        col = isI[i]
        if col.any():
            t = int(np.argmax(col))
            e = t + int(np.argmin(col[t:])) if not col[t:].all() else Dn
            if e - t >= minr:
                top[i], bot[i] = t * STEP, e * STEP
        eg = (L[i] == 2) & valid[i]
        if eg[:5].any():
            s0 = int(np.argmax(eg))
            e0 = s0 + int(np.argmin(eg[s0:])) if not eg[s0:].all() else Dn
            eb[i] = e0 * STEP
    return eb, top, bot


def _runs(m):
    m = np.r_[False, m, False]
    d = np.diff(m.astype(np.int8))
    return list(zip(np.flatnonzero(d == 1), np.flatnonzero(d == -1)))


def fit_line(v, ok, prm=PRM):
    """robust line through v (NaN = absent) inside each run of ok columns; returns fit, bridged."""
    N = len(v)
    fit = np.full(N, np.nan)
    bridged = np.zeros(N, bool)
    w = int(prm["med_win_um"] / STEP) | 1
    sig = prm["fit_sigma_um"] / STEP
    for a, b in _runs(ok):
        x = v[a:b]
        if np.isfinite(x).sum() < w * prm["med_frac"]:
            continue
        s = np.asarray(pd.Series(x).rolling(w, center=True, min_periods=max(3, int(w * prm["med_frac"]))).median())
        fin = np.isfinite(s)
        num = ndi.gaussian_filter1d(np.where(fin, s, 0.0), sig, mode="nearest")
        den = ndi.gaussian_filter1d(fin.astype(float), sig, mode="nearest")
        f = np.where(den > 0.25, num / np.maximum(den, 1e-6), np.nan)
        f[~fin] = np.nan
        for g0, g1 in _runs(~np.isfinite(f)):
            if g0 > 0 and g1 < len(f) and (g1 - g0) * STEP <= prm["bridge_um"]:
                f[g0:g1] = np.interp(np.arange(g0, g1), [g0 - 1, g1], [f[g0 - 1], f[g1]])
                bridged[a + g0:a + g1] = True
        fit[a:b] = f
    return fit, bridged


def _clean(m, prm=PRM):
    k = max(1, int(prm["min_flag_um"] / STEP))
    return morphology.remove_small_objects(m, max_size=k - 1) if m.any() else m


def igl_line(L, valid, ok, prm=PRM):
    eb, top, bot = bounds(L, valid, prm)
    t_ok = ok & valid[:, 3]
    ftop, br = fit_line(top, t_ok, prm)
    fbot, _ = fit_line(bot, t_ok, prm)
    present = np.isfinite(top)
    line = np.isfinite(ftop)
    ext = valid.sum(1) * STEP
    flags = dict(
        missing=_clean(line & ~present & t_ok & (ext >= ftop + 8), prm),
        off=_clean(present & line & ((np.abs(top - ftop) > prm["off_top_um"]) | (np.abs(bot - fbot) > prm["off_bot_um"])), prm),
        no_ml=_clean(present & np.isfinite(eb) & (top - eb <= prm["no_ml_um"]), prm))
    um = lambda m: int(m.sum() * STEP)
    segs = [(a, b) for a, b in _runs(line & ~br) if (b - a) * STEP >= 100]
    stats = dict(covered_um=um(t_ok), igl_present_um=um(present & t_ok), line_um=um(line), bridged_um=um(br),
                 segments=len(segs), short_columns_um=um(t_ok & ~present & (ext < 60)),
                 **{f"{k}_um": um(m) for k, m in flags.items()},
                 igl_top_med_um=round(float(np.nanmedian(top[t_ok])), 1) if np.isfinite(top[t_ok]).any() else None,
                 igl_thick_med_um=round(float(np.nanmedian((bot - top)[t_ok])), 1) if np.isfinite(top[t_ok]).any() else None)
    return dict(eb=eb, top=top, bot=bot, ftop=ftop, fbot=fbot, bridged=br, flags=flags, stats=stats)


def stains_reduced(sid):
    S, px = D.stains(sid)
    k = 3 if px > 0.6 else 4
    out = {name: downscale_local_mean(a.astype(np.float32), (k, k)) for name, a in S.items()}
    return out, k, next(iter(S.values())).shape


def reached(X, valid, ok, shape):
    hit = np.zeros(shape, bool)
    m = valid & ok[:, None]
    r = np.clip(np.round(X[..., 0][m]).astype(int), 0, shape[0] - 1)
    c = np.clip(np.round(X[..., 1][m]).astype(int), 0, shape[1] - 1)
    hit[r, c] = True
    return ndi.binary_dilation(hit)


def _columns(lab, tissue, mid, prm=PRM, full=False):
    refs, cut, mid = reference(lab, tissue, prm, mid=mid)
    Q = np.concatenate([r["Q"] for r in refs])
    n = np.concatenate([r["n"] for r in refs])
    X, valid = streamlines(Q, n, cut, prm)
    N, Dn = valid.shape
    L = ndi.map_coordinates(lab, X.reshape(-1, 2).T, order=0, mode="constant").reshape(N, Dn).astype(np.uint8)
    L[~valid] = 0
    if not full:
        return Q, valid, L
    ok = coverage(L, valid, prm)
    for a in np.cumsum([len(r["Q"]) for r in refs])[:-1]:
        ok[a - 1:a + 1] = False
    return dict(Q=Q, n=n, X=X, valid=valid, L=L, ok=ok)


def egl_depth(L, valid, gap=4):
    """per column (um): bottom of the EGL run starting in the first 20 um, across gaps of up to gap samples of
    unassigned / anomaly cells; NaN where no EGL starts there."""
    N, Dn = L.shape
    out = np.full(N, np.nan)
    for i in np.flatnonzero(((L[:, :10] == 2) & valid[:, :10]).any(1)):
        c, v = L[i], valid[i]
        k = int(np.argmax((c == 2) & v))
        last, g = k, 0
        while k < Dn and v[k]:
            if c[k] == 2:
                last, g = k, 0
            elif c[k] in (1, 6):
                g += 1
                if g > gap:
                    break
            else:
                break
            k += 1
        out[i] = (last + 1) * STEP
    return out


def egl_peaks(L, valid, ok, thr_um=40.0, win_um=200.0, sep_um=60.0):
    """sharp EGL peaks in the strips (the user's flatten v5 notes: a fissure cut that stops short of the furrow base
    leaves the columns at its end running down the rest of the fused EGL): columns on the EGL-covered line whose EGL
    runs thr_um or more deeper than its median over win_um around, one per sep_um. Returns (peaks [(index, excess
    um)], excess per column, the rolling median EGL depth)."""
    ed = egl_depth(L, valid)
    w = int(win_um / STEP) | 1
    base = pd.Series(ed).rolling(w, center=True, min_periods=20).median().to_numpy()
    exc = ed - base
    sep = int(sep_um / STEP / 2)
    out = []
    for i in np.flatnonzero(np.nan_to_num(exc, nan=-1) >= thr_um):
        a, b = max(0, i - sep), min(len(ed), i + sep + 1)
        if ok[i] and exc[i] >= np.nanmax(exc[a:b]) and (not out or i - out[-1][0] >= 2 * sep):
            out.append((int(i), float(exc[i])))
    return out, exc, base, ed


def furrow_walk(i, C, t_loc, ed_i, dcx, eglish, cortex, tissue, mid, prm=PRM):
    """from line point i (a sharp EGL peak), follow the fused EGL down to the furrow base: steps of one cell in the
    column's starting direction, turning at most pk_turn_deg per step and pk_max_turn_deg overall, through EGL-like
    cells, each step the one furthest from the cortex (smoothed, looking pk_look cells ahead: the middle of the two
    EGLs back to back, where the pia is). It stops where that distance falls below pk_stop_frac x its median along
    the walk (>= pk_stop_min_um; the two EGLs wrap round the base), where it gains less than pk_prog_frac cells of
    depth per step over pk_prog_k steps (running along the base), within 2 cells of another cut, or after the
    column's own EGL run + t_loc. A step counts as two EGLs back to back when, of pk_rays rays up to pk_side x t_loc
    long, one to its left and one to its right (within pk_side_deg of its normals) reach cortex before the background
    or a cut, the farther side at most pk_balance x the nearer (midway between two cortices, not along one), and the walk keeps moving into the tissue (its distance to the
    outside >= pk_into_frac x its length); it ends after pk_max_fail steps that are not. Returns (path cells,
    two-sided flags, why it stopped)."""
    H, W = tissue.shape
    Q, X, valid = C["Q"], C["X"], C["valid"]
    k = min(prm["pk_start_k"], int(valid[i].sum()) - 1)
    if k < 2:
        return [], [], "short column"
    d0 = X[i, k].astype(float) - X[i, 0].astype(float)
    d0 /= max(np.hypot(*d0), 1e-9)
    pos, d = Q[i].astype(float), d0.copy()
    bg = ~tissue
    rr, cc = np.indices(tissue.shape)
    far0 = np.hypot(rr - Q[i][0], cc - Q[i][1]) * G > prm["pk_start_um"]
    blk = (bg | mid) & far0 if prm["pk_start_um"] > 0 else (bg | mid)
    dmid = ndi.distance_transform_edt(~mid) if mid.any() else np.full(tissue.shape, np.inf)
    dout = ndi.distance_transform_edt(~(bg | mid)) * G
    dcs = ndi.gaussian_filter(dcx, 1.0)
    side = np.arange(1, int(prm["pk_side"] * t_loc / G) + 2)
    ang = np.linspace(0, 2 * np.pi, prm["pk_rays"], endpoint=False)
    U = np.stack([np.cos(ang), np.sin(ang)], 1)

    def two_sided(q, dq):
        pp = np.round(np.array(q, float)[None, None, :] + U[:, None, :] * side[None, :, None]).astype(int)
        pp[..., 0] = np.clip(pp[..., 0], 0, H - 1)
        pp[..., 1] = np.clip(pp[..., 1], 0, W - 1)
        hc, hb = cortex[pp[..., 0], pp[..., 1]], blk[pp[..., 0], pp[..., 1]]
        big = 10 ** 6
        fc = np.where(hc.any(1), np.argmax(hc, 1), big)
        fb = np.where(hb.any(1), np.argmax(hb, 1), big)
        hit = (fc < big) & (fc < fb)
        # cortex to the left and to the right of the walk (rays within pk_side_deg of its two normals), the nearer side
        # no more than pk_balance x closer than the other: the walk runs midway between two cortices, not along one
        rel = np.degrees(np.arctan2(U[:, 0] * dq[1] - U[:, 1] * dq[0], U[:, 0] * dq[0] + U[:, 1] * dq[1]))
        off = 90.0 - prm["pk_side_deg"]
        lm, rm = hit & (rel >= off) & (rel <= 180 - off), hit & (rel <= -off) & (rel >= off - 180)
        if not lm.any() or not rm.any():
            return False
        dl, dr = side[fc[lm]].min(), side[fc[rm]].min()
        return bool(max(dl, dr) <= prm["pk_balance"] * max(min(dl, dr), 2))
    rot = lambda v, a: np.array([v[0] * np.cos(a) - v[1] * np.sin(a), v[0] * np.sin(a) + v[1] * np.cos(a)])
    turns = [0.0] + [sg * np.radians(t) for t in (prm["pk_turn_deg"] / 2, prm["pk_turn_deg"]) for sg in (1, -1)]
    path, good, fails, why = [], [], 0, "length"
    stop = prm["pk_t_frac"] * t_loc
    for step in range(int((ed_i + t_loc) / G) + 2):
        best = None
        for a in turns:
            dd = rot(d, a)
            if np.degrees(np.arccos(np.clip((dd * d0).sum(), -1, 1))) > prm["pk_max_turn_deg"]:
                continue
            q = tuple(np.round(pos + dd).astype(int))
            if not (0 <= q[0] < H and 0 <= q[1] < W) or not tissue[q] or not eglish[q]:
                continue
            if step > 3 and dmid[q] <= 2:
                continue
            qa = tuple(np.clip(np.round(pos + prm["pk_look"] * dd).astype(int), 0, [H - 1, W - 1]))
            sc = dcs[q] + dcs[qa] + (0.5 if a == 0 else 0.0)
            if best is None or sc > best[0]:
                best = (sc, dd, q)
        if best is None:
            why = "no EGL ahead"
            break
        _, dd, q = best
        if len(path) >= 2:
            stop = max(prm["pk_stop_min_um"], prm["pk_stop_frac"] * float(np.median([dcx[x] for x in path])))
        if step >= 2 and dcx[q] < stop:
            why = "base"
            break
        kk = prm["pk_prog_k"]
        if len(path) >= kk and dout[q] - dout[path[-kk]] < prm["pk_prog_frac"] * kk * G:
            why = "along the base"
            break
        ts = bool(two_sided(q, dd) and dout[q] >= prm["pk_into_frac"] * (len(path) + 1) * G)
        d, pos = dd, pos + dd
        path.append(q)
        good.append(ts)
        fails = 0 if ts else fails + 1
        if fails > prm["pk_max_fail"]:
            why = "not two EGLs"
            break
    while good and not good[-1]:
        good.pop()
        path.pop()
    return path, good, why


def peak_cuts(lab, tissue, tis_ref, mid, egl_s, prm=PRM):
    """v7 (the user's flatten v6 / v7 notes: furrows the cuts still stop short of, s5_E at 10.6 mm, s5_F at 6 and 10):
    the error in the linearization finds them. Flatten with the current cuts; every sharp EGL peak (egl_peaks with
    pk_thr_um) is a column running down the rest of a fused EGL, so a cut is walked from its line point down the
    middle of that EGL (furrow_walk), on the tissue as the cuts see it (cut_tissue), with the cortex = ML / IGL /
    white matter / anomaly cells the stain does not call EGL (specks of <= pk_speck_cells cells inside the EGL
    ignored) and EGL-like = EGL label, unassigned, or stain EGL. A walk is kept when it is at least
    max(pk_min_len_um, the local EGL depth, half the peak's excess) long, pk_pass_frac of its steps look like two EGLs
    back to back, pk_egl_frac of it runs over EGL / unassigned labels (not stain-only EGL inside the ML) and it is
    straight (end-to-end >= pk_straight x its length; not along an EGL / IGL edge). Returns (new cut cells, report)."""
    from skimage import morphology as mo
    C = _columns(lab, tis_ref, mid, prm, full=True)
    pk, exc, base, ed = egl_peaks(C["L"], C["valid"], C["ok"], prm["pk_thr_um"], prm["pk_win_um"], prm["pk_sep_um"])
    cortex = tissue & np.isin(lab, (3, 4, 5, 6)) & ~egl_s
    cortex = mo.remove_small_objects(cortex, max_size=prm["pk_speck_cells"], connectivity=2)
    eglish = tissue & (np.isin(lab, (1, 2)) | egl_s)
    dcx = ndi.distance_transform_edt(~cortex) * G
    add = np.zeros_like(mid)
    rep = []
    H, W = mid.shape
    for i, e in pk:
        if not np.isfinite(base[i]):
            continue
        path, good, why = furrow_walk(i, C, float(base[i]), float(ed[i]), dcx, eglish, cortex, tissue, mid, prm)
        L_um = len(path) * G
        keep = bool(path) and L_um >= max(prm["pk_min_len_um"], base[i], 0.5 * e) and np.mean(good) >= prm["pk_pass_frac"]
        if keep:
            pa = np.array(path)
            ef = float(np.isin(lab[pa[:, 0], pa[:, 1]], (1, 2)).mean())
            straight = float(np.hypot(*(pa[-1] - C["Q"][i]))) / len(path)
            if ef < prm["pk_egl_frac"] or straight < prm["pk_straight"]:
                keep, why = False, f"{why}; EGL labels {ef:.2f}, straightness {straight:.2f}"
        rep.append(dict(i=int(i), excess_um=round(e), egl_um=round(float(ed[i])), t_um=round(float(base[i])), walk_um=L_um,
                        end=why, kept=bool(keep), at=[round(float(v), 1) for v in C["Q"][i]]))
        if keep:
            for q in path:
                add[q] = True
            for s_ in np.arange(0, 3.01, 0.5):                   # joined to the outside (or the cut) the line runs on
                q = tuple(np.round(C["Q"][i] - s_ * C["n"][i]).astype(int))
                if 0 <= q[0] < H and 0 <= q[1] < W:
                    add[q] = True
    add = ndi.binary_dilation(add, structure=np.ones((2, 2), bool)) & tissue
    if prm.get("pk_check") and add.any():
        add = _check_peak_cuts(lab, tis_ref, mid, add, C, rep, prm)
    return add, rep


def _short_um(C, near):
    """line length (um) among the points in near whose columns reach no IGL or white matter."""
    reach = (np.isin(C["L"], (4, 5)) & C["valid"]).any(1)
    return float((near & ~reach).sum() * STEP)


def _check_peak_cuts(lab, tis_ref, mid, add, C0, rep, prm=PRM):
    """v7 check of each kept walk by what it does to the strip: flattened with that cut alone, the line within
    pk_check_um of the walk must not gain more than pk_short_um of columns that reach no IGL or white matter (s5_C
    at 1.2 mm: a walk between the ML of a folium's thick EGL cap and the same folium below sliced the cap off, its
    columns ending in the EGL). Dropped walks get kept=False in the report."""
    from scipy.spatial import cKDTree
    comp, nc = ndi.label(add, structure=np.ones((3, 3), bool))
    out = np.zeros_like(add)
    r0 = prm["pk_check_um"] / G
    t0 = cKDTree(C0["Q"])
    for k in range(1, nc + 1):
        m = comp == k
        cells = np.argwhere(m)
        C1 = _columns(lab, tis_ref, mid | m, prm, full=True)
        tc_ = cKDTree(cells)
        near1 = tc_.query(C1["Q"])[0] <= r0
        near0 = tc_.query(C0["Q"])[0] <= r0
        gain = _short_um(C1, near1) - _short_um(C0, near0)
        hit = [r for r in rep if r["kept"] and tc_.query(r["at"])[0] <= 3]
        ok = gain <= prm["pk_short_um"]
        for r in hit:
            r["short_gain_um"] = round(gain)
            if not ok:
                r["kept"], r["end"] = False, r["end"] + f"; +{gain:.0f} um of columns reaching no IGL"
        if ok:
            out |= m
    return out


def prune_cuts(lab, tissue, mid, prm=PRM, egl_s=None):
    """drop the stretches of a fissure cut that have no folium on one side: the line points beside the cut
    whose columns end within dead_short_um with no ML, IGL or white matter (a strip of EGL between the cut and
    the outside, where the stain-based cut ran along thick EGL at an open surface; s2_2, s3_2) or (dead_noegl, v8)
    that have no EGL at their top (a cut through a folium's ML, s5_B), dead_frac or
    more of them within dead_win_um along the cut (one dead side = 0.5). Columns into a folium tip that has ML
    but no IGL in the section keep their cut. A deeper part of the cut left without a way out is joined to the
    background by the shortest path through EGL, unassigned or anomaly cells (ML 20 x dearer, IGL and white
    matter closed)."""
    if egl_s is not None and prm.get("dead_stain"):
        C = _columns(lab, tissue, mid, prm, full=True)
        Q, valid, L = C["Q"], C["valid"], C["L"]
        N, Dn = valid.shape
        Es = ndi.map_coordinates(egl_s.astype(np.uint8), C["X"].reshape(-1, 2).T, order=0, mode="constant").reshape(N, Dn).astype(bool) & valid
    else:
        Q, valid, L = _columns(lab, tissue, mid, prm)
        Es = np.zeros_like(valid)
    gh, gw = tissue.shape
    dead = ~(np.isin(L, (3, 4, 5)) & valid).any(1) & (valid.sum(1) * STEP < prm["dead_short_um"])
    if prm.get("dead_noegl"):
        # v8 (the user's flatten v5 notes, s5_B at 5.4 mm): a side of the cut whose columns have no EGL at the top
        # (as in coverage: within cover_um past leading unassigned / anomaly cells) is not a fissure wall; a cut with
        # ML or white matter right beside it runs through a folium. The repaired labels already have EGL added
        # wherever the stain shows it beside a cut.
        k = int(prm["cover_um"] / STEP)
        skip = np.isin(L, (1, 6)) & valid
        lead = np.minimum(np.argmin(skip, axis=1), int(prm["skip_um"] / STEP))
        idx = np.minimum(lead[:, None] + np.arange(k)[None, :], L.shape[1] - 1)
        rows = np.arange(L.shape[0])[:, None]
        # v7: or EGL in the stain (s4_A at 7.5 mm: the fused EGL at the furrow base is labelled white matter)
        top_egl = ((L[rows, idx] == 2) | (Es[rows, idx] & np.isin(L[rows, idx], (1, 5, 6)))) & valid[rows, idx]
        dead |= valid[:, 0] & ~top_egl.any(1)
    qi = np.clip(np.round(Q).astype(int), 0, [gh - 1, gw - 1])
    sliver, merged = np.zeros_like(mid), np.zeros_like(mid)
    if prm.get("drop_sliver"):
        # v6 (the user's flatten v6 notes, s5_E at 3 mm): a piece of the cut tissue under sliver_um2 with under
        # sliver_cortex_um2 of ML, IGL or white matter is a sliver, which reference() discards. A sliver whose edge is
        # sliver_open or more background is a strip sliced off an open surface (outer EGL beside a fissure the pia
        # trim opened): the cut beside it separates no folia and goes. A sliver closed in by cuts lies between two
        # parallel cuts down one fissure (the pia strip between two EGLs; s1_2 at 6, s5_C at 11): they become one cut
        # down its middle.
        rest = tissue & ~mid
        comp, nc = ndi.label(rest)
        k8 = np.ones((3, 3), bool)
        if nc:
            idx = np.arange(1, nc + 1)
            size = np.asarray(ndi.sum(rest, comp, index=idx)) * G * G
            cx = np.asarray(ndi.sum(np.isin(lab, (3, 4, 5)), comp, index=idx)) * G * G
            small = (size < prm["sliver_um2"]) & (cx < prm["sliver_cortex_um2"])
            for k, sl_ in enumerate(ndi.find_objects(comp), start=1):
                if not small[k - 1]:
                    continue
                sl_ = tuple(slice(max(a.start - 4, 0), a.stop + 4) for a in sl_)
                m = comp[sl_] == k
                ring = ndi.binary_dilation(m, k8) & ~m
                border = mid[sl_] & ndi.binary_dilation(m, k8, iterations=2)
                sliver[sl_] |= border
                if (ring & ~tissue[sl_]).sum() < prm["sliver_open"] * ring.sum():
                    merged[sl_] |= morphology.skeletonize(m | border) & tissue[sl_]
    dist, (ir, ic) = ndi.distance_transform_edt(~mid, return_indices=True)
    on = dist[qi[:, 0], qi[:, 1]] <= 1.5
    tot, bad = np.zeros(mid.shape), np.zeros(mid.shape)
    r, c = ir[qi[:, 0], qi[:, 1]][on], ic[qi[:, 0], qi[:, 1]][on]
    np.add.at(tot, (r, c), 1.0)
    np.add.at(bad, (r, c), dead[on].astype(float))
    k = int(prm["dead_win_um"] / G) | 1
    num, den = ndi.uniform_filter(bad * mid, k), ndi.uniform_filter(tot * mid, k)
    drop = mid & ((np.where(den > 0, num / np.maximum(den, 1e-9), 0.0) >= prm["dead_frac"]) | sliver)
    keep = (mid & ~drop) | merged
    drop &= ~merged
    bg = ~tissue
    cost = np.where(np.isin(lab, (1, 2, 6)), 1.0, np.where(lab == 3, 20.0, np.inf))
    cost[keep] = 0.05
    import egl_repair as ER
    way = ER._connect(keep, cost, bg)
    if prm.get("smooth_way") and way.any():
        way, gone = _smooth_ways(keep, way, cost, drop, bg, prm)
        keep &= ~gone
        drop |= gone
    keep = single_line(keep | way, bg, final_prm(prm)) & tissue
    keep = ER._prune_shallow(keep, tissue, prm)
    if prm.get("join_stranded"):
        keep = join_stranded(keep, tissue, lab)
    return keep, drop


def _axis(cells):
    """main direction (unit vector) of a set of cells."""
    c = cells - cells.mean(0)
    w, v = np.linalg.eigh(c.T @ c)
    return v[:, -1]


def _smooth_ways(keep, way, cost, drop, bg, prm=PRM):
    """v9 (the user's note on the walkthrough, 2026-11-24: the line turns into the ML; the EGL almost never makes sharp
    turns): a deeper piece of a cut left without a way out is joined to the background by the cheapest path, which can
    leave the piece at a right angle straight through the EGL (s2_2 at 11.2, s3_2 at 11.0). Where that path turns
    from the piece by way_turn_deg or more (main directions within way_um of where they meet), the piece goes out along
    its own dropped course instead (the dropped cells cost as little as the kept ones)."""
    k8 = np.ones((3, 3), bool)
    cost_b = cost.copy()
    cost_b[drop] = np.minimum(cost_b[drop], 0.05)
    import egl_repair as ER
    way_b = ER._connect(keep, cost_b, bg)
    comp, nc = ndi.label(keep, structure=k8)
    wcomp, _ = ndi.label(way, structure=k8)
    out = way.copy()
    dropped_pieces = np.zeros_like(way)
    r = prm["way_um"] / G
    for k in range(1, nc + 1):
        piece = comp == k
        near = ndi.binary_dilation(piece, k8, iterations=2)
        ids = np.unique(wcomp[near & way])
        ids = ids[ids > 0]
        for w in ids:
            path = (wcomp == w) & ~piece
            touch = np.argwhere(path & near)
            if not len(touch):
                continue
            j = touch.mean(0)
            pc = np.argwhere(piece)
            pc = pc[np.hypot(*(pc - j).T) <= r]
            wc = np.argwhere(path)
            wc = wc[np.hypot(*(wc - j).T) <= r]
            if len(pc) < 3 or len(wc) < 3:
                continue
            ang = np.degrees(np.arccos(min(abs(float(_axis(pc) @ _axis(wc))), 1.0)))
            if ang < prm["way_turn_deg"]:
                continue
            alt = way_b & ndi.binary_dilation(piece, k8, iterations=int(prm["way_reach_um"] / G))
            alt_c, _ = ndi.label(alt | piece, structure=k8)
            lab_p = np.unique(alt_c[piece])
            alt = np.isin(alt_c, lab_p[lab_p > 0]) & ~piece & way_b
            if alt.any() and (ndi.binary_dilation(alt, k8) & bg).any() and alt.sum() * G <= prm["way_alt_max_um"]:
                if prm.get("way_debug"):
                    print("smooth way", tuple(np.round(j).astype(int)), "turn", round(ang), "shortcut", int(path.sum() * G), "um -> own course", int(alt.sum() * G), "um")
                out &= ~(wcomp == w)
                out |= alt
            elif prm.get("way_drop"):
                # its own course runs on along thick EGL (dead columns, the reason it was dropped): the piece goes too
                if prm.get("way_debug"):
                    print("dropped piece", tuple(np.round(j).astype(int)), "turn", round(ang), "shortcut", int(path.sum() * G), "own course", int(alt.sum() * G), "piece", int(piece.sum() * G))
                out &= ~(wcomp == w)
                dropped_pieces |= piece
            elif prm.get("way_debug"):
                print("kept shortcut", tuple(np.round(j).astype(int)), "turn", round(ang), "shortcut", int(path.sum() * G), "own course", int(alt.sum() * G))
    return out, dropped_pieces


def join_stranded(mid, tissue, lab):
    """v11 (the user's notes on flattening Version 11, s5_E: "a discontinuity, this should come right after the 13"):
    every cut reaches the open surface. A cut piece left inside the tissue is a closed slit, and the outline round it is
    a loop of its own that the walk cannot reach, placed at the end of the strip. In s5_E prune_cuts left one: the way
    out it shared with a short piece near the surface went with that piece (_smooth_ways drops a way whole), stranding
    the 0.5 mm of cut behind it 24 um from the rest of its fissure. Such a piece is joined to the background again by
    the cheapest path (as prune_cuts: EGL, unassigned or anomaly cells 1, ML 20, IGL and white matter closed, other
    cuts almost free), widened as single_line widens; a piece with no such path is dropped."""
    import egl_repair as ER
    k8 = np.ones((3, 3), bool)
    bg = ~tissue
    touch = ndi.binary_dilation(bg, k8)
    comp, n = ndi.label(mid, structure=k8)
    if not n:
        return mid
    hit = np.asarray(ndi.maximum(touch, comp, index=np.arange(1, n + 1))).astype(bool)
    if hit.all():
        return mid
    cost = np.where(np.isin(lab, (1, 2, 6)), 1.0, np.where(lab == 3, 20.0, np.inf))
    cost[mid] = 0.05
    way = ER._connect(mid, cost, bg)
    out = mid | (ndi.binary_dilation(way, structure=np.ones((2, 2), bool)) & tissue)
    comp, n = ndi.label(out, structure=k8)
    if n:
        hit = np.asarray(ndi.maximum(touch, comp, index=np.arange(1, n + 1))).astype(bool)
        out &= np.r_[True, hit][comp] | ~out
    return out


OUT2 = ROOT / "outputs/flatten_v11"
PREV2 = ROOT / "outputs/flatten_v10"          # the first walk (flatten_v9: before the walk, flatten_v8: before the shape check)
SETS = {False: ("matlab", "current"), True: ("before", "repaired")}   # label sets sampled; the last one is the main one


def flatten(short, prm=PRM, save=True, repaired=False):
    """repaired=False: the v1 flattening (reference from the EGL labels; MATLAB and current labels) -> OUT.
    repaired=True: the reference line from the final EGL-repair cuts, less the stretches with no folium on one
    side (prune_cuts), and the repaired labels (egl_repair.repaired); the labels before the repair are sampled on
    the same columns for comparison -> OUT2."""
    t0 = time.time()
    sid = sid_of(short)
    auto, best, tissue = labels(short)
    if repaired:
        import egl_repair as ER
        lab, tis = ER.repaired(short)
        assert (tis == tissue).all()
        tis_p = ER.pia_tissue(short) if prm.get("use_pia", True) else tissue   # v6: inside the pia (arachnoid trimmed)
        ZR = np.load(ER.OUT / f"{short}.npz")
        mid_r, egl_s = ZR["mid"], ZR["egl_stain"]                   # the final EGL-repair fissure cuts, stain EGL
        tc = cut_tissue(tis_p, prm)
        mid_p, dropped = prune_cuts(lab, tc, mid_r, prm, egl_s=egl_s)
        pk_add, pk_rep = np.zeros_like(mid_p), []
        if prm.get("peak_cuts"):                                     # v7: on down to the furrow base at EGL peaks
            pk_add, pk_rep = peak_cuts(lab, tc, tis_p, mid_p, egl_s, prm)
            mid_p = mid_p | pk_add
        sh_fill, sh_trim, sh_rep = np.zeros_like(mid_p), np.zeros_like(mid_p), None
        if prm.get("shape_check"):                                   # v9: no sharp turns of the line off the EGL
            import shape_check as SC
            lab, tis_p, sh_fill, sh_trim, sh_rep = SC.check(short, lab, tissue, tis_p, mid_p)
            if prm.get("reprune") and sh_trim.any():
                # a cut beside meninges the check trimmed can be left with no folium on that side (a strip of EGL
                # between the cut and the trimmed space, s2_4 at 16.8): prune again, then check the line again
                mid_p, dropped2 = prune_cuts(lab, cut_tissue(tis_p, prm), mid_p, prm, egl_s=egl_s)
                dropped = dropped | dropped2
                lab, tis_p, f2, t2, rep2 = SC.check(short, lab, tissue, tis_p, mid_p)
                sh_fill, sh_trim = (sh_fill | f2) & tis_p, (sh_trim | t2) & ~tis_p
                sh_rep = dict(sh_rep, after=rep2["after"], fixes=sh_rep["fixes"] + [dict(f, round=f["round"] + 10) for f in rep2["fixes"]],
                              filled_um2=int(sh_fill.sum() * G * G), trimmed_um2=int(sh_trim.sum() * G * G),
                              open_after=rep2["open_after"], cut_after=rep2["cut_after"])
        wk_rep, tip_add = None, np.zeros_like(mid_p)
        tip_band = np.zeros_like(mid_p)
        if prm.get("walk"):                                          # v10: the line is one walk (walk.py)
            import walk as WK
            if WK.PRM["tip"]:                                        # cuts ending inside a fused EGL go on to the base
                tip_add, tip_rep = WK.tip_cuts(lab, cut_tissue(tis_p, prm), tis_p, mid_p, egl_s, fprm=prm)
                mid_p = mid_p | tip_add
                if short != WK.GOLD:                                 # the fused EGL it runs through is EGL
                    tip_band = WK.band(lab, tis_p, tip_add, tip_rep, egl_s) & ~mid_p
                    lab = lab.copy()
                    lab[tip_band] = 2
            if prm.get("join_stranded"):                             # v11: no closed slits (s5_E)
                mid_p = join_stranded(mid_p, cut_tissue(tis_p, prm), lab)
            refs, cut, cut_old, wk_rep = WK.reference(short, lab, tis_p, mid_p, egl_s, fprm=prm)
            lab, wfill, wtrim = WK.labels_after(short, lab, cut_old, cut, egl_s)
            tis_p = (tis_p | wfill) & ~wtrim
            sh_fill, sh_trim = (sh_fill | wfill) & tis_p, (sh_trim | wtrim) & ~tis_p
            wk_rep["tips"] = [r for r in tip_rep if r["kept"]] if WK.PRM["tip"] else []
            wk_rep["tip_walks"] = len(tip_rep) if WK.PRM["tip"] else 0
            wk_rep["tip_band_um2"] = int(tip_band.sum() * G * G)
            mid = mid_p
        else:
            refs, cut, mid = reference(lab, tis_p, prm, mid=mid_p)
        grids, out = {"before": best, "repaired": lab}, OUT2
    else:
        refs, cut, mid = reference(best, tissue, prm)
        grids, out = {"matlab": auto, "current": best}, OUT
    names = SETS[repaired]
    Q = np.concatenate([r["Q"] for r in refs])
    n = np.concatenate([r["n"] for r in refs])
    piece = np.concatenate([np.full(len(r["Q"]), i) for i, r in enumerate(refs)])
    X, valid = streamlines(Q, n, cut, prm)
    N, Dn = valid.shape
    cap_cells = np.zeros(cut.shape, bool)
    if repaired and prm.get("walk") and WK.PRM["cap"] and short != WK.GOLD:   # v10: the line is on the EGL
        cap_cells = WK.cap(X, valid, lab, egl_s)
        lab = lab.copy()
        lab[cap_cells] = 2
        grids["repaired"] = lab
        wk_rep["cap_um2"] = int(cap_cells.sum() * G * G)
    main = grids[names[-1]]
    flat = X.reshape(-1, 2).T
    L = {}
    for nm in names:
        L[nm] = ndi.map_coordinates(grids[nm], flat, order=0, mode="constant").reshape(N, Dn).astype(np.uint8)
        L[nm][~valid] = 0
    ok = coverage(L[names[-1]], valid, prm)
    starts = np.r_[0, np.cumsum([len(r["Q"]) for r in refs])]
    for a in starts[1:-1]:
        ok[a - 1:a + 1] = False
    lines = {nm: igl_line(L[nm], valid, ok, prm) for nm in names}
    hit = reached(X, valid, ok, tissue.shape)
    E = main == 2
    dcut = ndi.distance_transform_edt(cut) * G
    diag = dict(midline_um2=int(mid.sum() * G * G), egl_deep_before=round(float((ndi.distance_transform_edt(tissue)[E] * G > 60).mean()), 3),
                egl_deep_after=round(float((dcut[E] > 60).mean()), 3), pieces=len(refs), length_um=int(N * STEP),
                covered_frac=round(float(ok.mean()), 3),
                **{f"{nm}_reached": round(float(hit[main == c].mean()), 3) for c, nm in ((2, "egl"), (3, "ml"), (4, "igl"), (5, "wm"))})
    if repaired:
        diag["cut_dropped_um"] = int(dropped.sum() * G)
        if sh_rep is not None:
            diag["shape"] = {k: v for k, v in sh_rep.items() if k != "fixes"} | dict(fixes=[{k: v for k, v in f.items() if k != "labels"} | dict(labels=f["labels"]) for f in sh_rep["fixes"]])
        diag["peak_cuts"] = [r for r in pk_rep if r["kept"]]
        if wk_rep is not None:
            diag["walk"] = wk_rep
        diag["peak_walks"] = len(pk_rep)
    res = dict(section=short, group="control" if short in D.controls() else "cKO", sets=list(names), **diag,
               **{nm: lines[nm]["stats"] for nm in names})
    if save:
        R, k, (Hn, Wn) = stains_reduced(sid)
        gh, gw = tissue.shape
        red = (np.stack([flat[0] * (Hn - 1) / (gh - 1), flat[1] * (Wn - 1) / (gw - 1)]) + 0.5) / k - 0.5
        strips = {name: np.where(valid, ndi.map_coordinates(a, red, order=1, mode="nearest").reshape(N, Dn), 0) for name, a in R.items()}
        out.mkdir(parents=True, exist_ok=True)
        lab_keys = dict(Lm=L["matlab"], Lb=L["current"]) if not repaired else dict({f"L_{nm}": L[nm] for nm in names}, peak_add=pk_add,
                                                                                  lab_grid=lab.astype(np.uint8), pia=tis_p, shape_fill=sh_fill, shape_trim=sh_trim,
                                                                                  tip_add=tip_add, tip_band=tip_band, cap=cap_cells)
        np.savez_compressed(out / f"{short}.npz", Q=Q.astype(np.float32), n=n.astype(np.float32), piece=piece, valid=valid, ok=ok,
                            X=X.astype(np.float32), mid=mid, cut=cut, **lab_keys,
                            **{f"S_{k2}": np.round(v).astype(np.uint16) for k2, v in strips.items()},
                            **{f"{ls}_{k2}": v for ls, d in lines.items() for k2, v in d.items() if isinstance(v, np.ndarray)},
                            **{f"{ls}_flag_{k2}": v for ls, d in lines.items() for k2, v in d["flags"].items()})
    res["seconds"] = round(time.time() - t0, 1)
    print(json.dumps(res), flush=True)
    return res


def all_shorts():
    return sorted(GA.short(s) for s in GA.section_list())


if __name__ == "__main__":
    rep = "--repaired" in sys.argv
    args = [a for a in sys.argv[1:] if a != "--repaired"]
    shorts = args or all_shorts()
    rows = [flatten(s, repaired=rep) for s in shorts]
    if not args:
        (OUT2 if rep else OUT).joinpath("summary.json").write_text(json.dumps(dict(prm=PRM, repaired=rep, sections=rows), indent=1))
