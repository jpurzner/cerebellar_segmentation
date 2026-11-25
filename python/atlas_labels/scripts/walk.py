"""
The walk (2026-11-24, the user's notes on flattening Version 10: "establish a mandatory movement of the flattening.
Start on the left side of the cut surface and move clockwise following the boundary without jumping. It will help
enforce our shape restrictions. We need to fix these small jogs.").

The reference line is built as one walk round the tissue:

1. Start: the cut surface is the longest stretch of the outline (cs_min_um or more, gaps up to cs_gap_um closed) with
   no EGL (label or stain) within cs_egl_um inside and white matter within cs_wm_um: the white matter cut where the
   cerebellum was taken off the brainstem. The walk starts at its left end, the end where a clockwise walk (as the
   section is shown, image rows down) leaves it, and goes clockwise, so the cut surface comes last. Sections with no
   such stretch (the closed 10x controls) start where the white matter core is nearest the outline in the lowest
   fb_low_frac of the tissue.
2. Mandatory movement: the walk moves along the outline 2 um at a time, in order, and never jumps. Each step may turn
   only so far: the turn over +-turn_win_um stays under turn_deg (shape_check), except at a furrow base (a hairpin at
   the end of a cut) and at a fissure mouth. Where the outline turns sharper (a jog where a cut steps sideways, a
   crack, a flap), the walk cuts the corner with a smooth arc that leaves and rejoins the outline along its own
   direction, its points as far apart as the walk's (v2, after the user's notes on Version 11: v1 stiffened the line
   instead, and the stiffness jumping at the ends of each stiff stretch tore the walk there, steps of 5-14 um, at
   times across a cut). The arc is the shortest that removes the turn, moves the line at most disp_um, sweeps over
   no cut and takes in no furrow base or mouth; a turn no such arc removes is left and reported with the reason.
3. The tissue the columns run through is redrawn to the walk (the outline polygon, less the cuts): cells inside the
   walk that were outside are filled, tissue outside it is trimmed. Labels as in shape_check: trimmed cells become
   unassigned, filled cells EGL where the stain shows EGL, else unassigned (gold s1_2: only background cells change).
4. The walk turns round only at a furrow base. A cut that ends inside a fused EGL (the user's v10 marks: s5_E at
   3.5 mm, s2_2 at 7.8 mm; the two EGLs still run on back to back) is continued from its deep tip down the middle of
   the fused EGL (as flatten.furrow_walk: steps of one cell, at most 30 deg per step and tip_max_turn overall, through
   EGL-like cells, each the one furthest from the cortex) until the two EGLs part (distance to the cortex below
   tip_stop_frac of its median, at least tip_stop_min_um). Kept when it runs at least max(tip_min_um, the EGL
   thickness beside the cut's end), tip_pass of its steps have cortex on both sides at similar distances and it is
   straight; then checked by what it does to the strip (no more columns reaching no IGL, as for the peak cuts).
   The fused EGL it runs through (within one EGL thickness) becomes EGL where the stain shows EGL (band; s2_2 at 7.8 mm
   it was labelled white matter).
5. The line is on the EGL: a thin cap (at most cap_um) of ML / unassigned / anomaly cells between the line and the
   EGL that the stain calls EGL is relabelled EGL (s1_4 at 10.3 mm: ML on top of the EGL along the far wall of a cut).
"""
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi
from matplotlib.path import Path as MPath
from scipy.spatial import cKDTree
from skimage import measure, morphology
from skimage.draw import polygon as draw_polygon

sys.path.insert(0, str(Path(__file__).parent))
import flatten as F  # noqa: E402
import shape_check as SC  # noqa: E402

G, STEP = F.G, F.STEP
GOLD = "s1_2"
PRM = dict(cs_egl_um=40.0, cs_wm_um=60.0, cs_gap_um=100.0, cs_min_um=200.0, fb_low_frac=0.75,
           disp_um=12.0, arc_reach=60, arc_tan=4, arc_mid_tol=0.75, arc_rounds=300,
           tip=True, tip_max_um=300.0, tip_max_turn=60.0, tip_stop_frac=0.5, tip_stop_min_um=12.0, tip_min_um=32.0,
           tip_pass=0.7, tip_straight=0.75, tip_side=2.5, tip_balance=2.5, tip_rays=24, tip_side_deg=45.0, tip_back_cells=10,
           tip_check=True,
           cap=True, cap_um=20.0)


# --------------------------------------------------------------------------------------------------- the outline
def outline(cut):
    """the closed contours of the cut tissue (grid units, 2 um apart), longest first; each clockwise as shown."""
    out = []
    for C in measure.find_contours(np.pad(cut, 1).astype(np.float32), 0.5):
        C = C - 1.0
        if not np.allclose(C[0], C[-1]):
            continue
        L = np.hypot(*np.diff(C, axis=0).T).sum()
        if L * G < F.PRM["min_contour_um"]:
            continue
        P = F._resample(C, True, STEP / G)
        y, x = P[:, 0], P[:, 1]
        if 0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y) < 0:   # clockwise on screen (rows down)
            P = P[::-1].copy()
        out.append((L, P))
    return [P for _, P in sorted(out, key=lambda t: -t[0])]


def normals(P, cut):
    """inward normals as flatten.reference makes them."""
    phi = ndi.distance_transform_edt(cut) - ndi.distance_transform_edt(~cut)
    phi = ndi.gaussian_filter(phi.astype(np.float32), F.PRM["near_sigma_um"] / G)
    gy, gx = np.gradient(phi)
    t = ndi.gaussian_filter1d(P, F.PRM["dir_sigma_um"] / STEP, axis=0, order=1, mode="wrap")
    n = np.c_[-t[:, 1], t[:, 0]]
    n /= np.maximum(np.hypot(*n.T), 1e-6)[:, None]
    g = np.stack([ndi.map_coordinates(gy, P.T, order=1), ndi.map_coordinates(gx, P.T, order=1)], 1)
    if np.mean((n * g).sum(1)) < 0:
        n = -n
    return n


def info(P, cut, mid):
    """turns and cut contact of a closed line, as shape_check.line."""
    n = normals(P, cut)
    piece = np.zeros(len(P), int)
    tu = SC.turning(P, n, piece, int(round(SC.PRM["turn_win_um"] / STEP)))
    dm = ndi.distance_transform_edt(~mid) if mid.any() else np.full(mid.shape, 1e9)
    qi = np.clip(np.round(P).astype(int), 0, np.array(mid.shape) - 1)
    oncut = dm[qi[:, 0], qi[:, 1]] <= SC.PRM["oncut_cells"]
    return dict(Q=P, n=n, piece=piece, tu=tu, oncut=oncut)


# ---------------------------------------------------------------------------------------------------- the start
def ray_hits(P, n, lab, tis, codes, um, extra=None):
    """line points whose inward ray meets one of the label codes (or an extra mask) within um, inside the tissue."""
    H, W = lab.shape
    hit = np.zeros(len(P), bool)
    alive = np.ones(len(P), bool)
    for s in range(1, int(um / G) + 1):
        p = np.round(P + n * s).astype(int)
        inb = (p[:, 0] >= 0) & (p[:, 0] < H) & (p[:, 1] >= 0) & (p[:, 1] < W)
        p = np.clip(p, 0, [H - 1, W - 1])
        alive &= inb & (tis[p[:, 0], p[:, 1]] | (s <= 2))
        v = np.isin(lab[p[:, 0], p[:, 1]], codes)
        if extra is not None:
            v |= extra[p[:, 0], p[:, 1]]
        hit |= alive & v
    return hit


def runs_circ(m):
    """(start, length) of the True runs of a circular mask."""
    nn = len(m)
    if m.all():
        return [(0, nn)]
    if not m.any():
        return []
    k0 = int(np.argmin(m))
    mm = np.roll(m, -k0)
    out, i = [], 0
    while i < nn:
        if mm[i]:
            j = i
            while j < nn and mm[j]:
                j += 1
            out.append(((i + k0) % nn, j - i))
            i = j
        else:
            i += 1
    return out


def start(P, n, lab, tis, egl_s, prm=PRM):
    """index of the walk's start on the clockwise outline P, and how it was found."""
    noegl = ~ray_hits(P, n, lab, tis, (2,), prm["cs_egl_um"], egl_s)
    wm = ray_hits(P, n, lab, tis, (5,), prm["cs_wm_um"])
    cs = noegl & wm
    g = int(prm["cs_gap_um"] / STEP)
    for a, L in runs_circ(~cs):
        if L <= g:
            cs[(a + np.arange(L)) % len(cs)] = True
    rr = sorted([r for r in runs_circ(cs) if r[1] * STEP >= prm["cs_min_um"]], key=lambda r: -r[1])
    if rr:
        a, L = rr[0]
        k = (a + L) % len(P)                       # the first point past its end: the walk leaves the cut surface here
        return k, dict(kind="cut surface", length_um=int(L * STEP), surface_from=[round(float(v), 1) for v in P[a]],
                       at=[round(float(v), 1) for v in P[k]])
    big = lab == 5
    lw, nw = ndi.label(big)
    if nw:
        big = lw == (np.argmax(np.bincount(lw.ravel())[1:]) + 1)
    dw = ndi.distance_transform_edt(~big) * G if big.any() else np.zeros(lab.shape)
    qi = np.clip(np.round(P).astype(int), 0, np.array(lab.shape) - 1)
    d = dw[qi[:, 0], qi[:, 1]]
    rows = np.where(tis.any(1))[0]
    low = P[:, 0] >= rows.min() + prm["fb_low_frac"] * (rows.max() - rows.min())
    k = int(np.where(low)[0][np.argmin(d[low])]) if low.any() else int(np.argmax(P[:, 0]))
    return k, dict(kind="white matter nearest the lower outline", wm_um=round(float(d[k])), at=[round(float(v), 1) for v in P[k]])


# ----------------------------------------------------------------------------------------------------- the walk
def _bezier(p0, t0, p1, t1, step):
    """points of the cubic arc from p0 (leaving along t0) to p1 (arriving along t1), step apart by arc length,
    ends excluded."""
    c = float(np.hypot(*(p1 - p0)))
    c1, c2 = p0 + t0 * c / 3, p1 - t1 * c / 3
    s = np.linspace(0, 1, 400)[:, None]
    B = (1 - s) ** 3 * p0 + 3 * (1 - s) ** 2 * s * c1 + 3 * (1 - s) * s ** 2 * c2 + s ** 3 * p1
    L = np.r_[0, np.cumsum(np.hypot(*np.diff(B, axis=0).T))]
    m = max(int(round(L[-1] / step)), 1)
    u = np.linspace(0, L[-1], m + 1)[1:-1]
    return np.c_[np.interp(u, L, B[:, 0]), np.interp(u, L, B[:, 1])]


def _dense(P, k=4):
    """a polyline with k points per segment (for distances to it)."""
    if len(P) < 2:
        return P
    t = np.linspace(0, 1, k + 1)[:-1]
    D = (P[:-1, None, :] + t[None, :, None] * np.diff(P, axis=0)[:, None, :]).reshape(-1, 2)
    return np.r_[D, P[-1:]]


def _turn_open(P, w):
    """|turn| (deg) between the directions w points back and w ahead on an open stretch (0 within w of its ends)."""
    d = np.gradient(P, axis=0)
    th = np.unwrap(np.arctan2(d[:, 0], d[:, 1]))
    out = np.zeros(len(P))
    if len(P) > 2 * w:
        out[w:len(P) - w] = np.degrees(th[2 * w:] - th[:-2 * w])
    return np.abs(out)


def arc(P, e, allowed, midc, prm=PRM):
    """the corner cut for one sharp turn e of the closed walk P: the shortest stretch round it (at least the turn's own
    window) that can be replaced by a smooth arc leaving and joining the walk along its own direction, points as far
    apart as the walk's, with no turn of turn_deg or more on it or made at its ends, moving the line at most disp_um,
    sweeping over no cut (midc: cut cells) and taking in no furrow base or fissure mouth. Returns ((a, b, points),
    moved_um) or (None, why)."""
    n = len(P)
    w = int(round(SC.PRM["turn_win_um"] / STEP))
    k = prm["arc_tan"]
    keep = np.array([x["i"] for x in allowed], int)
    why = "no arc"
    for h in range(w, prm["arc_reach"] + 1, 2):
        a, b = e["a"] - h, e["b"] + h
        if b - a >= n - 4 * w:
            break
        if len(keep) and np.any(((keep - a) % n) <= (b - a)):
            why = "furrow base or fissure mouth in reach"
            break
        old = P[np.arange(a, b + 1) % n]
        p0, p1 = old[0], old[-1]
        t0 = P[a % n] - P[(a - k) % n]
        t1 = P[(b + k) % n] - P[b % n]
        t0, t1 = t0 / max(np.hypot(*t0), 1e-9), t1 / max(np.hypot(*t1), 1e-9)
        step = float(np.median(np.hypot(*np.diff(P[np.r_[np.arange(a - 10, a), np.arange(b, b + 10)] % n], axis=0).T)))
        C = _bezier(p0, t0, p1, t1, step)
        new = np.r_[p0[None], C, p1[None]]
        t_old, t_new = cKDTree(_dense(old)), cKDTree(_dense(new))
        moved = max(t_old.query(C)[0].max() if len(C) else 0.0, t_new.query(old)[0].max()) * G
        if moved > prm["disp_um"]:
            why = f"needs {moved:.0f} um"
            continue
        lo = np.floor(np.minimum(old.min(0), new.min(0))) - 1
        hi = np.ceil(np.maximum(old.max(0), new.max(0))) + 1
        mc = midc[(midc[:, 0] >= lo[0]) & (midc[:, 0] <= hi[0]) & (midc[:, 1] >= lo[1]) & (midc[:, 1] <= hi[1])]
        if len(mc) and len(C):
            inside = MPath(np.r_[old, C[::-1]]).contains_points(mc.astype(float))
            if inside.any() and (t_old.query(mc[inside])[0] > prm["arc_mid_tol"]).any():
                why = "would cross a cut"
                break
        m = 3 * w
        before, after = P[np.arange(a - m, a) % n], P[np.arange(b + 1, b + 1 + m) % n]
        tn = _turn_open(np.r_[before, new, after], w)
        to = _turn_open(P[np.arange(a - m, b + 1 + m) % n], w)
        if tn[m:m + len(new)].max() >= SC.PRM["turn_deg"]:
            why = f"still turns {tn[m:m + len(new)].max():.0f} deg"
            continue
        lim = SC.PRM["turn_deg"] - 1
        if (tn[w:m] > np.maximum(to[w:m], lim)).any() or (tn[m + len(new):-w] > np.maximum(to[m + len(old):-w], lim)).any():
            why = "turn moved to its end"
            continue
        return (a, b, new), moved
    return None, why


def smooth(P0, cut, mid, prm=PRM):
    """the walk: P0 (the outline, smoothed as before) with every sharp turn that is not a furrow base or a fissure
    mouth replaced by a smooth arc (arc), worst first, the whole walk checked again after each. v2 (the user's notes
    on flattening Version 11, s5_E: "failure to maintain the perimeter"): v1 stiffened the line (a smoothing spline
    whose stiffness jumped at the ends of each stiff stretch), which tore the walk at those ends (steps of 5-14 um in
    18 of 19 sections) and in places took it across a cut instead of round it. An arc keeps the walk's spacing and
    cannot cross a cut. Returns (P, report)."""
    midc = np.argwhere(mid)
    P, ids = P0.copy(), np.arange(len(P0))
    first, stuck, log = None, [], []
    for _ in range(prm["arc_rounds"]):
        ev = SC.classify(info(P, cut, mid))
        first = ev if first is None else first
        allowed = [x for x in ev if x["kind"] in ("base", "mouth")]
        todo = [x for x in ev if x["kind"] in ("open", "cut") and all(np.hypot(*(P[x["i"]] - s)) >= 3 for s, _ in stuck)]
        if not todo:
            break
        e = max(todo, key=lambda x: abs(x["deg"]))
        res, why = arc(P, e, allowed, midc, prm)
        if res is None:
            stuck.append((P[e["i"]].copy(), why))
            continue
        (a, b, new), moved = res, why
        n = len(P)
        rot = np.arange(a, a + n) % n
        R, Ri, L = P[rot], ids[rot], b - a + 1
        P = np.r_[new, R[L:]]
        ids = np.r_[Ri[:1], np.full(len(new) - 2, -1), Ri[L - 1:L], Ri[L:]]
        stuck = [(s, y) for s, y in stuck if np.hypot(*(new - s).T).min() > 30]   # the walk changed there: try again
        log.append(dict(kind=e["kind"], deg=e["deg"], moved_um=round(moved, 1), at=[round(float(v), 1) for v in P[len(new) // 2]]))
    ok = np.where(ids >= 0)[0]
    P = np.roll(P, -int(ok[np.argmin(ids[ok])]), axis=0)          # the walk still starts where it started
    ev = SC.classify(info(P, cut, mid))
    left = [x for x in ev if x["kind"] in ("open", "cut")]
    reasons = []
    for x in left:
        d = [(np.hypot(*(P[x["i"]] - s)), y) for s, y in stuck]
        reasons.append(min(d)[1] if d and min(d)[0] < 6 else "not tried")
    steps = np.hypot(*np.diff(np.r_[P, P[:1]], axis=0).T) * G
    rep = dict(sharp=sum(x["kind"] in ("open", "cut") for x in first), smoothed=len(log), max_step_um=round(float(steps.max()), 2),
               moved_max_um=max([r["moved_um"] for r in log], default=0.0), arcs=log,
               left=[dict(i=int(x["i"]), deg=x["deg"], kind=x["kind"], why=r) for x, r in zip(left, reasons)])
    return P, rep


def redraw(P, cut, mid):
    """the tissue the columns run through, redrawn to the walk: inside the walk's polygon, less the cuts, plus the
    tissue the walk does not go round (other pieces)."""
    H, W = cut.shape
    poly = np.zeros_like(cut)
    rr, cc = draw_polygon(P[:, 0], P[:, 1], shape=(H, W))
    poly[rr, cc] = True
    lab_, n_ = ndi.label(cut)
    qi = np.clip(np.round(P).astype(int), 0, [H - 1, W - 1])
    inside = np.zeros_like(cut)
    for k in range(1, n_ + 1):                      # the pieces the walk goes round
        m = lab_ == k
        if (poly & m).sum() > 0.5 * m.sum():
            inside |= m
    return (cut & ~inside) | (poly & ~mid)


# ------------------------------------------------------------------------------------------------- the cut tips
def tips(mid, tissue, back=10):
    """deep ends of the cuts: skeleton end points away from the outside, with the cut's direction there."""
    sk = morphology.skeletonize(mid)
    ends = F.skel_ends(sk)
    dout = ndi.distance_transform_edt(tissue | mid)
    nb = [(a, b) for a in (-1, 0, 1) for b in (-1, 0, 1) if a or b]
    H, W = mid.shape
    out = []
    for t in np.argwhere(ends):
        if dout[tuple(t)] <= 3:
            continue
        path, seen = [tuple(t)], {tuple(t)}
        for _ in range(back):
            r, c = path[-1]
            nxt = None
            for a, b in nb:
                q = (r + a, c + b)
                if 0 <= q[0] < H and 0 <= q[1] < W and sk[q] and q not in seen:
                    nxt = q
                    break
            if nxt is None:
                break
            path.append(nxt)
            seen.add(nxt)
        if len(path) < 4:
            continue
        d = np.array(path[0], float) - np.array(path[min(len(path) - 1, 6)], float)
        out.append((tuple(int(v) for v in t), d / np.hypot(*d), path))
    return out


def _two_sided(q, dq, cortex, blk, side, U, prm):
    """cortex to the left and to the right of a step, at similar distances (as flatten.furrow_walk)."""
    H, W = cortex.shape
    pp = np.round(np.array(q, float)[None, None, :] + U[:, None, :] * side[None, :, None]).astype(int)
    pp[..., 0] = np.clip(pp[..., 0], 0, H - 1)
    pp[..., 1] = np.clip(pp[..., 1], 0, W - 1)
    hc, hb = cortex[pp[..., 0], pp[..., 1]], blk[pp[..., 0], pp[..., 1]]
    big = 10 ** 6
    fc = np.where(hc.any(1), np.argmax(hc, 1), big)
    fb = np.where(hb.any(1), np.argmax(hb, 1), big)
    hit = (fc < big) & (fc < fb)
    rel = np.degrees(np.arctan2(U[:, 0] * dq[1] - U[:, 1] * dq[0], U[:, 0] * dq[0] + U[:, 1] * dq[1]))
    off = 90.0 - prm["tip_side_deg"]
    lm, rm = hit & (rel >= off) & (rel <= 180 - off), hit & (rel <= -off) & (rel >= off - 180)
    if not lm.any() or not rm.any():
        return False
    dl, dr = side[fc[lm]].min(), side[fc[rm]].min()
    return bool(max(dl, dr) <= prm["tip_balance"] * max(min(dl, dr), 2))


def extend(t, d0, t_um, tissue, mid, cortex, eglish, dcx, dcs, dmid, blk, prm=PRM):
    """from a cut's deep tip, on down the middle of the fused EGL to where the two EGLs part."""
    H, W = tissue.shape
    rot = lambda v, a: np.array([v[0] * np.cos(a) - v[1] * np.sin(a), v[0] * np.sin(a) + v[1] * np.cos(a)])
    turns = [0.0] + [sg * np.radians(x) for x in (15.0, 30.0) for sg in (1, -1)]
    side = np.arange(1, int(prm["tip_side"] * max(t_um, 16.0) / G) + 2)
    ang = np.linspace(0, 2 * np.pi, prm["tip_rays"], endpoint=False)
    U = np.stack([np.cos(ang), np.sin(ang)], 1)
    pos, d = np.array(t, float), np.array(d0, float)
    path, good, why = [], [], "length"
    for step in range(int(prm["tip_max_um"] / G)):
        best = None
        for a in turns:
            dd = rot(d, a)
            if (dd * d0).sum() < np.cos(np.radians(prm["tip_max_turn"])):
                continue
            q = tuple(np.round(pos + dd).astype(int))
            if not (0 <= q[0] < H and 0 <= q[1] < W) or not tissue[q] or not eglish[q]:
                continue
            if step > 2 and dmid[q] <= 1.5:
                continue
            qa = tuple(np.clip(np.round(pos + 3 * dd).astype(int), 0, [H - 1, W - 1]))
            sc = dcs[q] + dcs[qa] + (0.5 if a == 0 else 0.0)
            if best is None or sc > best[0]:
                best = (sc, dd, q)
        if best is None:
            why = "no EGL ahead"
            break
        _, dd, q = best
        if len(path) >= 2:
            stop = max(prm["tip_stop_min_um"], prm["tip_stop_frac"] * float(np.median([dcx[x] for x in path])))
            if dcx[q] < stop:
                why = "base"
                break
        d, pos = dd, pos + dd
        path.append(q)
        good.append(_two_sided(q, dd, cortex, blk, side, U, prm))
    while good and not good[-1]:
        good.pop()
        path.pop()
    return path, good, why


def tip_cuts(lab, tissue, tis_ref, mid, egl_s, prm=PRM, fprm=None):
    """cuts continued from their deep tips down the fused EGL (point 4). Returns (new cut cells, report)."""
    fprm = fprm or F.PRM
    cortex = tissue & np.isin(lab, (3, 4, 5, 6)) & ~egl_s
    cortex = morphology.remove_small_objects(cortex, max_size=fprm["pk_speck_cells"], connectivity=2)
    eglish = tissue & (np.isin(lab, (1, 2)) | egl_s)
    dcx = ndi.distance_transform_edt(~cortex) * G
    dcs = ndi.gaussian_filter(dcx, 1.0)
    dmid = ndi.distance_transform_edt(~mid)
    blk = ~tissue | mid
    add = np.zeros_like(mid)
    rep = []
    for t, d0, back in tips(mid, tissue, prm["tip_back_cells"]):
        b = np.array(back)
        t_um = float(np.median(dcx[b[:, 0], b[:, 1]]))           # one EGL's thickness beside the cut's end
        path, good, why = extend(t, d0, t_um, tissue, mid, cortex, eglish, dcx, dcs, dmid, blk, prm)
        L_um = len(path) * G
        keep = bool(path) and L_um >= max(prm["tip_min_um"], t_um) and np.mean(good) >= prm["tip_pass"]
        straight = float(np.hypot(*(np.array(path[-1]) - np.array(t)))) / len(path) if path else 0.0
        if keep and straight < prm["tip_straight"]:
            keep, why = False, f"{why}; straightness {straight:.2f}"
        r = dict(at=[int(v) for v in t], egl_um=round(t_um), walk_um=L_um, end=why, two_sided=round(float(np.mean(good)), 2) if good else 0.0,
                 kept=bool(keep))
        if path:
            pa = np.array(path)
            r["egl_label_frac"] = round(float(np.isin(lab[pa[:, 0], pa[:, 1]], (1, 2)).mean()), 2)
            r["to"] = [int(v) for v in path[-1]]
        rep.append(r)
        if keep:
            for q in path:
                add[q] = True
    add = ndi.binary_dilation(add, structure=np.ones((2, 2), bool)) & tissue & ~mid
    if prm.get("tip_check") and add.any():
        add = _check(lab, tis_ref, mid, add, rep, fprm)
    return add, rep


def _check(lab, tis_ref, mid, add, rep, fprm):
    """each kept continuation must not add columns that reach no IGL (flatten._check_peak_cuts)."""
    from scipy.spatial import cKDTree
    C0 = F._columns(lab, tis_ref, mid, fprm, full=True)
    comp, nc = ndi.label(add, structure=np.ones((3, 3), bool))
    out = np.zeros_like(add)
    r0 = fprm["pk_check_um"] / G
    for k in range(1, nc + 1):
        m = comp == k
        tc_ = cKDTree(np.argwhere(m))
        C1 = F._columns(lab, tis_ref, mid | m, fprm, full=True)
        gain = F._short_um(C1, tc_.query(C1["Q"])[0] <= r0) - F._short_um(C0, tc_.query(C0["Q"])[0] <= r0)
        ok = gain <= fprm["pk_short_um"]
        for r in rep:
            if r["kept"] and tc_.query(r["at"])[0] <= 3:
                r["short_gain_um"] = round(gain)
                if not ok:
                    r["kept"], r["end"] = False, r["end"] + f"; +{gain:.0f} um of columns reaching no IGL"
        if ok:
            out |= m
    return out


def band(lab, tis, add, rep, egl_s):
    """the fused EGL a continued cut runs through, relabelled EGL where the stain shows EGL (s2_2 at 7.8 mm: the band
    below the old cut end is labelled white matter): cells within one EGL thickness of a kept continuation, inside the
    tissue, stain EGL, labelled unassigned / ML / white matter / anomaly."""
    out = np.zeros_like(add)
    if not add.any():
        return out
    comp, nc = ndi.label(add, structure=np.ones((3, 3), bool))
    near = ndi.distance_transform_edt(comp == 0, return_indices=True)[1]
    cand = tis & egl_s & np.isin(lab, (1, 3, 5, 6))
    for r in rep:
        if not r["kept"]:
            continue
        k = comp[tuple(near[:, r["at"][0], r["at"][1]])]
        if k:
            out |= (ndi.distance_transform_edt(comp != k) * G <= r["egl_um"]) & cand
    return out


# ------------------------------------------------------------------------------------------------- the EGL cap
def cap(X, valid, lab, egl_s, prm=PRM):
    """thin caps on the EGL, found in the columns (so a cap is between the line and the EGL, never across a cut):
    from the top of a column (skipping up to two background samples), a run of ML / unassigned / anomaly samples, at
    most cap_um, that the stain calls EGL, directly followed by an EGL label. Returns the grid cells (point 5)."""
    H, W = lab.shape
    N, Dn = valid.shape
    qi = np.clip(np.round(X).astype(int), 0, [H - 1, W - 1])
    L = np.where(valid, lab[qi[..., 0], qi[..., 1]], 0)
    E = egl_s[qi[..., 0], qi[..., 1]] & valid
    kmax = int(prm["cap_um"] / STEP)
    out = np.zeros_like(lab, dtype=bool)
    for i in range(N):
        k = 0
        while k < 2 and k < Dn and L[i, k] == 0:
            k += 1
        a = k
        while k < Dn and k - a <= kmax and np.isin(L[i, k], (1, 3, 6)) and E[i, k]:
            k += 1
        if k > a and k - a <= kmax and k < Dn and L[i, k] == 2:
            c = qi[i, a:k]
            out[c[:, 0], c[:, 1]] = True
    return out & np.isin(lab, (1, 3, 6)) & egl_s


# --------------------------------------------------------------------------------------------- the whole line
def reference(short, lab, tis_p, mid, egl_s, prm=PRM, fprm=None):
    """the walked reference line: refs as flatten.reference (main walk first, then any other pieces), the tissue the
    columns run through, and a report."""
    fprm = fprm or F.PRM
    cut = F.cut_tissue(tis_p, fprm) & ~mid
    cut = morphology.remove_small_objects(cut, max_size=int(20000 / G ** 2))
    pieces = outline(cut)
    P_raw = pieces[0]
    P0 = ndi.gaussian_filter1d(P_raw, fprm["pos_sigma_um"] / STEP, axis=0, mode="wrap")
    n0 = normals(P0, cut)
    k, st = start(P0, n0, lab, cut, egl_s, prm)
    P0 = np.roll(P0, -k, axis=0)
    ev0 = SC.classify(info(P0, cut, mid))
    P, srep = smooth(P0, cut, mid, prm)
    cut2 = redraw(P, cut, mid)
    refs = []
    n = normals(P, cut2)
    refs.append(dict(Q=P, n=n, closed=True))
    for Pk in pieces[1:]:
        Qk = ndi.gaussian_filter1d(Pk, fprm["pos_sigma_um"] / STEP, axis=0, mode="wrap")
        refs.append(dict(Q=Qk, n=normals(Qk, cut2), closed=True))
    ev1 = SC.classify(info(P, cut2, mid))
    count = lambda ev: {k_: sum(e["kind"] == k_ for e in ev) for k_ in ("open", "cut", "base", "mouth")}
    rep = dict(start=st, clockwise=True, pieces=len(pieces), turns_before=count(ev0), turns_after=count(ev1), **srep,
               filled_um2=int((cut2 & ~cut).sum() * G * G), trimmed_um2=int((cut & ~cut2).sum() * G * G))
    return refs, cut2, cut, rep


def labels_after(short, lab, cut_old, cut_new, egl_s):
    """labels for the redrawn tissue: trimmed -> unassigned, filled -> EGL where the stain shows EGL, else unassigned
    (gold: only filled background cells change, to unassigned)."""
    lab = lab.copy()
    fill = cut_new & ~cut_old
    trim = cut_old & ~cut_new
    if short == GOLD:
        lab[fill & (lab == 0)] = 1
        return lab, fill, trim
    lab[trim & (lab > 0)] = 1
    lab[fill & egl_s] = 2
    lab[fill & ~egl_s & (lab == 0)] = 1
    lab[fill & ~egl_s & np.isin(lab, (3, 4, 5, 6))] = 1
    return lab, fill, trim
