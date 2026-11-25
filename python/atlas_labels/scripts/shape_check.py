"""
Shape check of the reference line (2026-11-24, the user's note on the flattening walkthrough: "the EGL label and the
line that denotes the flattening take a turn into the arachnoid or into the ML ... the EGL almost never makes sharp
turns. You should be able to detect shape violation and then correct").

The line is the outer face of the EGL. The EGL is a smooth ribbon, so the line only turns sharply where the ribbon
folds: at the base of a fissure (the line runs down one side of a cut and back up the other, a hairpin at the cut's
end) and at a fissure mouth (the open surface turns into a cut). Any other sharp turn is a shape violation.

1. Turn: the angle between the line's direction turn_win_um back and turn_win_um ahead (2 um samples), signed + where
   the line bends toward the tissue (convex), - away from it (concave). Sharp: |turn| >= turn_deg.
2. Not a violation:
   - furrow base: a hairpin (turn <= -base_deg) with base_cut or more of the turn on the sides of a cut;
   - fissure mouth: a convex turn within mouth_um of where the line enters or leaves a cut.
   Line points within oncut_cells of a cut are on a cut.
3. Fixes on the open surface (the cuts are fixed in flatten.prune_cuts):
   a. Pockets: a stretch of the line whose ends are at most pk_chord_um apart while the stretch is at least pk_ratio
      times longer (and pk_extra_um more) is a detour. With the stretch on the tissue side of its chord it is a notch
      (a crack the meninges trim cut into EGL, a gap between the EGL and a strip of EGL lying on it): filled. On the
      outside it is a bump (a flap of meninges or a loose strip of EGL): cut off. Not touched: notches deeper than
      pk_fissure_um whose two walls both have cortex behind them (a narrow open fissure), detours deeper than
      pk_max_depth_um, bumps holding pk_igl_frac or more stain IGL, and anything on a cut.
   b. Peel: a violation the pockets leave is the line running round meninges the trim had to keep because they carry
      a cortex label (s2_2 at 8.4: meninges labelled ML on the outer face of the EGL). Within peel_um of it, tissue
      with the meninges signature of egl_repair.pia_trim (low p27 + NeuN, not bright on DAPI), not stain EGL or IGL,
      in pieces touching the outside, is trimmed whatever its label.
   Repeated for up to `rounds` rounds on the new line. A fix is kept only when the sharp turns it was made for are gone
   and it adds none within fix_check_um.
Labels (not on the gold s1_2): cells trimmed become unassigned; filled cells become EGL where the stain shows EGL, else
unassigned. On the gold only background cells that are filled change (to unassigned), so the columns can start there.
"""
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi
from skimage.draw import polygon as draw_polygon

sys.path.insert(0, str(Path(__file__).parent))
import flatten as F  # noqa: E402

G, STEP = F.G, F.STEP
GOLD = "s1_2"
PRM = dict(turn_win_um=20.0, turn_deg=60.0, base_deg=140.0, base_cut=0.8, mouth_um=30.0, oncut_cells=2.0,
           pk_chord_um=48.0, pk_ratio=2.0, pk_extra_um=16.0, pk_max_arc_um=300.0, pk_max_depth_um=80.0, pk_fissure_um=40.0,
           pk_wall_um=120.0, pk_wall_frac=0.5, pk_igl_frac=0.3,
           pk_menin_frac=0.5, pk_menin_egl=0.3, peel_um=40.0, peel_on_egl=0.5, peel_egl_gap=2, peel_grow_um=300.0, snap_min_um=8.0, snap_max_um=48.0, snap_egl_um=12.0, snap_gap_um=8.0, snap_min_run_um=8.0, snap_run_um=40.0, snap_sig_frac=0.5, snap_bad_cells=1, snap_bad_frac=0.2, snap_taper=2.0, fix_check_um=60.0, fix_keep=0.5, fix_slack_um=24.0, fix_egl_gain=0.15, fix_at_cells=3, rounds=3)


def turning(Q, n, piece, w):
    """signed turn (deg) between the directions w samples back and w ahead; + = toward the tissue (convex)."""
    out = np.zeros(len(Q))
    for p in np.unique(piece):
        ix = np.where(piece == p)[0]
        P = Q[ix]
        closed = np.hypot(*(P[0] - P[-1])) < 2
        if closed:
            d = (np.roll(P, -1, 0) - np.roll(P, 1, 0)) / 2
            dd = np.r_[d[-w:], d, d[:w]]
            th = np.unwrap(np.arctan2(dd[:, 0], dd[:, 1]))
            dt = th[2 * w:] - th[:-2 * w]
        else:
            d = np.gradient(P, axis=0)
            th = np.unwrap(np.arctan2(d[:, 0], d[:, 1]))
            k = np.arange(len(P))
            dt = th[np.clip(k + w, 0, len(P) - 1)] - th[np.clip(k - w, 0, len(P) - 1)]
        t = d / np.maximum(np.hypot(*d.T), 1e-9)[:, None]
        side = np.sign(t[:, 0] * n[ix, 1] - t[:, 1] * n[ix, 0])
        out[ix] = -np.degrees(dt) * side
    return out


def line(lab, tis, mid, prm=PRM, egl=None):
    """the reference line as flatten builds it, with turns and the points on a cut (egl: the EGL map for egl_top)."""
    refs, cut, mid2 = F.reference(lab, tis, F.PRM, mid=mid)
    Q = np.concatenate([r["Q"] for r in refs])
    n = np.concatenate([r["n"] for r in refs])
    piece = np.concatenate([np.full(len(r["Q"]), i) for i, r in enumerate(refs)])
    tu = turning(Q, n, piece, int(round(prm["turn_win_um"] / STEP)))
    dm = ndi.distance_transform_edt(~mid2)
    qi = np.clip(np.round(Q).astype(int), 0, np.array(mid.shape) - 1)
    oncut = dm[qi[:, 0], qi[:, 1]] <= prm["oncut_cells"]
    return dict(Q=Q, n=n, piece=piece, tu=tu, oncut=oncut, cut=cut, mid=mid2, egl=egl if egl is not None else lab == 2)


def classify(Ln, prm=PRM):
    """events (runs of |turn| >= turn_deg of one sign) with their kind: base, mouth or violation."""
    tu, oncut, piece = Ln["tu"], Ln["oncut"], Ln["piece"]
    jn = np.where((np.diff(oncut.astype(int)) != 0) & (piece[1:] == piece[:-1]))[0]
    ev = []
    for sg in (1, -1):
        lab_, _ = ndi.label(sg * tu >= prm["turn_deg"])
        for sl in ndi.find_objects(lab_):
            a, b = sl[0].start, sl[0].stop
            k = a + int(np.argmax(sg * tu[a:b]))
            dj = float(np.min(np.abs(jn - k))) * STEP if len(jn) else 1e9
            cutf = float(oncut[a:b].mean())
            if sg < 0 and tu[k] <= -prm["base_deg"] and cutf >= prm["base_cut"]:
                kind = "base"
            elif sg > 0 and dj <= prm["mouth_um"]:
                kind = "mouth"
            elif cutf >= 0.5:
                kind = "cut"
            else:
                kind = "open"
            ev.append(dict(i=k, a=a, b=b, deg=int(round(tu[k])), cut=round(cutf, 2), kind=kind))
    return sorted(ev, key=lambda e: e["i"])


def _pockets(Ln, prm=PRM):
    """detours: line stretches [a, b] (indices into Q) whose ends are close, maximal and merged."""
    Q, piece, oncut = Ln["Q"], Ln["piece"], Ln["oncut"]
    kmin, kmax = 4, int(prm["pk_max_arc_um"] / STEP)
    out = []
    for p in np.unique(piece):
        ix = np.where(piece == p)[0]
        P = Q[ix]
        nn = len(P)
        closed = np.hypot(*(P[0] - P[-1])) < 2
        best = np.full(nn, -1)
        for k in range(kmin, min(kmax, nn - 1)):
            b = (np.arange(nn) + k) % nn if closed else np.minimum(np.arange(nn) + k, nn - 1)
            ch = np.hypot(*(P[b] - P).T) * G
            arc = k * STEP
            ok = (ch <= prm["pk_chord_um"]) & (arc >= prm["pk_ratio"] * ch) & (arc - ch >= prm["pk_extra_um"])
            if not closed:
                ok &= np.arange(nn) + k < nn
            best[ok] = k
        iv = sorted((a, a + best[a]) for a in range(nn) if best[a] > 0)
        merged = []
        for a, b in iv:
            if merged and a <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], b)
            else:
                merged.append([a, b])
        for a, b in merged:
            loc = np.arange(a, b + 1) % nn if closed else np.arange(a, min(b, nn - 1) + 1)
            out.append(ix[loc])
    return out


def _wall_cortex(Ln, gi, lab, tis, prm=PRM):
    """share of the line points gi whose inward ray reaches ML / IGL / white matter within pk_wall_um inside the tissue."""
    Q, n = Ln["Q"][gi], Ln["n"][gi]
    steps = np.arange(1, int(prm["pk_wall_um"] / G) + 1)
    hit = np.zeros(len(gi), bool)
    alive = np.ones(len(gi), bool)
    H, W = lab.shape
    for s in steps:
        p = np.round(Q + n * s).astype(int)
        inb = (p[:, 0] >= 0) & (p[:, 0] < H) & (p[:, 1] >= 0) & (p[:, 1] < W)
        p = np.clip(p, 0, [H - 1, W - 1])
        alive &= inb & tis[p[:, 0], p[:, 1]] if s > 2 else inb
        hit |= alive & np.isin(lab[p[:, 0], p[:, 1]], (3, 4, 5))
    return hit


def pocket_fixes(Ln, viol_idx, lab, tis, igl_s, sig=None, egl_s=None, tissue=None, prm=PRM):
    """candidate fills and trims for the detours that hold a violation."""
    fixes = []
    vset = np.zeros(len(Ln["Q"]), bool)
    vset[viol_idx] = True
    for gi in _pockets(Ln, prm):
        if not vset[gi].any() or Ln["oncut"][gi].mean() > 0.0:
            continue
        P = Ln["Q"][gi]
        A, B = P[0], P[-1]
        v = B - A
        L = max(np.hypot(*v), 1e-6)
        depth = float((np.abs((P[:, 0] - A[0]) * v[1] - (P[:, 1] - A[1]) * v[0]) / L).max() * G)
        if depth > prm["pk_max_depth_um"]:
            continue
        rr, cc = draw_polygon(P[:, 0], P[:, 1], shape=tis.shape)
        R = np.zeros_like(tis)
        R[rr, cc] = True
        m = (A + B) / 2
        notch = not Ln["cut"][int(round(m[0])), int(round(m[1]))]
        if notch:
            if depth >= prm["pk_fissure_um"]:
                h = len(gi) // 2
                w1 = _wall_cortex(Ln, gi[:h], lab, tis, prm).mean()
                w2 = _wall_cortex(Ln, gi[h:], lab, tis, prm).mean()
                if min(w1, w2) >= prm["pk_wall_frac"]:
                    continue                                   # a narrow open fissure: cortex behind both walls
            cells = R & ~tis
            if sig is not None:
                inner = cells & tissue                         # tissue the meninges trim took
                if inner.sum() >= 4 and sig[inner].mean() >= prm["pk_menin_frac"] and egl_s[inner].mean() < prm["pk_menin_egl"]:
                    continue                                   # meninges fill the notch: a real gap, not a crack in the EGL
            kind = "fill"
        else:
            cells = R & tis
            if cells.any() and igl_s[cells].mean() >= prm["pk_igl_frac"]:
                continue
            kind = "trim"
        if cells.any():
            fixes.append(dict(kind=kind, cells=cells, idx=gi, pts=Ln["Q"][gi], depth_um=round(depth, 1), chord_um=round(float(L * G), 1)))
    return fixes


def peel_fixes(Ln, viol_idx, tis, lab, sig, egl_s, igl_s, prm=PRM):
    """meninges the line runs round: pieces of meninges-signature tissue (not stain EGL or IGL) touching the outside
    within peel_um of a violation. A piece lying on the EGL (peel_on_egl or more of its border with the rest of the
    tissue within peel_egl_gap cells of EGL, label or stain) goes as far as it reaches (up to peel_grow_um from the violation): it is the
    meninges strip the line follows outside the EGL. Other pieces only within peel_um."""
    fixes = []
    if not len(viol_idx):
        return fixes
    k8 = np.ones((3, 3), bool)
    r = int(round(prm["peel_um"] / G))
    outside = ndi.binary_dilation(~tis, k8)
    cand = tis & sig & ~egl_s & ~igl_s
    comp, nc = ndi.label(cand, structure=k8)
    if not nc:
        return fixes
    touch = np.unique(comp[outside & cand])
    touch = set(touch[touch > 0].tolist())
    slices = ndi.find_objects(comp)
    eglish = ndi.binary_dilation(egl_s | (lab == 2), iterations=prm["peel_egl_gap"])   # a cell of mixed stain between
    on_egl = {}
    lab_v, _ = ndi.label(np.isin(np.arange(len(Ln["Q"])), viol_idx))
    for sl in ndi.find_objects(lab_v):
        gi = np.arange(sl[0].start, sl[0].stop)
        W = np.zeros_like(tis)
        q = np.clip(np.round(Ln["Q"][gi]).astype(int), 0, np.array(tis.shape) - 1)
        W[q[:, 0], q[:, 1]] = True
        Wg = ndi.distance_transform_edt(~W) * G
        W = Wg <= prm["peel_um"]
        cells = np.zeros_like(tis)
        for k in set(np.unique(comp[W & cand]).tolist()) & touch:
            if k not in on_egl:
                s_ = tuple(slice(max(a_.start - 2, 0), a_.stop + 2) for a_ in slices[k - 1])
                m = comp[s_] == k
                ring = ndi.binary_dilation(m, k8) & ~m & tis[s_] & ~cand[s_]
                on_egl[k] = bool(ring.any() and eglish[s_][ring].mean() >= prm["peel_on_egl"])
            reach = prm["peel_grow_um"] if on_egl[k] else prm["peel_um"]
            cells |= (comp == k) & (Wg <= reach)
        if cells.any():
            fixes.append(dict(kind="peel", cells=cells, idx=gi, pts=Ln["Q"][gi]))
    return fixes


def snap_fixes(Ln, viol_idx, tis, lab, sig, egl_s, igl_s, ml_s, prm=PRM):
    """the line running outside the EGL: along the inward normal the first snap_min_um or more is not EGL but
    meninges (signature or no stain call; not stain IGL, not stain ML without the signature) and then the stain
    EGL starts, at most snap_max_um down and at least snap_egl_um thick. Such columns are followed along the line
    from a violation (up to peel_grow_um, allowing snap_gap_um of other columns) and their cells above the EGL are
    trimmed, whatever their label: nothing of the cortex lies outside the EGL."""
    fixes = []
    Q, n, piece = Ln["Q"], Ln["n"], Ln["piece"]
    H, W = tis.shape
    K = int(prm["snap_max_um"] / G) + int(prm["snap_egl_um"] / G) + 1
    s = np.arange(K) * 1.0
    P = Q[:, None, :] + n[:, None, :] * s[None, :, None] + 0.5 * n[:, None, :]
    pi = np.clip(np.round(P).astype(int), 0, [H - 1, W - 1])
    inb = tis[pi[..., 0], pi[..., 1]]
    E_ = egl_s[pi[..., 0], pi[..., 1]] & inb
    bad = (igl_s[pi[..., 0], pi[..., 1]] | (ml_s & ~sig)[pi[..., 0], pi[..., 1]]) & inb
    first = np.where(E_.any(1), np.argmax(E_, 1), K)
    kt = int(prm["snap_egl_um"] / G)
    idx = np.minimum(first[:, None] + np.arange(kt)[None, :], K - 1)
    thick = E_[np.arange(len(Q))[:, None], idx].all(1)
    above = np.arange(K)[None, :] < first[:, None]
    nbad = (bad & above).sum(1)
    clean = nbad <= np.maximum(prm["snap_bad_cells"], prm["snap_bad_frac"] * first)   # a cell of mixed stain at the EGL edge
    ok = (first * G >= prm["snap_min_um"]) & (first * G <= prm["snap_max_um"]) & thick & clean
    sg = sig[pi[..., 0], pi[..., 1]] & inb
    sigf = (sg & above).sum(1) / np.maximum(above.sum(1), 1)
    grow = int(prm["peel_grow_um"] / STEP)
    gap = int(prm["snap_gap_um"] / STEP)
    lab_v, _ = ndi.label(np.isin(np.arange(len(Q)), viol_idx))
    for sl in ndi.find_objects(lab_v):
        a, b = sl[0].start, sl[0].stop
        p = piece[a]
        sel = np.zeros(len(Q), bool)
        for direction, start in ((-1, a - 1), (1, b)):
            j, miss = start, 0
            for _ in range(grow):
                if j < 0 or j >= len(Q) or piece[j] != p or miss > gap:
                    break
                if ok[j]:
                    sel[j], miss = True, 0
                else:
                    miss += 1
                j += direction
        sel[a:b] |= ok[a:b]
        if sel.sum() * STEP < prm["snap_min_run_um"]:
            continue
        cells = np.zeros_like(tis)
        for j in np.where(sel)[0]:
            q = pi[j, :first[j]]
            cells[q[:, 0], q[:, 1]] = True
        cells = ndi.binary_closing(cells, iterations=1) & tis & ~egl_s
        if cells.any():
            fixes.append(dict(kind="snap", cells=cells, idx=np.arange(a, b), pts=Q[a:b]))
    # the same offset without a sharp turn: runs of snap_run_um or more columns whose cells above the EGL are mostly
    # meninges (signature on snap_sig_frac or more of them)
    run_ok = ok & (sigf >= prm["snap_sig_frac"])
    lab_r, _ = ndi.label(run_ok)
    for sl in ndi.find_objects(lab_r):
        a, b = sl[0].start, sl[0].stop
        if (b - a) * STEP < prm["snap_run_um"] or len(np.unique(piece[a:b])) > 1:
            continue
        cells = np.zeros_like(tis)
        for j in range(a, b):
            q = pi[j, :first[j]]
            cells[q[:, 0], q[:, 1]] = True
        # taper the ends (2 um along the line per 1 um of depth) through cells above the EGL that are not cortex by
        # stain, so the line does not step in at a corner
        for j0, dj in ((a, -1), (b - 1, 1)):
            T = int(prm["snap_taper"] * first[j0] * G / STEP)
            for t in range(1, T + 1):
                j = j0 + dj * t
                if j < 0 or j >= len(Q) or piece[j] != piece[j0]:
                    break
                d = int(min(first[j], np.floor(first[j0] * (1 - t / (T + 1)))))
                badj = np.where((igl_s[pi[j, :d, 0], pi[j, :d, 1]] | E_[j, :d]))[0]
                d = badj[0] if len(badj) else d
                if d <= 0:
                    break
                q = pi[j, :d]
                cells[q[:, 0], q[:, 1]] = True
        cells = ndi.binary_closing(cells, iterations=1) & tis & ~egl_s
        if cells.any():
            fixes.append(dict(kind="snap_run", cells=cells, idx=np.arange(a, b), pts=Q[a:b]))
    return fixes


def meninges_signature(S, lab, prm):
    """egl_repair.pia_trim's arachnoid cells: p27 + NeuN (log, against the section's ML) below arach_thr, nuclei not
    bright (below arach_bright of the way from the ML to the EGL level)."""
    M, E = lab == 3, lab == 2
    nm, pm = float(np.median(S["neun_a"][M])), float(np.median(S["p27a"][M]))
    um, ue = float(np.median(S["nuc_a"][M])), float(np.median(S["nuc_a"][E]))
    return ((S["neun_a"] - nm) + (S["p27a"] - pm) < prm["arach_thr"]) & (S["nuc_a"] < um + prm["arach_bright"] * (ue - um))


def egl_top(Ln, sel):
    """share of the line points sel whose first cell inward (either of the first two) is EGL, label or stain."""
    if not sel.any():
        return 0.0
    q = Ln["Q"][sel]
    nn = Ln["n"][sel]
    H, W = Ln["egl"].shape
    hit = np.zeros(len(q), bool)
    for k in (0.5, 1.5):
        p = np.clip(np.round(q + nn * k).astype(int), 0, [H - 1, W - 1])
        hit |= Ln["egl"][p[:, 0], p[:, 1]]
    return float(hit.mean())


def _viol(Ln, prm):
    ev = classify(Ln, prm)
    idx = np.concatenate([np.arange(e["a"], e["b"]) for e in ev if e["kind"] == "open"] or [np.zeros(0, int)])
    return ev, idx


def check(short, lab, tissue, tis_p, mid, prm=PRM, log=None):
    """shape-check the line of one section. Returns the corrected labels, pia tissue, filled and trimmed cells and a
    report (violations before / after, the fixes kept)."""
    import egl_repair as ER
    gold = short == GOLD
    E = np.load(ER.OUT / f"{short}.npz")
    S = ER.stain_maps(short, E["lab1"], tissue)
    sig = meninges_signature(S, E["lab1"], ER.PRM) & tissue
    egl_s, igl_s = S["egl"], S["igl"]
    lab2, tis2 = lab.copy(), tis_p.copy()
    filled, trimmed = np.zeros_like(tissue), np.zeros_like(tissue)
    eglish = egl_s | (lab == 2)
    Ln = line(lab2, tis2, mid, prm, eglish)
    ev0, vidx = _viol(Ln, prm)
    kept = []
    for rnd in range(prm["rounds"]):
        if not len(vidx) and rnd > 0:
            break
        # meninges first (the line runs round them), then the detours left
        fx = (peel_fixes(Ln, vidx, tis2, lab2, sig, egl_s, igl_s, prm) + snap_fixes(Ln, vidx, tis2, lab2, sig, egl_s, igl_s, S["ml"], prm)
              + pocket_fixes(Ln, vidx, lab2, tis2, igl_s, sig, egl_s, tissue, prm))
        if not fx:
            break
        changed = False
        for f in fx:
            t2 = tis2.copy()
            if f["kind"] == "fill":
                t2 |= f["cells"]
            else:
                t2 &= ~f["cells"]
            L2 = line(lab2, t2, mid, prm, eglish)
            # the violations near this fix before and after
            q0 = f["pts"]
            near = lambda L_: np.min(np.hypot(*(L_["Q"][:, None, :] - q0[None, ::4, :]).transpose(2, 0, 1)), axis=1) * G <= prm["fix_check_um"]
            ev_a, idx_a = _viol(L2, prm)
            n_before = int(np.isin(np.where(near(Ln))[0], vidx).sum())
            n_after = int(np.isin(np.where(near(L2))[0], idx_a).sum())
            if f["kind"] in ("fill", "trim"):                  # shape only: the turns must go
                keep = n_before and n_after <= prm["fix_keep"] * n_before
            else:                                              # stains: the line must come onto the EGL, and may add
                dc = ndi.distance_transform_edt(~f["cells"])                    # at most one sharp corner
                at = lambda L_: dc[tuple(np.clip(np.round(L_["Q"]).astype(int), 0, np.array(dc.shape) - 1).T)] <= prm["fix_at_cells"]
                e0, e1 = egl_top(Ln, at(Ln)), egl_top(L2, at(L2))
                keep = (n_after <= n_before + prm["fix_slack_um"] / STEP) and e1 >= e0 + prm["fix_egl_gain"]
            if keep:
                tis2 = t2
                if f["kind"] == "fill":
                    filled |= f["cells"]
                else:
                    trimmed |= f["cells"]
                c = f["cells"]
                rc = np.argwhere(c).mean(0)
                kept.append(dict(kind=f["kind"], round=rnd, at=[round(float(rc[0]), 1), round(float(rc[1]), 1)],
                                 um2=int(c.sum() * G * G), q=[round(float(x), 1) for x in q0[len(q0) // 2]],
                                 labels={int(k): int(v) for k, v in zip(*np.unique(lab[c], return_counts=True))},
                                 before=n_before, after=n_after, **{k: f[k] for k in ("depth_um", "chord_um") if k in f}))
                Ln, vidx = L2, idx_a
                changed = True
        if not changed:
            break
    ev1, _ = _viol(Ln, prm)
    filled &= tis2
    trimmed &= ~tis2
    if not gold:
        lab2[trimmed & np.isin(lab2, (2, 3, 4, 5, 6))] = 1
        lab2[filled & egl_s & np.isin(lab2, (0, 1))] = 2
    lab2[filled & (lab2 == 0)] = 1
    count = lambda ev: {k: sum(e["kind"] == k for e in ev) for k in ("open", "cut", "base", "mouth")}
    rep = dict(before=count(ev0), after=count(ev1), fixes=kept,
               filled_um2=int(filled.sum() * G * G), trimmed_um2=int(trimmed.sum() * G * G),
               open_after=[dict(i=e["i"], deg=e["deg"]) for e in ev1 if e["kind"] == "open"],
               cut_after=[dict(i=e["i"], deg=e["deg"]) for e in ev1 if e["kind"] == "cut"])
    if log:
        log(f"{short} shape: {rep['before']} -> {rep['after']}, fixes {len(kept)}, filled {rep['filled_um2']} um2, trimmed {rep['trimmed_um2']} um2")
    return lab2, tis2, filled, trimmed, rep


def events_saved(Z, prm=PRM):
    """the sharp turns of a saved flattening (flatten npz: Q, n, piece, mid), classified as in check()."""
    Q, n, piece, mid = Z["Q"], Z["n"], Z["piece"], Z["mid"]
    tu = turning(Q, n, piece, int(round(prm["turn_win_um"] / STEP)))
    dm = ndi.distance_transform_edt(~mid)
    qi = np.clip(np.round(Q).astype(int), 0, np.array(mid.shape) - 1)
    oncut = dm[qi[:, 0], qi[:, 1]] <= prm["oncut_cells"]
    return classify(dict(Q=Q, n=n, piece=piece, tu=tu, oncut=oncut), prm)
