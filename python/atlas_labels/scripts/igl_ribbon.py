"""
The IGL as a ribbon, described from the IGL itself (step 1 of the IGL repair).

    from igl_ribbon import grid_labels, describe
    g = grid_labels(native_labels, px_um)          # 4 um grid: tissue, igl, egl, wm masks
    rib = describe(g)

describe():
  1. IGL mask cleaned: holes < 2000 um2 filled, closed over 8 um, opened over 6 um,
     pieces < 2000 um2 dropped.
  2. Centre line = skeleton of the cleaned mask; end branches shorter than
     max(2 x local width, 40 um) pruned, repeated until stable.
  3. Branches traced into ordered pixel paths. For every centre-line pixel: width
     (2 x distance to the IGL edge) and depth below the EGL (distance to the nearest
     EGL pixel, so closed fissures count as surface).
  4. Ends: position, direction (from the last 60 um of the branch), and width and depth
     as the median over the ribbon 40-200 um behind the end (the tip itself tapers).
  5. Flags along each branch: thinning where the width drops below `thin` x the median
     width of the same branch 150-450 um to either side, bulging above `bulge` x, each
     over at least 20 um.
"""
import numpy as np
from scipy import ndimage as ndi
from skimage import morphology

GRID = 4.0
ST = np.ones((3, 3), bool)
P = dict(hole_um2=2000.0, close_um=8.0, open_um=6.0, min_piece_um2=2000.0, prune_factor=2.0, prune_min_um=40.0,
         dir_um=60.0, win_um=200.0, thin=0.5, bulge=1.8, min_flag_um=20.0, flank_in_um=150.0, flank_out_um=450.0,
         flank_from_um=40.0, flank_to_um=200.0)


def grid_labels(lab, px, grid=GRID):
    z = px / grid
    out = {}
    for name, codes in (("tissue", (1, 2, 3, 4, 5, 6, 7, 8)), ("igl", (4,)), ("egl", (2, 3)), ("wm", (6, 8)),
                        ("ml", (5, 7)), ("anom", (9,)), ("unass", (1,))):
        out[name] = ndi.zoom(np.isin(lab, codes).astype(np.float32), z, order=1) >= 0.5
    return out


def clean_igl(igl, prm=P):
    r = lambda um: max(1, int(round(um / GRID)))
    m = igl.copy()
    holes = ndi.binary_fill_holes(m) & ~m
    hl, n = ndi.label(holes, structure=ST)
    if n:
        a = np.bincount(hl.ravel()) * GRID ** 2
        m |= np.isin(hl, 1 + np.nonzero(a[1:] < prm["hole_um2"])[0])
    k = r(prm["close_um"])
    m = ndi.binary_closing(np.pad(m, k + 1), structure=morphology.disk(k))[k + 1:-k - 1, k + 1:-k - 1]
    m = ndi.binary_opening(m, structure=morphology.disk(r(prm["open_um"])))
    return morphology.remove_small_objects(m, max_size=int(prm["min_piece_um2"] / GRID ** 2))


def _nbr_count(sk):
    return ndi.convolve(sk.astype(np.uint8), np.ones((3, 3), np.uint8), mode="constant") - sk


def _segments(sk):
    nb = _nbr_count(sk) * sk
    junc = sk & (nb >= 3)
    seg = sk & ~junc
    sl, ns = ndi.label(seg, structure=ST)
    jl, nj = ndi.label(junc, structure=ST)
    return sl, ns, jl, nj, nb


def prune(sk, width_um, prm=P):
    """drop end branches shorter than max(prune_factor x width, prune_min_um)."""
    sk = sk.copy()
    for _ in range(12):
        sl, ns, jl, nj, nb = _segments(sk)
        if ns == 0:
            break
        ends = sk & (nb == 1)
        removed = False
        objs = ndi.find_objects(sl)
        jd = ndi.binary_dilation(jl > 0, structure=ST)
        for i, so in enumerate(objs, 1):
            if so is None:
                continue
            m = sl[so] == i
            n_end = int((ends[so] & m).sum())
            touches = jd[so] & m
            if n_end == 0 or not touches.any():
                continue                               # interior branch or isolated piece
            length = m.sum() * GRID * 1.1
            w = np.median(width_um[so][m])
            if length < max(prm["prune_factor"] * w, prm["prune_min_um"]):
                sk[so][m] = False
                removed = True
        if not removed:
            break
        sk = morphology.skeletonize(sk)
    return sk


def _order(mask_pts_set, start):
    """walk a branch from start through 8-neighbours, returning an ordered list."""
    path, seen, cur = [start], {start}, start
    while True:
        nxt = None
        for dy in (-1, 0, 1):
            for dx in (-1, 0, 1):
                q = (cur[0] + dy, cur[1] + dx)
                if (dy or dx) and q in mask_pts_set and q not in seen:
                    nxt = q
                    break
            if nxt:
                break
        if nxt is None:
            return path
        seen.add(nxt)
        path.append(nxt)
        cur = nxt


def describe(g, prm=P):
    igl = clean_igl(g["igl"], prm)
    edt = ndi.distance_transform_edt(igl) * GRID
    width = 2 * edt
    dE = ndi.distance_transform_edt(~g["egl"]) * GRID
    sk = prune(morphology.skeletonize(igl), width, prm)
    sl, ns, jl, nj, nb = _segments(sk)
    ends_mask = sk & (nb == 1)
    branches, flags = [], []
    win = int(round(prm["win_um"] / GRID))
    minf = max(1, int(round(prm["min_flag_um"] / GRID)))
    for i, so in enumerate(ndi.find_objects(sl), 1):
        if so is None:
            continue
        pts = np.argwhere(sl[so] == i) + [so[0].start, so[1].start]
        S = set(map(tuple, pts))
        # start at a branch end: a pixel with one neighbour inside the branch
        start = None
        for p in map(tuple, pts):
            k = sum((p[0] + dy, p[1] + dx) in S for dy in (-1, 0, 1) for dx in (-1, 0, 1) if dy or dx)
            if k <= 1:
                start = p
                break
        path = np.array(_order(S, start or tuple(pts[0])))
        w = width[path[:, 0], path[:, 1]]
        d = dE[path[:, 0], path[:, 1]]
        branches.append(dict(path=path, width=w, depth=d, length_um=float(len(path) * GRID * 1.1)))
        if len(path) >= 2 * win // 3:
            # reference = the ribbon on either side (flank_in..flank_out um away), so a
            # long defect does not pull its own baseline toward itself
            fi = int(round(prm["flank_in_um"] / (GRID * 1.1)))
            fo = int(round(prm["flank_out_um"] / (GRID * 1.1)))
            ref = np.empty(len(w))
            for j in range(len(w)):
                fl = np.concatenate([w[max(0, j - fo):max(0, j - fi)], w[j + fi:j + fo]])
                ref[j] = np.median(fl) if fl.size >= 5 else np.median(w)
            for kind, test in (("thin", w < prm["thin"] * ref), ("bulge", w > prm["bulge"] * ref)):
                lab_, n_ = ndi.label(test)
                for r_ in range(1, n_ + 1):
                    idx = np.nonzero(lab_ == r_)[0]
                    if len(idx) >= minf:
                        flags.append(dict(kind=kind, branch=len(branches) - 1, pts=path[idx],
                                          length_um=float(len(idx) * GRID * 1.1),
                                          ratio=float(np.median(w[idx] / ref[idx]))))
    ends = []
    nd = int(round(prm["dir_um"] / GRID))
    for e in map(tuple, np.argwhere(ends_mask)):
        b = next((br for br in branches if tuple(br["path"][0]) == e or tuple(br["path"][-1]) == e), None)
        if b is None:
            continue
        p = b["path"] if tuple(b["path"][0]) == e else b["path"][::-1]
        back = p[min(nd, len(p) - 1)]
        v = np.array(e, float) - back
        nv = np.hypot(*v)
        # width and depth of the ribbon behind the end (the tip itself tapers)
        a, b_ = int(round(prm["flank_from_um"] / (GRID * 1.1))), int(round(prm["flank_to_um"] / (GRID * 1.1)))
        seg = p[min(a, len(p) - 1):min(b_, len(p))]
        wf = float(np.median(width[seg[:, 0], seg[:, 1]])) if len(seg) else float(width[e])
        df = float(np.median(dE[seg[:, 0], seg[:, 1]])) if len(seg) else float(dE[e])
        ends.append(dict(pos=np.array(e), dir=v / nv if nv else v, width=wf, depth=df,
                         tip_width=float(width[e]), tip_depth=float(dE[e])))
    return dict(igl=igl, skeleton=sk, width=width, depth_egl=dE, branches=branches, ends=ends, flags=flags,
                length_mm=float(sum(b["length_um"] for b in branches) / 1000),
                median_width_um=float(np.median(np.concatenate([b["width"] for b in branches]))) if branches else 0.0)
