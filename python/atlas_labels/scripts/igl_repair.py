"""
IGL gap repair from the ribbon itself (step 2 + 4 of the IGL repair).

    from igl_repair import repair_gaps
    new_igl, bridges = repair_gaps(g)        # g = igl_ribbon.grid_labels(labels, px)

For every end of the IGL ribbon (igl_ribbon.describe):
  1. Walk forward from the end, steering along the ribbon's own depth below the EGL:
     direction = mix of the current direction, the tangent of the iso-depth line,
     and a pull back toward the end's depth. Step 4 um, at most max_gap_um.
  2. The walk must leave the IGL, cross non-IGL tissue (ML, PCL, unassigned, or
     anomaly = bleach), and enter IGL again. It fails on white matter, EGL or
     background, so natural ends (the base, the core between the anterior and
     posterior halves) stay open.
  3. Accepted only if the walk stays within 50 % (or 25 um) of the ribbon's depth, never
     crosses a tear (MATLAB tear mask), and - when it runs mostly through ML / PCL - its
     NeuN is IGL-like (median above the geometric mean of the section's IGL and ML NeuN):
     NeuN-dim ML between IGL pieces is real anatomy, NeuN-bright "ML" is mislabelled IGL.
  4. Extra start points: ends found with lighter pruning (1x width), which recovers the
     stub side of a gap; a walk from one of these is kept only if the two sides are not
     already joined through the IGL by a route shorter than 5x the bridge.
  5. Bridge painted as the band the ribbon occupies: between the 5th and 95th percentile
     depth below the EGL of the IGL on each flank (interpolated along the bridge), within
     reach of the walk, plus the tube of the flank width; never over EGL or white matter.
"""
import numpy as np
from scipy import ndimage as ndi

import igl_ribbon as RB

GRID = RB.GRID
DEF = dict(max_gap_um=500.0, extra_prune_factor=1.0, route_factor=5.0, band_near_um=60.0, band_slack_um=8.0, mix_dir=0.5, mix_tan=0.5, pull=0.6, through_anomaly=True, min_width_um=20.0,
           max_depth_dev=0.5, min_depth_dev_um=25.0, max_ml_frac_without_neun=0.5)


def _bilin(a, y, x):
    return ndi.map_coordinates(a, [[y], [x]], order=1, mode="nearest")[0]


def walk(start, v, depth, dE, gy, gx, ok, igl_c, egl, wm, prm, bleach=None):
    """returns (status, points) with status 'joined', 'blocked:<why>' or 'too_long'."""
    H, W = dE.shape
    p = start.astype(float)
    v = v / max(np.hypot(*v), 1e-9)
    pts, left = [], False
    n_max = int(prm["max_gap_um"] / GRID)
    for k in range(n_max):
        g = np.array([_bilin(gy, *p), _bilin(gx, *p)])
        ng = np.hypot(*g)
        if ng > 1e-6:
            g = g / ng
            t = np.array([-g[1], g[0]])
            if t @ v < 0:
                t = -t
            err = (depth - _bilin(dE, *p)) / max(depth, 8.0)
            nd = prm["mix_dir"] * v + prm["mix_tan"] * t + prm["pull"] * err * g
            v = nd / max(np.hypot(*nd), 1e-9)
        p = p + v
        iy, ix = int(round(p[0])), int(round(p[1]))
        if not (0 <= iy < H and 0 <= ix < W):
            return "blocked:edge", pts
        if igl_c[iy, ix]:
            if left:
                return "joined", pts + [("v", float(v[0]), float(v[1]))]
            continue
        left = True
        if wm[iy, ix] and not (bleach is not None and bleach[iy, ix]):
            return "blocked:white matter", pts
        if egl[iy, ix]:
            return "blocked:EGL", pts
        if not ok[iy, ix]:
            return "blocked:background", pts
        pts.append((iy, ix))
    return "too_long", pts


def _band_depths(igl_c, dE, pts, near_px):
    """5th / 95th percentile depth below the EGL of IGL pixels within near_px of pts."""
    m = np.zeros(igl_c.shape, bool)
    pts = np.asarray(pts, int)
    m[pts[:, 0], pts[:, 1]] = True
    zone = (ndi.distance_transform_edt(~m) <= near_px) & igl_c
    if zone.sum() < 5:
        return None
    v = dE[zone]
    return float(np.percentile(v, 5)), float(np.percentile(v, 95))


def _route_ok(igl_c, a, b, bridge_um, factor):
    """the two sides must not already be close through the IGL (a spur looping back)."""
    from skimage.graph import MCP_Geometric
    cost = np.where(igl_c, 1.0, np.inf)
    ent = None
    y, x = b
    win = igl_c[max(y - 2, 0):y + 3, max(x - 2, 0):x + 3]
    ys, xs = np.nonzero(win)
    if not len(ys):
        return True, float("inf")
    ent = (int(ys[0] + max(y - 2, 0)), int(xs[0] + max(x - 2, 0)))
    cum, _ = MCP_Geometric(cost).find_costs([tuple(int(v) for v in a)], [ent])
    route = float(cum[ent] * GRID)
    return (not np.isfinite(route)) or route > factor * bridge_um, route


def repair_gaps(g, prm=None, rib=None, neun=None, tear=None, bleach=None):
    """g: grid masks; neun: NeuN image on the same grid (optional); tear: tear mask
    (walks never cross tears; other anomaly pixels, e.g. bleach, may be crossed);
    bleach: user-confirmed bleached regions or breaches - inside them walks may cross
    white-matter labels and the repair may paint over them; a walk lying at least half
    inside them skips the depth, NeuN and loop-back checks."""
    prm = {**DEF, **(prm or {})}
    rib = rib or RB.describe(g)
    t_neun = None
    if neun is not None:
        t_neun = float(np.sqrt(np.median(neun[g["igl"]]) * np.median(neun[g["ml"]])))
    igl_c = rib["igl"]
    dE = ndi.gaussian_filter(rib["depth_egl"], 1.0)
    gy, gx = np.gradient(dE)
    anom_ok = g["anom"] & ~tear if tear is not None else g["anom"]
    ok = g["tissue"] | (anom_ok if prm["through_anomaly"] else False)
    width = rib["width"]
    new = g["igl"].copy()
    paintable = ok & ~g["egl"] & (~g["wm"] | (bleach if bleach is not None else False))
    ok = ok | (bleach & g["wm"] if bleach is not None else False)
    # start points: ends of the ribbon, plus ends found with lighter pruning (gap stubs)
    starts = [(e, False) for e in rib["ends"]]
    if prm["extra_prune_factor"]:
        aux = RB.describe(g, {**RB.P, "prune_factor": prm["extra_prune_factor"]})
        for e in aux["ends"]:
            if all(np.hypot(*(e["pos"] - q["pos"])) * GRID > 20 for q in rib["ends"]):
                starts.append((e, True))
    bridges = []
    for e, extra in starts:
        status, pts = walk(np.array(e["pos"], float), e["dir"].astype(float), max(e["depth"], 8.0), dE, gy, gx, ok,
                           igl_c, g["egl"], g["wm"], prm, bleach)
        if status == "joined":
            pts = pts[:-1]
        b = dict(end=e["pos"].tolist(), status=status, extra_start=extra, length_um=round(len(pts) * GRID, 1))
        if status == "joined" and pts:
            P = np.array(pts)
            tgt = P[-1]
            y0, x0 = max(tgt[0] - 4, 0), max(tgt[1] - 4, 0)
            wt = max(float(np.max(width[y0:tgt[0] + 5, x0:tgt[1] + 5])), prm["min_width_um"])
            w0 = max(e["width"], prm["min_width_um"])
            dpath = dE[P[:, 0], P[:, 1]]
            dev = float(np.max(np.abs(dpath - e["depth"])))
            ml_frac = float(g["ml"][P[:, 0], P[:, 1]].mean())
            b.update(depth_dev_um=round(dev, 1), ml_frac=round(ml_frac, 2))
            reasons = []
            in_bleach = float(bleach[P[:, 0], P[:, 1]].mean()) if bleach is not None else 0.0
            b["bleach_frac"] = round(in_bleach, 2)
            confirmed = in_bleach >= 0.5      # user-confirmed breach: trust the user over the checks
            if dev > max(prm["max_depth_dev"] * e["depth"], prm["min_depth_dev_um"]) and not confirmed:
                reasons.append("left the ribbon depth")
            if ml_frac > prm["max_ml_frac_without_neun"] and in_bleach < 0.5:
                if neun is None:
                    reasons.append("through ML, no NeuN image to check")
                else:
                    nm = float(np.median(neun[P[:, 0], P[:, 1]]))
                    b["neun_ratio"] = round(nm / t_neun, 2)
                    if nm < t_neun:
                        reasons.append("through NeuN-dim ML (real ML)")
            if extra and not reasons and not confirmed:
                okr, route = _route_ok(igl_c, e["pos"], tgt, len(P) * GRID, prm["route_factor"])
                b["igl_route_um"] = round(route, 1) if np.isfinite(route) else None
                if not okr:
                    reasons.append("loops back onto its own ribbon")
            if reasons:
                b["status"] = "rejected: " + ", ".join(reasons)
                bridges.append(b)
                continue
            # band: depths of the ribbon on both flanks, interpolated along the bridge
            near = prm["band_near_um"] / GRID
            d_start = _band_depths(igl_c, dE, [e["pos"]], near)
            d_end = _band_depths(igl_c, dE, [tgt], near)
            d_start = d_start or (e["depth"] - w0 / 2, e["depth"] + w0 / 2)
            d_end = d_end or (dpath[-1] - wt / 2, dpath[-1] + wt / 2)
            f = np.linspace(0, 1, len(P))
            lo = (1 - f) * d_start[0] + f * d_end[0]
            hi = (1 - f) * d_start[1] + f * d_end[1]
            ws = (1 - f) * w0 + f * wt
            line = np.zeros(new.shape, bool)
            line[P[:, 0], P[:, 1]] = True
            d, (iy, ix) = ndi.distance_transform_edt(~line, return_indices=True)
            idx = np.full(new.shape, -1)
            idx[P[:, 0], P[:, 1]] = np.arange(len(P))
            k = idx[iy, ix]
            reach = (d * GRID <= ws[k] / 2 + prm["band_slack_um"]) & (k >= 0)
            tol = prm["band_slack_um"]
            inband = reach & (dE >= lo[k] - tol) & (dE <= hi[k] + tol)
            tube = reach & (d * GRID <= ws[k] / 2)
            paint = (inband | tube) & paintable
            new |= paint
            b.update(painted_um2=round(float(paint.sum() * GRID ** 2), 1), width_um=[round(w0), round(wt)],
                     band_um=[[round(v, 1) for v in d_start], [round(v, 1) for v in d_end]])
        bridges.append(b)
    return new, bridges, rib
