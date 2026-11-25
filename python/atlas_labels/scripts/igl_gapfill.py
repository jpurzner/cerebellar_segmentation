"""
IGL gap filling for the repair pipeline: the combined outline / ribbon-model filler
(benchmarked on synthetic gaps in igl_outline_eval.py) in place of the walk (igl_repair.py).

    from igl_gapfill import fill_gaps
    new_igl, bridges, rib = fill_gaps(g, neun=neun, tear=tear, bleach=conf)
    # same call and return as igl_repair.repair_gaps; accepted bridges have status "joined"

1. Ends: the ribbon's ends (igl_ribbon.describe) plus ends found with lighter pruning
   (1x width), which recovers the stub side of a gap.
2. Far side (search="hybrid", 2026-11-19, after the user's review: the bridges they marked
   wrong jumped to the nearest IGL and "turn at right angles"; the right one "connects to the
   deeper part"):
   a) the predicted ribbon (predict): from where the IGL actually stops in front of the end
      (the skeleton end lies up to ~100 um short of it), the centre line is carried on at the
      depth below the EGL the ribbon ran at, turning as the outline turns; white-matter labels
      in the shallow half of the IGL band do not stop it, deeper white matter, EGL, background
      and tears do. The IGL it runs into is the far side only if it continues the ribbon
      (_continues: an end facing back or at the predicted depth; a stretch of ribbon at the
      predicted depth running the same way within 50 degrees, or any direction over a gap
      shorter than the ribbon is wide); other IGL it passes over. The fill follows the course.
   b) for ends the prediction cannot join, the earlier searches: a geodesic search outward
      through the expected IGL depth band (m - 8 to m + 1.5 w + 8 um below the surface, outline
      core excluded, white-matter labels blocking only deeper than m + w/2) to the first IGL
      within 500 um that is not the end's own ribbon (IGL route from the end longer than 3x
      the distance, or none), and the walk of igl_repair (steered along the ribbon's depth,
      stopped by white matter). They find the short crossings of white-matter labels inside a
      wide ribbon that a forward course misses (user: s4_C-1, s4_C-3, s5_C-1 correct).
   Far side = the end there if one is within 60 um, else the nearest centre-line point.
   c) main_axis / axis_target (user, marker D4: "find the main axis of the IGL above and below,
      penalize changes in angle from this axis"): the IGL ahead within 60 degrees of the axis
      of the last 30-200 um of ribbon, turning least away from both axes (each at most 45
      degrees). Off by default (axis_search): on the 30 judged fills it changes nothing, on
      synthetic gaps it adds false fills (307 -> 481 um2 per section). A 45-degree axis check
      on the old searches' long bridges (axis_check) is off too: it cost 12 points of
      synthetic restoration (80 -> 68 %) and changed no judged fill. Both are kept for the
      marker review (igl_finalize draws each end's axis and best candidates).
   Predicted pairs are taken before the others (step 5). Scored against the user's verdicts
   (igl_gapfill_score.py): prediction alone 16/21, earlier searches + the checks below 18/21,
   hybrid 18/21 (the version the user judged: 12/21); with the second review's 9 bridges
   added, hybrid 26/30, earlier searches + checks 25/30, prediction alone 22/30.
3. Per pair, with m = ML thickness at the two ends (depth below the EGL minus half the
   ribbon width, median over 40-200 um behind each end) and w = mean flank width:
   outline test - the band m <= depth <= m + w joining the ends, drawn without white
                  matter and minus the outline core; wm_shallow = share of its shallow half
                  (depth <= m + w/2) labelled white matter. Real white matter lies deeper
                  than the IGL band, so above 0.35 the labels are taken to be wrong
                  (synthetic test: fires on 96 % of gaps whose IGL was labelled white matter,
                  4 % of plain gaps).
   labels trusted - IGL rule: tissue at least m below the surface (EGL, or background near
                  the EGL) that is not white matter (labels; inside user-confirmed bleach
                  regions white-matter labels are ignored). Fill = the rule's IGL in the
                  geodesic corridor between the ends; if the ends do not connect inside it,
                  the rule's IGL in a tube around the ribbon model's centre line.
   labels wrong - fill = the ribbon model's band (Hermite width between the flanks) with
                  white-matter labels ignored (EGL, background and tears still barriers).
4. Checks; a pair is rejected if any fails:
   - the ribbon model's centre line is at most 600 um long
   - the fill and the IGL together join the two ends, and the fill covers at least 30 %
     of chord x flank width
   - tube / white-matter-ignored fills: the centre line stays within 50 % (or 25 um) of
     the ends' depth below the EGL
   - white-matter-ignored fills: at most 30 % inside the outline's white-matter core (the
     IGL cannot lie where the shape puts the core), and the two ends' depths below the EGL
     differ by at most 50 um (the IGL's depth varies by about +-18 um); on the s1_2
     automatic labels these stop a bridge across the real white matter of the base
   - fills that are mostly ML / PCL: NeuN must be IGL-like (median above the geometric
     mean of the section's IGL and ML NeuN), as in igl_repair: NeuN-dim ML between IGL
     pieces is real anatomy. Skipped inside user-confirmed bleach.
   - (user's cKO notes, 2026-11-19) every fill: NeuN at least half that mean ("the IGL is
     super diffuse"), and neither side's IGL may touch the EGL (ML under 15 um: "the fuzzy
     inner EGL causes these bridges"). Both skipped inside user-confirmed bleach.
5. Accepted pairs are taken best first (predicted pairs first, then ribbon-model energy +
   length / 300 um, as in igl_model_apply); a pair is dropped if one of its ends is already used or less than 30 %
   of its fill is new (the same gap reached from its other side); a pair that overlaps a
   taken fill but extends it (e.g. a short bridge to a fragment inside a longer gap) is kept.
6. Layout check (user review 2026-11-19: s4_A bridge 1 "cuts off the deep white matter and
   creates a circular piece"), against the IGL and white matter as they stand when the pair
   comes up, so it applies everywhere, user-confirmed bleach included. The white matter is a
   tree (every folium's core joins the deep white matter through its stalk) and the IGL an
   open ribbon along it, so a real repair never
   - cuts white matter off the rest of its white-matter piece: rejected if the pieces cut off
     add up to more than 4000 um2 (thin slivers of mislabelled IGL left beside a fill stay
     below that; the user's correct bridges cut off at most 900 um2, the wrong one 45,800),
   - closes the IGL into a ring around more than 4000 um2 of tissue that is mostly not white
     matter (a fill across a fissure closes off its bottom), or around white matter that was
     joined to white matter outside the ring (a fill across a folium's stalk, or across the
     base where the white matter leaves the cerebellum: it closes the IGL around the whole
     arbor). A ring around a white-matter island (a folium cut obliquely) is allowed.
   With the predicted course (2026-11-19) only cut-off white matter deeper than the ribbon
   band counts (shallower is mislabelled IGL; user: s2_2 "should be joined"), and cut-off
   white matter inside user-confirmed regions does not (user: s4_A "right to fix").
"""
import numpy as np
from scipy import ndimage as ndi
from skimage.graph import MCP_Geometric

import igl_model as M
import igl_model_eval as ME
import igl_outline_fill as OF
import igl_repair as IR
import igl_ribbon as RB

G = RB.GRID
DEF = dict(max_gap_um=500.0, max_curve_um=600.0, min_gap_um=20.0, facing=-0.17, route_factor=3.0, extra_prune_factor=1.0,
           wm_shallow=0.35, min_cover=0.3, min_width_um=20.0, max_depth_dev=0.5, min_depth_dev_um=25.0,
           max_ml_frac_without_neun=0.5, len_scale_um=300.0, max_core_frac=0.3, max_depth_gap_um=50.0,
           max_wm_cut_um2=4000.0, max_ring_um2=4000.0, island_wm=0.5,
           search="hybrid", mix_dir=0.5, mix_tan=0.5, pull=0.6, min_ml_um=15.0, min_neun_ratio=0.5,
           contact_depth_um=30.0, contact_cos=0.64, fill="corridor", axis_cone=60.0, max_axis_turn=45.0,
           axis_search=False, axis_check=False)
MODEL = dict(sigma_depth_um=15.0, bend_scale=0.3)     # best ribbon-model setting (igl_model_eval)


def _route_from(igl_c, a, r_px):
    """geodesic distance (px) through the IGL from end a, inside a window of radius r_px."""
    H, W = igl_c.shape
    y0, x0 = max(int(a[0]) - r_px, 0), max(int(a[1]) - r_px, 0)
    y1, x1 = min(int(a[0]) + r_px + 1, H), min(int(a[1]) + r_px + 1, W)
    sub = igl_c[y0:y1, x0:x1]
    out = np.full(igl_c.shape, np.inf)
    ys, xs = np.nonzero(sub)
    if not len(ys):
        return out
    k = np.argmin((ys + y0 - a[0]) ** 2 + (xs + x0 - a[1]) ** 2)
    cum, _ = MCP_Geometric(np.where(sub, 1.0, np.inf)).find_costs([(int(ys[k]), int(xs[k]))])
    out[y0:y1, x0:x1] = cum
    return out


def _joins(igl_c, fill, a, b, pad):
    y0 = max(int(min(a[0], b[0])) - pad, 0)
    x0 = max(int(min(a[1], b[1])) - pad, 0)
    y1 = min(int(max(a[0], b[0])) + pad + 1, fill.shape[0])
    x1 = min(int(max(a[1], b[1])) + pad + 1, fill.shape[1])
    lab, _ = ndi.label((igl_c | fill)[y0:y1, x0:x1], structure=np.ones((3, 3), bool))
    la = lab[int(round(a[0])) - y0, int(round(a[1])) - x0]
    lb = lab[int(round(b[0])) - y0, int(round(b[1])) - x0]
    return bool(la and la == lb)


S8 = np.ones((3, 3), bool)


def _split(mask, fill, conn=None):
    """pieces of every connected piece of mask that fill touches, once fill is taken out:
    returns (piece labels, piece sizes px, list of piece-index arrays, one per touched piece)."""
    lab, _ = ndi.label(mask, structure=conn)
    hit = np.unique(lab[fill & mask])
    hit = hit[hit > 0]
    if not hit.size:
        return None, None, []
    lab2, n2 = ndi.label(mask & ~fill & np.isin(lab, hit), structure=conn)
    if not n2:
        return lab2, np.zeros(0, int), []
    sizes = np.bincount(lab2.ravel(), minlength=n2 + 1)[1:]
    parent = np.asarray(ndi.maximum(lab, lab2, index=np.arange(1, n2 + 1)), int)
    return lab2, sizes, [np.nonzero(parent == h)[0] for h in hit]


def topology(igl, wm, fill, deep=None):
    """what a fill does to the layout: white matter it cuts off the rest of its piece (um2,
    every resulting piece but the largest; with deep given, only the cut-off pixels inside
    deep count, as wm_cut_deep_um2), and the largest region of non-IGL tissue it closes off
    from the rest (um2) with that region's white-matter share."""
    out = dict(wm_cut_um2=0.0, ring_um2=0.0, ring_wm=0.0)
    lab2, sizes, groups = _split(wm, fill, S8)
    cut = [j for i in groups if i.size > 1 for j in i[np.argsort(sizes[i])[::-1][1:]]]
    out["wm_cut_um2"] = round(float(sum(sizes[j] for j in cut)) * G * G)
    if deep is not None:
        cm = np.isin(lab2, np.asarray(cut) + 1) if cut else np.zeros(wm.shape, bool)
        out["wm_cut_deep_um2"] = round(float((cm & deep).sum()) * G * G)
    lab2, sizes, groups = _split(~igl, fill)               # 4-connected: the IGL is the barrier
    best = None
    for i in groups:
        if i.size < 2:
            continue
        for j in i[np.argsort(sizes[i])[::-1][1:]]:
            if sizes[j] * G * G > out["ring_um2"]:
                out["ring_um2"] = round(float(sizes[j] * G * G))
                best = j
    if best is not None:
        P_ = lab2 == best + 1
        out["ring_wm"] = round(float(wm[P_].mean()), 2)
        # does the ring's white matter join white matter outside it (through the fill)?
        wl, _ = ndi.label(wm, structure=S8)
        inside = set(np.unique(wl[P_ & wm])) - {0}
        outside = set(np.unique(wl[~P_ & ~fill & wm])) - {0}
        out["ring_severs_wm"] = bool(inside & outside)
    return out


def layout_reasons(t, prm, confirmed=False):
    r = []
    cut = t.get("wm_cut_deep_um2", t["wm_cut_um2"])
    if cut > prm["max_wm_cut_um2"] and not confirmed:
        r.append("cuts white matter off (%d um2)" % cut)
    if t["ring_um2"] > prm["max_ring_um2"] and (t["ring_wm"] < prm["island_wm"] or t.get("ring_severs_wm")):
        r.append("closes the IGL into a ring (%d um2)" % t["ring_um2"])
    return r


def _crop(mask, pad=2):
    ys, xs = np.nonzero(mask)
    if not len(ys):
        return None
    y0, x0 = max(ys.min() - pad, 0), max(xs.min() - pad, 0)
    y1, x1 = ys.max() + pad + 1, xs.max() + pad + 1
    return [int(y0), int(y1), int(x0), int(x1)], mask[y0:y1, x0:x1].copy()


def section_setup(g, tear=None, bleach=None, rib=None, prm=None):
    """everything per section that the pairs share (fits, geometry, fields)."""
    prm = {**DEF, **(prm or {})}
    rib = rib or RB.describe(g)
    shape = g["igl"].shape
    tear = np.zeros(shape, bool) if tear is None else tear
    bleach = np.zeros(shape, bool) if bleach is None else bleach
    ok = (g["tissue"] | (g["anom"] & ~tear)) & ~tear
    wm_eff = g["wm"] & ~bleach
    geo = OF.geometry(g)
    allowed = geo["inner"] & ~tear
    igl_c = rib["igl"]
    fit_label = OF.fit(geo, igl_c, wm_eff, "label")
    fit_outline = OF.fit(geo, igl_c, wm_eff, "outline")
    dE = rib["depth_egl"]
    return dict(g=g, rib=rib, prm=prm, tear=tear, bleach=bleach, ok=ok, wm_eff=wm_eff, geo=geo, allowed=allowed,
                igl_c=igl_c, dE=dE, fit_label=fit_label, fit_outline=fit_outline,
                core=OF.wm_outline(geo, fit_outline),
                field=M.Field(dE, g["egl"] | wm_eff | ~ok),
                field_nowm=M.Field(dE, g["egl"] | ~ok))


def ends_of(S):
    prm, g, rib = S["prm"], S["g"], S["rib"]
    ends = [dict(e, extra=False) for e in rib["ends"]]
    if prm["extra_prune_factor"]:
        aux = RB.describe(g, {**RB.P, "prune_factor": prm["extra_prune_factor"]})
        for e in aux["ends"]:
            if all(np.hypot(*(e["pos"] - q["pos"])) * G > 20 for q in rib["ends"]):
                ends.append(dict(e, extra=True))
    for i, e in enumerate(ends):
        e["id"] = i
        e["width"] = max(float(e["width"]), prm["min_width_um"])
    return ends


def skeleton_points(S):
    """centre-line points with local width, depth and tangent (for targets inside a ribbon)."""
    rib = S["rib"]
    pts, W, D, T = [], [], [], []
    for b in rib["branches"]:
        P_ = b["path"].astype(float)
        if len(P_) < 3:
            continue
        Sm = np.stack([ndi.gaussian_filter1d(P_[:, k], 2.0, mode="nearest") for k in range(2)], 1)
        Tg = np.gradient(Sm, axis=0)
        Tg /= np.maximum(np.hypot(*Tg.T), 1e-9)[:, None]
        pts.append(b["path"])
        W.append(b["width"])
        D.append(b["depth"])
        T.append(Tg)
    if not pts:
        return None
    return dict(pos=np.vstack(pts), width=np.concatenate(W), depth=np.concatenate(D), tan=np.vstack(T))


def tip(S, e, max_um=150.0):
    """where the IGL actually stops in front of end e: march from the skeleton end along its
    direction while inside the IGL (pruning leaves the skeleton end up to ~100 um short of it)."""
    igl_c = S["igl_c"]
    p = np.asarray(e["pos"], float)
    v = np.asarray(e["dir"], float)
    v = v / max(np.hypot(*v), 1e-9)
    last = p.copy()
    for _ in range(int(max_um / G * 2)):
        p = p + 0.5 * v
        iy, ix = int(round(p[0])), int(round(p[1]))
        if not (0 <= iy < igl_c.shape[0] and 0 <= ix < igl_c.shape[1]) or not igl_c[iy, ix]:
            break
        last = np.array([iy, ix], float)
    return last


def discover(S, e, ends):
    """the far side of the gap in front of end e: search outward through the expected IGL
    depth band (no white-matter labels, outline core excluded) and take the first IGL that is
    not e's own ribbon (IGL route from e longer than route_factor x the distance, or none).
    Returns a target end-like dict or None."""
    prm, geo = S["prm"], S["geo"]
    igl_c = S["igl_c"]
    a = np.asarray(e["fill_pos"], float)
    m = max(e["depth"] - e["width"] / 2, 8.0)
    R = int((prm["max_gap_um"] + 100) / G)
    H, W = igl_c.shape
    y0, x0 = max(int(a[0]) - R, 0), max(int(a[1]) - R, 0)
    y1, x1 = min(int(a[0]) + R + 1, H), min(int(a[1]) + R + 1, W)
    win = (slice(y0, y1), slice(x0, x1))
    dS = geo["dS"][win]
    band = S["allowed"][win] & ~S["core"][win] & (dS >= m - 8) & (dS <= m + 1.5 * e["width"] + 8)
    # white-matter labels block the search except in the shallow half of the IGL band, where
    # real white matter almost never is (mislabelled IGL is passable, a folium's core is not)
    band &= ~(S["wm_eff"][win] & (dS > m + e["width"] / 2))
    iglw = igl_c[win]
    cand = band & ~ndi.binary_dilation(iglw, iterations=1)
    seeds = OF._seeds(cand, a - [y0, x0], (e["width"] / 2 + 12.0) / G)
    if len(seeds) == 0:
        return None
    D = OF._geodesic(cand, seeds)
    Dadj = ndi.minimum_filter(np.where(np.isfinite(D), D, 1e9), size=5)
    reach = iglw & (Dadj * G <= prm["max_gap_um"])
    if not reach.any():
        return None
    route = S["_route"].get(e["id"])
    if route is None:
        route = _route_from(igl_c, np.asarray(e["pos"], float), int(prm["route_factor"] * prm["max_gap_um"] / G) + 4)
        S["_route"][e["id"]] = route
    ys, xs = np.nonzero(reach)
    d_e = np.hypot(ys + y0 - a[0], xs + x0 - a[1])
    rr = route[ys + y0, xs + x0]
    valid = ~np.isfinite(rr) | (rr >= prm["route_factor"] * np.maximum(Dadj[ys, xs], d_e))
    valid &= d_e * G >= prm["min_gap_um"]
    if not valid.any():
        return None
    k = np.flatnonzero(valid)[np.argmin(Dadj[ys[valid], xs[valid]])]
    q = np.array([ys[k] + y0, xs[k] + x0], float)
    return target_at(S, q, a, e, ends, float(Dadj[ys[k], xs[k]] * G), "band search")


def target_at(S, q, a, e, ends, corridor_um, found_by):
    """far-side description at IGL pixel q: the end there if one is within 60 um, else the
    nearest centre-line point (width, depth and direction from the ribbon around it)."""
    prm = S["prm"]
    near_end = [f for f in ends if f["id"] != e["id"] and np.hypot(*(f["pos"] - q)) * G <= 60]
    if near_end:
        f = min(near_end, key=lambda f: np.hypot(*(f["pos"] - q)))
        return dict(pos=np.asarray(f["pos"]), fill_pos=q, dir=np.asarray(f["dir"], float), width=f["width"], depth=f["depth"],
                    id=f["id"], extra=f["extra"], kind="end", corridor_um=corridor_um, found_by=found_by)
    sk = S["_sk"]
    if sk is None:
        return None
    dd = np.hypot(*(sk["pos"] - q).T)
    j = int(np.argmin(dd))
    if dd[j] * G > 80:
        return None
    near = dd <= max(60.0 / G, dd[j] + 1.0)
    t = sk["tan"][j]
    if t @ (sk["pos"][j] - a) < 0:
        t = -t
    return dict(pos=sk["pos"][j], fill_pos=q, dir=-t, width=max(float(np.median(sk["width"][near])), prm["min_width_um"]),
                depth=float(np.median(sk["depth"][near])), id=None, extra=False, kind="ribbon",
                corridor_um=corridor_um, found_by=found_by)


def walk_target(S, e, ends):
    """the far side as the walk of igl_repair finds it (steered along the ribbon's depth,
    stopped by white matter, EGL and background): where it re-enters the IGL."""
    g, igl_c = S["g"], S["igl_c"]
    if "_walkf" not in S:
        dEs = ndi.gaussian_filter(S["dE"], 1.0)
        gy, gx = np.gradient(dEs)
        S["_walkf"] = (dEs, gy, gx)
    dEs, gy, gx = S["_walkf"]
    status, pts = IR.walk(np.array(e["pos"], float), np.asarray(e["dir"], float), max(e["depth"], 8.0), dEs, gy, gx,
                          S["ok"] | (S["bleach"] & g["wm"]), igl_c, g["egl"], g["wm"], IR.DEF, S["bleach"])
    if status != "joined":
        return None
    P = [p for p in pts if p and p[0] != "v"]
    if not P:
        return None
    last = np.array(P[-1], float)
    v = np.array(pts[-1][1:], float)
    q = np.round(last + v).astype(int)
    q = np.clip(q, 0, np.array(igl_c.shape) - 1)
    if not igl_c[q[0], q[1]]:
        ys, xs = np.nonzero(igl_c[max(q[0] - 2, 0):q[0] + 3, max(q[1] - 2, 0):q[1] + 3])
        if not len(ys):
            return None
        q = np.array([ys[0] + max(q[0] - 2, 0), xs[0] + max(q[1] - 2, 0)])
    return target_at(S, q.astype(float), np.asarray(e["fill_pos"], float), e, ends, len(P) * G, "walk")


def predict(S, e, ends):
    """the ribbon's predicted continuation from end e (2026-11-19, after the user's review:
    bridges to the nearest IGL "turn at right angles"; the IGL "should connect to the deeper
    part"). From where the IGL stops, carry the centre line forward at the depth below the EGL
    the ribbon ran at, turning as the outline turns (the depth field's level line) and not by
    more than that; white-matter labels do not stop it (mislabelled or bleached IGL lies on
    the course), EGL, background and tears do. The far side is the first IGL of another
    stretch (IGL route from e longer than route_factor x the distance, or none) it runs into.
    Returns (target, path) or None."""
    prm, g, igl_c = S["prm"], S["g"], S["igl_c"]
    if "_walkf" not in S:
        dEs = ndi.gaussian_filter(S["dE"], 1.0)
        gy, gx = np.gradient(dEs)
        S["_walkf"] = (dEs, gy, gx)
    dEs, gy, gx = S["_walkf"]
    if "_stop" not in S:
        S["_stop"] = g["egl"] | ~S["ok"] | S["tear"]
    stop = S["_stop"]
    route = S["_route"].get(e["id"])
    if route is None:
        route = _route_from(igl_c, np.asarray(e["pos"], float), int(prm["route_factor"] * prm["max_gap_um"] / G) + 4)
        S["_route"][e["id"]] = route
    H, W = igl_c.shape
    a = np.asarray(e["fill_pos"], float)
    p = a.copy()
    v = np.asarray(e["dir"], float)
    v = v / max(np.hypot(*v), 1e-9)
    depth = max(float(e["depth"]), 8.0)
    # real white matter (deeper than the shallow half of the IGL band) stops the course, as in
    # discover: without this it runs out of the cerebellum at the base and joins the two
    # natural ends of the ribbon there (synthetic test, 2026-11-19)
    m_e = max(depth - e["width"] / 2, 8.0)
    wm_block = S["wm_eff"] & (S["geo"]["dS"] > m_e + e["width"] / 2)
    if "_cc" not in S:
        S["_cc"] = ndi.label(igl_c, structure=S8)[0]
    pts, left, passed = [], False, set()
    for _ in range(int(prm["max_gap_um"] / G)):
        gr = np.array([IR._bilin(gy, *p), IR._bilin(gx, *p)])
        ng = np.hypot(*gr)
        if ng > 1e-6:
            gr = gr / ng
            t = np.array([-gr[1], gr[0]])
            if t @ v < 0:
                t = -t
            err = (depth - IR._bilin(dEs, *p)) / depth
            nd = prm["mix_dir"] * v + prm["mix_tan"] * t + prm["pull"] * err * gr
            v = nd / max(np.hypot(*nd), 1e-9)
        p = p + v
        iy, ix = int(round(p[0])), int(round(p[1]))
        if not (0 <= iy < H and 0 <= ix < W):
            return None
        if igl_c[iy, ix]:
            if not left:
                continue
            d_e = float(np.hypot(iy - a[0], ix - a[1]))
            rr = route[iy, ix]
            if np.isfinite(rr) and rr < prm["route_factor"] * max(len(pts), d_e):
                continue                          # its own ribbon again: keep going
            if d_e * G < prm["min_gap_um"] and len(pts) * G < prm["min_gap_um"]:
                return None
            lab = S["_cc"][iy, ix]
            if lab in passed:
                continue
            tg = target_at(S, np.array([iy, ix], float), a, e, ends, len(pts) * G, "prediction")
            if "_trace" in S and tg is not None:
                d_ = np.asarray(tg["dir"], float) / max(np.hypot(*tg["dir"]), 1e-9)
                S["_trace"].append(dict(end=e["id"], at=[iy, ix], kind=tg["kind"], depth=round(float(tg["depth"])), own_depth=round(depth),
                                        width=round(float(tg["width"])), own_width=round(float(e["width"])), cos=round(float(d_ @ v), 2),
                                        path_um=round(len(pts) * G), ok=_continues(tg, v, depth, e["width"], len(pts) * G, prm),
                                        _path=np.array(pts, float).reshape(-1, 2)))
            if tg is not None and _continues(tg, v, depth, e["width"], len(pts) * G, prm):
                return tg, np.array(pts, float).reshape(-1, 2)
            passed.add(lab)                       # skims another ribbon: carry on past it
            continue
        left = True
        if stop[iy, ix] or wm_block[iy, ix]:
            if "_trace" in S:
                S["_trace"].append(dict(end=e["id"], stopped_at=[iy, ix], by="white matter" if wm_block[iy, ix] else "EGL" if g["egl"][iy, ix] else "background / tear",
                                        path_um=round(len(pts) * G), _path=np.array(pts, float).reshape(-1, 2)))
            return None
        pts.append(p.copy())
    if "_trace" in S:
        S["_trace"].append(dict(end=e["id"], stopped_at=[int(p[0]), int(p[1])], by="length", path_um=round(len(pts) * G),
                                _path=np.array(pts, float).reshape(-1, 2)))
    return None


def main_axis(S, p, lo_um=30.0, hi_um=200.0, both=False):
    """the ribbon's main axis near IGL pixel p (user, s4_C marker D4: "find the main axis of the
    IGL above and below, penalize changes in angle from this axis"): principal direction of the
    IGL pixels lo..hi um from p along the IGL (0..hi on both sides with both=True), pointing
    from them towards p. Returns (unit vector, elongation) or (None, 0) if the IGL there is not
    elongated (elongation < 1.5)."""
    r = _route_from(S["igl_c"], np.asarray(p, float), int(hi_um / G) + 4) * G
    sel = np.isfinite(r) & (r >= (0.0 if both else lo_um)) & (r <= hi_um)
    ys, xs = np.nonzero(sel)
    if len(ys) < 15:
        return None, 0.0
    Q = np.stack([ys, xs], 1).astype(float)
    mu = Q.mean(0)
    _, sv, vt = np.linalg.svd(Q - mu, full_matrices=False)
    el = float(sv[0] / max(sv[1], 1e-9))
    if el < 1.5:
        return None, el
    u = vt[0]
    if u @ (np.asarray(p, float) - mu) < 0:
        u = -u
    return u, el


def _hermite_course(a, ua, b, ub, step=0.5):
    """cubic Hermite centre line from a (leaving along ua) to b (arriving along ub), tangent
    length = chord: the straightest course that keeps both main axes."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    L = float(np.hypot(*(b - a)))
    t = np.linspace(0, 1, max(int(L / step), 2))[:, None]
    h00, h10, h01, h11 = 2 * t ** 3 - 3 * t ** 2 + 1, t ** 3 - 2 * t ** 2 + t, -2 * t ** 3 + 3 * t ** 2, t ** 3 - t ** 2
    return h00 * a + h10 * L * np.asarray(ua) + h01 * b + h11 * L * np.asarray(ub)


def _angle(u, v, undirected=False):
    c = float(u @ v / max(np.hypot(*u) * np.hypot(*v), 1e-9))
    return float(np.degrees(np.arccos(np.clip(abs(c) if undirected else c, -1, 1))))


def axis_angles(S, a, b):
    """how far a straight bridge a -> b turns away from the main axis at each side (degrees)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    ua, _ = main_axis(S, a)
    ub, _ = main_axis(S, b, 0.0, 150.0, both=True)
    ch = b - a
    return (None if ua is None else _angle(ua, ch)), (None if ub is None else _angle(ub, ch, undirected=True)), ua, ub


def axis_target(S, e, ends):
    """far side along the ribbon's main axis: the IGL of another stretch within reach and within
    axis_cone degrees of the axis behind end e whose own axis the bridge turns least away from
    (cost = angle at e + angle at the far side + length / 10 um); both angles at most
    max_axis_turn. Returns a target dict (directions set from the two axes) or None."""
    prm, igl_c = S["prm"], S["igl_c"]
    a = np.asarray(e["fill_pos"], float)
    ua, _ = main_axis(S, a)
    if ua is None:
        return None
    if "_cc" not in S:
        S["_cc"] = ndi.label(igl_c, structure=S8)[0]
    route = S["_route"].get(e["id"])
    if route is None:
        route = _route_from(igl_c, np.asarray(e["pos"], float), int(prm["route_factor"] * prm["max_gap_um"] / G) + 4)
        S["_route"][e["id"]] = route
    R = int(prm["max_gap_um"] / G) + 1
    y0, x0 = max(int(a[0]) - R, 0), max(int(a[1]) - R, 0)
    ys, xs = np.nonzero(igl_c[y0:int(a[0]) + R + 1, x0:int(a[1]) + R + 1])
    ys, xs = ys + y0, xs + x0
    ch = np.stack([ys - a[0], xs - a[1]], 1)
    d = np.hypot(ch[:, 0], ch[:, 1])
    cosang = (ch @ ua) / np.maximum(d, 1e-9)
    rr = route[ys, xs]
    ok = (~np.isfinite(rr) | (rr >= prm["route_factor"] * np.maximum(d, 1))) & (d * G >= prm["min_gap_um"]) \
        & (d * G <= prm["max_gap_um"]) & (cosang >= np.cos(np.radians(prm["axis_cone"])))
    best, S.setdefault("_axis_cands", {})[e["id"]] = None, []
    for lab in np.unique(S["_cc"][ys[ok], xs[ok]]):
        m = ok & (S["_cc"][ys, xs] == lab)
        th = np.degrees(np.arccos(np.clip(cosang[m], -1, 1)))
        k = int(np.argmin(d[m] * (1 + th / 45)))
        q = np.array([ys[m][k], xs[m][k]], float)
        ub, _ = main_axis(S, q, 0.0, 150.0, both=True)
        if ub is None:
            continue
        c = q - a
        tha, thb = _angle(ua, c), _angle(ub, c, undirected=True)
        S["_axis_cands"][e["id"]].append(dict(q=[float(q[0]), float(q[1])], dist_um=round(float(np.hypot(*c) * G)), deg=[round(tha), round(thb)],
                                              axis_far=[round(float(ub[0]), 3), round(float(ub[1]), 3)]))
        if tha > prm["max_axis_turn"] or thb > prm["max_axis_turn"]:
            continue
        cost = tha + thb + float(np.hypot(*c) * G) / 10
        if best is None or cost < best[0]:
            best = (cost, q, ub if ub @ c > 0 else -ub, tha, thb)
    S["_axis_cands"][e["id"]].sort(key=lambda c: c["deg"][0] + c["deg"][1] + c["dist_um"] / 10)
    S.setdefault("_axis", {})[e["id"]] = [round(float(ua[0]), 3), round(float(ua[1]), 3)]
    if best is None:
        return None
    _, q, ub, tha, thb = best
    tg = target_at(S, q, a, e, ends, float(np.hypot(*(q - a)) * G), "main axis")
    if tg is None:
        # no centre line of the main description near q (a stretch pruned as small): take
        # width and depth from the IGL within 40 um of q
        yy, xx = np.ogrid[:igl_c.shape[0], :igl_c.shape[1]]
        near = igl_c & ((yy - q[0]) ** 2 + (xx - q[1]) ** 2 <= (40 / G) ** 2)
        if "_edt" not in S:
            S["_edt"] = ndi.distance_transform_edt(igl_c) * G
        tg = dict(pos=q, fill_pos=q, dir=-ub, width=max(float(2 * np.percentile(S["_edt"][near], 90)), prm["min_width_um"]),
                  depth=float(np.median(S["dE"][near])), id=None, extra=False, kind="ribbon",
                  corridor_um=float(np.hypot(*(q - a)) * G), found_by="main axis")
    tg["dir"] = -ub                         # evaluate_pair arrives along -dir
    tg["axis_deg"] = [round(tha), round(thb)]
    return tg, ua


def _continues(tg, v, depth, width, path_um, prm):
    """does the IGL the prediction ran into continue the ribbon? (user review 2026-11-19: the
    bridges to IGL that is close but "turns at right angles" / lies shallower were wrong.)
    End there: it lies at the predicted depth (within max(30 um, half the width)) or faces
    back along the course (within 50 degrees). A stretch of ribbon with no end there: the
    predicted depth is inside it (within max(30 um, half its width)) and it runs along the
    course (within 50 degrees); over a gap shorter than the ribbon is wide the far side's
    direction cannot be measured, so only the depth counts."""
    d = np.asarray(tg["dir"], float)
    d = d / max(np.hypot(*d), 1e-9)
    if tg["kind"] == "end":
        same_depth = abs(float(tg["depth"]) - depth) <= max(prm["contact_depth_um"], width / 2)
        return bool(same_depth or d @ v <= -prm["contact_cos"])
    same_depth = abs(float(tg["depth"]) - depth) <= max(prm["contact_depth_um"], float(tg["width"]) / 2)
    short = path_um <= max(width, float(tg["width"]))
    return bool(same_depth and (short or abs(d @ v) >= prm["contact_cos"]))


def evaluate_path(S, ea, eb, path, neun=None, t_neun=None):
    """fill and checks for a predicted bridge: the band of Hermite width along the predicted
    course; white-matter labels in it are kept (not filled) unless the outline test says they
    are wrong (white matter in the shallow half of the expected IGL band)."""
    prm, g, geo = S["prm"], S["g"], S["geo"]
    igl_c = S["igl_c"]
    a, b = np.asarray(ea["fill_pos"], float), np.asarray(eb["fill_pos"], float)
    wa, wb = ea["width"], eb["width"]
    da, db = float(ea["depth"]), float(eb["depth"])
    m = max(((da - wa / 2) + (db - wb / 2)) / 2, 8.0)
    w = (wa + wb) / 2
    C = np.vstack([a[None], path, b[None]])
    if len(C) >= 7:
        Cs = np.stack([ndi.gaussian_filter1d(C[:, k], 2.0, mode="nearest") for k in range(2)], 1)
        Cs[0], Cs[-1] = a, b
        C = Cs
    L = float(np.hypot(*np.diff(C, axis=0).T).sum() * G)
    chord_um = float(np.hypot(*(b - a)) * G)
    d = np.diff(C, axis=0)
    th = np.unwrap(np.arctan2(d[:, 0], d[:, 1])) if len(d) > 1 else np.zeros(1)
    out = dict(ends=[int(ea["id"]), None if eb["id"] is None else int(eb["id"])], a=[round(float(x), 1) for x in a], b=[round(float(x), 1) for x in b],
               extra=[bool(ea["extra"]), bool(eb["extra"])], chord_um=round(chord_um, 1), curve_um=round(L, 1),
               width_um=[round(wa), round(wb)], depth_um=[round(da), round(db)], m_local=round(m, 1),
               ml_um=[round(da - wa / 2), round(db - wb / 2)], turn_deg=round(float(np.degrees(np.abs(np.diff(th)).sum()))),
               energy=0.0, score=round(L / prm["len_scale_um"], 3))
    keep = S["ok"] & ~g["egl"] & ~igl_c
    band_all = ME.band(C, wa, wb, igl_c.shape, keep)
    shallow = band_all & (geo["dS"] <= m + w / 2)
    wm_shallow = float(np.count_nonzero(shallow & S["wm_eff"]) / max(np.count_nonzero(shallow), 1))
    out["wm_shallow"] = round(wm_shallow, 3)
    if wm_shallow > prm["wm_shallow"]:
        fill, method = band_all, "ribbon model, white-matter labels ignored"
    else:
        fill, method = band_all & ~S["wm_eff"], "predicted course"
    fill = fill & ~g["igl"]
    out.update(method=method, wm_ignored=method.startswith("ribbon model, white"), fill_um2=round(float(fill.sum() * G * G)),
               core_frac=round(float(S["core"][fill].mean()), 3) if fill.any() else 0.0, depth_gap_um=round(abs(da - db), 1))
    reasons = []
    pad = int(max(chord_um / G, L / G / 2, 60))
    if not fill.any() or not _joins(igl_c, fill, a, b, pad):
        reasons.append("does not join the two ends")
    if fill.sum() * G * G < prm["min_cover"] * chord_um * w:
        reasons.append("covers too little of the gap")
    if out["wm_ignored"] and out["core_frac"] > prm["max_core_frac"]:
        reasons.append("lies in the outline's white-matter core")
    in_bleach = float(S["bleach"][fill].mean()) if fill.any() else 0.0
    out["bleach_frac"] = round(in_bleach, 2)
    confirmed = in_bleach >= 0.5
    if not confirmed:
        if min(out["ml_um"]) < prm["min_ml_um"]:
            reasons.append("IGL touches the EGL on one side")
        if out["wm_ignored"] and out["depth_gap_um"] > prm["max_depth_gap_um"]:
            reasons.append("ends at different depths")
    if fill.any():
        ml_frac = float(g["ml"][fill].mean())
        out["ml_frac"] = round(ml_frac, 2)
        out["wm_frac"] = round(float(g["wm"][fill].mean()), 2)
        if neun is not None:
            nm = float(np.median(neun[fill]))
            out["neun_ratio"] = round(nm / t_neun, 2)
            if not confirmed and ml_frac > prm["max_ml_frac_without_neun"] and nm < t_neun:
                reasons.append("through NeuN-dim ML (real ML)")
            elif not confirmed and nm < prm["min_neun_ratio"] * t_neun:
                reasons.append("NeuN too dim for IGL")
        elif ml_frac > prm["max_ml_frac_without_neun"] and not confirmed:
            reasons.append("through ML, no NeuN image to check")
    out["status"] = "rejected: " + ", ".join(reasons) if reasons else "ok"
    c = _crop(fill)
    if c is not None:
        out["_box"], out["_mask"] = c
    out["_curve"] = C
    return out, fill


def _path_curve(a, path, b):
    C = np.vstack([np.asarray(a, float)[None], np.asarray(path, float).reshape(-1, 2), np.asarray(b, float)[None]])
    if len(C) >= 7:
        Cs = np.stack([ndi.gaussian_filter1d(C[:, k], 2.0, mode="nearest") for k in range(2)], 1)
        Cs[0], Cs[-1] = C[0], C[-1]
        C = Cs
    return C


def evaluate_pair(S, ea, eb, neun=None, t_neun=None, course=None):
    """course: the predicted course (predict); None = the ribbon model's curve between the ends."""
    prm, g, geo = S["prm"], S["g"], S["geo"]
    # fill / curve positions: where the IGL stops on each side (not the skeleton ends)
    a, b = np.asarray(ea.get("fill_pos", ea["pos"]), float), np.asarray(eb.get("fill_pos", eb["pos"]), float)
    ta, tb_in = np.asarray(ea["dir"], float), -np.asarray(eb["dir"], float)
    wa, wb = ea["width"], eb["width"]
    da, db = float(ea["depth"]), float(eb["depth"])
    m = max(((da - wa / 2) + (db - wb / 2)) / 2, 8.0)
    w = (wa + wb) / 2
    chord_um = float(np.hypot(*(b - a)) * G)
    igl_c = S["igl_c"]
    yy, xx = np.ogrid[:igl_c.shape[0], :igl_c.shape[1]]
    avoid = igl_c.copy()
    for p in (a, b):
        avoid &= (yy - p[0]) ** 2 + (xx - p[1]) ** 2 > 16
    if course is None:
        C, info = M.complete(a, ta, b, tb_in, S["field"], da, db, MODEL, avoid)
    else:
        C, info = np.asarray(course, float), dict(energy=0.0)
    L = float(np.hypot(*np.diff(C, axis=0).T).sum() * G)
    out = dict(ends=[int(ea["id"]), None if eb["id"] is None else int(eb["id"])], a=[round(float(v), 1) for v in a], b=[round(float(v), 1) for v in b],
               extra=[bool(ea["extra"]), bool(eb["extra"])], chord_um=round(chord_um, 1), curve_um=round(L, 1),
               width_um=[round(wa), round(wb)], depth_um=[round(da), round(db)], m_local=round(m, 1),
               energy=round(float(info["energy"]), 3), score=round(float(info["energy"]) + L / prm["len_scale_um"], 3))
    if L > prm["max_curve_um"]:
        out.update(status="rejected: longer than %d um" % prm["max_curve_um"])
        return out, None
    # outline test: white matter in the shallow half of the expected IGL band
    Pd = OF.predict_depth(geo, m, w) & ~S["core"]
    Fd, st_d = OF.fill_between(Pd, igl_c, S["allowed"], a, b, wa, wb)
    if st_d != "ok":
        Fd = OF.band(Pd, igl_c, S["allowed"], C, wa, wb)
    shallow = Fd & (geo["dS"] <= m + w / 2)
    wm_shallow = float(np.count_nonzero(shallow & S["wm_eff"]) / max(np.count_nonzero(shallow), 1))
    out["wm_shallow"] = round(wm_shallow, 3)
    curve = C
    if wm_shallow > prm["wm_shallow"]:
        C2 = C if course is not None else M.complete(a, ta, b, tb_in, S["field_nowm"], da, db, MODEL, avoid)[0]
        curve = C2
        fill = ME.band(C2, wa, wb, igl_c.shape, S["ok"] & ~g["egl"] & ~igl_c)
        method = "ribbon model, white-matter labels ignored"
    else:
        P = OF.predict(geo, S["fit_label"], S["wm_eff"], "label", m=m)
        fill, st = OF.fill_between(P, igl_c, S["allowed"], a, b, wa, wb)
        method = "outline corridor"
        if st != "ok":
            fill = OF.band(P, igl_c, S["allowed"], C, wa, wb)
            method = "outline rule along the ribbon model"
    fill = fill & ~g["igl"]
    out["method"] = method
    out["wm_ignored"] = method.startswith("ribbon model, white")
    out["ml_um"] = [round(da - wa / 2), round(db - wb / 2)]
    out["fill_um2"] = round(float(fill.sum() * G * G))
    out["core_frac"] = round(float(S["core"][fill].mean()), 3) if fill.any() else 0.0
    out["depth_gap_um"] = round(abs(da - db), 1)
    reasons = []
    pad = int(max(chord_um / G, 60))
    if not fill.any() or not _joins(igl_c, fill, a, b, pad):
        reasons.append("does not join the two ends")
    if fill.sum() * G * G < prm["min_cover"] * chord_um * w:
        reasons.append("covers too little of the gap")
    if method.startswith("ribbon model, white") and out["core_frac"] > prm["max_core_frac"]:
        reasons.append("lies in the outline's white-matter core")
    in_bleach = float(S["bleach"][fill].mean()) if fill.any() else 0.0
    out["bleach_frac"] = round(in_bleach, 2)
    confirmed = in_bleach >= 0.5
    if method != "outline corridor" and not confirmed:
        Cr = ME.resample(curve, 1.0)
        iy = np.clip(np.round(Cr[:, 0]).astype(int), 0, igl_c.shape[0] - 1)
        ix = np.clip(np.round(Cr[:, 1]).astype(int), 0, igl_c.shape[1] - 1)
        f = np.linspace(0, 1, len(Cr))
        dev = float(np.max(np.abs(S["dE"][iy, ix] - ((1 - f) * da + f * db))))
        out["depth_dev_um"] = round(dev, 1)
        if dev > max(prm["max_depth_dev"] * (da + db) / 2, prm["min_depth_dev_um"]):
            reasons.append("left the ribbon depth")
        if method.startswith("ribbon model, white") and out["depth_gap_um"] > prm["max_depth_gap_um"]:
            reasons.append("ends at different depths")
    if not confirmed and min(out["ml_um"]) < prm["min_ml_um"]:
        reasons.append("IGL touches the EGL on one side")         # cKO: fuzzy inner EGL labelled IGL
    if fill.any():
        ml_frac = float(g["ml"][fill].mean())
        out["ml_frac"] = round(ml_frac, 2)
        out["wm_frac"] = round(float(g["wm"][fill].mean()), 2)
        if neun is not None:
            nm = float(np.median(neun[fill]))
            out["neun_ratio"] = round(nm / t_neun, 2)
        if ml_frac > prm["max_ml_frac_without_neun"] and not confirmed:
            if neun is None:
                reasons.append("through ML, no NeuN image to check")
            elif nm < t_neun:
                reasons.append("through NeuN-dim ML (real ML)")
        elif neun is not None and not confirmed and nm < prm["min_neun_ratio"] * t_neun:
            reasons.append("NeuN too dim for IGL")                  # cKO: diffuse IGL
    out["status"] = "rejected: " + ", ".join(reasons) if reasons else "ok"
    c = _crop(fill)
    if c is not None:
        out["_box"], out["_mask"] = c
    out["_curve"] = curve
    return out, fill


def fill_gaps(g, prm=None, rib=None, neun=None, tear=None, bleach=None, verbose=False, trace=False):
    """g: 4 um label grid (igl_ribbon.grid_labels); neun: NeuN on the same grid; tear: tear
    mask (never crossed); bleach: user-confirmed bleach / breach regions (white-matter labels
    ignored there, checks relaxed). Returns (new IGL mask, bridges, ribbon description)."""
    S = section_setup(g, tear, bleach, rib, prm)
    if trace:
        S["_trace"] = []                          # every prediction: contacts, stops, courses
    t_neun = None
    if neun is not None:
        t_neun = float(np.sqrt(np.median(neun[g["igl"]]) * np.median(neun[g["ml"]])))
    ends = ends_of(S)
    S["_route"], S["_sk"] = {}, skeleton_points(S)
    rows, fills, seen = [], {}, set()
    for e in ends:
        e["fill_pos"] = tip(S, e)
    for e in ends:
        if S["prm"]["search"] in ("predict", "hybrid"):
            res = predict(S, e, ends)
            if res is not None:
                tg, path = res
                if S["prm"]["fill"] == "corridor":
                    r, f = evaluate_pair(S, e, tg, neun, t_neun, course=_path_curve(e["fill_pos"], path, tg["fill_pos"]))
                else:
                    r, f = evaluate_path(S, e, tg, path, neun, t_neun)
                r.update(target=tg["kind"], found_by=tg["found_by"], corridor_um=round(tg["corridor_um"], 1))
                if f is not None and r["status"] == "ok":
                    fills[len(rows)] = f
                rows.append(r)
            if S["prm"]["search"] == "predict":
                continue
            ax = axis_target(S, e, ends) if S["prm"]["axis_search"] else None
            if ax is not None:
                tg, ua = ax
                ea = dict(e, dir=ua)
                r, f = evaluate_pair(S, ea, tg, neun, t_neun, course=_hermite_course(e["fill_pos"], ua, tg["fill_pos"], -np.asarray(tg["dir"])))
                r.update(target=tg["kind"], found_by="main axis", corridor_um=round(tg["corridor_um"], 1), axis_deg=tg["axis_deg"])
                if f is not None and r["status"] == "ok":
                    fills[len(rows)] = f
                rows.append(r)
        tgs = [t for t in (discover(S, e, ends), walk_target(S, e, ends)) if t is not None]
        if len(tgs) == 2 and np.hypot(*(tgs[0]["fill_pos"] - tgs[1]["fill_pos"])) * G <= 40:
            tgs = tgs[:1]                 # both searches reached the same place
        for tg in tgs:
            key = (e["id"], tg["id"], tg["found_by"]) if tg["id"] is not None else None
            if key is not None and (tg["id"], e["id"], tg["found_by"]) in seen:
                continue                  # the same pair found from the other side
            if key is not None:
                seen.add(key)
            r, f = evaluate_pair(S, e, tg, neun, t_neun)
            r.update(target=tg["kind"], found_by=tg["found_by"], corridor_um=round(tg["corridor_um"], 1))
            if S["prm"]["search"] == "hybrid" and S["prm"]["axis_check"] and r["chord_um"] > max(r["width_um"]) and r.get("bleach_frac", 0) < 0.5:
                tha, thb, _, _ = axis_angles(S, e["fill_pos"], tg["fill_pos"])
                r["axis_deg"] = [None if tha is None else round(tha), None if thb is None else round(thb)]
                if (tha is not None and tha > S["prm"]["max_axis_turn"]) or (thb is not None and thb > S["prm"]["max_axis_turn"]):
                    r["status"] = (r["status"] + ", " if r["status"].startswith("rejected") else "rejected: ") + "turns away from the ribbon's main axis"
            if f is not None and r["status"] == "ok":
                fills[len(rows)] = f
            rows.append(r)
    if verbose:
        print(f"    {len(ends)} ends ({sum(e['extra'] for e in ends)} from light pruning), {len(rows)} with a far side in reach", flush=True)
    used, new = set(), g["igl"].copy()
    taken = np.zeros_like(new)
    rank = {"prediction": 0, "main axis": 1}
    first = (lambda r: rank.get(r["found_by"], 2)) if S["prm"]["search"] == "hybrid" else (lambda r: 0)
    order = sorted([k for k, r in enumerate(rows) if r["status"] == "ok"], key=lambda k: (first(rows[k]), rows[k]["score"]))
    for k in order:
        r, f = rows[k], fills[k]
        ids = {i for i in r["ends"] if i is not None}
        if used & ids or (f & ~taken).sum() < 0.3 * f.sum():
            r["status"] = "superseded: the same gap was filled from a better pair"
            continue
        deep = None
        if r.get("found_by") == "prediction":
            # white matter cut off on the shallow side of a predicted course is mislabelled IGL
            # (user: s2_2 bridge "should be joined"); only white matter deeper than the
            # ribbon band (m + 1.5 w below the EGL) counts
            deep = S["dE"] > r["m_local"] + 1.5 * float(np.mean(r["width_um"]))
        t = topology(new, g["wm"] & ~new, f, deep)
        r.update(t)
        bad = layout_reasons(t, S["prm"], confirmed=r.get("bleach_frac", 0.0) >= 0.5)
        if bad:
            r["status"] = "rejected: " + ", ".join(bad)
            continue
        used |= ids
        r["status"] = "joined"
        new |= f
        taken |= ndi.binary_dilation(f, iterations=2)
    for r in rows:
        r["length_um"] = r["chord_um"]
    S["ends"] = ends
    S["rib"]["fill_setup"] = dict(fit_label=S["fit_label"], fit_outline=S["fit_outline"], n_ends=len(ends),
                                  ends=[dict(id=e["id"], tip=[float(v) for v in e["fill_pos"]], depth=float(e["depth"]), width=float(e["width"])) for e in ends],
                                  trace=S.get("_trace"), axis=S.get("_axis", {}), axis_candidates=S.get("_axis_cands", {}))
    return new, rows, S["rib"]
