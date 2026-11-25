"""
Final IGL-repaired label maps for s4_A and s4_C from the user's review (2026-11-18).

    python scripts/igl_finalize.py

Second review: s4_A marker 1 still disconnected, marker 2 still absorbed the EGL;
s4_C "failed to bridge". On top of the ribbon repair:
  - the accepted box fixes B (s4_A B5, s4_C B2: IGL shell around the white matter);
  - s4_A folium 6 smear from the image: DAPI-dense tissue within 120 um of the surface
    = EGL (the image shows the EGL ~50 um deeper than in neighbouring folia, under a
    dim outer layer), DAPI-poor tissue under it = ML, the rest stays IGL.
Gap filling: igl_gapfill.fill_gaps with the user's bleach box (+ 100 um) and markers (150 um)
as confirmed breaches (white-matter labels ignored there, depth / NeuN checks relaxed).
Third review (2026-11-19): the shortest-path joins forced at the markers and the E4 band fill
were wrong ("connecting to the wrong part of the IGL, it should connect to the deeper part,
it shouldn't turn at right angles"), so they are gone; the filler's predicted ribbon
continuation now decides where each break joins. Per marker, the bridges that reach it and
the prediction from every ribbon end within 150 um are logged (and their courses saved).
Fourth review (2026-11-19): at markers D4 and E4 the user traced the IGL by clicking on the
NeuN image (outputs/igl_fixed/user_traces.json, grid coordinates); trace_band paints it.
Fifth review: E4's end carried on to the IGL; then "you are missing this whole segment below
the bleached area": further strokes (round5 in user_traces.json) painted in click order, over
EGL labels too, ends joined to IGL / earlier traces within 60 um.
Earlier outputs: outputs/igl_fixed_walk_backup (walk), outputs/igl_fixed_prelayout_backup,
outputs/igl_fixed_layout_backup, outputs/igl_fixed_round3_backup.
Grid (4 um) decisions are stamped onto native pixels with an allowed source label.
Output: outputs/igl_fixed/<slide>_igl_fixed_labels.tif, igl_fixed.json,
        <slide>_final_display.npz (before / after in the review-page frame),
        <slide>_marker_courses.npz (predicted courses near the markers, 4 um grid),
        <slide>_bridges.npz (fill of every bridge, in igl_fixed.json order, 4 um grid)
"""
import json
import sys
from pathlib import Path

import numpy as np
import tifffile
from scipy import ndimage as ndi

sys.path.insert(0, str(Path(__file__).parent))
import folium_repair as FR  # noqa: E402
import igl_check_page as P  # noqa: E402
import igl_ribbon as RB  # noqa: E402
import igl_gapfill as GF  # noqa: E402
import igl_repair_eval as EV  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/igl_fixed"
RIB = ROOT / "outputs/igl_ribbon"
PX = 0.65
MARKERS_ALL = {"s4_A_p27": {"m1": (0.6168, 0.4863)},
               "s4_C_p27": {"C3": (0.3321, 0.2319), "D4": (0.4805, 0.3616), "E4": (0.6192, 0.3807)}}
TRACES = json.loads((OUT / "user_traces.json").read_text()) if (OUT / "user_traces.json").exists() else {}
ACCEPTED_BOX_FIX = {"s4_A_p27": "B5", "s4_C_p27": "B2"}       # fix B (shape) accepted by the user


def smear_fix_dapi(labr, terr, folium, dapi):
    """folium 6 smear, from the image: in the folium, DAPI-dense tissue within 120 um of
    the surface is EGL (oEGL in its outer third), DAPI-poor tissue under it is ML, the
    rest stays IGL. Thresholds from the section's own layers outside the folium."""
    out = labr.copy()
    R = (terr == folium) & (labr > 0)
    outside = (terr != folium) & (labr > 0)
    d = ndi.gaussian_filter(dapi, 4 / (0.65 * P.DS))
    med = {k: float(np.median(d[outside & np.isin(labr, c)])) for k, c in (("egl", (2, 3)), ("igl", (4,)), ("ml", (5,)))}
    t_egl = (med["egl"] + med["igl"]) / 2
    t_ml = (med["ml"] + med["igl"]) / 2
    um = 0.65 * P.DS
    filled = ndi.binary_fill_holes(labr > 0)
    dO = ndi.distance_transform_edt(filled) * um
    src = R & np.isin(labr, (1, 4, 7, 9))
    egl = src & (d >= t_egl) & (dO <= 120)
    egl = ndi.binary_opening(egl, structure=np.ones((3, 3), bool))
    # EGL band: from the surface down to the deepest DAPI-dense pixel along each folium column
    deep = float(np.percentile(dO[egl], 95)) if egl.any() else 0.0
    egl_band = src & (dO <= deep) & ((d >= t_ml) | egl)
    out[egl_band] = np.where(dO[egl_band] <= deep / 3, 3, 2)
    ml_zone = src & ~egl_band & (dO <= deep + 50) & (d < t_ml)
    out[ml_zone] = 5
    return out, dict(egl_depth_um=round(deep, 1), dapi_medians={k: round(v) for k, v in med.items()},
                     egl_um2=round(float(egl_band.sum() * um * um)), ml_um2=round(float(ml_zone.sum() * um * um)))


def disp_to_grid(tr, x, y):
    h, w = tr["shape"]
    q = np.array(tr["A"]) @ (np.array([y * h, x * w]) + np.array(tr["crop"])) + np.array(tr["offset"])
    return q * P.DS * PX / RB.GRID


def confirmed_regions(sid, g, tr=None):
    """the user's bleach box (+ 100 um) and break markers (150 um) on the 4 um grid."""
    if tr is None:
        tr = {t["sid"]: t for t in json.loads((P.OUT / "transforms.json").read_text())}[sid]
    box = np.load(RIB / f"{sid}_ribbon_repair_grid.npz")["bleach"]
    conf = ndi.distance_transform_edt(~box) * RB.GRID <= 100
    yy, xx = np.indices(g["igl"].shape)
    for name, (x, y) in MARKERS_ALL.get(sid, {}).items():
        c = disp_to_grid(tr, x, y)
        conf |= ((yy - c[0]) ** 2 + (xx - c[1]) ** 2) * RB.GRID ** 2 <= 150 ** 2
    return conf


def order_clicks(pts, igl):
    """the user's clicks as one path: start at the click nearest the IGL, then always the
    nearest remaining click (the page stores clicks in the order made, and the user adds
    fill-in clicks out of order)."""
    P = [np.asarray(q, float) for q in pts]
    dI = ndi.distance_transform_edt(~igl)
    k = int(np.argmin([dI[int(round(q[0])), int(round(q[1]))] for q in P]))
    path = [P.pop(k)]
    while P:
        j = int(np.argmin([np.hypot(*(q - path[-1])) for q in P]))
        path.append(P.pop(j))
    return np.array(path)


def trace_band(g, clicks, tear, max_width_um=70.0, connect_um=0.0, ordered=False, egl_full=False, connect_start=False, also=None, width_um=None):
    """IGL along the user's traced course (third review, markers D4 / E4: "follow the clicks"):
    a smooth centre line through the ordered clicks (Catmull-Rom), painted at the width of the
    IGL where the trace starts (2 x 90th percentile of its distance transform, 30-70 um) over
    every tissue label but EGL; EGL labels
    only within 8 um of the centre line (the user clicked on them)."""
    igl = g["igl"]
    P = np.asarray(clicks, float) if ordered else order_clicks(clicks, igl)
    if len(P) == 1:
        P = np.vstack([P, P + 0.01])
    Q = np.vstack([P[:1], P, P[-1:]])
    C = []
    for i in range(1, len(Q) - 2):
        p0, p1, p2, p3 = Q[i - 1], Q[i], Q[i + 1], Q[i + 2]
        n = max(int(np.hypot(*(p2 - p1)) * 2), 2)
        for t in np.linspace(0, 1, n, endpoint=False):
            C.append(0.5 * ((2 * p1) + (-p0 + p2) * t + (2 * p0 - 5 * p1 + 4 * p2 - p3) * t ** 2 + (-p0 + 3 * p1 - 3 * p2 + p3) * t ** 3))
    C.append(P[-1])
    C = np.array(C)
    if connect_um:
        # the user (E4): "extend it more to the rest of the igl region": carry the open end
        # straight on to the nearest IGL (or, with also, an earlier trace) within connect_um
        target = igl | (also if also is not None else False)
        dI, (iy_, ix_) = ndi.distance_transform_edt(~target, return_indices=True)
        for at_start in ((False, True) if connect_start else (False,)):
            e = np.round(C[0] if at_start else C[-1]).astype(int)
            e = np.clip(e, 0, np.array(igl.shape) - 1)
            if 0 < dI[e[0], e[1]] * RB.GRID <= connect_um:
                q = np.array([iy_[e[0], e[1]], ix_[e[0], e[1]]], float)
                n = max(int(np.hypot(*(q - e)) * 2), 2)
                seg = e + (q - e) * np.linspace(0, 1, n)[1:, None]
                C = np.vstack([seg[::-1], C]) if at_start else np.vstack([C, seg])
    line = np.zeros(igl.shape, bool)
    iy = np.clip(np.round(C[:, 0]).astype(int), 0, igl.shape[0] - 1)
    ix = np.clip(np.round(C[:, 1]).astype(int), 0, igl.shape[1] - 1)
    line[iy, ix] = True
    d = ndi.distance_transform_edt(~line) * RB.GRID
    yy, xx = np.indices(igl.shape)
    near0 = igl & ((yy - P[0][0]) ** 2 + (xx - P[0][1]) ** 2 <= (60 / RB.GRID) ** 2)
    width = width_um or float(np.clip(np.percentile(2 * ndi.distance_transform_edt(igl)[near0] * RB.GRID, 90) if near0.any() else 50.0, 30.0, max_width_um))
    tissue = (g["tissue"] | g["anom"]) & ~tear
    band = (d <= width / 2) & tissue & ~g["egl"] & ~igl
    egl = ((d <= width / 2) if egl_full else (d <= 8)) & g["egl"] & ~tear
    L = float(np.hypot(*np.diff(C, axis=0).T).sum() * RB.GRID)
    info = dict(length_um=round(L), width_um=round(width), clicks=len(P), joined_to_igl=bool(connect_um) and bool(ndi.distance_transform_edt(~igl)[int(round(C[-1][0])), int(round(C[-1][1]))] * RB.GRID <= 4),
                painted_um2={k: round(float(((band | egl) & g[k]).sum() * RB.GRID ** 2)) for k in ("ml", "wm", "egl", "unass", "anom") if k in g},
                end_to_igl_um=round(float(ndi.distance_transform_edt(~igl)[int(round(C[-1][0])), int(round(C[-1][1]))] * RB.GRID)))
    return band, egl, C, info


def bright_next_to(g, band, sid, tear, reach_um=120.0):
    """tissue within reach_um of a traced band that is as bright as the band in channels 0 and 1
    of the native image (there the bleached IGL still stands out: band ~ IGL level, the tissue
    around ~ ML level), connected to it. Shown to the user as an option, not painted: about
    half of it carries EGL labels in s4_C."""
    raw = tifffile.imread(EV.SB10 / f"{sid.replace('_p27', '')}_p27_native.tif") if not (EV.SB10 / f"{sid}_native.tif").exists() else tifffile.imread(EV.SB10 / f"{sid}_native.tif")
    f = RB.GRID / PX
    H, W = g["igl"].shape
    ys = np.minimum((np.arange(H) * f).astype(int), raw.shape[1] - 1)
    xs = np.minimum((np.arange(W) * f).astype(int), raw.shape[2] - 1)
    img = [ndi.gaussian_filter(ndi.uniform_filter(raw[c].astype(np.float32), size=int(f))[np.ix_(ys, xs)], 1.0) for c in (0, 1)]
    tissue = (g["tissue"] | g["anom"]) & ~tear
    d = ndi.distance_transform_edt(~band) * RB.GRID
    ring = (d <= 80) & ~band & tissue & ~g["egl"] & ~g["igl"]
    t = [(np.median(img[c][band]) + np.median(img[c][ring])) / 2 for c in (0, 1)]
    egl_hi = (np.median(img[1][g["igl"]]) + np.median(img[1][g["egl"]])) / 2
    bright = tissue & (img[0] >= t[0]) & (img[1] >= t[1]) & ~(g["egl"] & (img[1] >= egl_hi)) & (d <= reach_um)
    bright = ndi.binary_opening(bright, iterations=1)
    lab, _ = ndi.label(bright | band, structure=np.ones((3, 3), bool))
    keep = np.isin(lab, np.unique(lab[band]))
    return ndi.binary_closing(keep | band, iterations=2) & tissue & ~band & ~g["igl"]


def stamp(native, grid_paint, grid_allow_wm, grid_allow_egl=None, px=PX):
    f = RB.GRID / px
    H, W = native.shape
    ys = np.minimum((np.arange(H) / f).astype(int), grid_paint.shape[0] - 1)
    xs = np.minimum((np.arange(W) / f).astype(int), grid_paint.shape[1] - 1)
    p_up = grid_paint[np.ix_(ys, xs)]
    w_up = grid_allow_wm[np.ix_(ys, xs)]
    allowed = np.isin(native, (1, 5, 7, 9)) | (w_up & np.isin(native, (6, 8)))
    if grid_allow_egl is not None:
        allowed |= grid_allow_egl[np.ix_(ys, xs)] & np.isin(native, (2, 3))
    out = native.copy()
    out[p_up & allowed] = 4
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    trs = {t["sid"]: t for t in json.loads((P.OUT / "transforms.json").read_text())}
    log = {}
    for sid in ("s4_A_p27", "s4_C_p27"):
        tr = trs[sid]
        native = tifffile.imread(EV.RE.FX / f"{sid}_fixed_labels.tif")
        g = RB.grid_labels(native, PX)
        neun, tear = EV.images(sid.replace("_p27", ""), g["igl"].shape)
        conf = confirmed_regions(sid, g, tr)
        new, br, rib = GF.fill_gaps(g, neun=neun, tear=tear, bleach=conf, trace=True)
        paint = new & ~g["igl"]
        # white matter may be painted in the confirmed regions and where a bridge found the
        # white-matter labels wrong (outline test)
        allow_wm = conf.copy()
        for b in br:
            if b["status"] == "joined" and b.get("wm_ignored"):
                y0, y1, x0, x1 = b["_box"]
                allow_wm[y0:y1, x0:x1] |= b["_mask"]
        entry = dict(bridges=[{k: v for k, v in b.items() if not k.startswith("_")} for b in br if b["status"] == "joined"])
        # the user's traced courses (third review): painted as they are, labels overridden
        allow_egl = np.zeros_like(paint)
        traces, tmasks = {}, {}
        for name, tr_ in TRACES.get(sid, {}).items():
            band, egl, C, info = trace_band(g, tr_["clicks_grid"], tear, connect_um=60.0 if tr_.get("connect_end") else 0.0)
            tmasks[f"{name}_bright"] = bright_next_to(g, band | egl, sid, tear)
            # fifth review: "you are missing this whole segment below the bleached area": the
            # user's further strokes, painted in click order over every label (EGL included)
            done = band | egl
            w_sec = float(np.median([np.median(b["width"]) for b in RB.describe(g)["branches"]]))   # the section's ribbon width
            for k, stroke in enumerate(tr_.get("round5", {}).get("strokes", [])):
                b2, e2, C2, i2 = trace_band(g, stroke, tear, connect_um=60.0, ordered=True, egl_full=True, connect_start=True, also=done, width_um=w_sec)
                tmasks[f"{name}_stroke{k}_line"] = C2
                band |= b2
                egl |= e2
                done |= b2 | e2
                info.setdefault("strokes", []).append(i2)
            paint |= band | egl
            allow_wm |= band
            allow_egl |= egl
            traces[name] = info
            tmasks[f"{name}_band"], tmasks[f"{name}_egl"], tmasks[f"{name}_line"] = band, egl, C
        entry["user_traces"] = traces
        if tmasks:
            np.savez_compressed(OUT / f"{sid}_traces.npz", **tmasks)
        # per break marker: the bridges that reach it and what the prediction did from every
        # ribbon end within 150 um (courses kept for the review page)
        fs = rib["fill_setup"]
        tips = {e["id"]: np.array(e["tip"]) for e in fs["ends"]}
        # main axes and the best far sides along them (axis search on for this report only)
        fs_ax = GF.fill_gaps(g, prm={"axis_search": True}, neun=neun, tear=tear, bleach=conf)[2]["fill_setup"]
        courses, marks = {}, {}
        for name, (x, y) in MARKERS_ALL.get(sid, {}).items():
            c = disp_to_grid(tr, x, y)
            near_b = []
            for b in br:
                if b["status"] == "joined":
                    by0, by1, bx0, bx1 = b["_box"]
                    ys, xs = np.nonzero(b["_mask"])
                    if (np.hypot(ys + by0 - c[0], xs + bx0 - c[1]) * RB.GRID <= 150).any():
                        near_b.append(dict(a=b["a"], b=b["b"], chord_um=b["chord_um"], method=b["method"]))
            preds = []
            for t in fs["trace"]:
                if np.hypot(*(tips[t["end"]] - c)) * RB.GRID > 150:
                    continue
                k = f"{name}_{len(preds)}"
                courses[k] = np.vstack([tips[t["end"]][None], t["_path"]]) if len(t["_path"]) else tips[t["end"]][None]
                preds.append(dict({kk: v for kk, v in t.items() if not kk.startswith("_")}, key=k))
            axes = []
            for eid, tp in tips.items():
                if np.hypot(*(tp - c)) * RB.GRID <= 150 and eid in fs_ax["axis"]:
                    axes.append(dict(end=eid, tip=[float(v) for v in tp], axis=fs_ax["axis"][eid], candidates=fs_ax["axis_candidates"].get(eid, [])[:3]))
            marks[name] = dict(centre=[round(float(v), 1) for v in c], bridges=near_b, predictions=preds, axes=axes)
        entry["markers"] = marks
        np.savez_compressed(OUT / f"{sid}_marker_courses.npz", **courses)
        fills = {}
        for i, b in enumerate([b for b in br if b["status"] == "joined"]):
            f = np.zeros(g["igl"].shape, bool)
            y0, y1, x0, x1 = b["_box"]
            f[y0:y1, x0:x1] = b["_mask"]
            fills[f"fill_{i}"] = f
        np.savez_compressed(OUT / f"{sid}_bridges.npz", **fills)
        out = stamp(native, paint, allow_wm, allow_egl)
        # display-frame fixes: accepted box fix B (both), image-guided smear fix (s4_A folium 6)
        y0, x0 = tr["crop"]
        h, w = tr["shape"]
        small_now = out[::P.DS, ::P.DS]
        labr_now, _ = P.rotate(small_now, np.array(tr["A"]), tr["canvas"], 0)
        labr_now = labr_now[y0:y0 + h, x0:x0 + w]
        terr = FR.territories(labr_now, sid)
        img3 = np.stack([FR.display_image(sid, tr, c) for c in range(3)])
        fixed = labr_now.copy()
        filled = ndi.binary_fill_holes(labr_now > 0)
        dO = ndi.distance_transform_edt(filled) * FR.PX
        dW = ndi.distance_transform_edt(~np.isin(labr_now, (6, 8))) * FR.PX
        t = FR.thickness(labr_now, terr, [k for k in range(1, 40) if k not in (6,)], dO, dW)
        bid = ACCEPTED_BOX_FIX[sid]
        box = [b for b in tr["bleach_candidates"] if b["id"] == bid][0]
        Rb = FR.box_mask_display(box, tr, labr_now.shape)
        fixed, _ = FR.bleach_fix(fixed, Rb, "B", img3, t)
        entry["box_fix_B"] = dict(box=bid, changed_um2=round(float((fixed != labr_now).sum() * FR.PX ** 2)))
        if sid == "s4_A_p27":
            before_smear = fixed.copy()
            fixed, info = smear_fix_dapi(fixed, terr, 6, img3[1])
            info["changed_um2"] = round(float((fixed != before_smear).sum() * FR.PX ** 2))
            entry["smear_fix_folium6_dapi"] = info
        out = FR.apply_native(out, P.DS, FR.back_to_section(fixed, labr_now, small_now, tr))
        tifffile.imwrite(OUT / f"{sid}_igl_fixed_labels.tif", out, compression="zlib")
        before_d, _ = P.rotate(native[::P.DS, ::P.DS], np.array(tr["A"]), tr["canvas"], 0)
        after_d, _ = P.rotate(out[::P.DS, ::P.DS], np.array(tr["A"]), tr["canvas"], 0)
        y0, x0 = tr["crop"]
        h, w = tr["shape"]
        np.savez_compressed(OUT / f"{sid}_final_display.npz", before=before_d[y0:y0 + h, x0:x0 + w],
                            after=after_d[y0:y0 + h, x0:x0 + w])
        g_out = RB.grid_labels(out, PX)
        entry["layout_gap_paint"] = GF.topology(g["igl"], g["wm"], paint)
        entry["layout_written"] = GF.topology(g["igl"], g["wm"], g_out["igl"] & ~g["igl"])
        ch = out != native
        entry.update(changed_mm2=round(float(ch.sum() * PX * PX / 1e6), 4),
                     igl_added_mm2=round(float(((out == 4) & (native != 4)).sum() * PX * PX / 1e6), 4),
                     igl_removed_mm2=round(float(((native == 4) & (out != 4)).sum() * PX * PX / 1e6), 4))
        log[sid] = entry
        print(sid, {k: v for k, v in entry.items() if k not in ("bridges", "markers")},
              "| bridges:", [(round(b["length_um"]), round(b.get("fill_um2", 0)), b["method"]) for b in entry["bridges"]], flush=True)
    (OUT / "igl_fixed.json").write_text(json.dumps(log, indent=1, default=float))


if __name__ == "__main__":
    main()
