"""
IGL Gap Fill review page: every bridge the pipeline gap filler (igl_gapfill_apply.py) made on
the real sections, with the NeuN image, the labels before and after, and a verdict per bridge;
plus the end-to-end synthetic test (igl_gapfill_eval.py), the user's earlier verdicts against
each filler version (igl_gapfill_score.py) and, per break marker, the predicted ribbon courses.

    python scripts/igl_gapfill_page.py      # -> outputs/igl_gapfill/igl_gap_fill.html
"""
import base64
import io
import json
import sys
from pathlib import Path

import numpy as np
import tifffile
from scipy import ndimage as ndi
from PIL import Image, ImageDraw, ImageFont

sys.path.insert(0, str(Path(__file__).parent))
import igl_gapfill_apply as GA  # noqa: E402
import igl_gapfill as GF  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/igl_gapfill"
TEMPLATE = Path(__file__).parent / "igl_gapfill_template.html"
G = 4.0
PAD = 30            # grid px (120 um) around each bridge
DISP_UM = 2.0       # crop panels: um per display pixel
# display codes from igl_gapfill_apply.grid_code: 0 bg, 1 tissue, 2 EGL, 3 ML/PCL, 4 IGL, 5 WM, 6 anomaly
PAL = np.array([(11, 14, 18), (58, 66, 77), (122, 90, 166), (78, 90, 104), (170, 180, 192), (140, 115, 72), (96, 70, 84)], np.uint8)
FILL = (59, 224, 255)
WALK = (255, 176, 0)
CURVE = (255, 255, 255)
CUT = (255, 84, 84)
COURSE = (255, 214, 10)
COURSE_OK = (80, 220, 120)
AXIS = (255, 80, 220)
FIXED = ROOT / "outputs/igl_fixed"


def png(rgb):
    b = io.BytesIO()
    Image.fromarray(rgb.astype(np.uint8)).save(b, "PNG", optimize=True)
    return "data:image/png;base64," + base64.b64encode(b.getvalue()).decode()


def jpg(gray):
    b = io.BytesIO()
    Image.fromarray(gray.astype(np.uint8)).save(b, "JPEG", quality=82)
    return "data:image/jpeg;base64," + base64.b64encode(b.getvalue()).decode()


def up(rgb, k):
    return np.repeat(np.repeat(rgb, k, 0), k, 1)


def overview(code, fills, walk, scale=2):
    rgb = PAL[code].copy()
    rgb[walk] = WALK
    for f in fills.values():
        rgb[f] = FILL
    im = Image.fromarray(rgb[::scale, ::scale])
    d = ImageDraw.Draw(im)
    try:
        font = ImageFont.load_default(size=14)
    except TypeError:
        font = ImageFont.load_default()
    for n, f in fills.items():
        ys, xs = np.nonzero(f)
        cy, cx = ys.mean() / scale, xs.mean() / scale
        r = 9
        d.ellipse([cx - r, cy - r, cx + r, cy + r], outline=FILL, width=2)
        d.text((cx + r + 3, cy - 8), str(n), fill=(255, 255, 255), font=font)
    b = io.BytesIO()
    im.save(b, "PNG", optimize=True)
    return "data:image/png;base64," + base64.b64encode(b.getvalue()).decode()


def native_crop(raw, px, box, ch, disp_um=DISP_UM):
    """the native image over grid box (y0, y1, x0, x1), resampled (area average) to exactly
    disp_um per pixel, i.e. (y1 - y0) * G / disp_um rows: the same frame as the label panels.
    (Until 2026-11-19 this block-averaged by a whole number of native pixels, 1.95 um per
    pixel at 10x and 2.05 at 20x, so images were 2.5 % off the labels' scale.)"""
    y0, y1, x0, x1 = box
    f = G / px
    a = raw[ch][int(round(y0 * f)):int(round(y1 * f)), int(round(x0 * f)):int(round(x1 * f))].astype(np.float32)
    hh, ww = int(round((y1 - y0) * G / disp_um)), int(round((x1 - x0) * G / disp_um))
    a = np.asarray(Image.fromarray(a, mode="F").resize((ww, hh), Image.BOX))
    lo, hi = np.percentile(a, [1, 99.5])
    return np.clip((a - lo) / max(hi - lo, 1) * 255, 0, 255)


def cut_mask(code, f):
    """white matter the fill f cuts off the rest of its piece (pieces of 2000 um2 or more)."""
    wm = code == 5
    lab2, sizes, groups = GF._split(wm, f, GF.S8)
    out = np.zeros(code.shape, bool)
    for i in groups:
        if i.size > 1:
            for j in i[np.argsort(sizes[i])[::-1][1:]]:
                if sizes[j] * G * G >= 2000:
                    out |= lab2 == j + 1
    return out


def crop_panels(code, raw, px, ch, f, pts, cut=None):
    ys, xs = np.nonzero(f | (cut if cut is not None else False))
    pts = np.asarray(pts, float).reshape(-1, 2)
    by0 = int(max(min(ys.min(), pts[:, 0].min()) - PAD, 0))
    by1 = int(min(max(ys.max(), pts[:, 0].max()) + PAD + 1, code.shape[0]))
    bx0 = int(max(min(xs.min(), pts[:, 1].min()) - PAD, 0))
    bx1 = int(min(max(xs.max(), pts[:, 1].max()) + PAD + 1, code.shape[1]))
    sl = (slice(by0, by1), slice(bx0, bx1))
    kk = int(G / DISP_UM)
    before = PAL[code[sl]]
    after = before.copy()
    if cut is not None:
        after[cut[sl]] = CUT
    after[f[sl]] = FILL
    img = native_crop(raw, px, (by0, by1, bx0, bx1), ch)
    hh, ww = before.shape[0] * kk, before.shape[1] * kk
    img = np.pad(img, ((0, max(hh - img.shape[0], 0)), (0, max(ww - img.shape[1], 0))))[:hh, :ww]
    lab = "After (fill in cyan)" if cut is None else "Fill (cyan) · white matter it cuts off (red)"
    panels = [dict(key="image", label="NeuN", img=jpg(img)),
              dict(key="before", label="Labels before", img=png(up(before, kk))),
              dict(key="after", label=lab, img=png(up(after, kk)))]
    return panels, round((bx1 - bx0) * G), round((by1 - by0) * G)


def course_panel(code, sl, kk, courses, marker=None, fills=None, axes=None):
    """labels with the predicted courses (yellow), their stop / contact points, the marker
    (white ring) and, with axes, each end's main axis (magenta arrow) and its lettered best
    far sides (thin magenta line to a letter)."""
    rgb = PAL[code[sl]].copy()
    if fills is not None:
        rgb[fills[sl]] = FILL
    im = Image.fromarray(up(rgb, kk))
    d = ImageDraw.Draw(im)
    y0, x0 = sl[0].start, sl[1].start
    try:
        font = ImageFont.load_default(size=15)
    except TypeError:
        font = ImageFont.load_default()
    for ax in axes or []:
        ty, tx = ax["tip"]
        P0 = ((tx - x0 + 0.5) * kk, (ty - y0 + 0.5) * kk)
        uy, ux = ax["axis"]
        L = 90 / G * kk
        P1 = (P0[0] + ux * L, P0[1] + uy * L)
        d.line([P0, P1], fill=AXIS, width=3)
        hx, hy = -ux, -uy
        for s_ in (1, -1):
            d.line([P1, (P1[0] + (hx * 0.9 - s_ * hy * 0.5) * 12, P1[1] + (hy * 0.9 + s_ * hx * 0.5) * 12)], fill=AXIS, width=3)
        for c in ax["candidates"]:
            qy, qx = c["q"]
            Q = ((qx - x0 + 0.5) * kk, (qy - y0 + 0.5) * kk)
            d.line([P0, Q], fill=AXIS, width=1)
            d.ellipse([Q[0] - 11, Q[1] - 11, Q[0] + 11, Q[1] + 11], fill=(20, 20, 20), outline=AXIS, width=2)
            d.text((Q[0] - 5, Q[1] - 9), c["letter"], fill=(255, 255, 255), font=font)
    for c, ok in courses:
        pts = [((x - x0 + 0.5) * kk, (y - y0 + 0.5) * kk) for y, x in c]
        if len(pts) > 1:
            d.line(pts, fill=COURSE_OK if ok else COURSE, width=2)
        if pts:
            ex, ey = pts[-1]
            d.ellipse([ex - 4, ey - 4, ex + 4, ey + 4], outline=COURSE_OK if ok else COURSE, width=2)
            sx, sy = pts[0]
            d.rectangle([sx - 3, sy - 3, sx + 3, sy + 3], fill=COURSE_OK if ok else COURSE)
    if marker is not None:
        my, mx = (marker[0] - y0 + 0.5) * kk, (marker[1] - x0 + 0.5) * kk
        d.ellipse([mx - 9, my - 9, mx + 9, my + 9], outline=(255, 255, 255), width=2)
    b = io.BytesIO()
    im.save(b, "PNG", optimize=True)
    return "data:image/png;base64," + base64.b64encode(b.getvalue()).decode()


def judged():
    tag = "v2" if (OUT / "judged_v2.json").exists() else "v1"
    meta = json.loads((OUT / f"judged_{tag}.json").read_text())
    Mk = np.load(OUT / f"judged_{tag}.npz")
    return meta, {m["key"]: Mk[m["key"]] for m in meta}


def prior_verdicts(sid, f, meta, masks):
    """earlier verdicts on this fill: judged fills of which this one paints at least 30 % of
    the part no fill with the opposite verdict shares (as igl_gapfill_score scores them)."""
    out = {}
    good = {"correct": True, "join": True, "wrong": False, "refuse": False}
    for m in meta:
        if m["sid"] == sid and m["verdict"]:
            jm = masks[m["key"]].copy()
            for o in meta:
                if o["sid"] == sid and o["verdict"] and good[o["verdict"]] != good[m["verdict"]]:
                    jm &= ~masks[o["key"]]
            ov = float((f & jm).sum() / max(jm.sum(), 1))
            if ov >= 0.3:
                out[m["id"]] = dict(id=m["id"], verdict=m["verdict"], note=m.get("note"))
    return list(out.values())


def bridge_items(A, meta, masks):
    """bridges of the current run: s4_A / s4_C from igl_finalize (bleach boxes + markers, what is
    written), every other section from igl_gapfill_apply."""
    J = json.loads((FIXED / "igl_fixed.json").read_text())
    items, overviews = [], {}
    for sid, e in A.items():
        fin = sid in J
        Z = np.load(OUT / f"{sid}_gapfill.npz")
        code = Z["code"]
        if fin:
            BZ = np.load(FIXED / f"{sid}_bridges.npz")
            rows = [(b, BZ[f"fill_{i}"]) for i, b in enumerate(J[sid]["bridges"])]
        else:
            rows = []
            for b in e["candidates"]:
                if b["status"] == "joined":
                    y0, y1, x0, x1 = [int(v) for v in Z[f"box_{b['id']}"]]
                    f = np.zeros(code.shape, bool)
                    f[y0:y1, x0:x1] = Z[f"fill_{b['id']}"]
                    rows.append((b, f))
        if not rows:
            continue
        px, _, big = GA.info(sid)
        raw = tifffile.imread((GA.EV.SB20 if big else GA.EV.SB10) / f"{sid}_native.tif")
        ch = 1 if big else 2
        fills = {}
        for n, (b, f) in enumerate(rows, 1):
            fills[n] = f
            panels, w_um, h_um = crop_panels(code, raw, px, ch, f, [b["a"], b["b"]])
            items.append(dict(id=f"p-{e['short']}-{n}", sid=sid, short=e["short"], n=n, group=e["group"], method=b["method"],
                              found_by=b.get("found_by"), target=b.get("target"), chord_um=round(b["chord_um"]), curve_um=round(b.get("curve_um", 0)),
                              fill_um2=b["fill_um2"], wm_shallow=b.get("wm_shallow"), wm_frac=b.get("wm_frac"), neun_ratio=b.get("neun_ratio"),
                              width=b["width_um"], depth=b["depth_um"], written=bool(e["group"] == "control"),
                              source="run with your bleach boxes and markers (written)" if fin else "plain run",
                              prior=prior_verdicts(sid, f, meta, masks), panels=panels, w_um=w_um, h_um=h_um))
        overviews[sid] = overview(code, fills, Z["walk_paint"])
    order = {"control": 0, "cKO": 1}
    items.sort(key=lambda i: (order[i["group"]], i["short"], i["n"]))
    return items, overviews


def marker_items():
    J = json.loads((FIXED / "igl_fixed.json").read_text())
    T = json.loads((FIXED / "user_traces.json").read_text()) if (FIXED / "user_traces.json").exists() else {}
    out = []
    kk = int(G / DISP_UM)
    for sid in ("s4_A_p27", "s4_C_p27"):
        Z = np.load(OUT / f"{sid}_gapfill.npz")
        code = Z["code"]
        CZ = np.load(FIXED / f"{sid}_marker_courses.npz")
        BZ = np.load(FIXED / f"{sid}_bridges.npz")
        fills = np.zeros(code.shape, bool)
        for k in BZ.files:
            fills |= BZ[k]
        if (FIXED / f"{sid}_traces.npz").exists():
            TZ = np.load(FIXED / f"{sid}_traces.npz")
            for k in TZ.files:
                if k.endswith("_band") or k.endswith("_egl"):
                    fills |= TZ[k]
        px, _, big = GA.info(sid)
        raw = tifffile.imread((GA.EV.SB20 if big else GA.EV.SB10) / f"{sid}_native.tif")
        for name, mk in J[sid]["markers"].items():
            c = np.array(mk["centre"])
            traced = T.get(sid, {}).get(name)
            ask = not mk["bridges"] and not traced    # not joined: ask the user where it goes
            courses = [(CZ[p["key"]], bool(p.get("ok"))) for p in mk["predictions"] if p["key"] in CZ.files]
            pts = [c[None] + d for d in (np.array([[-37, -37]]), np.array([[37, 37]]))]
            for cc, _ in courses:
                pts.append(cc)
            axes, choices, n_l = [], [], 0
            if ask:
                for ax in mk.get("axes", []):
                    cands = []
                    for cd in ax["candidates"][:2]:
                        letter = "ABCDEFGH"[n_l]
                        n_l += 1
                        cands.append(dict(cd, letter=letter))
                        choices.append([letter, f"{letter}: end {ax['end']} to the IGL {cd['dist_um']} µm away "
                                                f"({cd['deg'][0]}° off its own axis, {cd['deg'][1]}° off the far IGL's)"])
                        pts.append(np.array([cd["q"]]))
                    axes.append(dict(ax, candidates=cands))
            P_ = np.vstack(pts)
            y0, x0 = np.maximum(P_.min(0).astype(int) - 12, 0)
            y1, x1 = np.minimum(P_.max(0).astype(int) + 13, code.shape)
            lim = 190 if ask else 160
            if y1 - y0 > lim:
                y0, y1 = int(max(c[0] - lim // 2, 0)), int(min(c[0] + lim // 2, code.shape[0]))
            if x1 - x0 > lim:
                x0, x1 = int(max(c[1] - lim // 2, 0)), int(min(c[1] + lim // 2, code.shape[1]))
            if traced:                               # same crop as when the user clicked
                y0, x0, y1, x1 = traced["crop"]
            sl = (slice(y0, y1), slice(x0, x1))
            img = native_crop(raw, px, (y0, y1, x0, x1), 1 if big else 2)
            hh, ww = (y1 - y0) * kk, (x1 - x0) * kk
            img = np.pad(img, ((0, max(hh - img.shape[0], 0)), (0, max(ww - img.shape[1], 0))))[:hh, :ww]
            preds = []
            for p in mk["predictions"]:
                if "kind" in p:
                    what = (f"ran {p['path_um']} µm into IGL at {p['depth']} µm below the EGL (it started at {p['own_depth']}), "
                            + ("an end facing back" if p["kind"] == "end" else "the side of a stretch of ribbon")
                            + (": continues the ribbon, joined" if p["ok"] else f": does not continue it (angle {round(np.degrees(np.arccos(min(abs(p['cos']), 1))))}°), passed over"))
                else:
                    what = f"stopped after {p['path_um']} µm ({p['by']})"
                preds.append(f"end {p['end']}: {what}")
            if traced:
                ti = J[sid]["user_traces"][name]
                TZ = np.load(FIXED / f"{sid}_traces.npz")
                tb = TZ[f"{name}_band"] | TZ[f"{name}_egl"]
                br_ = TZ[f"{name}_bright"] if f"{name}_bright" in TZ.files else np.zeros_like(tb)
                r5 = T[sid][name].get("round5")
                yy_, xx_ = np.nonzero(tb | (br_ if not r5 else False))
                y0, y1 = int(max(min(yy_.min(), c[0]) - 60, 0)), int(min(max(yy_.max(), c[0]) + 61, code.shape[0]))
                x0, x1 = int(max(min(xx_.min(), c[1]) - 60, 0)), int(min(max(xx_.max(), c[1]) + 61, code.shape[1]))
                sl = (slice(y0, y1), slice(x0, x1))
                imgs = []
                for chn, lab_ in ((1, "Channel 1 (nuclei; click to add points)"), (2, "NeuN")):
                    im_ = native_crop(raw, px, (y0, y1, x0, x1), chn)
                    hh, ww = (y1 - y0) * kk, (x1 - x0) * kk
                    imgs.append(dict(key=f"ch{chn}", label=lab_, img=jpg(np.pad(im_, ((0, max(hh - im_.shape[0], 0)), (0, max(ww - im_.shape[1], 0))))[:hh, :ww]),
                                     click=(chn == 1)))
                after = PAL[code[sl]].copy()
                after[fills[sl]] = FILL
                if not r5:
                    edge = br_ & ~ndi.binary_erosion(br_)
                    after[edge[sl]] = (255, 214, 10)
                prev = [f"you said: {T[sid][name]['round4']['verdict'] or 'no verdict'}, “{T[sid][name]['round4']['note']}”"] if T[sid][name].get("round4") else []
                before = PAL[code[sl]].copy()
                bim = Image.fromarray(up(before, kk))
                dr = ImageDraw.Draw(bim)
                for cy_, cx_ in (r5["clicks_grid"] if r5 else []):
                    X, Y = (cx_ - x0 + 0.5) * kk, (cy_ - y0 + 0.5) * kk
                    dr.ellipse([X - 5, Y - 5, X + 5, Y + 5], fill=(255, 47, 214), outline=(255, 255, 255))
                bb = io.BytesIO()
                bim.save(bb, "PNG", optimize=True)
                before_img = "data:image/png;base64," + base64.b64encode(bb.getvalue()).decode()
                if r5:
                    prev = [f"you said: needs more, \u201c{r5['note']}\u201d"]
                    st = ti.get("strokes", [])
                    prev.append("painted your strokes in click order over every label: " + "; ".join(
                        f"{x['length_um']} µm ({', '.join(f'{v:,} µm² {k.upper()}' for k, v in x['painted_um2'].items() if v)})" for x in st))
                out.append(dict(id=f"mk6-{sid[:4]}-{name}", short=sid[:4], name=name, group="marker", ask=True, traced=True,
                                choices=[["right", "Right as painted"], ["extend", "Needs more: I clicked along it"], ["wrong", "Wrong"]],
                                crop=[y0, x0, y1, x1], bridges=[], predictions=prev + [
                                    f"your trace: {ti['clicks']} clicks, {ti['length_um']} µm long, painted {ti['width_um']} µm wide"
                                    + (", carried on to the IGL at its end" if ti.get("joined_to_igl") else f", its end {ti['end_to_igl_um']} µm from the nearest IGL")
                                    ],
                                panels=imgs + [dict(key="before", label="Labels before, with your last clicks (magenta)", img=before_img),
                                               dict(key="after", label="After: your traces as IGL (cyan)", img=png(up(after, kk)))],
                                w_um=round((x1 - x0) * G), h_um=round((y1 - y0) * G)))
                continue
                ti = J[sid]["user_traces"][name]
                out.append(dict(id=f"mk3-{sid[:4]}-{name}", short=sid[:4], name=name, group="marker", ask=True, traced=True,
                                choices=[["right", "Painted right"], ["extend", "Needs more: I added clicks"], ["wrong", "Wrong"]],
                                crop=[int(y0), int(x0), int(y1), int(x1)], bridges=[], predictions=[
                                    f"your trace: {ti['clicks']} clicks, {ti['length_um']} µm long, painted {ti['width_um']} µm wide over "
                                    + ", ".join(f"{v:,} µm² {k.upper() if k in ('ml', 'wm', 'egl') else k}" for k, v in ti["painted_um2"].items() if v)
                                    + f"; its end is {ti['end_to_igl_um']} µm from the nearest IGL"],
                                panels=[dict(key="image", label="NeuN with your clicks (click to add more)", img=jpg(img), click=True),
                                        dict(key="before", label="Labels before", img=course_panel(code, sl, kk, [], c)),
                                        dict(key="after", label="After: your trace painted as IGL (cyan)", img=course_panel(code, sl, kk, [], c, fills))],
                                w_um=round((x1 - x0) * G), h_um=round((y1 - y0) * G)))
                continue
            panels = [dict(key="image", label="NeuN" + (" (click where the IGL continues)" if ask else ""), img=jpg(img), click=ask),
                      dict(key="course", label="Predicted courses (yellow; green = joined), marker (white ring)",
                           img=course_panel(code, sl, kk, courses, c))]
            if ask:
                panels.append(dict(key="axes", label="Main axis of each end (magenta) and its best far sides",
                                   img=course_panel(code, sl, kk, [], c, None, axes)))
            else:
                panels.append(dict(key="after", label="After (fill in cyan)", img=course_panel(code, sl, kk, [], c, fills)))
            out.append(dict(id=(f"mk3-{sid[:4]}-{name}" if ask else f"mk-{sid[:4]}-{name}"), short=sid[:4], name=name, group="marker", ask=ask,
                            choices=choices, crop=[int(y0), int(x0), int(y1), int(x1)],
                            bridges=[f"{round(b['chord_um'])} µm ({b['method']})" for b in mk["bridges"]], predictions=preds,
                            panels=panels, w_um=round((x1 - x0) * G), h_um=round((y1 - y0) * G)))
    return out


def verdict_table(meta):
    S = {t: json.loads((OUT / f"score_{t}.json").read_text()) for t in ("hybrid", "band", "predict") if (OUT / f"score_{t}.json").exists()}
    notes = {}
    for m in meta:
        notes.setdefault(m["id"], m.get("note"))
    rows = []
    for r in S["hybrid"]["rows"]:
        if r["igl_there"] is None:
            continue
        shown = not r["id"].startswith("ref-")      # first-run bridges were painted, refused ones not
        row = dict(id=r["id"], verdict=r["verdict"], note=notes.get(r["id"]), igl_there=r["igl_there"], shown=shown)
        for t, sc in S.items():
            row[t] = next(x["filled"] for x in sc["rows"] if x["id"] == r["id"])
        rows.append(row)
    agree = {t: sum(r[t] == r["igl_there"] for r in rows) for t in list(S) + ["shown"]}
    return rows, agree


def main():
    A = json.loads((OUT / "apply.json").read_text())
    E = json.loads((OUT / "eval.json").read_text())
    P0 = json.loads((OUT / "eval_prelayout.json").read_text())["summary"]
    meta, masks = judged()
    items, overviews = bridge_items(A, meta, masks)
    markers = marker_items()
    vt, agree = verdict_table(meta)
    J = json.loads((FIXED / "igl_fixed.json").read_text())
    sections = []
    for sid, e in A.items():
        fin = sid in J
        nb = sum(1 for i in items if i["sid"] == sid)
        sec = dict(sid=sid, short=e["short"], group=e["group"], ends=e["fits"]["n_ends"], joined=nb,
                   painted_mm2=round(sum(i["fill_um2"] for i in items if i["sid"] == sid) / 1e6, 4),
                   walk=e["walk"], stamped=e["stamped"] or fin, source="with your boxes and markers" if fin else "plain run")
        if sid in overviews:
            sec["overview"] = overviews[sid]
            sec["overview_um_per_px"] = 2 * G
        sections.append(sec)
    order = {"control": 0, "cKO": 1}
    sections.sort(key=lambda s: (order[s["group"]], s["short"]))
    summ = E["summary"]
    ev = []
    for kind, lab in (("ML", "Gap's IGL labelled ML (NeuN still IGL-bright)"), ("anomaly", "Gap's IGL labelled anomaly (bleach-like)"),
                      ("WM", "Gap's IGL labelled white matter")):
        ev.append(dict(scenario=lab, method="Walk (original pipeline)", restored=summ[f"{kind} | walk"]["restored"], dice=summ[f"{kind} | walk"]["mean_dice"],
                       false_um2=summ[f"{kind} | walk"]["false_fill_um2_per_section"]))
        ev.append(dict(scenario=lab, method="Gap filler you judged", restored=P0[f"{kind} | gap filler"]["restored"], dice=P0[f"{kind} | gap filler"]["mean_dice"],
                       false_um2=P0[f"{kind} | gap filler"]["false_fill_um2_per_section"]))
        ev.append(dict(scenario=lab, method="Predicted ribbon + checks (now)", restored=summ[f"{kind} | gap filler"]["restored"], dice=summ[f"{kind} | gap filler"]["mean_dice"],
                       false_um2=summ[f"{kind} | gap filler"]["false_fill_um2_per_section"]))
    clean = {k: {m: v["painted_um2"] for m, v in c.items()} for k, c in E["clean"].items()}
    s_ = lambda k: summ[k]
    n_items = len(items)
    EA = {k: json.loads((OUT / f"{k}.json").read_text())["summary"] for k in ("eval_axis_both", "eval_no_axis_check", "eval_no_axis_search")
          if (OUT / f"{k}.json").exists()}
    ml = lambda d: f"{round(100 * d['ML | gap filler']['restored'])}%"
    Jt = J.get("s4_C_p27", {}).get("user_traces", {})
    Jt = J.get("s4_C_p27", {}).get("user_traces", {})
    added = J.get("s4_C_p27", {}).get("igl_added_mm2")
    findings = [
        "<b>The segment below the bleached area is painted.</b> Your two new paths trace the same U from both sides: from the bottom of the D4 trace down around the bend and up "
        "the E4 trace to the IGL (they run 13 µm apart, median, where both cover it). I painted every stroke in the order you clicked it, at the section's typical IGL width "
        "(48 µm), over every label, EGL included, since the labels there call it EGL and white matter; stroke ends within 60 µm of IGL or a trace are joined to it.",
        "<b>So D4 and E4 now form one ribbon</b> from the IGL on the D4 side, down the D4 trace, around the U and up the E4 trace into the IGL at the top. "
        f"Written into the s4_C label file (IGL added in s4_C by all repairs: {added} mm²).",
        "<b>Please check the two cards at the bottom</b>, zoomed out to take in the whole U: channel 1 (nuclei) and NeuN, the labels before with your last clicks, and the painted IGL.",
        f"<b>Unchanged:</b> the {len(items)} bridges you judged, {agree['hybrid']} of your 30 verdicts agreed with.",
    ]
    data = dict(sections=sections, items=items, markers=markers, verdicts=vt, agree=agree, eval=ev, clean=clean, gold=E["gold_s1_2_auto"],
                n_gaps=summ["ML | walk"]["n"], findings=findings)
    page = TEMPLATE.read_text().replace("__DATA__", json.dumps(data, default=float))
    dst = OUT / "igl_gap_fill.html"
    dst.write_text(page)
    print(dst, round(dst.stat().st_size / 1e6, 2), "MB", len(items), "bridges,", len(markers), "markers,", sum("overview" in s for s in sections), "overviews", agree)


if __name__ == "__main__":
    main()
