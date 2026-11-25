"""
Page for the error regions in the flattened cerebellum and their correction (flat_fix.py), version 12 of the
flattening page: the flattening strips (flatten v11: the line is one walk, clockwise from the left end of the cut
surface, with no sharp turns except at furrow bases and fissure mouths, corners cut by arcs, no jumps, no closed
slits; walk.py) with three label sets (before the
EGL repair, repaired, corrected), bars for the error regions (green: label errors, corrected; grey: geometry, columns
end before the IGL; violet: the line moved; red: a sharp turn left), and the corrected cells marked on the overview.

    python scripts/flat_fix_page.py      # -> outputs/flat_fix/page/
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import flatten as F  # noqa: E402
import flatten_page as FP  # noqa: E402
import flat_fix as FF  # noqa: E402

PAGE = FF.OUT / "page"
NOTES = json.loads((FF.OUT / "notes_v11.json").read_text()) if (FF.OUT / "notes_v11.json").exists() else {}
PREV_FIX = FF.ROOT / "outputs/flat_fix_v8"        # version 11 of the page (flat_fix on flatten_v10)
ARACH_RGB = (150, 120, 255)                       # outline of the meninges the shape check trimmed (v10)
FILL_RGB = (255, 255, 255)                        # outline of the cracks the shape check filled (v10)
TURN_RGB = (255, 64, 64)                          # sharp turns left (v10)
PREV_CUT = (255, 160, 0)                          # the version 7 line, drawn where the new line is not
PEAK_RGB = (255, 235, 0)                          # the new cuts down to the furrow base (v8; v11: also cuts continued from their tips)
CAP_RGB = (255, 120, 200)                         # thin ML caps relabelled EGL (v11)
LABSETS = [["before", "Before repair"], ["repaired", "Repaired"], ["corrected", "Corrected"]]
FIX_RGB = (60, 255, 120)


def km(x):
    return f"{x / 1000:.1f}"


def section(short, row, fx):
    Z = dict(np.load(F.OUT2 / f"{short}.npz"))
    R = np.load(FF.OUT / f"{short}.npz")
    Z["L_corrected"] = R["L_new"]
    if "ok_new" in R.files:
        Z["ok"] = R["ok_new"]                                # not-cortex stretches leave the EGL-covered line
    line = F.igl_line(R["L_new"], Z["valid"], Z["ok"], F.PRM)
    for k, v in line.items():
        if isinstance(v, np.ndarray):
            Z[f"corrected_{k}"] = v
    for k, v in line["flags"].items():
        Z[f"corrected_flag_{k}"] = v
    regs = json.loads((FF.OUT / f"{short}_regions.json").read_text())
    regions = [[r["a"], r["b"], 2 if r.get("noncortex") else int(bool(r.get("geometry"))),
                r["type"] + ("" if r.get("geometry") else f" ({r['um']} µm, {round(r['changed_um2'] / 1000, 1)}k µm² changed)")] for r in regs]
    row = dict(row, corrected=line["stats"], extra=dict(err=km(fx["error_um"]) + " mm", chg=f"{round(fx['changed_um2'] / 1000)}k µm²",
                                                         geo=km(fx["geometry_um"]) + " mm", ncx=km(fx.get("noncortex_um", 0)) + " mm"))
    marks = [(R["changed"], FIX_RGB)]
    from scipy import ndimage as ndi
    import shape_check as SC
    top = [(Z["peak_add"].astype(bool) & Z["mid"], PEAK_RGB)] if "peak_add" in Z else []
    if "tip_add" in Z:
        top.append((Z["tip_add"].astype(bool) & Z["mid"], PEAK_RGB))
    for key in ("cap", "tip_band"):
        if key in Z:
            marks.append((Z[key].astype(bool), CAP_RGB))
    H, W = Z["mid"].shape
    qi = np.clip(np.round(Z["Q"]).astype(int), 0, [H - 1, W - 1])
    for key, rgb, what in (("shape_trim", ARACH_RGB, "tissue trimmed"), ("shape_fill", FILL_RGB, "gap filled")):
        m = Z[key].astype(bool)
        marks.append((m & ~ndi.binary_erosion(m), rgb))
        near = ndi.distance_transform_edt(~m)[qi[:, 0], qi[:, 1]] <= 3 if m.any() else np.zeros(len(qi), bool)
        for a, b in F._runs(near):
            regions.append([int(a), int(b), 3, f"line moved (shape check or walk), {what}"])
    turns = np.zeros((H, W), bool)
    for e in SC.events_saved(Z):
        if e["kind"] in ("open", "cut"):
            regions.append([int(e["a"]), int(e["b"]), 4, f"sharp turn {e['deg']:+d}° {'along a cut' if e['kind'] == 'cut' else 'on the surface'}, not fixed".replace("-", "−")])
            r_, c_ = qi[e["i"]]
            turns[max(r_ - 1, 0):r_ + 2, max(c_ - 1, 0):c_ + 2] = True
    top.append((turns, TURN_RGB))
    if (F.PREV2 / f"{short}.npz").exists():                  # the last version's line, where the new one is not
        H, W = Z["mid"].shape
        raster = lambda Q: (lambda m, q: (m.__setitem__((np.clip(q[:, 0], 0, H - 1), np.clip(q[:, 1], 0, W - 1)), True), m)[1])(
            np.zeros((H, W), bool), np.round(Q).astype(int))
        old = raster(np.load(F.PREV2 / f"{short}.npz")["Q"])
        marks.append((old & (ndi.distance_transform_edt(~raster(Z["Q"])) > 2), PREV_CUT))
    out = FP.build(short, row, Z=Z, regions=regions, over_marks=marks, top_marks=top)
    out["prev_note"] = NOTES.get(short, {}).get("note", "")
    return out


def _peaks(run_dir, fix_dir, short):
    Z, R = np.load(run_dir / f"{short}.npz"), np.load(fix_dir / f"{short}.npz")
    pk, exc = FF.egl_peaks(R["L_new"], Z["valid"], R["ok_new"])
    return Z, pk, exc


def _mark_stats(Z, L, i, half=100):
    """within +-200 um of line point i: line length (um) whose columns start in unassigned / anomaly cells, and whose
    columns reach no IGL or white matter."""
    a, b = max(0, i - half), min(len(L), i + half)
    v, LL = Z["valid"][a:b], L[a:b]
    unas = (np.isin(LL[:, :4], (1, 6))).all(1) & v[:, 0]
    reach = (np.isin(LL, (4, 5)) & v).any(1)
    return int(unas.sum() * F.STEP), int((~reach).sum() * F.STEP)


def _top_egl(Z, L, i, half=50):
    """within +-100 um of line point i: share of the columns with an EGL label in their top 8 um."""
    a, b = max(0, i - half), min(len(L), i + half)
    v, LL = Z["valid"][a:b], L[a:b]
    return float(((LL[:, :4] == 2) & v[:, :4]).any(1).mean())


def _perimeter(Z):
    """steps of the line longer than 4 um (the points are 2 um apart) within one loop, those that pass over a cut cell,
    and the number of loops."""
    Q, piece, mid = Z["Q"].astype(float), Z["piece"], Z["mid"]
    d = np.hypot(*np.diff(Q, axis=0).T) * F.G
    j = np.where((d > 4.0) & (piece[1:] == piece[:-1]))[0]
    t = np.linspace(0, 1, 12)[:, None]
    hop = [int(mid[tuple(np.round(Q[i] + t * (Q[i + 1] - Q[i])).astype(int).T)].any()) for i in j]
    return dict(jumps=len(j), hops=int(sum(hop)), loops=len(np.unique(piece)), max_um=float(d[piece[1:] == piece[:-1]].max()),
                hop_at=[int(i) for i, h in zip(j, hop) if h], jump_um=[float(d[i]) for i in j])


def findings(rows, fx, ev):
    """version 12: the user's note on version 11 (s5_E: a discontinuity, the piece should come right after 13 mm, "a
    failure to maintain the perimeter"; find other examples and fix the problem), what was wrong, the fix, tests."""
    import shape_check as SC
    from collections import Counter
    ctl = ["s1_2"] + [s for s, r in rows.items() if r["group"] == "control" and s != "s1_2"]
    cko = [s for s, r in rows.items() if r["group"] == "cKO"]
    pc = lambda x: f"{round(x * 100)}"
    sg = lambda x: f"{round(x * 100):+d}".replace("-", "−")
    pres = lambda r, ls: r[ls]["igl_present_um"] / max(r[ls]["covered_um"], 1)
    wk = {s: rows[s]["walk"] for s in rows}
    old = {s: _perimeter(np.load(F.PREV2 / f"{s}.npz")) for s in rows}
    new = {s: _perimeter(np.load(F.OUT2 / f"{s}.npz")) for s in rows}
    nj = sum(old[s]["jumps"] for s in rows)
    nh = sum(old[s]["hops"] for s in rows)
    with_j = [s for s in rows if old[s]["jumps"]]
    cnt = lambda Z: {kk: sum(e["kind"] == kk for e in SC.events_saved(Z)) for kk in ("open", "cut")}
    t11 = {s: cnt(np.load(F.PREV2 / f"{s}.npz")) for s in rows}
    t12 = {s: cnt(np.load(F.OUT2 / f"{s}.npz")) for s in rows}
    t0 = {s: wk[s]["turns_before"] for s in rows}
    ts = lambda T, grp, kk: sum(T[s][kk] for s in grp)
    tri = lambda grp, kk: f"{ts(t0, grp, kk)} → {ts(t11, grp, kk)} → {ts(t12, grp, kk)}"
    why = Counter()
    for s in rows:
        for x in wk[s]["left"]:
            w = x["why"]
            why["need more than 12 µm" if w.startswith("needs") else "have a furrow base or fissure mouth within reach" if w.startswith("furrow")
                else "would cross a cut" if w.startswith("would") else "have no arc that removes them"] += 1
    arcs = sum(wk[s]["smoothed"] for s in rows)
    moved = max(wk[s]["moved_max_um"] for s in rows)
    Zo, Zn = np.load(F.PREV2 / "s5_E.npz"), np.load(F.OUT2 / "s5_E.npz")
    base = np.array([309, 292])
    from scipy.spatial import cKDTree
    jb = int(cKDTree(Zn["Q"]).query(base)[1])
    tail = int((Zo["piece"] == 1).sum() * F.STEP)
    old_ev = {r["section"]: r for r in json.loads((PREV_FIX / "eval.json").read_text())["sections"]}
    fx11 = {r["section"]: r for r in json.loads((PREV_FIX / "summary.json").read_text())["sections"]}

    def prev_pres(s):                                 # IGL present on the version 11 page (corrected labels)
        Z, R = np.load(F.PREV2 / f"{s}.npz"), np.load(PREV_FIX / f"{s}.npz")
        st = F.igl_line(R["L_new"], Z["valid"], R["ok_new"], F.PRM)["stats"]
        return st["igl_present_um"] / max(st["covered_um"], 1)
    dp = {s: pres(rows[s], "corrected") - prev_pres(s) for s in rows}
    evs = lambda s: f"{pc(old_ev[s]['recall'])} / {pc(old_ev[s]['precision'])}% → {pc(ev[s]['recall'])} / {pc(ev[s]['precision'])}%"
    s32 = [min(old["s3_2"]["hop_at"], key=lambda i: abs(i - 4153))]   # the one checked: 8.3 mm on the version 11 strip
    return [
        "<b>Your notes on version 11:</b> s1_2 and s5_F right; in s5_E a discontinuity at 17.7 mm, a piece that should come right after "
        "13 mm. <b>Cause:</b> one fissure cut in s5_E had a 24 µm gap. When the cuts are pruned, a cut piece left without a way out is "
        "joined to the surface by the cheapest path; here two pieces shared one such path, and when the short piece near the surface was "
        "dropped (its path turned sharply off it), the shared path went with it. The deeper 0.5 mm of the cut was left as a closed slit "
        f"inside the tissue. The outline round a slit is a loop of its own that the walk cannot reach, so it was put at the end of the strip "
        f"(the last {tail / 1000:.1f} mm). <b>Fix:</b> every cut must reach the open surface; a stranded piece is joined back by the cheapest "
        "path through EGL (as before), or dropped if there is none. This runs after the pruning and again just before the walk. In s5_E "
        f"the line now goes down that cut after 13.5 mm, rounds its base at {jb * F.STEP / 1000:.1f} mm and is back at the surface by 14.4 mm; "
        f"every section is one loop ({sum(new[s]['loops'] for s in rows)} loops in 19 sections).",
        "<b>Other failures to keep to the perimeter</b> (looked for in all 19 sections): (1) closed slits: only the one in s5_E. "
        f"(2) <b>Jumps</b>: the version 11 walk stepped 5–{max(max(old[s]['jump_um'], default=0) for s in rows):.0f} µm between neighbouring "
        f"points (they are 2 µm apart) in {len(with_j)} of 19 sections, {nj} places in all. To cut a corner it stiffened the line along a "
        "stretch, and the stiffness fell to nothing in one step at each end of the stretch, so the last stiff point was pulled off the "
        "outline while its neighbour stayed on it: the walk tore there. "
        f"(3) <b>Hops over cuts</b>: {nh} of those jumps crossed a cut, the walk leaving one wall of a fissure for the other without going "
        f"round its base (s3_2 at {s32[0] * F.STEP / 1000:.1f} mm on the version 11 strip skipped the last 50 µm of a cut). The version 11 "
        "counts of turns fixed included these: a tear or a hop spreads or skips a turn.",
        "<b>The walk now cuts corners with arcs.</b> Each sharp turn, worst first, is replaced by the shortest smooth arc that leaves and "
        "rejoins the outline along its own direction, with points 2 µm apart like the rest of the walk. The arc may move the line at most "
        "12 µm, may not sweep over a cut, and may not take in a furrow base or a fissure mouth; the whole walk is checked again after each "
        f"arc. {arcs} arcs in all, moving the line up to {moved:.0f} µm. Largest step now "
        f"{max(new[s]['max_um'] for s in rows):.1f} µm in every section, no hops. "
        f"Sharp turns, outline before the walk → version 11 → now: open surface controls {tri(ctl, 'open')}, cKO {tri(cko, 'open')}; "
        f"along cuts controls {tri(ctl, 'cut')}, cKO {tri(cko, 'cut')}. "
        "More are left than version 11 showed, because no arc may cross a cut or tear the walk. Of those left, "
        + ", ".join(f"{v} {k}" for k, v in why.most_common()) + ". "
        "The ones along cuts are mostly corners a line round the cuts has to make: the tip of a wedge of tissue where two cuts meet, the "
        "bottom of a furrow the pia trim opened (the cut runs down each wall), a cut meeting the surface at a slant. The ones on the "
        "surface needing more than 12 µm are knobs and notches 15–35 µm deep. Red marks show them all.",
        "<b>Still open from version 10</b> (the line is not the problem there; the EGL labels are): s3_2 at 8.2 mm, s3_5 at 11.8 and "
        "13.3 mm, s2_4 at 14.1 mm, s2_2 at 9.6 mm (positions as on the version 11 page).",
        f"<b>Tested.</b> Eval against your corrections, recall / precision (version 11 → 12): gold s1_2 {evs('s1_2')}, s4_C {evs('s4_C')}, "
        f"s5_E {evs('s5_E')}. IGL present (corrected labels): controls {sg(min(dp[s] for s in ctl))} to "
        f"{sg(max(dp[s] for s in ctl))} points; cKO " + ", ".join(f"{s} {sg(dp[s])}" for s in cko) + ". "
        f"Geometry errors (columns end before the IGL): cKO " + ", ".join(f"{s} {km(fx11[s]['geometry_um'])} → {km(fx[s]['geometry_um'])}" for s in cko)
        + " mm. Please mark anywhere the line still leaves the perimeter, and any red turn that is a real jog.",
    ]


def main():
    (PAGE / "img").mkdir(parents=True, exist_ok=True)
    (PAGE / "data").mkdir(parents=True, exist_ok=True)
    FP.PAGE, FP.LABSETS, FP.MAIN, FP.REP = PAGE, LABSETS, "corrected", False
    FP.NPZ_KEY = dict(FP.NPZ_KEY, corrected="L_corrected")
    summ = json.loads((F.OUT2 / "summary.json").read_text())
    rows = {r["section"]: r for r in summ["sections"]}
    fx = {r["section"]: r for r in json.loads((FF.OUT / "summary.json").read_text())["sections"]}
    ev = {r["section"]: r for r in json.loads((FF.OUT / "eval.json").read_text())["sections"]}
    order = ["s1_2"] + sorted(s for s, r in rows.items() if r["group"] == "control" and s != "s1_2") + \
        sorted(s for s, r in rows.items() if r["group"] == "cKO")
    secs = []
    for s in order:
        secs.append(section(s, rows[s], fx[s]))
        rows[s] = secs[-1]["stats"]
        print("built", s, flush=True)
    data = dict(collection="flatten_check_v12", shape_bars=True, sections=secs, prm=dict(flatten=summ["prm"], flat_fix=FF.PRM),
                findings=findings(rows, fx, ev), labsets=LABSETS, main="corrected", region_bars=True,
                extra_cols=[["Label errors", "err"], ["Changed", "chg"], ["Geometry", "geo"], ["Not cortex", "ncx"]])
    data["lede"] = ("Twelfth version, after your note on version 11: the walk keeps to the perimeter. It is one loop in every section "
                    "(the s5_E cut that was left as a closed slit reaches the surface again), it never steps more than 2 µm, and it "
                    "cuts sharp corners with smooth arcs that never cross a cut. "
                    "Each section is unrolled along the outer face of the EGL; columns run inward at right angles, so the layers come out as horizontal bands. "
                    "Where the bands break, the labels are compared with a layer-by-layer reading of the stains and the neighbouring columns, and corrected.")
    data["table_note"] = ("Numbers use the corrected labels, over the stretch of reference line with EGL under it. &ldquo;Label errors&rdquo;: length of line in "
                          "corrected error regions; &ldquo;Changed&rdquo;: area relabelled; &ldquo;Geometry&rdquo;: length where the columns end before the IGL "
                          "(flagged, not corrected). Click a row to open it.")
    data["over_note"] = ("Nuclei grey, NeuN green; <b style=\"color:#2fbf5f\">green</b>: cells the correction relabelled; magenta: fissure cuts; yellow: cuts down to the furrow base; "
                         "violet: outline of tissue trimmed by the shape check or the walk; white: gaps filled; pink: thin ML caps relabelled EGL; red dots: sharp turns left; orange: the last version\'s line where the new one runs elsewhere. The reference line is cyan "
                         "where it has EGL under it, grey where not; numbers are millimetres along it. Hover a strip to see that column here; click here to jump the strips there.")
    data["legend_extra"] = ('<span><i class="line" style="background:#4dff88"></i>label error region, corrected (hover for the type)</span>'
                            '<span><i class="line" style="background:#c8d0da"></i>geometry error: columns end before the IGL (not corrected)</span>'
                            '<span><i class="line" style="background:#ff9a5a"></i>not cortex (deep nuclei boundary, cut edge): EGL label removed</span>'
                            '<span><i class="line" style="background:#9678ff"></i>the line moved (shape check or walk: tissue trimmed or gap filled)</span>'
                            '<span><i class="line" style="background:#ff4040"></i>sharp turn left (hover for the angle)</span>'
                            '<span><i style="background:#ffa000"></i>overview: the last version\'s line, where the new one runs elsewhere</span>'
                            '<span><i style="background:#ffeb00"></i>overview: cut down to the furrow base or continued from its tip</span>'
                            '<span><i style="background:#9678ff"></i>overview: tissue trimmed (outline)</span>'
                            '<span><i style="background:#ff78c8"></i>overview: thin ML cap relabelled EGL</span>')
    (PAGE / "index.html").write_text(FP.TEMPLATE.read_text().replace("__DATA__", json.dumps(data, separators=(",", ":"))))
    tot = sum(f.stat().st_size for f in PAGE.rglob("*") if f.is_file())
    print(f"page: {len(list(PAGE.rglob('*')))} files, {tot / 1e6:.1f} MB")
    for f in data["findings"]:
        print("-", f[:400])


if __name__ == "__main__":
    main()
