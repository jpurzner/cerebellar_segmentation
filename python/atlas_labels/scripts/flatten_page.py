"""
Page for the flattened cerebellum (2026-11-21): per section, the section overview with the
reference line, and the unrolled strips (nuclei / NeuN / p27) with the MATLAB or current label
call and the line fitted through the IGL. The user checks the flattening and flags spots.

    python scripts/flatten_page.py              # -> outputs/flatten/page/ (index.html, img/, data/)
    python scripts/flatten_page.py --repaired   # version 3, on EGL repair v6 (cuts on the pia) -> outputs/flatten_v3/page/
"""
import json
import sys
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage as ndi

sys.path.insert(0, str(Path(__file__).parent))
import flatten as F  # noqa: E402
import egl_repair as ER  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
REP = "--repaired" in sys.argv
OUTD = F.OUT2 if REP else F.OUT
PAGE = OUTD / "page"
TEMPLATE = Path(__file__).parent / "flatten_template.html"
LABSETS = [["before", "Before repair"], ["repaired", "Repaired"]] if REP else [["matlab", "MATLAB"], ["current", "Current"]]
MAIN = LABSETS[-1][0]
NPZ_KEY = {"matlab": "Lm", "current": "Lb", "before": "L_before", "repaired": "L_repaired"}
NOTES = {}                     # no notes were left on version 2 (db flatten_check_v2 empty); v1 notes: outputs/flatten/notes_v1.json
PREV_LINE = (255, 160, 0)      # the version 2 cut (EGL repair v5), drawn where the new cut is not


def stretch(v, mask, lo=1.0, hi=99.5):
    a, b = np.percentile(v[mask], [lo, hi]) if mask.any() else (0, 1)
    return np.clip((v - a) / max(b - a, 1e-6) * 255, 0, 255).astype(np.uint8)


def runs(m):
    return [[int(a), int(b)] for a, b in F._runs(np.asarray(m, bool))]


def ints(v, scale=1.0):
    return [int(round(x * scale)) if np.isfinite(x) else -1 for x in v]


def build(short, row, Z=None, regions=None, over_marks=(), top_marks=()):
    """one section's images and data; Z (default: the flatten npz) may be a dict with extra label sets;
    regions = [[a, b, geometry, type], ...] for error bars; over_marks = [(mask, rgb)] painted on the overview."""
    Z = Z if Z is not None else np.load(OUTD / f"{short}.npz")
    N, Dn = Z["valid"].shape
    valid, ok = Z["valid"], Z["ok"]
    m = valid & ok[:, None]
    files = {}
    for name, key in (("dapi", "S_nuclei"), ("neun", "S_neun"), ("p27", "S_p27")):
        img = stretch(Z[key].astype(np.float32), m).T
        img[~valid.T] = 0
        f = PAGE / "img" / f"{short}_{name}.jpg"
        Image.fromarray(img).save(f, quality=74, optimize=True)
        files[name] = f"img/{short}_{name}.jpg"
    for ls, _ in LABSETS:
        f = PAGE / "img" / f"{short}_lab_{ls}.png"
        Image.fromarray((Z[NPZ_KEY[ls]].T.astype(np.uint8) * 40)).save(f, optimize=True)
        files[f"lab_{ls}"] = f"img/{short}_lab_{ls}.png"
    # overview: nuclei (grey) + NeuN (green) at the label-grid resolution
    sid = F.sid_of(short)
    R, k, (Hn, Wn) = F.stains_reduced(sid)
    cut = Z["cut"]
    gh, gw = cut.shape
    tissue = cut | Z["mid"]
    zf = lambda a: ndi.zoom(a, (gh / a.shape[0], gw / a.shape[1]), order=1)
    nuc, neun = zf(R["nuclei"]), zf(R["neun"])
    nuc, neun = stretch(nuc, tissue), stretch(neun, tissue)
    rgb = np.stack([nuc * 0.85, np.maximum(nuc * 0.85, neun), nuc * 0.85], -1).astype(np.uint8)
    if REP and (F.PREV2 / f"{short}.npz").exists():
        old = np.load(F.PREV2 / f"{short}.npz")["mid"]
        rgb[old & ~Z["mid"]] = PREV_LINE
    for m_, c_ in over_marks:
        rgb[m_] = c_
    rgb[Z["mid"]] = (255, 60, 200)
    for m_, c_ in top_marks:                                   # drawn over the cuts (flat_fix_page v8: the new cuts)
        rgb[m_] = c_
    f = PAGE / "img" / f"{short}_over.jpg"
    Image.fromarray(rgb).save(f, quality=82, optimize=True)
    files["over"] = f"img/{short}_over.jpg"
    Q, X = Z["Q"], Z["X"]
    every, sub = 8, 5
    cols = []
    for i in range(0, N, every):
        L = int(valid[i].sum())
        pts = np.vstack([Q[i:i + 1], X[i, :L:sub]]) * F.G if L else Q[i:i + 1] * F.G
        cols.append([int(round(v)) for v in pts.ravel()])
    data = dict(short=short, group=row["group"], N=int(N), Dn=int(Dn), step=F.STEP, grid=[int(gh), int(gw)], grid_um=F.G,
                ok=runs(ok), pieces=[[int(a), int(b)] for a, b in zip(np.r_[0, np.flatnonzero(np.diff(Z["piece"])) + 1],
                                                                      np.r_[np.flatnonzero(np.diff(Z["piece"])) + 1, N])],
                Q=[int(round(v)) for v in (Q[::4] * F.G).ravel()], col_every=every, col_sub_um=sub * F.STEP, cols=cols,
                fits={}, flags={}, files=files, stats=row)
    if regions is not None:
        data["regions"] = regions
    for ls, _ in LABSETS:
        data["fits"][ls] = dict(top=ints(Z[f"{ls}_ftop"][::2]), bot=ints(Z[f"{ls}_fbot"][::2]), bridged=runs(Z[f"{ls}_bridged"]))
        data["flags"][ls] = {k2: runs(Z[f"{ls}_flag_{k2}"]) for k2 in ("missing", "off", "no_ml")}
    (PAGE / "data" / f"{short}.json").write_text(json.dumps(data, separators=(",", ":")))
    return dict(short=short, group=row["group"], N=int(N), data=f"data/{short}.json", stats=row, prev_note=NOTES.get(short, {}).get("note", ""))


def findings(rows):
    ctl = [r for s, r in rows.items() if r["group"] == "control" and s != "s1_2"]
    cko = [r for r in rows.values() if r["group"] == "cKO"]
    g = rows["s1_2"]
    pres = lambda r: r["current"]["igl_present_um"] / max(r["current"]["covered_um"], 1)
    rng = lambda v, f=lambda x: x: f"{f(min(v))}–{f(max(v))}"
    pc = lambda x: f"{round(x * 100)}"
    km = lambda x: f"{x / 1000:.1f}"
    return [
        f"Built for all {len(rows)} sections. The fissures are unrolled: before, about half of the EGL lay more than 60 µm inside the tissue "
        f"(the fissures are closed in the tissue mask); after cutting along the fissure midlines, {rng([r['egl_deep_after'] for r in rows.values()], pc)}% does.",
        f"<b>Controls:</b> IGL present along {rng([pres(r) for r in ctl], pc)}% of the reference line that has EGL under it, the fitted line in "
        f"{rng([r['current']['segments'] for r in ctl])} pieces (s1_2: {pc(pres(g))}%, {g['current']['segments']} pieces). "
        f"<b>cKO:</b> {rng([pres(r) for r in cko], pc)}%, {rng([r['current']['segments'] for r in cko])} pieces, and "
        f"{rng([r['current']['no_ml_um'] for r in cko], km)} mm per section where the IGL touches the EGL with no ML between "
        f"(controls {rng([r['current']['no_ml_um'] for r in ctl], km)} mm): the inner EGL being called IGL shows up as exactly this.",
        f"<b>Not trustworthy everywhere yet:</b> in the cKO sections {rng([r['current']['short_columns_um'] for r in cko], km)} mm of columns end "
        f"within 60 µm with no IGL (controls {rng([r['current']['short_columns_um'] for r in ctl], km)} mm). Most sit where the cKO EGL is a thick, "
        f"ragged band and the reference line is probably in the wrong place. Please flag those spots.",
        "The fitted line goes through the existing IGL calls (MATLAB or current), so for now it shows where those calls break. "
        "The next step fits the boundaries from the stains themselves, using the line as the expected course.",
        f"IGL top below the reference (median): {rng([r['current']['igl_top_med_um'] for r in ctl], lambda x: round(x))} µm in the 10x controls, "
        f"{round(g['current']['igl_top_med_um'])} µm in s1_2, {rng([r['current']['igl_top_med_um'] for r in cko], lambda x: round(x))} µm in the cKO "
        f"sections. s1_2 is the only control from the cKO batch, so any depth comparison has to lean on it.",
    ]


def findings_v2(rows):
    v1 = {r["section"]: r for r in json.loads((F.OUT / "summary.json").read_text())["sections"]}
    ctl = [s for s, r in rows.items() if r["group"] == "control" and s != "s1_2"]
    cko = [s for s, r in rows.items() if r["group"] == "cKO"]
    pres = lambda r, ls: r[ls]["igl_present_um"] / max(r[ls]["covered_um"], 1)
    pc = lambda x: f"{round(x * 100)}"
    km = lambda x: f"{x / 1000:.1f}"
    rng = lambda v, f=lambda x: x: f"{f(min(v))}–{f(max(v))}"
    new = lambda ss, f: rng([f(rows[s]) for s in ss])
    old = lambda ss, f: rng([f(v1[s]) for s in ss])
    drop = sorted(((rows[s]["cut_dropped_um"], s) for s in rows), reverse=True)[:3]
    worse = [s for s in cko if pres(rows[s], "repaired") < pres(v1[s], "current")]
    return [
        "<b>What is new.</b> The labels are the final EGL repair, and the reference line uses its fissure cuts: down the middle of the fused EGL the stains show, "
        "continued into the furrows. Your two notes on the first version (s1_2 missing the tip of the EGL, s4_A missing part of the EGL) led to that repair.",
        "<b>One fix to those cuts here.</b> In the cKO sections some ran along thick EGL at an open surface, leaving a strip of EGL with only background beyond it; "
        "columns into that strip stopped after about 40 µm with nothing under them. A stretch of cut whose columns on one side end within 60 µm with no ML, IGL or "
        "white matter is dropped, and the deeper part of that cut is joined to the surface by the shortest way out. Most dropped: "
        + ", ".join(f"{s} {km(u)} mm" for u, s in drop) + ". Cuts into a folium tip that has ML but no IGL in the section stay.",
        f"<b>Controls:</b> IGL present along {new(ctl, lambda r: pc(pres(r, 'repaired')))}% of the line with EGL under it (first version {old(ctl, lambda r: pc(pres(r, 'current')))}%); "
        f"s4_A {pc(pres(v1['s4_A'], 'current'))} → {pc(pres(rows['s4_A'], 'repaired'))}%, s4_F {pc(pres(v1['s4_F'], 'current'))} → {pc(pres(rows['s4_F'], 'repaired'))}%. "
        f"<b>cKO:</b> {new(cko, lambda r: pc(pres(r, 'repaired')))}% (first version {old(cko, lambda r: pc(pres(r, 'current')))}%); short columns "
        f"{new(cko, lambda r: km(r['repaired']['short_columns_um']))} mm per section (first version {old(cko, lambda r: km(r['current']['short_columns_um']))} mm).",
        f"<b>The repair itself changes little here.</b> On the same columns the labels before and after the repair give IGL present within a few points: "
        f"s4_C {pc(pres(rows['s4_C'], 'before'))} → {pc(pres(rows['s4_C'], 'repaired'))}% (IGL that was inner EGL or lay over the ML taken out), "
        f"s4_F {pc(pres(rows['s4_F'], 'before'))} → {pc(pres(rows['s4_F'], 'repaired'))}% (EGL that was IGL). Most of the gain comes from the new reference line. "
        "Switch Labels between Before repair and Repaired to compare.",
        f"<b>Still weak:</b> {', '.join(worse)} are below the first version. There, long columns run through thick EGL, unassigned and anomaly cells with no IGL "
        "under them (EGL masses), and short columns remain at fissure mouths in thick EGL and at folium tips. Please flag spots where the line itself is wrong.",
        f"IGL top below the line (median): {new(ctl, lambda r: round(r['repaired']['igl_top_med_um']))} µm in the 10x controls, "
        f"{round(rows['s1_2']['repaired']['igl_top_med_um'])} µm in s1_2, {new(cko, lambda r: round(r['repaired']['igl_top_med_um']))} µm in the cKO sections "
        f"(first version {old(cko, lambda r: round(r['current']['igl_top_med_um']))}): the repair added EGL, and the line now runs along the outer face of thick EGL "
        "instead of through it. "
        "s1_2 is the only control from the cKO batch.",
    ]


def findings_v3(rows):
    v2 = {r["section"]: r for r in json.loads((F.PREV2 / "summary.json").read_text())["sections"]}
    er = {r["section"]: r for r in json.loads((ER.OUT / "summary.json").read_text())["sections"]}
    er5 = {r["section"]: r for r in json.loads((ER.OUT / "npz_v5" / "summary.json").read_text())["sections"]}
    ctl = [s for s, r in rows.items() if r["group"] == "control" and s != "s1_2"]
    cko = [s for s, r in rows.items() if r["group"] == "cKO"]
    pres = lambda r: r["repaired"]["igl_present_um"] / max(r["repaired"]["covered_um"], 1)
    pc = lambda x: f"{round(x * 100)}"
    km = lambda x: f"{x / 1000:.1f}"
    rng = lambda v, f=lambda x: x: f"{f(min(v))}–{f(max(v))}"
    pia = {s: er[s]["cuts"]["pia"] for s in rows}
    ch = lambda s: f"{s} {pc(pres(v2[s]))} → {pc(pres(rows[s]))}%"
    fell = sum(er[s]["egl_added_um2"] < er5[s]["egl_added_um2"] for s in rows)
    return [
        "<b>What is new: the fissure cuts sit on the pia.</b> In a closed fissure two EGLs face each other across a thin layer of pia and arachnoid. "
        "As layer_quant splits its fissure profiles at the p27 minimum, each point of a closed cut now looks across the fissure for the lowest p27 + NeuN "
        "(the pia has no p27 or NeuN background) nearest the cut, with EGL on both sides and below the section's ML level, and the cut moves there, smoothed along it. "
        "Before, the cut ran down the middle of the fused EGL, which is the pia only where both EGLs are equally thick. "
        "<b>Orange</b> on the overview is the previous line where it differs.",
        f"The pia shows along {rng([pia[s]['found_frac'] for s in cko], pc)}% of the closed cuts in the cKO sections and {rng([pia[s]['found_frac'] for s in ctl + ['s1_2']], pc)}% in the controls; "
        f"elsewhere the cut stays where it was. Where it shows, the old cut was a median 6–10 µm off it and is now 3–4 µm off (the cut is 8 µm wide). "
        f"Cut moved by more than 4 µm: {rng([pia[s]['moved_over_4um_um'] for s in cko], km)} mm per cKO section, {rng([pia[s]['moved_over_4um_um'] for s in ctl + ['s1_2']], km)} mm per control.",
        "<b>Against your s1_2 tracing:</b> the unassigned strips you left between two EGLs (the pia) have their centreline within 4 µm of the cut 75% of the time, "
        "against 54% before. In places the old cut ran down both edges of such a strip; now it runs inside it.",
        f"<b>Flattening:</b> cKO IGL present {rng([pres(rows[s]) for s in cko], pc)}% (version 2: {rng([pres(v2[s]) for s in cko], pc)}%), "
        f"{ch('s3_5')}, {ch('s3_2')}, {ch('s2_2')}; short columns {rng([rows[s]['repaired']['short_columns_um'] for s in cko], km)} mm per cKO section "
        f"(version 2: {rng([v2[s]['repaired']['short_columns_um'] for s in cko], km)} mm). Controls {rng([pres(rows[s]) for s in ctl], pc)}% ({ch('s4_C')}). "
        f"The line is {rng([(v2[s]['length_um'] - rows[s]['length_um']) for s in rows], km)} mm shorter per section, with the same share of it over EGL.",
        f"<b>EGL repair:</b> the repair reruns on the new cuts. The EGL it adds fell in {fell} of {len(rows)} sections "
        f"(s1_4 {round(er5['s1_4']['egl_added_um2'] / 1000)} → {round(er['s1_4']['egl_added_um2'] / 1000)}k µm², s2_4 {round(er5['s2_4']['egl_added_um2'] / 1000)} → "
        f"{round(er['s2_4']['egl_added_um2'] / 1000)}k), probably because columns from the cut no longer start partway into one EGL, where it looked too thin; "
        f"s3_2 rose ({round(er5['s3_2']['egl_added_um2'] / 1000)} → {round(er['s3_2']['egl_added_um2'] / 1000)}k). The Labels switch compares the labels before the repair with the new repaired labels.",
        f"IGL top below the line (median): controls {rng([rows[s]['repaired']['igl_top_med_um'] for s in ctl], lambda x: round(x))} µm, "
        f"s1_2 {round(rows['s1_2']['repaired']['igl_top_med_um'])} µm, cKO {rng([rows[s]['repaired']['igl_top_med_um'] for s in cko], lambda x: round(x))} µm "
        f"(version 2: {rng([v2[s]['repaired']['igl_top_med_um'] for s in cko], lambda x: round(x))}). Please flag furrows where the line still runs inside one EGL.",
    ]


def main():
    (PAGE / "img").mkdir(parents=True, exist_ok=True)
    (PAGE / "data").mkdir(parents=True, exist_ok=True)
    summ = json.loads((OUTD / "summary.json").read_text())
    rows = {r["section"]: r for r in summ["sections"]}
    order = ["s1_2"] + sorted(s for s, r in rows.items() if r["group"] == "control" and s != "s1_2") + \
        sorted(s for s, r in rows.items() if r["group"] == "cKO")
    secs = []
    for s in order:
        secs.append(build(s, rows[s]))
        print("built", s, flush=True)
    data = dict(collection="flatten_check_v3" if REP else "flatten_check", sections=secs, prm=summ["prm"],
                findings=findings_v3(rows) if REP else findings(rows), labsets=LABSETS, main=MAIN)
    if REP:
        data["over_note"] = ("Nuclei grey, NeuN green; magenta: fissure cuts the line was cut along, now on the pia between the two EGLs; "
                             "<b style=\"color:#e69500\">orange</b>: the version 2 cut where it differs. The reference line is cyan where it has EGL under it, "
                             "grey where not; numbers are millimetres along it. Hover a strip to see that column here; click here to jump the strips to that spot.")
        data["lede"] = ("Third version: the fissure cuts now follow the pia between the two EGLs, as in layer_quant. Each section is unrolled along the outer face of the EGL: the reference line runs along "
                        "the surface, down one side of every fissure and back up the other. Columns run inward from it at right angles, so the layers "
                        "should come out as horizontal bands and a wrong IGL call shows up as a break in the band.")
        data["table_note"] = ("Numbers use the repaired labels, over the stretch of reference line with EGL under it. &ldquo;IGL present&rdquo;: share of "
                              "that stretch with IGL in the column. &ldquo;Line pieces&rdquo;: separate stretches of the fitted IGL line over 100 µm. "
                              "&ldquo;Short columns&rdquo;: columns that end within 60 µm with no IGL. Click a row to open it.")
    page = TEMPLATE.read_text().replace("__DATA__", json.dumps(data, separators=(",", ":")))
    (PAGE / "index.html").write_text(page)
    tot = sum(f.stat().st_size for f in PAGE.rglob("*") if f.is_file())
    print(f"page: {len(list(PAGE.rglob('*')))} files, {tot / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
