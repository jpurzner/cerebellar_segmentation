"""
EGL repair page (2026-11-21): per section, the IGL (current labels), the new reference line (fissure
cuts placed in the stain-based EGL, see egl_repair.py), the EGL labels and the EGL the repair adds
from the nuclei stain where the IGL shows cortex below and the EGL label is missing or thin.

    python scripts/egl_repair_page.py      # -> outputs/egl_repair/page/
"""
import json
import sys
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage as ndi

sys.path.insert(0, str(Path(__file__).parent))
import egl_repair as R  # noqa: E402
import flatten as F  # noqa: E402
from igl_reference_page import overlay, stretch  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
PAGE = R.OUT / "page"
TEMPLATE = Path(__file__).parent / "egl_repair_template.html"
NOTES_FILE = R.OUT / "notes_v4.json"          # the user's verdicts and notes on the last version of the page
NOTES = json.loads(NOTES_FILE.read_text()) if NOTES_FILE.exists() else {}


def build(short, row):
    Z = np.load(R.OUT / f"{short}.npz")
    auto, best, tissue = F.labels(short)
    gh, gw = tissue.shape
    Rr, k, _ = F.stains_reduced(F.sid_of(short))
    nuc = ndi.zoom(Rr["nuclei"], (2 * gh / Rr["nuclei"].shape[0], 2 * gw / Rr["nuclei"].shape[1]), order=1)
    t2 = np.repeat(np.repeat(tissue, 2, 0), 2, 1)[:nuc.shape[0], :nuc.shape[1]]
    base = (stretch(nuc, t2) * 0.8).astype(np.uint8)
    files = {}
    Image.fromarray(base).save(PAGE / "img" / f"{short}_base.jpg", quality=80, optimize=True)
    files["base"] = f"img/{short}_base.jpg"
    nl = Z["new_lab"]
    edge = lambda m: m & ~ndi.binary_erosion(m)
    for name, m, rgb, fa, ea in (("igl", nl == 4, (255, 214, 10), 110, 255), ("egl", (nl == 2) & ~Z["add"], (59, 224, 255), 95, 200),
                                 ("add", Z["add"], (60, 230, 90), 200, 255), ("cut", Z["mid"], (255, 60, 200), 255, 255),
                                 ("old", Z["mid_egl"], (170, 120, 255), 230, 230),
                                 ("i2e", Z["swap_igl_to_egl"] | Z["swap_igl_above_ml"], (255, 79, 160), 70, 255), ("e2i", Z["swap_egl_to_igl"], (255, 255, 255), 50, 255),
                                 ("i2m", Z["swap_igl_to_ml"] | Z["swap_igl_over_ml"], (205, 120, 60), 90, 255)):
        Image.fromarray(overlay(m, rgb, fa, ea), "RGBA").save(PAGE / "img" / f"{short}_{name}.png", optimize=True)
        files[name] = f"img/{short}_{name}.png"
    Q, X, valid, piece = Z["Q"], Z["X"], Z["valid"], Z["piece"]
    lab_ok, fill, bridge, gap, cortex = Z["lab_ok"], Z["fill"], Z["bridge"], Z["gap"], Z["cortex"]
    cls = np.full(len(Q), 4, np.int8)
    cls[lab_ok] = 0
    cls[fill] = 1
    cls[bridge] = 2
    cls[gap & cortex & ~fill & ~bridge] = 3
    N = len(Q)
    um = lambda m: int(m.sum() * F.STEP)
    stats = dict(row, length_um=int(N * F.STEP), labelled_um=um(cls == 0), filled_um=um(cls == 1), bridged_um=um(cls == 2),
                 left_um=um(cls == 3), none_um=um(cls == 4), note_v1=NOTES.get(short, {}).get("note", ""),
                 verdict_v1=NOTES.get(short, {}).get("verdict", ""))
    line = []
    for i in range(0, N, 2):
        line += [int(round(Q[i, 0] * F.G)), int(round(Q[i, 1] * F.G)), int(cls[i]), int(piece[i])]
    cols = []
    for i in range(0, N, 25):
        L = int(valid[i].sum())
        pts = np.vstack([Q[i:i + 1], X[i, :L:5]]) * F.G
        cols.append([int(round(v)) for v in pts.ravel()])
    data = dict(short=short, group=row["group"], base=[int(base.shape[0]), int(base.shape[1])], grid=[int(gh), int(gw)],
                line=line, line_step=2, cols=cols, files=files, stats=stats)
    (PAGE / "data" / f"{short}.json").write_text(json.dumps(data, separators=(",", ":")))
    return dict(short=short, group=row["group"], data=f"data/{short}.json", stats=stats)


def findings(secs):
    st = {s["short"]: s["stats"] for s in secs}
    k = lambda x: f"{x / 1000:.0f}k"
    def rng(grp, f, fmt):
        v = [f(s["stats"]) for s in secs if s["group"] == grp]
        return f"{fmt(min(v))}–{fmt(max(v))}"
    mm = lambda x: f"{x / 1000:.1f}"
    sw = lambda x, key: x["swaps_um2"][key]
    old = {r["section"]: r for r in json.loads((R.OUT / "npz_v4" / "summary.json").read_text())["sections"]}
    ch = lambda sec, key: f"{k(old[sec][key])} → {k(st[sec][key])}"
    return [
        f"<b>Arachnoid and pia are not EGL</b> (s1_2 at 11 and 15, s2_4 at 15, s3_2 between 10 and 21). Two rules on the outer surface. "
        f"A column with no EGL near the line but labelled EGL further down, above its IGL, has the EGL under something else (a flap, arachnoid), so nothing is added on top. "
        f"And EGL is not added where p27 and NeuN are both well below the ML (p27 under 0.6 ×, NeuN under 0.7 × the section's ML level): the EGL has 1.3–1.5 × the ML's p27, "
        f"the spots you marked 0.3–0.4 ×. Added EGL: s1_2 {ch('s1_2', 'egl_added_um2')} µm², s2_4 {ch('s2_4', 'egl_added_um2')}, s3_2 {ch('s3_2', 'egl_added_um2')}.",
        f"<b>EGL called IGL in s3_2</b> (between 10 and 21): those cells have EGL-level NeuN (far below this section's IGL) but look dim next to the very bright EGL around them. "
        f"IGL → EGL now needs nuclei only brighter than the IGL / ML split, with low NeuN, within 60 µm of the line: {k(sw(st['s3_2'], 'igl_to_egl'))} µm² in s3_2 "
        f"(was {k(old['s3_2']['swaps_um2']['igl_to_egl'])}).",
        f"<b>EGL deep to the IGL in s1_5</b> (12 and 13): the cut there was a fallback cut halfway between two facing IGL bands, and the white-matter core of a folium looks the same. "
        f"Fallback cuts now need stain EGL beside them (at least 20 %) and no white matter (under 10 %). Without the cut that EGL lies deep: deep EGL with low p27 and NeuN becomes "
        f"white matter, other deep EGL takes the nearest label ({k(st['s1_5']['deep_egl_relabelled_um2'])} µm² in s1_5).",
        f"<b>EGL down the furrow in s4_A</b> (deep to 5): the fused EGL at the furrow base is labelled white matter, and the cut stopped there. The cut now continues through "
        f"white-matter labels where the stain shows EGL (bright nuclei, low NeuN), still stopping 30 µm from the IGL, and those cells become EGL. "
        f"Cut extensions in s4_A {k(old['s4_A']['cuts']['extension_um2'])} → {k(st['s4_A']['cuts']['extension_um2'])} µm². "
        f"The same happens in s2_2, along an EGL streak in white matter (your cKO corrections had EGL there; the last version turned it into white matter as deep EGL).",
        f"EGL gaps over IGL left open (orange): {rng('control', lambda x: x['left_um'], mm)} mm per control section, {rng('cKO', lambda x: x['left_um'], mm)} mm per cKO section. "
        "The mm marks have moved again; your note from the last version is shown on each section.",
    ]


def main():
    (PAGE / "img").mkdir(parents=True, exist_ok=True)
    (PAGE / "data").mkdir(parents=True, exist_ok=True)
    rows = {r["section"]: r for r in json.loads((R.OUT / "summary.json").read_text())["sections"]}
    order = ["s1_2"] + sorted(s for s, r in rows.items() if r["group"] == "control" and s != "s1_2") + \
        sorted(s for s, r in rows.items() if r["group"] == "cKO")
    secs = [build(s, rows[s]) for s in order]
    (PAGE / "index.html").write_text(TEMPLATE.read_text().replace(
        "__DATA__", json.dumps(dict(collection="egl_repair_check_v5", sections=secs, findings=findings(secs)), separators=(",", ":"))))
    tot = sum(f.stat().st_size for f in PAGE.rglob("*") if f.is_file())
    print(f"page: {len(list(PAGE.rglob('*.*')))} files, {tot / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
