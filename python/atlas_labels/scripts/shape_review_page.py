"""
Tick-box review page for the shape profile: orientation, IGL and white matter of every
IF section, IGL and white matter of the Allen P56 sections within +-0.6 mm.

    python scripts/shape_review_page.py  -> outputs/shape_review/shape_review.html

Answers are stored in the page's db (collection "reviews", one document per item id)
and read back with ArtifactData. Label maps are drawn in atlas orientation
(anterior left, dorsal up) at ~8 um/px, using the rotations in
outputs/wm_library/library.json.
"""
import base64
import html
import io
import json
import sys
from pathlib import Path

import numpy as np
import tifffile
from PIL import Image
from scipy import ndimage as ndi

sys.path.insert(0, str(Path(__file__).parent))
import wm_library as WL  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/shape_review"
DRAW_UM = 8.0
# 9-class IF codes; Allen codes 3 EGL, 4 IGL, 5 ML, 6 WM map onto the same colours
PALETTE = {0: (0, 0, 0, 0), 1: (150, 152, 160, 255), 2: (110, 175, 235, 255), 3: (0, 196, 214, 255),
           4: (232, 186, 28, 255), 5: (146, 212, 140, 255), 6: (226, 122, 50, 255), 7: (214, 64, 64, 255),
           8: (196, 84, 196, 255), 9: (80, 80, 86, 255)}


def rotate_labels(lab, rot_deg, mirror):
    th = np.deg2rad(rot_deg)
    c, s = np.cos(th), np.sin(th)
    A = (np.diag([1.0, -1.0]) if mirror else np.eye(2)) @ np.array([[c, -s], [s, c]]).T
    H, W = lab.shape
    n = int(np.ceil(np.hypot(H, W))) + 4
    o = np.array([H, W]) / 2.0 - A @ (np.array([n, n]) / 2.0)
    return ndi.affine_transform(lab, A, o, output_shape=(n, n), order=0)


def crop_labels(lab, pad=3):
    ys, xs = np.nonzero(lab)
    return lab[max(ys.min() - pad, 0):ys.max() + pad + 1, max(xs.min() - pad, 0):xs.max() + pad + 1]


def png_uri(lab):
    rgba = np.zeros(lab.shape + (4,), np.uint8)
    for k, col in PALETTE.items():
        rgba[lab == k] = col
    buf = io.BytesIO()
    Image.fromarray(rgba, "RGBA").save(buf, "PNG", optimize=True)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode(), lab.shape


def if_labels(sid, orient):
    lab = tifffile.imread(WL.FIXED / f"{sid}_fixed_labels.tif")
    px = 0.5119 if sid.startswith("2018_") else WL.PX_10X
    step = max(1, int(round(DRAW_UM / px)))
    return crop_labels(rotate_labels(lab[::step, ::step], orient["rot_deg"], orient["mirror"]))


def allen_labels(age, sec, px):
    lab = tifffile.imread(ROOT / f"outputs/allen_{age.lower()}/sections/sec{sec:03d}_labels.tif")
    z = px / DRAW_UM
    return crop_labels(ndi.zoom(lab, z, order=0))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    wl = {e["id"]: e for e in json.loads((WL.OUT / "library.json").read_text())["entries"]}
    il = {e["id"]: e for e in json.loads((ROOT / "outputs/igl_library/library.json").read_text())["entries"]}
    items = []
    for sid in WL.SLIDES_20X + WL.SLIDES_10X:
        eid = f"if_{sid}"
        o = wl[eid]["orientation"]
        uri, shp = png_uri(if_labels(sid, o))
        items.append(dict(id=eid, group="if", name=sid.replace("2018_05_22_", ""),
                          sub="20x" if sid.startswith("2018_") else "10x", img=uri, w=shp[1], h=shp[0],
                          facts=f"orientation fit margin {o['mirror_margin']:.3f}"
                                + (" (uncertain)" if o["mirror_margin"] < 0.02 else "")
                                + f" · IGL {il[eid]['igl_mm2']:.2f} mm² · WM {wl[eid]['wm_mm2']:.2f} mm²"))
    st = {r["sec"]: r["offset_um"] for r in json.loads((ROOT / "outputs/allen_p56/p56_nissl_stats.json").read_text())}
    for sec in sorted(st, key=lambda s: st[s]):
        if abs(st[sec]) > 600:
            continue
        eid = f"allen_P56_sec{sec:03d}"
        uri, shp = png_uri(allen_labels("P56", sec, 4.188))
        items.append(dict(id=eid, group="p56", name="midline" if st[sec] == 0 else f"{st[sec]:+.0f} µm",
                          sub=f"sec {sec}", img=uri, w=shp[1], h=shp[0],
                          facts=f"IGL {il[eid]['igl_mm2']:.2f} mm² · WM {wl[eid]['wm_mm2']:.2f} mm²"))
    gold = tifffile.imread(WL.SEG / "labelled" / "2018_05_22_s1_2_p27-0021_labelled.tiff").astype(np.int16)
    gl = gold.copy()
    gl[gold == 4], gl[gold == 5] = 5, 4
    og = wl["gold_s1_2"]["orientation"]
    ref1 = png_uri(crop_labels(rotate_labels(gl.astype(np.uint8)[::16, ::16], og["rot_deg"], og["mirror"])))
    ref2 = png_uri(allen_labels("P7", 218, 3.96))
    n_if = sum(i["group"] == "if" for i in items)
    n_p56 = len(items) - n_if
    data = json.dumps(items)

    page = """<title>Cerebellum Shape Review</title>
<meta name="claude:viewport-fit-sides" content="cover">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans+Condensed:wght@500;600&family=IBM+Plex+Sans:wght@400;500&display=swap">
<style>
/* layout: sticky progress bar, reference plates, then one review card per section */
:root {
  --bg: #f2f4f7; --surface: #ffffff; --fg: #18202b; --muted: #5a6575; --line: #d6dce4;
  --accent: #1d5b86; --accent-fg: #ffffff; --good: #1f7a4d; --ok: #9a6a00; --bad: #b3362f;
  --chip-bg: #eef1f5; --plate: #e9edf2;
  --display: "IBM Plex Sans Condensed", "Arial Narrow", system-ui, sans-serif;
  --text: "IBM Plex Sans", system-ui, -apple-system, sans-serif;
  --mono: "IBM Plex Mono", ui-monospace, Menlo, monospace;
}
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) {
  --bg: #11151b; --surface: #19202a; --fg: #e5e9ef; --muted: #9aa5b3; --line: #2c3542;
  --accent: #82b8e2; --accent-fg: #0d1520; --good: #6fd29c; --ok: #f0c060; --bad: #ff8b80;
  --chip-bg: #232b36; --plate: #141a22; color-scheme: dark; } }
:root[data-theme="dark"] {
  --bg: #11151b; --surface: #19202a; --fg: #e5e9ef; --muted: #9aa5b3; --line: #2c3542;
  --accent: #82b8e2; --accent-fg: #0d1520; --good: #6fd29c; --ok: #f0c060; --bad: #ff8b80;
  --chip-bg: #232b36; --plate: #141a22; color-scheme: dark; }
body { background: var(--bg); color: var(--fg); font-family: var(--text); font-size: 15px; line-height: 1.5; }
body::before { content: ""; position: fixed; inset: 0 0 auto; height: env(safe-area-inset-top, 0px); z-index: 10; background: var(--surface); }
.wrap { max-width: 1180px; margin: 0 auto; padding-block: 28px 64px;
  padding-left: max(16px, env(safe-area-inset-left, 0px)); padding-right: max(16px, env(safe-area-inset-right, 0px));
  display: flex; flex-direction: column; gap: 28px; }
h1, h2 { font-family: var(--display); font-weight: 600; margin: 0; text-wrap: balance; }
h1 { font-size: clamp(28px, 4vw, 40px); line-height: 1.1; }
h2 { font-size: 22px; }
.lede { color: var(--muted); max-width: 70ch; margin: 8px 0 0; }
.bar { position: sticky; top: env(safe-area-inset-top, 0px); z-index: 5; background: var(--surface); border: 1px solid var(--line);
  border-radius: 6px; padding: 10px 14px; display: flex; flex-wrap: wrap; gap: 6px 20px; align-items: center; font-size: 14px; }
.bar .count { font-family: var(--mono); font-variant-numeric: tabular-nums; }
.bar .status { color: var(--muted); }
.bar .status.err { color: var(--bad); }
.refs { display: grid; grid-template-columns: repeat(auto-fit, minmax(260px, 1fr)); gap: 12px; }
.ref, .card { background: var(--surface); border: 1px solid var(--line); border-radius: 6px; overflow: hidden; display: flex; flex-direction: column; min-width: 0; }
.plate { background: var(--plate); display: flex; align-items: center; justify-content: center; padding: 8px; position: relative; }
.plate img { max-width: 100%; height: auto; image-rendering: pixelated; }
.axes { position: absolute; left: 8px; top: 6px; font: 11px var(--mono); color: var(--muted); }
.cap { padding: 10px 12px; display: flex; flex-direction: column; gap: 8px; border-top: 1px solid var(--line); }
.cap-top { display: flex; justify-content: space-between; gap: 8px; align-items: baseline; }
.name { font-weight: 500; overflow-wrap: anywhere; }
.sub { font: 12px var(--mono); color: var(--muted); white-space: nowrap; }
.facts { font: 12px var(--mono); color: var(--muted); }
.legend { display: flex; flex-wrap: wrap; gap: 6px 14px; font-size: 13px; color: var(--muted); }
.legend i { display: inline-block; width: 12px; height: 12px; border-radius: 2px; margin-right: 6px; vertical-align: -1px; }
.grid { display: grid; grid-template-columns: repeat(auto-fill, minmax(300px, 1fr)); gap: 12px; }
.grid.small { grid-template-columns: repeat(auto-fill, minmax(220px, 1fr)); }
.q { display: grid; grid-template-columns: 7.5em 1fr; gap: 6px; align-items: center; }
.q .lab { font-size: 13px; color: var(--muted); }
.seg { display: flex; flex-wrap: wrap; gap: 4px; }
.seg button { font: 13px var(--text); padding: 5px 10px; border-radius: 4px; border: 1px solid var(--line); background: var(--chip-bg);
  color: var(--fg); cursor: pointer; min-height: 32px; }
.seg button:hover { border-color: var(--accent); }
.seg button:focus-visible, .note:focus-visible { outline: 2px solid var(--accent); outline-offset: 1px; }
.seg button[aria-pressed="true"] { background: var(--accent); border-color: var(--accent); color: var(--accent-fg); }
.seg button[aria-pressed="true"].v-good { background: var(--good); border-color: var(--good); }
.seg button[aria-pressed="true"].v-ok { background: var(--ok); border-color: var(--ok); }
.seg button[aria-pressed="true"].v-bad { background: var(--bad); border-color: var(--bad); }
.seg button:disabled { opacity: 0.5; cursor: default; }
.note { width: 100%; box-sizing: border-box; font: 13px var(--text); padding: 6px 8px; border: 1px solid var(--line); border-radius: 4px; background: var(--bg); color: var(--fg); }
.card.done { border-color: var(--good); }
.tag { font: 11px var(--mono); color: var(--good); }
@media (prefers-reduced-motion: no-preference) { .seg button { transition: background-color .12s, border-color .12s; } }
</style>
<div class="wrap">
<header>
  <h1>Cerebellum Shape Review</h1>
  <p class="lede">Tick what you see. Every map is drawn the way the shape profile will use it: anterior to the left, dorsal up, about 8 µm per pixel. Your answers save as you click and stay for next time. Untouched questions count as not reviewed.</p>
</header>
<div class="bar" role="status"><span class="count" id="count">0 of __N__ reviewed</span><span class="status" id="status">Connecting to saved answers…</span></div>
<section>
  <div class="legend"><span><i style="background:rgb(0,196,214)"></i>oEGL</span><span><i style="background:rgb(110,175,235)"></i>iEGL</span><span><i style="background:rgb(146,212,140)"></i>ML</span><span><i style="background:rgb(214,64,64)"></i>PCL</span><span><i style="background:rgb(232,186,28)"></i>IGL</span><span><i style="background:rgb(226,122,50)"></i>white matter</span><span><i style="background:rgb(196,84,196)"></i>DCN</span><span><i style="background:rgb(150,152,160)"></i>unassigned</span><span><i style="background:rgb(80,80,86)"></i>tear / bleach</span></div>
</section>
<section>
  <h2>Orientation references</h2>
  <div class="refs">
    <div class="ref"><div class="plate"><span class="axes">← anterior · dorsal ↑</span><img src="__REF1__" width="__R1W__" height="__R1H__" alt="s1_2 gold standard in atlas orientation"></div><div class="cap"><span class="name">s1_2 gold standard</span><span class="facts">your hand-painted labels, same rotation as its automatic labels</span></div></div>
    <div class="ref"><div class="plate"><span class="axes">← anterior · dorsal ↑</span><img src="__REF2__" width="__R2W__" height="__R2H__" alt="Allen P7 midline Nissl labels"></div><div class="cap"><span class="name">Allen P7 midline (section 218)</span><span class="facts">template the orientations were matched to</span></div></div>
  </div>
</section>
<section><h2>Your sections (__NIF__)</h2><div class="grid" id="if"></div></section>
<section><h2>Allen P56 within ±0.6 mm (__NP56__)</h2><p class="lede">Adult reference. Judge whether the IGL and the white matter look complete enough to use as shape.</p><div class="grid small" id="p56"></div></section>
</div>
<script>
const ITEMS = __DATA__;
const QUESTIONS = {
  if: [["orientation", "Orientation", [["correct", "Correct", "good"], ["flipped", "Flipped front–back", "ok"], ["wrong", "Wrong", "bad"]]],
       ["igl", "IGL", [["good", "Good", "good"], ["usable", "Usable", "ok"], ["bad", "Bad", "bad"]]],
       ["wm", "White matter", [["good", "Good", "good"], ["usable", "Usable", "ok"], ["bad", "Bad", "bad"]]]],
  p56: [["igl", "IGL", [["good", "Good", "good"], ["bad", "Bad", "bad"]]],
        ["wm", "White matter", [["good", "Good", "good"], ["bad", "Bad", "bad"]]]]
};
const reviews = {};
const chains = {};
let db = null;
const $ = (id) => document.getElementById(id);
const statusEl = $("status");
function setStatus(t, err) { statusEl.textContent = t; statusEl.classList.toggle("err", !!err); }

function card(it) {
  const el = document.createElement("article");
  el.className = "card"; el.id = "c-" + it.id;
  const plate = document.createElement("div"); plate.className = "plate";
  const img = document.createElement("img"); img.src = it.img; img.width = it.w; img.height = it.h; img.alt = "Label map of " + it.name; img.loading = "lazy";
  plate.appendChild(img);
  const cap = document.createElement("div"); cap.className = "cap";
  const top = document.createElement("div"); top.className = "cap-top";
  const nm = document.createElement("span"); nm.className = "name"; nm.textContent = it.name;
  const sb = document.createElement("span"); sb.className = "sub"; sb.textContent = it.sub;
  top.append(nm, sb); cap.appendChild(top);
  const f = document.createElement("span"); f.className = "facts"; f.textContent = it.facts; cap.appendChild(f);
  for (const [field, label, opts] of QUESTIONS[it.group]) {
    const row = document.createElement("div"); row.className = "q";
    const l = document.createElement("span"); l.className = "lab"; l.textContent = label;
    const seg = document.createElement("div"); seg.className = "seg"; seg.setAttribute("role", "group"); seg.setAttribute("aria-label", label + " for " + it.name);
    for (const [val, text, tone] of opts) {
      const b = document.createElement("button"); b.type = "button"; b.textContent = text; b.className = "v-" + tone;
      b.dataset.id = it.id; b.dataset.field = field; b.dataset.value = val; b.setAttribute("aria-pressed", "false"); b.disabled = true;
      b.addEventListener("click", () => choose(it.id, field, val));
      seg.appendChild(b);
    }
    row.append(l, seg); cap.appendChild(row);
  }
  if (it.group === "if") {
    const n = document.createElement("input"); n.type = "text"; n.className = "note"; n.id = "note-" + it.id; n.placeholder = "Note (optional)";
    n.disabled = true; n.setAttribute("aria-label", "Note for " + it.name);
    n.addEventListener("change", () => choose(it.id, "note", n.value.trim()));
    cap.appendChild(n);
  }
  el.append(plate, cap);
  return el;
}

function isDone(it) {
  const r = reviews[it.id] || {};
  return QUESTIONS[it.group].every(([field]) => r[field]);
}

function paint() {
  let done = 0;
  for (const it of ITEMS) {
    const r = reviews[it.id] || {};
    const el = $("c-" + it.id);
    el.querySelectorAll("button").forEach((b) => b.setAttribute("aria-pressed", String(r[b.dataset.field] === b.dataset.value)));
    const n = $("note-" + it.id);
    if (n && document.activeElement !== n) n.value = r.note || "";
    const d = isDone(it); el.classList.toggle("done", d); if (d) done++;
  }
  $("count").textContent = done + " of " + ITEMS.length + " reviewed";
}

function choose(id, field, value) {
  if (!db) return;
  const cur = Object.assign({}, reviews[id] || {});
  if (field !== "note" && cur[field] === value) delete cur[field]; else cur[field] = value;
  if (field === "note" && !value) delete cur.note;
  cur.updatedAt = new Date().toISOString();
  reviews[id] = cur;
  paint();
  setStatus("Saving…");
  chains[id] = (chains[id] || Promise.resolve()).then(() => db.doc("reviews/" + id).set(cur))
    .then(() => setStatus("Saved"))
    .catch((e) => setStatus(e && e.code === "invalid_argument" ? "You can view this page but not save answers." : "Not saved: " + (e && e.message ? e.message : "connection problem") + ". Click again to retry.", true));
}

for (const it of ITEMS) $(it.group).appendChild(card(it));
paint();

(async () => {
  const use = window.claude && window.claude.use;
  db = use ? await window.claude.use("db") : null;
  if (!db) { setStatus("Answers can't be saved in this view. Open the page on claude.ai to review.", true); return; }
  document.querySelectorAll(".seg button, .note").forEach((b) => { b.disabled = false; });
  const timer = setTimeout(() => setStatus("Can't reach saved answers yet. Clicks are kept on this page and saved when the connection returns.", true), 10000);
  db.collection("reviews").onSnapshot((snap) => {
    for (const d of snap.docs) if (!chains[d.id] || !snap.metadata.hasPendingWrites) reviews[d.id] = d.data();
    if (!snap.metadata.fromCache) { clearTimeout(timer); if (statusEl.textContent.startsWith("Connecting") || statusEl.textContent.startsWith("Can't reach")) setStatus("Answers load and save automatically"); }
    paint();
  }, (e) => { clearTimeout(timer); setStatus("Saved answers can't be loaded right now (" + e.code + ").", true); });
})();
</script>
"""
    page = (page.replace("__DATA__", data).replace("__N__", str(len(items)))
            .replace("__NIF__", str(n_if)).replace("__NP56__", str(n_p56))
            .replace("__REF1__", ref1[0]).replace("__R1W__", str(ref1[1][1])).replace("__R1H__", str(ref1[1][0]))
            .replace("__REF2__", ref2[0]).replace("__R2W__", str(ref2[1][1])).replace("__R2H__", str(ref2[1][0])))
    p = OUT / "shape_review.html"
    p.write_text(page)
    print(p, f"{p.stat().st_size / 1e6:.2f} MB", f"{n_if} IF + {n_p56} P56 items")


if __name__ == "__main__":
    main()
