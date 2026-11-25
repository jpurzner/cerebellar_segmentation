"""
HTML gallery of the white-matter library (run wm_library.py first).

    python scripts/wm_library_page.py   -> outputs/wm_library/wm_library.html
"""
import base64
import html
import json
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / "outputs/wm_library"


def b64(p):
    return "data:image/png;base64," + base64.b64encode(p.read_bytes()).decode()


def pieces(wm, grid):
    cc, n = ndi.label(wm, structure=np.ones((3, 3), bool))
    a = np.bincount(cc.ravel())[1:] * grid * grid
    return int((a >= 1000).sum())


def main():
    lib = json.loads((D / "library.json").read_text())
    z = np.load(D / "library.npz")
    grid, thumb_um = lib["grid_um"], lib["thumb_um"]
    H, W = lib["canvas_px"]
    bar_pct = 100 * (1000 / thumb_um) / W          # 1 mm as % of thumbnail width
    E = lib["entries"]
    cur = json.loads((D / "curation.json").read_text())
    RANK = {"really good": 0, "pretty good": 1, "usable, in several pieces": 2}
    for e in E:
        e["pieces"] = pieces(z[f"wm_{e['id']}"], grid)
        e["wm_pct"] = round(100 * e["wm_mm2"] / max(e["tissue_mm2"], 1e-6))
        e["rating"] = cur["ratings"].get(e["id"])
        if e["source"].startswith("Allen"):
            e["shape"] = e["age"] in cur["allen_shape_ages"] and abs(e["offset_um"]) <= cur["allen_shape_max_offset_um"]
        else:
            e["shape"] = e["rating"] is not None
    shape_set = [e for e in E if e["shape"]]

    user = sorted([e for e in E if e["source"].startswith("IF")],
                  key=lambda e: (RANK.get(e["rating"], 3), e["source"] != "IF gold", e["name"]))
    allen = {a: sorted([e for e in E if e["source"] == f"Allen {a}"], key=lambda e: e["offset_um"])
             for a in ("P4", "P7", "P56")}

    def pick(rows, step, lim):
        out = []
        for t in np.arange(-lim, lim + 1, step):
            r = min(rows, key=lambda e: abs(e["offset_um"] - t))
            if abs(r["offset_um"] - t) <= step / 2 and r not in out:
                out.append(r)
        return out

    gallery = {"P7": pick(allen["P7"], 200, 2000), "P4": pick(allen["P4"], 200, 2000),
               "P56": pick(allen["P56"], 150, 1200)}

    def card(e, label, sub):
        chips = ""
        if e.get("rating"):
            chips += f'<span class="chip good">your grade: {html.escape(e["rating"])}</span>'
        if e["source"].startswith("Allen"):
            chips += '<span class="chip good">used for shape</span>' if e["shape"] else '<span class="chip off">not used for shape</span>'
        chips += "".join(f'<span class="chip warn">{html.escape(n)}</span>' for n in e["notes"])
        return f"""<figure class="card">
  <div class="plate"><img src="{b64(D / 'thumbs' / (e['id'] + '.png'))}" alt="White matter of {html.escape(label)}" width="{W}" height="{H}" loading="lazy"><span class="bar" style="width:{bar_pct:.1f}%"></span><span class="barlab">1 mm</span></div>
  <figcaption><div class="cap-top"><span class="name">{html.escape(label)}</span><span class="sub">{html.escape(sub)}</span></div>
  <div class="stats"><span>{e['wm_mm2']:.2f}<small> mm²</small></span><span>{e['wm_pct']}<small>% of tissue</small></span><span>{e['pieces']}<small> {'piece' if e['pieces'] == 1 else 'pieces'}</small></span></div>{('<div class="chips">' + chips + '</div>') if chips else ''}</figcaption>
</figure>"""

    def off(e):
        return "midline" if e["offset_um"] == 0 else f"{e['offset_um']:+.0f} µm"

    user_cards = "\n".join(card(e, e["name"], e["source"].replace("IF ", "")) for e in user)
    allen_cards = {a: "\n".join(card(e, off(e), e["name"]) for e in gallery[a]) for a in gallery}

    # chart: WM area vs medial-lateral offset
    cw, ch, ml, mr, mt, mb = 760, 340, 56, 20, 16, 46
    xmin, xmax, ymin, ymax = -2000, 2000, 0.0, 1.2
    X = lambda v: ml + (v - xmin) / (xmax - xmin) * (cw - ml - mr)
    Y = lambda v: mt + (1 - (v - ymin) / (ymax - ymin)) * (ch - mt - mb)
    u20 = [e["wm_mm2"] for e in user if e["source"] == "IF 20x"]
    svg = [f'<svg viewBox="0 0 {cw} {ch}" role="img" aria-label="White-matter area against medial-lateral position">']
    svg.append(f'<rect x="{ml}" y="{Y(max(u20)):.1f}" width="{cw - ml - mr}" height="{Y(min(u20)) - Y(max(u20)):.1f}" class="band"/>')
    svg.append(f'<text x="{cw - mr - 6}" y="{Y(max(u20)) - 6:.1f}" text-anchor="end" class="bandlab">your 20x sections: {min(u20):.2f}–{max(u20):.2f} mm²</text>')
    for v in np.arange(0, 1.21, 0.2):
        svg.append(f'<line x1="{ml}" x2="{cw - mr}" y1="{Y(v):.1f}" y2="{Y(v):.1f}" class="grid"/>'
                   f'<text x="{ml - 8}" y="{Y(v) + 4:.1f}" text-anchor="end" class="tick">{v:.1f}</text>')
    for v in range(-2000, 2001, 500):
        svg.append(f'<text x="{X(v):.1f}" y="{ch - mb + 18}" text-anchor="middle" class="tick">{v:+d}</text>'.replace("+0<", "0<"))
    svg.append(f'<text x="{(ml + cw - mr) / 2}" y="{ch - 8}" text-anchor="middle" class="axis">distance from midline (µm)</text>')
    svg.append(f'<text x="14" y="{(mt + ch - mb) / 2}" text-anchor="middle" class="axis" transform="rotate(-90 14 {(mt + ch - mb) / 2})">white matter (mm²)</text>')
    for a, cls in (("P4", "p4"), ("P7", "p7"), ("P56", "p56")):
        rows = [e for e in allen[a] if xmin <= e["offset_um"] <= xmax]
        xs = np.array([e["offset_um"] for e in rows])
        ys = np.array([e["wm_mm2"] for e in rows])
        k = 5
        ys_s = np.convolve(ys, np.ones(k) / k, mode="same") if len(ys) > k else ys
        pts = " ".join(f"{X(x):.1f},{Y(y):.1f}" for x, y in zip(xs, ys_s))
        svg.append(f'<polyline points="{pts}" class="line {cls}"/>')
        if a == "P56":
            for s in (-600, 600):
                svg.append(f'<line x1="{X(s):.1f}" x2="{X(s):.1f}" y1="{mt}" y2="{ch - mb}" class="limit"/>')
            svg.append(f'<text x="{X(600) + 6:.1f}" y="{mt + 12}" class="tick">P56 WM unreliable beyond ±600</text>')
    svg.append("</svg>")
    chart = "\n".join(svg)

    n_allen = {a: len(v) for a, v in allen.items()}
    page = f"""<title>Cerebellar White Matter Library</title>
<meta name="claude:viewport-fit-sides" content="cover">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans+Condensed:wght@500;600&family=IBM+Plex+Sans:wght@400;500&display=swap">
<style>
/* layout: specimen plates on one scale, grouped by source; summary and chart first */
:root {{
  --bg: #f2f4f7; --surface: #ffffff; --fg: #18202b; --muted: #5a6575; --line: #d6dce4;
  --accent: #1d5b86; --wm: #e27a32; --warn-fg: #8a4d00; --warn-bg: #fbe9d4;
  --good-fg: #1f6b45; --good-bg: #dcf0e4; --off-fg: #5a6575; --off-bg: #e8ecf1;
  --p4: #8fa3b8; --p7: #1d5b86; --p56: #6b4f9a; --band: rgba(226,122,50,0.16);
  --display: "IBM Plex Sans Condensed", "Arial Narrow", system-ui, sans-serif;
  --text: "IBM Plex Sans", system-ui, -apple-system, sans-serif;
  --mono: "IBM Plex Mono", ui-monospace, Menlo, monospace;
}}
@media (prefers-color-scheme: dark) {{ :root:not([data-theme="light"]) {{
  --bg: #11151b; --surface: #19202a; --fg: #e5e9ef; --muted: #9aa5b3; --line: #2c3542;
  --accent: #82b8e2; --warn-fg: #f3bb78; --warn-bg: #3a2a17; --p4: #6d8299; --p7: #82b8e2; --p56: #b39ddb;
  --good-fg: #8fd9b0; --good-bg: #183528; --off-fg: #9aa5b3; --off-bg: #232b36;
  --band: rgba(226,122,50,0.22); color-scheme: dark; }} }}
:root[data-theme="dark"] {{
  --bg: #11151b; --surface: #19202a; --fg: #e5e9ef; --muted: #9aa5b3; --line: #2c3542;
  --accent: #82b8e2; --warn-fg: #f3bb78; --warn-bg: #3a2a17; --p4: #6d8299; --p7: #82b8e2; --p56: #b39ddb;
  --good-fg: #8fd9b0; --good-bg: #183528; --off-fg: #9aa5b3; --off-bg: #232b36;
  --band: rgba(226,122,50,0.22); color-scheme: dark; }}
body {{ background: var(--bg); color: var(--fg); font-family: var(--text); font-size: 15px; line-height: 1.55; }}
.wrap {{ max-width: 1180px; margin: 0 auto; padding-block: 40px 64px;
  padding-left: max(16px, env(safe-area-inset-left, 0px)); padding-right: max(16px, env(safe-area-inset-right, 0px));
  display: flex; flex-direction: column; gap: 40px; }}
h1, h2 {{ font-family: var(--display); font-weight: 600; text-wrap: balance; margin: 0; letter-spacing: 0.01em; }}
h1 {{ font-size: clamp(30px, 4.5vw, 44px); line-height: 1.1; }}
h2 {{ font-size: 24px; }}
.lede {{ max-width: 68ch; color: var(--muted); margin: 10px 0 0; }}
.eyebrow {{ font-family: var(--mono); font-size: 12px; letter-spacing: 0.08em; text-transform: uppercase; color: var(--accent); }}
.legend {{ display: flex; flex-wrap: wrap; gap: 18px; margin-top: 16px; font-size: 13px; color: var(--muted); }}
.legend span {{ display: inline-flex; align-items: center; gap: 8px; }}
.sw {{ width: 18px; height: 12px; border-radius: 2px; display: inline-block; }}
.sw.wm {{ background: var(--wm); }} .sw.sk {{ background: #5a2606; }} .sw.ti {{ background: rgba(128,132,140,0.35); border: 1px solid rgba(128,132,140,0.8); }}
section {{ display: flex; flex-direction: column; gap: 14px; }}
.sec-head {{ display: flex; flex-wrap: wrap; align-items: baseline; justify-content: space-between; gap: 8px 24px; }}
.sec-head p {{ margin: 0; color: var(--muted); max-width: 70ch; font-size: 14px; }}
.tablewrap {{ overflow-x: auto; border: 1px solid var(--line); border-radius: 6px; background: var(--surface); }}
table {{ border-collapse: collapse; width: 100%; font-size: 14px; }}
th, td {{ text-align: left; padding: 9px 14px; border-bottom: 1px solid var(--line); vertical-align: top; }}
th {{ font-family: var(--mono); font-weight: 500; font-size: 12px; letter-spacing: 0.05em; text-transform: uppercase; color: var(--muted); }}
tr:last-child td {{ border-bottom: 0; }}
td.num {{ font-family: var(--mono); font-variant-numeric: tabular-nums; white-space: nowrap; }}
.chart {{ background: var(--surface); border: 1px solid var(--line); border-radius: 6px; padding: 14px; overflow-x: auto; }}
.chart svg {{ width: 100%; min-width: 520px; height: auto; display: block; }}
.chart .grid {{ stroke: var(--line); stroke-width: 1; }}
.chart .tick {{ fill: var(--muted); font: 11px var(--mono); }}
.chart .axis {{ fill: var(--muted); font: 12px var(--text); }}
.chart .band {{ fill: var(--band); }}
.chart .bandlab {{ fill: var(--warn-fg); font: 12px var(--text); }}
.chart .line {{ fill: none; stroke-width: 2.2; }}
.chart .p4 {{ stroke: var(--p4); }} .chart .p7 {{ stroke: var(--p7); }} .chart .p56 {{ stroke: var(--p56); stroke-dasharray: 5 4; }}
.chart .limit {{ stroke: var(--p56); stroke-width: 1; stroke-dasharray: 2 3; opacity: 0.7; }}
.keys {{ display: flex; flex-wrap: wrap; gap: 16px; font-size: 13px; color: var(--muted); margin-top: 8px; }}
.keys i {{ display: inline-block; width: 22px; height: 0; border-top: 3px solid; vertical-align: middle; margin-right: 6px; }}
.grid-cards {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(176px, 1fr)); gap: 12px; }}
.card {{ margin: 0; background: var(--surface); border: 1px solid var(--line); border-radius: 6px; overflow: hidden; display: flex; flex-direction: column; }}
.plate {{ position: relative; aspect-ratio: {W} / {H}; max-width: 100%; }}
.plate img {{ width: 100%; height: 100%; display: block; }}
.bar {{ position: absolute; left: 8px; bottom: 9px; height: 3px; background: var(--muted); border-radius: 1px; }}
.barlab {{ position: absolute; left: 8px; bottom: 14px; font: 10px var(--mono); color: var(--muted); }}
figcaption {{ padding: 8px 10px 10px; border-top: 1px solid var(--line); display: flex; flex-direction: column; gap: 5px; min-width: 0; }}
.cap-top {{ display: flex; justify-content: space-between; gap: 8px; align-items: baseline; min-width: 0; }}
.name {{ font-weight: 500; font-size: 13px; overflow-wrap: anywhere; }}
.sub {{ font: 11px var(--mono); color: var(--muted); white-space: nowrap; }}
.stats {{ display: flex; flex-wrap: wrap; gap: 4px 10px; font-family: var(--mono); font-variant-numeric: tabular-nums; font-size: 13px; }}
.stats small {{ font-size: 11px; color: var(--muted); }}
.chips {{ display: flex; flex-wrap: wrap; gap: 4px; }}
.chip {{ font-size: 11px; padding: 1px 7px; border-radius: 999px; }}
.chip.warn {{ color: var(--warn-fg); background: var(--warn-bg); }}
.chip.good {{ color: var(--good-fg); background: var(--good-bg); }}
.chip.off {{ color: var(--off-fg); background: var(--off-bg); }}
.shapeset {{ background: var(--surface); border: 1px solid var(--line); border-radius: 6px; padding: 14px 16px; display: flex; flex-direction: column; gap: 6px; max-width: 90ch; }}
.shapeset strong {{ font-family: var(--display); font-weight: 600; font-size: 17px; }}
.note {{ background: var(--surface); border: 1px dashed var(--line); border-radius: 6px; padding: 14px 16px; color: var(--muted); font-size: 14px; max-width: 80ch; }}
nav {{ display: flex; flex-wrap: wrap; gap: 8px 18px; font-size: 14px; }}
nav a {{ color: var(--accent); text-decoration: none; border-bottom: 1px solid transparent; }}
nav a:hover, nav a:focus-visible {{ border-bottom-color: var(--accent); outline: none; }}
</style>
<div class="wrap">
<header>
  <div class="eyebrow">Mouse cerebellum · sagittal · deep white matter</div>
  <h1>Cerebellar White Matter Library</h1>
  <p class="lede">Every white-matter shape we have: your P7 immunofluorescence sections (the automatic labels after the anatomy fixes, plus the hand-painted s1_2) and the Allen reference Nissl series at P4, P7 and P56. All plates share one scale (15 µm per pixel) and one orientation: anterior to the left, dorsal up. {len(E)} shapes are in the library file; this page shows all of yours and a sample of the Allen series every 150–200 µm across both hemispheres.</p>
  <div class="legend"><span><i class="sw wm"></i>white matter (DWL + DCN)</span><span><i class="sw sk"></i>white-matter centre line</span><span><i class="sw ti"></i>rest of the cerebellum</span></div>
  <nav style="margin-top:16px" aria-label="Sections"><a href="#sources">Sources</a><a href="#area">Area by position</a><a href="#yours">Your sections</a><a href="#p7">Allen P7</a><a href="#p4">Allen P4</a><a href="#p56">Allen P56</a></nav>
</header>

<div class="shapeset" id="shape-set"><strong>Shape set: {sum(1 for e in shape_set if not e['source'].startswith('Allen'))} of your sections + {sum(1 for e in shape_set if e['source'].startswith('Allen'))} Allen P56 sections</strong>
  <span>Graded by you: s5_A and s1_2 really good; s4_C and s5_G pretty good; s4_F usable but in several pieces. From the Allen series only the adult (P56) white matter is usable for the overall shape, and only within ±0.6 mm of the midline, where it was checked against the atlas. Allen P4 and P7 white matter is shown for reference but not used.</span></div>

<section id="sources">
  <div class="sec-head"><h2>Sources</h2></div>
  <div class="tablewrap"><table>
    <thead><tr><th>Source</th><th>Shapes</th><th>Resolution</th><th>White matter from</th><th>Reliability</th></tr></thead>
    <tbody>
      <tr><td>Your sections, 20x</td><td class="num">8</td><td class="num">0.51 µm/px</td><td>MATLAB labels after anatomy fixes R1–R4</td><td>s1_2 against gold: DWL IoU 0.87. Rotated into atlas orientation by template matching; 6 of 20 orientations uncertain (mirror image scored within 0.02).</td></tr>
      <tr><td>Your sections, 10x</td><td class="num">11</td><td class="num">0.65 µm/px</td><td>same (the MATLAB 10x run assumed 1.02 µm/px)</td><td>No gold standard at 10x. Controls from the EZH2 experiment.</td></tr>
      <tr><td>Gold standard s1_2</td><td class="num">1</td><td class="num">0.51 µm/px</td><td>hand-painted DWL + DCN</td><td>Reference.</td></tr>
      <tr><td>Allen P7 Nissl</td><td class="num">{n_allen['P7']}</td><td class="num">3.96 µm/px</td><td>palest dense Nissl class inside the cerebellum</td><td>Patchy (many small holes); not used for shape.</td></tr>
      <tr><td>Allen P4 Nissl</td><td class="num">{n_allen['P4']}</td><td class="num">3.96 µm/px</td><td>same method</td><td>Patchy, as P7; not used for shape.</td></tr>
      <tr><td>Allen P56 Nissl</td><td class="num">{n_allen['P56']}</td><td class="num">4.19 µm/px</td><td>pale tissue enclosed by the IGL</td><td>Checked against 13 atlas plates: Dice 0.82–0.88 within ±0.6 mm of the midline, under 0.3 beyond. Used for shape within ±0.6 mm.</td></tr>
      <tr><td>Allen P14 Nissl</td><td class="num">0</td><td class="num">4.19 µm/px</td><td>not separable</td><td>White matter is as dark as the IGL in P14 Nissl, so it is labelled with the IGL.</td></tr>
    </tbody></table></div>
</section>

<section id="area">
  <div class="sec-head"><h2>Area by position</h2><p>White-matter area of each Allen section against its distance from the midline (5-section running mean). The shaded band is the range of your 20x sections; the Allen P7 curve reaches it only between about −250 and +500 µm. Across all positions Allen P7 white matter stays at 18–22 % of the cerebellum, as in your 20x sections (20–25 %). Your 10x sections range from {min(e['wm_mm2'] for e in user if e['source'] == 'IF 10x'):.2f} to {max(e['wm_mm2'] for e in user if e['source'] == 'IF 10x'):.2f} mm².</p></div>
  <div class="chart">{chart}
    <div class="keys"><span><i style="border-color:var(--p4)"></i>P4</span><span><i style="border-color:var(--p7)"></i>P7</span><span><i style="border-color:var(--p56);border-top-style:dashed"></i>P56</span></div>
  </div>
</section>

<section id="yours">
  <div class="sec-head"><h2>Your P7 sections</h2><p>Eight 20x and eleven 10x slides plus the gold standard, your graded sections first. "Orientation uncertain" means the best rotation and its mirror image fit the template almost equally well, so this plate may be flipped front to back. 10x sizes use 0.65 µm per pixel, the camera's true scale (the MATLAB 10x port assumed 1.02). At that scale the 10x cerebella cover 2.1–3.1 mm², and their white matter is a smaller share (7–18 %) than at 20x (20–25 %), possibly because the 10x pipeline's size filters ran at the wrong scale. All 10x slides are controls.</p></div>
  <div class="grid-cards">
{user_cards}
  </div>
</section>

<section id="p7">
  <div class="sec-head"><h2>Allen P7</h2><p>Every 200 µm from −2 mm to +2 mm. The template section that best matched s1_2 sits at −160 µm.</p></div>
  <div class="grid-cards">
{allen_cards['P7']}
  </div>
</section>

<section id="p4">
  <div class="sec-head"><h2>Allen P4</h2><p>Every 200 µm from −2 mm to +2 mm.</p></div>
  <div class="grid-cards">
{allen_cards['P4']}
  </div>
</section>

<section id="p56">
  <div class="sec-head"><h2>Allen P56</h2><p>Every 150 µm from −1.2 mm to +1.2 mm. Beyond ±0.6 mm the white matter opens toward the peduncles and this method misses most of it.</p></div>
  <div class="grid-cards">
{allen_cards['P56']}
  </div>
</section>

<p class="note">Library files: outputs/wm_library/library.npz holds every white-matter and tissue mask at 4 µm per pixel, and outputs/wm_library/library.json holds their metadata (source, position, area, orientation). Built by scripts/wm_library.py; this page by scripts/wm_library_page.py.</p>
</div>
"""
    p = D / "wm_library.html"
    p.write_text(page)
    print(p, f"{p.stat().st_size / 1e6:.2f} MB", {a: len(v) for a, v in gallery.items()}, "user", len(user))


if __name__ == "__main__":
    main()
