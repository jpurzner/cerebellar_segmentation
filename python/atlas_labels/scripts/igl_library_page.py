"""
HTML gallery of the IGL library (run igl_library.py first).

    python scripts/igl_library_page.py   -> outputs/igl_library/igl_library.html

Every number in the page text is computed here from outputs/igl_library/library.json.
"""
import base64
import html
import json
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / "outputs/igl_library"
AGES = ("P4", "P7", "P14", "P56")
GALLERY = {"P4": (200, 2000), "P7": (200, 2000), "P14": (200, 2000), "P56": (150, 1200)}   # step, +-limit (um)
SMOOTH = 5            # running mean over this many Allen sections in the chart
SCALE_MM = 1          # scale bar on the plates
NOTE_CLS = {"orientation uncertain": "warn", "IGL + white matter merged": "warn",
            "hand-painted": "good", "IGL not validated": "off"}


def b64(p):
    return "data:image/png;base64," + base64.b64encode(p.read_bytes()).decode()


def um(v):
    """Signed micrometre value with a real minus sign."""
    return "midline" if v == 0 else f"{v:+.0f}".replace("-", "−")


def rng(vals, fmt="{:.2f}"):
    return f"{fmt.format(min(vals))}–{fmt.format(max(vals))}"


def med(vals, fmt="{:.0f}"):
    return fmt.format(float(np.median(vals)))


def main():
    lib = json.loads((D / "library.json").read_text())
    z = np.load(D / "library.npz")
    grid, thumb_um, defs, chk = lib["grid_um"], lib["thumb_um"], lib["definitions"], lib["checks"]
    H, W = lib["canvas_px"]
    bar_pct = 100 * (SCALE_MM * 1000 / thumb_um) / W
    E = lib["entries"]

    gold = [e for e in E if e["source"] == "IF gold"]
    if20 = sorted([e for e in E if e["source"] == "IF 20x"], key=lambda e: e["name"])
    if10 = sorted([e for e in E if e["source"] == "IF 10x"], key=lambda e: e["name"])
    user = gold + if20 + if10
    allen = {a: sorted([e for e in E if e["source"] == f"Allen {a}"], key=lambda e: e["offset_um"]) for a in AGES}

    def pick(rows, step, lim):
        out = []
        for t in np.arange(-lim, lim + 1, step):
            r = min(rows, key=lambda e: abs(e["offset_um"] - t))
            if abs(r["offset_um"] - t) <= step / 2 and r not in out:
                out.append(r)
        return out

    gallery = {a: pick(allen[a], *GALLERY[a]) for a in AGES}

    def card(e, label, sub):
        chips = "".join(f'<span class="chip {NOTE_CLS.get(n, "warn")}">{html.escape(n)}</span>' for n in e["notes"])
        p = e["pieces"]
        return f"""<figure class="card">
  <div class="plate"><img src="{b64(D / 'thumbs' / (e['id'] + '.png'))}" alt="Internal granule layer of {html.escape(label)}" width="{W}" height="{H}" loading="lazy"><span class="bar" style="width:{bar_pct:.1f}%"></span><span class="barlab">{SCALE_MM} mm</span></div>
  <figcaption><div class="cap-top"><span class="name">{html.escape(label)}</span><span class="sub">{html.escape(sub)}</span></div>
  <div class="stats"><span>{e['igl_mm2']:.2f}<small> mm²</small></span><span>{e['igl_pct']:.0f}<small>% of tissue</small></span><span>{p}<small> {'piece' if p == 1 else 'pieces'}</small></span><span>{e['thickness_um']:.0f}<small> µm thick</small></span></div>{('<div class="chips">' + chips + '</div>') if chips else ''}</figcaption>
</figure>"""

    user_cards = "\n".join(card(e, e["name"], e["source"].replace("IF ", "")) for e in user)
    allen_cards = {a: "\n".join(card(e, um(e["offset_um"]), e["name"]) for e in gallery[a]) for a in AGES}

    # ---- chart: IGL area vs medial-lateral offset -------------------------------------------------
    lines = {}
    for a in AGES:
        xs = np.array([e["offset_um"] for e in allen[a]])
        ys = np.array([e["igl_mm2"] for e in allen[a]])
        lines[a] = (xs, ndi.uniform_filter1d(ys, SMOOTH, mode="nearest") if len(ys) > SMOOTH else ys)
    a20 = [e["igl_mm2"] for e in if20]
    a10 = [e["igl_mm2"] for e in if10]
    lim = 500 * np.ceil(max(abs(x) for xs, _ in lines.values() for x in xs) / 500)
    xmin, xmax = -lim, lim
    ytop = max(max(ys.max() for _, ys in lines.values()), max(a20))
    ystep = 0.5 if ytop <= 4 else 1.0
    ymin, ymax = 0.0, ystep * np.ceil(ytop / ystep)
    cw, ch, ml, mr, mt, mb = 760, 340, 56, 20, 16, 46
    X = lambda v: ml + (v - xmin) / (xmax - xmin) * (cw - ml - mr)     # noqa: E731
    Y = lambda v: mt + (1 - (v - ymin) / (ymax - ymin)) * (ch - mt - mb)  # noqa: E731
    svg = [f'<svg viewBox="0 0 {cw} {ch}" role="img" aria-label="IGL area of each Allen section against distance from the midline, P4, P7, P14 (IGL plus white matter) and P56, with the range of your 20x sections shaded">']
    svg.append(f'<rect x="{ml}" y="{Y(max(a20)):.1f}" width="{cw - ml - mr}" height="{Y(min(a20)) - Y(max(a20)):.1f}" class="band"/>')
    svg.append(f'<text x="{cw - mr - 6}" y="{Y(max(a20)) - 6:.1f}" text-anchor="end" class="bandlab">your 20x sections: {rng(a20)} mm²</text>')
    for v in np.arange(ymin, ymax + 1e-9, ystep):
        svg.append(f'<line x1="{ml}" x2="{cw - mr}" y1="{Y(v):.1f}" y2="{Y(v):.1f}" class="grid"/>'
                   f'<text x="{ml - 8}" y="{Y(v) + 4:.1f}" text-anchor="end" class="tick">{v:.1f}</text>')
    for v in np.arange(xmin, xmax + 1, 500):
        svg.append(f'<text x="{X(v):.1f}" y="{ch - mb + 18}" text-anchor="middle" class="tick">{"0" if v == 0 else um(v)}</text>')
    svg.append(f'<text x="{(ml + cw - mr) / 2}" y="{ch - 8}" text-anchor="middle" class="axis">distance from midline (µm)</text>')
    svg.append(f'<text x="14" y="{(mt + ch - mb) / 2}" text-anchor="middle" class="axis" transform="rotate(-90 14 {(mt + ch - mb) / 2})">IGL area (mm²)</text>')
    for a in AGES:
        xs, ys = lines[a]
        pts = " ".join(f"{X(x):.1f},{Y(y):.1f}" for x, y in zip(xs, ys))
        svg.append(f'<polyline points="{pts}" class="line {a.lower()}"/>')
    svg.append("</svg>")
    chart = "\n".join(svg)

    # ---- numbers for the prose ----------------------------------------------------------------------
    xs7, ys7 = lines["P7"]

    def runs(mask):
        """Contiguous stretches of sections where mask holds, as 'a to b µm' text."""
        out, i = [], 0
        while i < len(mask):
            if mask[i]:
                j = i
                while j + 1 < len(mask) and mask[j + 1]:
                    j += 1
                out.append(f"{um(xs7[i])} to {um(xs7[j])} µm" if j > i else f"{um(xs7[i])} µm")
                i = j + 1
            else:
                i += 1
        return out

    def join(parts):
        return parts[0] if len(parts) == 1 else ", ".join(parts[:-1]) + " and " + parts[-1]

    inside, above = runs((ys7 >= min(a20)) & (ys7 <= max(a20))), runs(ys7 > max(a20))
    k7 = int(np.argmax(ys7))
    p7_band = ("the smoothed Allen P7 curve " + (f"lies inside it at {join(inside)}" if inside else "never enters it")
               + (f", and above it at {join(above)} (peak {ys7[k7]:.2f} mm² at {um(xs7[k7])} µm)" if above else ""))
    within = lambda a, l: [e for e in allen[a] if abs(e["offset_um"]) <= l]   # noqa: E731
    pct = {a: [e["igl_pct"] for e in within(a, GALLERY[a][1])] for a in AGES}
    n_unc = sum("orientation uncertain" in e["notes"] for e in user)
    s12a = next(e for e in if20 if e["name"].startswith("s1_2"))
    pcl_pct = [100 * e["pcl_mm2"] / e["tissue_mm2"] for e in user]
    p56_lo, p56_hi = chk["p56_plate_offset_um"]
    thick = {a: [e["thickness_um"] for e in within(a, GALLERY[a][1])] for a in AGES}
    pieces = {a: [e["pieces"] for e in within(a, GALLERY[a][1])] for a in AGES}
    lp = {a: [e["largest_piece_pct"] for e in within(a, GALLERY[a][1])] for a in AGES}
    t_if = {"20x": [e["thickness_um"] for e in if20], "10x": [e["thickness_um"] for e in if10]}
    n_allen = {a: len(allen[a]) for a in AGES}
    px = {s: next(e["px"] for e in E if e["source"] == s) for s in ("IF 20x", "IF 10x", "IF gold", *(f"Allen {a}" for a in AGES))}
    steps = sorted({s for s, _ in GALLERY.values()})

    holes = {a: [e["holes"] for e in within(a, GALLERY[a][1])] for a in AGES}
    tfill = {a: [e["thickness_filled_um"] for e in within(a, GALLERY[a][1])] for a in AGES}
    h20, tf20 = [e["holes"] for e in if20], [e["thickness_filled_um"] for e in if20]
    lp20, lp10 = [e["largest_piece_pct"] for e in if20], [e["largest_piece_pct"] for e in if10]
    top2_20 = []
    for e in if20:
        cc, _ = ndi.label(z[f"igl_{e['id']}"], structure=np.ones((3, 3), bool))
        s = np.sort(np.bincount(cc.ravel())[1:])[::-1]
        top2_20.append(100 * s[:2].sum() / s.sum())

    def allen_head(a):
        step, l = GALLERY[a]
        what = "IGL + white matter" if a == "P14" else "IGL"
        mp = float(np.median(pieces[a]))
        return (f"Every {step} µm from {um(-l)} to {um(l)} µm ({len(gallery[a])} of {n_allen[a]} sections). "
                f"Within that range: {what} {rng(pct[a], '{:.0f}')} % of tissue, median thickness {med(thick[a])} µm, "
                f"median {mp:.0f} {'piece' if mp == 1 else 'pieces'}, median {med(lp[a])} % of the {what} in the largest piece.")

    def holes_note(a):
        return (f" Small holes in the IGL (median {med(holes[a])} per section; your 20x sections {med(h20)}) lower the "
                f"measured thickness: with holes under {defs['hole_fill_um2']:.0f} µm² filled the median would be "
                f"{med(tfill[a])} µm (your 20x sections {med(tf20)} µm).")

    page = f"""<title>Cerebellar IGL Library</title>
<meta name="claude:viewport-fit-sides" content="cover">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans+Condensed:wght@500;600&family=IBM+Plex+Sans:wght@400;500&display=swap">
<style>
/* layout: specimen plates on one scale, grouped by source; sources and chart first */
:root {{
  --bg: #f2f4f7; --surface: #ffffff; --fg: #18202b; --muted: #5a6575; --line: #d6dce4;
  --accent: #1d5b86; --igl: #e2b41e; --warn-fg: #8a4d00; --warn-bg: #fbe9d4;
  --good-fg: #1f6b45; --good-bg: #dcf0e4; --off-fg: #5a6575; --off-bg: #e8ecf1;
  --p4: #7b8fa5; --p7: #1d5b86; --p14: #2e7d5b; --p56: #6b4f9a;
  --band: rgba(226,180,30,0.20); --band-fg: #7a5a00;
  --tissue: rgba(128,132,140,0.35); --tissue-edge: rgba(128,132,140,0.8);
  --display: "IBM Plex Sans Condensed", "Arial Narrow", system-ui, sans-serif;
  --text: "IBM Plex Sans", system-ui, -apple-system, sans-serif;
  --mono: "IBM Plex Mono", ui-monospace, Menlo, monospace;
}}
@media (prefers-color-scheme: dark) {{ :root:not([data-theme="light"]) {{
  --bg: #11151b; --surface: #19202a; --fg: #e5e9ef; --muted: #9aa5b3; --line: #2c3542;
  --accent: #82b8e2; --warn-fg: #f3bb78; --warn-bg: #3a2a17;
  --good-fg: #8fd9b0; --good-bg: #183528; --off-fg: #9aa5b3; --off-bg: #232b36;
  --p4: #8496ab; --p7: #82b8e2; --p14: #74c9a0; --p56: #b39ddb;
  --band: rgba(226,180,30,0.20); --band-fg: #f0cd5a; color-scheme: dark; }} }}
:root[data-theme="dark"] {{
  --bg: #11151b; --surface: #19202a; --fg: #e5e9ef; --muted: #9aa5b3; --line: #2c3542;
  --accent: #82b8e2; --warn-fg: #f3bb78; --warn-bg: #3a2a17;
  --good-fg: #8fd9b0; --good-bg: #183528; --off-fg: #9aa5b3; --off-bg: #232b36;
  --p4: #8496ab; --p7: #82b8e2; --p14: #74c9a0; --p56: #b39ddb;
  --band: rgba(226,180,30,0.20); --band-fg: #f0cd5a; color-scheme: dark; }}
body {{ background: var(--bg); color: var(--fg); font-family: var(--text); font-size: 15px; line-height: 1.55; margin: 0; }}
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
.sw.igl {{ background: var(--igl); }} .sw.ti {{ background: var(--tissue); border: 1px solid var(--tissue-edge); }}
section {{ display: flex; flex-direction: column; gap: 14px; }}
.sec-head {{ display: flex; flex-wrap: wrap; align-items: baseline; justify-content: space-between; gap: 8px 24px; }}
.sec-head p {{ margin: 0; color: var(--muted); max-width: 70ch; font-size: 14px; }}
.tablewrap {{ overflow-x: auto; border: 1px solid var(--line); border-radius: 6px; background: var(--surface); }}
table {{ border-collapse: collapse; width: 100%; font-size: 14px; }}
th, td {{ text-align: left; padding: 9px 14px; border-bottom: 1px solid var(--line); vertical-align: top; }}
th {{ font-family: var(--mono); font-weight: 500; font-size: 12px; letter-spacing: 0.05em; text-transform: uppercase; color: var(--muted); }}
tr:last-child td {{ border-bottom: 0; }}
td.num {{ font-family: var(--mono); font-variant-numeric: tabular-nums; white-space: nowrap; }}
.flag {{ color: var(--warn-fg); font-weight: 500; }}
.chart {{ background: var(--surface); border: 1px solid var(--line); border-radius: 6px; padding: 14px; overflow-x: auto; }}
.chart svg {{ width: 100%; min-width: 520px; height: auto; display: block; }}
.chart .grid {{ stroke: var(--line); stroke-width: 1; }}
.chart .tick {{ fill: var(--muted); font: 11px var(--mono); }}
.chart .axis {{ fill: var(--muted); font: 12px var(--text); }}
.chart .band {{ fill: var(--band); }}
.chart .bandlab {{ fill: var(--band-fg); font: 12px var(--text); }}
.chart .line {{ fill: none; stroke-width: 2.2; stroke-linecap: round; }}
.chart .p4 {{ stroke: var(--p4); }} .chart .p7 {{ stroke: var(--p7); }}
.chart .p14 {{ stroke: var(--p14); stroke-dasharray: 1 5; stroke-width: 3; }}
.chart .p56 {{ stroke: var(--p56); stroke-dasharray: 6 4; }}
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
.note {{ background: var(--surface); border: 1px dashed var(--line); border-radius: 6px; padding: 14px 16px; color: var(--muted); font-size: 14px; max-width: 80ch; }}
nav {{ display: flex; flex-wrap: wrap; gap: 8px 18px; font-size: 14px; }}
nav a {{ color: var(--accent); text-decoration: none; border-bottom: 1px solid transparent; }}
nav a:hover, nav a:focus-visible {{ border-bottom-color: var(--accent); outline: none; }}
</style>
<div class="wrap">
<header>
  <div class="eyebrow">Mouse cerebellum · sagittal · internal granule layer</div>
  <h1>Cerebellar IGL Library</h1>
  <p class="lede">Every internal granule layer (IGL) shape we have: your P7 immunofluorescence sections (the automatic labels after the anatomy fixes, plus the hand-painted s1_2) and the Allen reference Nissl series at P4, P7, P14 and P56. All plates share one scale ({thumb_um:.0f} µm per pixel) and one orientation: anterior to the left, dorsal up. {len(E)} shapes are in the library file; this page shows all {len(user)} of yours and a sample of the Allen series every {steps[0]}–{steps[-1]} µm.</p>
  <div class="legend"><span><i class="sw igl"></i>internal granule layer</span><span><i class="sw ti"></i>rest of the cerebellum</span></div>
  <nav style="margin-top:16px" aria-label="Sections"><a href="#sources">Sources</a><a href="#area">Area by position</a><a href="#yours">Your sections</a>{''.join(f'<a href="#{a.lower()}">Allen {a}</a>' for a in AGES)}</nav>
</header>

<section id="sources">
  <div class="sec-head"><h2>Sources</h2><p>Card numbers: IGL area, IGL as % of the cerebellum on the plate, pieces (8-connected IGL parts of at least {defs['piece_min_um2']:.0f} µm²), and thickness (twice the median distance from the IGL centre line to its edge, measured on the {grid:.0f} µm grid).</p></div>
  <div class="tablewrap"><table>
    <thead><tr><th>Source</th><th>Shapes</th><th>Resolution</th><th>IGL from</th><th>Reliability</th></tr></thead>
    <tbody>
      <tr><td>Your sections, 20x</td><td class="num">{len(if20)}</td><td class="num">{px['IF 20x']:.2f} µm/px</td><td>MATLAB labels after anatomy fixes; IGL label only (Purkinje layer kept apart)</td><td>s1_2 against the gold standard: IGL IoU {chk['s1_2_igl_iou_native_gold_assigned']:.3f} at native resolution, leaving out tissue the gold standard left unassigned ({chk['s1_2_igl_iou_native']:.3f} over all pixels). Rotated into atlas orientation with the white-matter library's rotations; {n_unc} of {len(user)} orientations uncertain (mirror image scored within {defs['uncertain_mirror_margin']}).</td></tr>
      <tr><td>Your sections, 10x</td><td class="num">{len(if10)}</td><td class="num">{px['IF 10x']:.2f} µm/px</td><td>same (camera scale; the MATLAB 10x scale is wrong)</td><td>No gold standard at 10x. Controls from the EZH2 experiment.</td></tr>
      <tr><td>Gold standard s1_2</td><td class="num">{len(gold)}</td><td class="num">{px['IF gold']:.2f} µm/px</td><td>hand-painted IGL (stored under the ML code; swapped back)</td><td>Reference.</td></tr>
      <tr><td>Allen P4 Nissl</td><td class="num">{n_allen['P4']}</td><td class="num">{px['Allen P4']:.2f} µm/px</td><td>IGL grey-level class of the Nissl stain</td><td>Not checked against any annotation.</td></tr>
      <tr><td>Allen P7 Nissl</td><td class="num">{n_allen['P7']}</td><td class="num">{px['Allen P7']:.2f} µm/px</td><td>same method</td><td>Not checked against any annotation.</td></tr>
      <tr><td>Allen P14 Nissl</td><td class="num">{n_allen['P14']}</td><td class="num">{px['Allen P14']:.2f} µm/px</td><td>dark Nissl core</td><td><span class="flag">IGL + white matter, not the IGL alone.</span> White matter is as dark as the IGL in P14 Nissl and cannot be separated, so every P14 area and thickness includes the white matter. Flagged on every P14 plate.</td></tr>
      <tr><td>Allen P56 Nissl</td><td class="num">{n_allen['P56']}</td><td class="num">{px['Allen P56']:.2f} µm/px</td><td>dark Nissl class</td><td>Checked against {chk['p56_plates']} annotated atlas plates ({um(p56_lo)} to {um(p56_hi)} µm): IGL Dice {chk['p56_igl_dice_min']:.2f}–{chk['p56_igl_dice_max']:.2f}.</td></tr>
    </tbody></table></div>
</section>

<section id="area">
  <div class="sec-head"><h2>Area by position</h2><p>IGL area of each Allen section against its distance from the midline ({SMOOTH}-section running mean). P14 (dotted) is the IGL + white-matter core, so it is not comparable with the other ages. The shaded band is the range of your 20x sections ({rng(a20)} mm²); {p7_band}. Your 10x sections range from {min(a10):.2f} to {max(a10):.2f} mm². As a share of the cerebellum, the IGL is {rng(pct['P7'], '{:.0f}')} % in Allen P7 within {um(-GALLERY['P7'][1])} to {um(GALLERY['P7'][1])} µm, against {rng([e['igl_pct'] for e in if20], '{:.0f}')} % in your 20x and {rng([e['igl_pct'] for e in if10], '{:.0f}')} % in your 10x sections.</p></div>
  <div class="chart">{chart}
    <div class="keys"><span><i style="border-color:var(--p4)"></i>P4</span><span><i style="border-color:var(--p7)"></i>P7</span><span><i style="border-color:var(--p14);border-top-style:dotted"></i>P14 (IGL + white matter)</span><span><i style="border-color:var(--p56);border-top-style:dashed"></i>P56</span></div>
  </div>
</section>

<section id="yours">
  <div class="sec-head"><h2>Your P7 sections</h2><p>All {len(user)}: {len(if20)} slides at 20x, {len(if10)} at 10x, and the gold standard. The automatic s1_2 has {s12a['igl_mm2']:.2f} mm² of IGL against {gold[0]['igl_mm2']:.2f} mm² in the gold standard. "Orientation uncertain" means the best rotation and its mirror image fit the template almost equally well, so this plate may be flipped front to back. The Purkinje cell layer is not counted as IGL; it is a further {rng(pcl_pct, '{:.1f}')} % of the tissue. Median IGL thickness is {med(t_if['20x'])} µm at 20x ({med(tf20)} µm with holes under {defs['hole_fill_um2']:.0f} µm² filled) and {med(t_if['10x'])} µm at 10x. At 20x the IGL usually falls into two large parts: the largest piece holds {rng(lp20, '{:.0f}')} % of the IGL and the two largest together {rng(top2_20, '{:.0f}')} %, while at 10x the median largest piece holds {med(lp10)} %. The hand-painted s1_2 is split the same way (largest piece {gold[0]['largest_piece_pct']:.0f} %), so this is in the sections, not only in the automatic labels.</p></div>
  <div class="grid-cards">
{user_cards}
  </div>
</section>

<section id="p4">
  <div class="sec-head"><h2>Allen P4</h2><p>{allen_head('P4')} IGL from Nissl grey-level classes, not checked against any annotation.{holes_note('P4')}</p></div>
  <div class="grid-cards">
{allen_cards['P4']}
  </div>
</section>

<section id="p7">
  <div class="sec-head"><h2>Allen P7</h2><p>{allen_head('P7')} IGL from Nissl grey-level classes, not checked against any annotation.{holes_note('P7')}</p></div>
  <div class="grid-cards">
{allen_cards['P7']}
  </div>
</section>

<section id="p14">
  <div class="sec-head"><h2>Allen P14</h2><p><span class="flag">These plates show the IGL + white-matter core, not the IGL alone.</span> {allen_head('P14')}</p></div>
  <div class="grid-cards">
{allen_cards['P14']}
  </div>
</section>

<section id="p56">
  <div class="sec-head"><h2>Allen P56</h2><p>{allen_head('P56')} IGL Dice against the atlas {chk['p56_igl_dice_min']:.2f}–{chk['p56_igl_dice_max']:.2f} on {chk['p56_plates']} plates.</p></div>
  <div class="grid-cards">
{allen_cards['P56']}
  </div>
</section>

<p class="note">Library files: outputs/igl_library/library.npz holds every IGL and tissue mask at {grid:.0f} µm per pixel, and outputs/igl_library/library.json holds their metadata (source, position, orientation, areas, pieces, thickness, skeleton length, Purkinje layer area for your sections). Built by scripts/igl_library.py; this page by scripts/igl_library_page.py.</p>
</div>
"""
    p = D / "igl_library.html"
    p.write_text(page)
    print(p, f"{p.stat().st_size / 1e6:.2f} MB", {a: len(v) for a, v in gallery.items()}, "user", len(user))


if __name__ == "__main__":
    main()
