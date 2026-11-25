"""
Page of the viable IGL / white-matter shapes and the shapes generated from them
(run shape_profile.py stage0 + features first).

    python scripts/shape_profile_page.py  -> outputs/shape_profile/shape_profiles.html
"""
import base64
import html
import io
import json
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage as ndi

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / "outputs/shape_profile"
GRID = 4.0
FEAT = {"igl": dict(name="IGL", rgb=(240, 190, 40)), "wm": dict(name="White matter", rgb=(236, 120, 48))}


def uri(rgba):
    buf = io.BytesIO()
    Image.fromarray(rgba, "RGBA").save(buf, "PNG", optimize=True)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def edge(m):
    return m & ~ndi.binary_erosion(m)


def plate(feature_mask, tissue, outline, rgb, step, prob=None):
    """RGBA plate: this shape's tissue faint, feature filled (or probability ramp), mean outline."""
    H, W = outline.shape
    img = np.zeros((H, W, 4), np.float32)
    if tissue is not None:
        img[tissue] = (150, 156, 166, 70)
    if prob is not None:
        a = np.clip(prob, 0, 1)
        lo = np.array([60, 40, 10], np.float32)
        hi = np.array(rgb, np.float32) * 0.55 + 255 * 0.45
        col = lo[None, None] * (1 - a[..., None]) + hi[None, None] * a[..., None]
        on = a > 0.02
        img[on, :3] = col[on]
        img[on, 3] = 60 + 195 * a[on]
    if feature_mask is not None:
        img[feature_mask] = (*rgb, 255)
    img[edge(outline)] = (235, 240, 245, 255)
    img = img[::step, ::step]
    return uri(img.astype(np.uint8)), img.shape[1], img.shape[0]


def main():
    z = np.load(D / "stage1_features.npz")
    st0 = np.load(D / "stage0_shape_set.npz")
    meta = {m["id"]: m for m in json.loads((D / "stage0.json").read_text())["shapes"]}
    res = json.loads((D / "stage1_features.json").read_text())
    outline = st0["ref_tissue"]
    A = st0["aligned"]
    order = list(st0["ids"])
    secs = []
    for key, f in FEAT.items():
        r = res[key]
        gen = []
        u, w, h = plate(None, None, outline, f["rgb"], 2, prob=z[f"{key}_prob"])
        gen.append(("Probability, all viable shapes", f"your {r['n_if']} sections and {r['n_p56']} P56 weighted equally", u, w, h))
        u, w, h = plate(None, None, outline, f["rgb"], 2, prob=z[f"{key}_prob_if"])
        gen.append(("Probability, your sections only", f"{r['n_if']} control sections", u, w, h))
        u, w, h = plate(None, None, outline, f["rgb"], 2, prob=z[f"{key}_prob_p56"])
        gen.append(("Probability, Allen P56 only", f"{r['n_p56']} sections within ±0.6 mm, scaled to P7", u, w, h))
        u, w, h = plate(z[f"{key}_mean"], None, outline, f["rgb"], 2)
        gen.append(("Mean shape", f"{r['area_mean_mm2']:.2f} mm² (sections: {min(r['area_shapes_mm2']):.2f}–{max(r['area_shapes_mm2']):.2f} mm²)", u, w, h))
        modes = z[f"{key}_modes"]
        for j in range(2):
            for k, sign in ((0, "−2 SD"), (1, "+2 SD")):
                u, w, h = plate(modes[j, k], None, outline, f["rgb"], 3)
                gen.append((f"Mode {j + 1}, {sign}", f"mode {j + 1} explains {100 * (r['variance_explained'][j] - (r['variance_explained'][j - 1] if j else 0)):.0f} % of variation", u, w, h))
        # viable shapes gallery
        idx = z[f"{key}_idx"]
        shapes = z[f"{key}_shapes"]
        gal_if, gal_p56 = [], []
        for n, i in enumerate(idx):
            sid = order[i]
            m = meta[sid]
            is_if = m["source"].startswith("IF")
            u, w, h = plate(shapes[n], A[i, 0], outline, f["rgb"], 2 if is_if else 3)
            area = shapes[n].sum() * GRID ** 2 / 1e6
            if is_if:
                grade = m["igl_grade"] if key == "igl" else m["wm_grade"]
                nm = "s1_2 (gold labels)" if sid == "gold_s1_2" else sid.replace("if_", "").replace("_p27", "")
                gal_if.append((nm, f"{grade} · {area:.2f} mm²", u, w, h))
            else:
                gal_p56.append(("midline" if m["offset_um"] == 0 else f"{m['offset_um']:+.0f} µm", f"{area:.2f} mm²", u, w, h))
        # LOO chart data
        loo = r["loo"]
        rows = [(("s1_2" if k == "gold_s1_2" else k.replace("if_", "").replace("_p27", "")), v["mean"], v["other_median"], v["other_best"])
                for k, v in loo.items()]
        secs.append((key, f, r, gen, gal_if, gal_p56, rows))

    def cards(lst, cls):
        out = []
        for t, sub, u, w, h in lst:
            out.append(f'<figure class="{cls}"><div class="plate"><img src="{u}" width="{w}" height="{h}" alt="{html.escape(t)}" loading="lazy"></div>'
                       f'<figcaption><span class="name">{html.escape(t)}</span><span class="sub">{html.escape(sub)}</span></figcaption></figure>')
        return "\n".join(out)

    def chart(rows, title):
        n = len(rows)
        W, rowh, ml, mr, mt = 640, 28, 110, 24, 30
        H = mt + n * rowh + 34
        X = lambda v: ml + v * (W - ml - mr)
        g = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="{html.escape(title)}">']
        for v in np.arange(0, 1.01, 0.2):
            g.append(f'<line x1="{X(v):.1f}" x2="{X(v):.1f}" y1="{mt - 8}" y2="{H - 30}" class="grid"/>'
                     f'<text x="{X(v):.1f}" y="{H - 12}" text-anchor="middle" class="tick">{v:.1f}</text>')
        for k, (nm, mean, med, best) in enumerate(rows):
            y = mt + k * rowh + rowh / 2
            g.append(f'<text x="{ml - 10}" y="{y + 4:.1f}" text-anchor="end" class="lab">{html.escape(nm)}</text>')
            g.append(f'<line x1="{X(min(mean, med)):.1f}" x2="{X(max(mean, med)):.1f}" y1="{y:.1f}" y2="{y:.1f}" class="link"/>')
            g.append(f'<circle cx="{X(med):.1f}" cy="{y:.1f}" r="5" class="s2"><title>{nm}: another section, median Dice {med:.2f} (best {best:.2f})</title></circle>')
            g.append(f'<circle cx="{X(mean):.1f}" cy="{y:.1f}" r="5" class="s1"><title>{nm}: generated mean shape, Dice {mean:.2f}</title></circle>')
        g.append("</svg>")
        return "\n".join(g)

    body = []
    for key, f, r, gen, gal_if, gal_p56, rows in secs:
        better = sum(1 for _, m_, md, _b in rows if m_ > md)
        body.append(f"""
<section id="{key}">
  <div class="sec-head"><h2>{f['name']}</h2>
  <p>{r['n_if']} of your control sections and {r['n_p56']} Allen P56 sections, all aligned on the outline (white line: mean control outline). The probability map is the share of shapes that have {f['name'].lower()} at each point.</p></div>
  <h3>Generated shapes</h3>
  <div class="gen">{cards(gen[:4], 'card big')}</div>
  <div class="gen small">{cards(gen[4:], 'card')}</div>
  <div class="scale"><span>0 %</span><span class="ramp ramp-{key}"></span><span>100 % of shapes</span></div>
  <h3>How well the mean shape predicts a held-out section</h3>
  <p class="note">Each of your sections was left out, the model rebuilt from the rest, and the held-out {f['name'].lower()} compared with the generated mean shape (blue). Orange is the typical overlap with another of your sections: a generated shape is only useful if it beats that. The mean beats it for {better} of {len(rows)} sections.</p>
  <div class="chart">{chart(rows, f['name'] + ' leave-one-out overlap')}
    <div class="keys"><span><i class="k1"></i>generated mean shape</span><span><i class="k2"></i>another of your sections (median)</span><span>Dice overlap, 0 = none, 1 = identical</span></div>
  </div>
  <h3>Viable shapes: your sections ({len(gal_if)})</h3>
  <div class="grid-if">{cards(gal_if, 'card')}</div>
  <h3>Viable shapes: Allen P56 ({len(gal_p56)})</h3>
  <div class="grid-p56">{cards(gal_p56, 'card tiny')}</div>
</section>""")

    igl, wm = res["igl"], res["wm"]
    page = f"""<title>Cerebellar Shape Profiles</title>
<meta name="claude:viewport-fit-sides" content="cover">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans+Condensed:wght@500;600&family=IBM+Plex+Sans:wght@400;500&display=swap">
<style>
/* layout: summary, then per layer: generated shapes, leave-one-out chart, viable-shape galleries */
:root {{
  --bg: #f2f4f7; --surface: #ffffff; --fg: #18202b; --muted: #5a6575; --line: #d6dce4; --plate: #0e1218;
  --accent: #1d5b86; --s1: #2a78d6; --s2: #eb6834;
  --display: "IBM Plex Sans Condensed", "Arial Narrow", system-ui, sans-serif;
  --text: "IBM Plex Sans", system-ui, -apple-system, sans-serif;
  --mono: "IBM Plex Mono", ui-monospace, Menlo, monospace;
}}
@media (prefers-color-scheme: dark) {{ :root:not([data-theme="light"]) {{
  --bg: #11151b; --surface: #19202a; --fg: #e5e9ef; --muted: #9aa5b3; --line: #2c3542; --accent: #82b8e2;
  --s1: #3987e5; --s2: #d95926; color-scheme: dark; }} }}
:root[data-theme="dark"] {{
  --bg: #11151b; --surface: #19202a; --fg: #e5e9ef; --muted: #9aa5b3; --line: #2c3542; --accent: #82b8e2;
  --s1: #3987e5; --s2: #d95926; color-scheme: dark; }}
body {{ background: var(--bg); color: var(--fg); font-family: var(--text); font-size: 15px; line-height: 1.55; }}
.wrap {{ max-width: 1180px; margin: 0 auto; padding-block: 36px 64px;
  padding-left: max(16px, env(safe-area-inset-left, 0px)); padding-right: max(16px, env(safe-area-inset-right, 0px));
  display: flex; flex-direction: column; gap: 36px; }}
h1, h2, h3 {{ font-family: var(--display); font-weight: 600; margin: 0; text-wrap: balance; }}
h1 {{ font-size: clamp(30px, 4.5vw, 44px); line-height: 1.1; }}
h2 {{ font-size: 26px; }}
h3 {{ font-size: 18px; margin-top: 8px; }}
.lede {{ max-width: 70ch; color: var(--muted); margin: 10px 0 0; }}
.eyebrow {{ font-family: var(--mono); font-size: 12px; letter-spacing: 0.08em; text-transform: uppercase; color: var(--accent); }}
section {{ display: flex; flex-direction: column; gap: 14px; }}
.sec-head p, .note {{ margin: 0; color: var(--muted); max-width: 72ch; font-size: 14px; }}
.tablewrap {{ overflow-x: auto; border: 1px solid var(--line); border-radius: 6px; background: var(--surface); }}
table {{ border-collapse: collapse; width: 100%; font-size: 14px; }}
th, td {{ text-align: left; padding: 9px 14px; border-bottom: 1px solid var(--line); }}
th {{ font-family: var(--mono); font-weight: 500; font-size: 12px; letter-spacing: 0.05em; text-transform: uppercase; color: var(--muted); }}
tr:last-child td {{ border-bottom: 0; }}
td.num {{ font-family: var(--mono); font-variant-numeric: tabular-nums; }}
.card {{ margin: 0; background: var(--surface); border: 1px solid var(--line); border-radius: 6px; overflow: hidden; display: flex; flex-direction: column; min-width: 0; }}
.plate {{ background: var(--plate); display: flex; align-items: center; justify-content: center; padding: 6px; }}
.plate img {{ max-width: 100%; height: auto; image-rendering: pixelated; }}
figcaption {{ padding: 8px 10px; display: flex; flex-direction: column; gap: 2px; border-top: 1px solid var(--line); }}
.name {{ font-weight: 500; font-size: 14px; }}
.sub {{ font: 12px var(--mono); color: var(--muted); }}
.gen {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(250px, 1fr)); gap: 12px; }}
.gen.small {{ grid-template-columns: repeat(auto-fill, minmax(180px, 1fr)); }}
.grid-if {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(180px, 1fr)); gap: 10px; }}
.grid-p56 {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(118px, 1fr)); gap: 8px; }}
.tiny figcaption {{ padding: 5px 7px; }} .tiny .name {{ font-size: 12px; }}
.scale {{ display: flex; align-items: center; gap: 10px; font: 12px var(--mono); color: var(--muted); }}
.ramp {{ display: inline-block; width: 180px; height: 10px; border-radius: 2px; }}
.ramp-igl {{ background: linear-gradient(90deg, rgb(60,40,10), rgb(247,219,137)); }}
.ramp-wm {{ background: linear-gradient(90deg, rgb(60,40,10), rgb(245,181,141)); }}
.chart {{ background: var(--surface); border: 1px solid var(--line); border-radius: 6px; padding: 12px 14px; overflow-x: auto; }}
.chart svg {{ width: 100%; max-width: 720px; min-width: 460px; height: auto; display: block; }}
.chart .grid {{ stroke: var(--line); stroke-width: 1; }}
.chart .tick {{ fill: var(--muted); font: 11px var(--mono); }}
.chart .lab {{ fill: var(--fg); font: 12px var(--mono); }}
.chart .link {{ stroke: var(--muted); stroke-width: 1; opacity: 0.5; }}
.chart .s1 {{ fill: var(--s1); stroke: var(--surface); stroke-width: 2; }}
.chart .s2 {{ fill: var(--s2); stroke: var(--surface); stroke-width: 2; }}
.chart circle:hover {{ r: 7; }}
.keys {{ display: flex; flex-wrap: wrap; gap: 8px 18px; font-size: 13px; color: var(--muted); margin-top: 8px; }}
.keys i {{ display: inline-block; width: 10px; height: 10px; border-radius: 50%; margin-right: 6px; vertical-align: -1px; }}
.k1 {{ background: var(--s1); }} .k2 {{ background: var(--s2); }}
nav {{ display: flex; flex-wrap: wrap; gap: 8px 18px; font-size: 14px; }}
nav a {{ color: var(--accent); text-decoration: none; border-bottom: 1px solid transparent; }}
nav a:hover, nav a:focus-visible {{ border-bottom-color: var(--accent); outline: none; }}
</style>
<div class="wrap">
<header>
  <div class="eyebrow">P7 control profile · IGL and deep white matter</div>
  <h1>Cerebellar Shape Profiles</h1>
  <p class="lede">The viable shapes are the ones you graded good or usable on your control sections, plus the Allen P56 sections you graded good. All are aligned on their outline in one frame (anterior left, dorsal up; adult P56 scaled to P7 size). From them I generate a probability map, a mean shape and the main modes of variation for each layer.</p>
  <nav style="margin-top:14px" aria-label="Sections"><a href="#summary">Summary</a><a href="#igl">IGL</a><a href="#wm">White matter</a></nav>
</header>
<section id="summary">
  <div class="tablewrap"><table>
    <thead><tr><th>Layer</th><th>Your sections</th><th>Allen P56</th><th>Mean shape vs held-out section (Dice)</th><th>Another section vs held-out (median Dice)</th></tr></thead>
    <tbody>
      <tr><td>IGL</td><td class="num">{igl['n_if']}</td><td class="num">{igl['n_p56']}</td><td class="num">{min(v['mean'] for v in igl['loo'].values()):.2f}–{max(v['mean'] for v in igl['loo'].values()):.2f}</td><td class="num">{min(v['other_median'] for v in igl['loo'].values()):.2f}–{max(v['other_median'] for v in igl['loo'].values()):.2f}</td></tr>
      <tr><td>White matter</td><td class="num">{wm['n_if']}</td><td class="num">{wm['n_p56']}</td><td class="num">{min(v['mean'] for v in wm['loo'].values()):.2f}–{max(v['mean'] for v in wm['loo'].values()):.2f}</td><td class="num">{min(v['other_median'] for v in wm['loo'].values()):.2f}–{max(v['other_median'] for v in wm['loo'].values()):.2f}</td></tr>
    </tbody></table></div>
  <p class="note">Aligned on the outline alone, folia sit in different places from section to section, so the layers overlap only partly. The IGL mean shape is slightly better than using another section's IGL. The white-matter mean collapses to the trunk ({wm['area_mean_mm2']:.2f} mm² against {min(wm['area_shapes_mm2']):.2f}–{max(wm['area_shapes_mm2']):.2f} mm² in the sections) because the branches never line up; the probability map is the usable white-matter profile until the branches are matched to each other.</p>
</section>
{''.join(body)}
</div>
"""
    p = D / "shape_profiles.html"
    p.write_text(page)
    print(p, f"{p.stat().st_size / 1e6:.2f} MB")


if __name__ == "__main__":
    main()
