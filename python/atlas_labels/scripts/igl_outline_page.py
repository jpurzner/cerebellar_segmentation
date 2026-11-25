"""
Gap Filler Test page: benchmark summary of the outline-driven IGL gap filler on synthetic
gaps (igl_outline_eval.py) and an example gallery with a pick-the-best answer per gap.

    python scripts/igl_outline_page.py     # -> outputs/igl_outline/gap_filler_test.html
"""
import base64
import io
import json
import sys
from pathlib import Path

import numpy as np
from PIL import Image

sys.path.insert(0, str(Path(__file__).parent))
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/igl_outline"
TEMPLATE = Path(__file__).parent / "igl_outline_template.html"

# methods shown in the gallery (key in examples -> label)
COMB = "combined"
SWITCH = 0.35        # white-matter share of the shallow IGL band above which the labels are distrusted
SHOW = [("ribbon model", "Ribbon model"), ("walk (current)", "Walk (current)"), ("label local geo|band", "Outline + WM labels"),
        ("ribbon model, WM ignored", "Ribbon model, WM ignored"), (COMB, "Combined (recommended)")]
# methods shown in the summary table
TABLE = [("walk (current)", "Walk (current pipeline)"),
         ("ribbon model", "Ribbon gap model (white-matter labels are a barrier)"),
         ("ribbon model, WM ignored", "Ribbon gap model, white-matter labels ignored"),
         ("label local geo|band", "Outline + white-matter labels (m from the flanks)"),
         ("label global geo|band", "  m fitted on the whole section"),
         ("label local band", "  outline edges along the ribbon model's route"),
         ("gated local geo|band", "Outline + white-matter labels gated by the outline core"),
         ("depth+core local geo|band", "Outline, no white-matter labels (depth band minus outline core)"),
         ("outline local geo|band", "Outline core as the only white matter"),
         (COMB, "Combined: outline + white-matter labels; ribbon model ignoring the labels when the outline test fires")]


def score(r, k):
    if k == COMB:
        k = "ribbon model, WM ignored" if r["wm_shallow"] > SWITCH else "label local geo|band"
    return r["scores"].get(k)


def stats(gaps, sc, k):
    rr = [score(r, k) for r in gaps if r["scenario"] == sc]
    if not rr or any(v is None for v in rr):
        return None
    d = np.array([v["dice"] for v in rr])
    return dict(mean=round(float(d.mean()), 3), median=round(float(np.median(d)), 3), ge08=round(float((d >= 0.8).mean()), 2),
                lt05=round(float((d < 0.5).mean()), 2), prec=round(float(np.mean([v["prec"] for v in rr])), 3),
                rec=round(float(np.mean([v["rec"] for v in rr])), 3))


PLATE = (11, 14, 18)
TISSUE = (58, 66, 77)
EGL = (122, 90, 166)
WM = (140, 115, 72)
IGL = (170, 180, 192)
TP = (59, 224, 255)
FP = (255, 176, 0)
FN = (255, 59, 138)


def png(rgb, scale):
    im = Image.fromarray(rgb.astype(np.uint8)).resize((rgb.shape[1] * scale, rgb.shape[0] * scale), Image.NEAREST)
    b = io.BytesIO()
    im.save(b, "PNG", optimize=True)
    return "data:image/png;base64," + base64.b64encode(b.getvalue()).decode()


def base(ex):
    rgb = np.zeros(ex["tissue"].shape + (3,), np.uint8)
    rgb[:] = PLATE
    rgb[ex["tissue"]] = TISSUE
    rgb[ex["wm"] & ex["tissue"]] = WM
    rgb[ex["egl"]] = EGL
    rgb[ex["igl"]] = IGL
    return rgb


def panel_truth(ex):
    rgb = base(ex)
    rgb[ex["gap"]] = TP
    return rgb


def panel_fill(ex, fill):
    rgb = base(ex)
    gap = ex["gap"]
    rgb[fill & gap] = TP
    rgb[fill & ~gap] = FP
    rgb[gap & ~fill] = FN
    return rgb


def main():
    B = json.loads((OUT / "benchmark_s1_2.json").read_text())
    E = np.load(OUT / "examples_s1_2.npz")
    others = []
    for t in ("s5_F", "s5_G", "s6_C", "s6_G"):
        p = OUT / f"benchmark_{t}.json"
        if p.exists():
            others.append(json.loads(p.read_text()))
    og = [r for x in others for r in x["gaps"]]
    rows_by = {}
    for sc in ("erased", "as_wm"):
        rows = []
        for k, lab in TABLE:
            v = stats(B["gaps"], sc, k)
            if v is None:
                continue
            o = stats(og, sc, k) if og else None
            rows.append(dict(key=k, label=lab, **v, others=o["mean"] if o else None))
        rows_by[sc] = rows
    gaps = {(r["scenario"], r["gap"]): r for r in B["gaps"]}
    items = []
    keys = sorted({tuple(k.split("|")[:2]) for k in E.files}, key=lambda t: (t[0] != "erased", int(t[1])))
    for sc, gi in keys:
        ex = {k.split("|", 2)[2]: E[k] for k in E.files if k.startswith(f"{sc}|{gi}|")}
        r = gaps[(sc, int(gi))]
        h = ex["tissue"].shape[0]
        scale = 3 if max(ex["tissue"].shape) <= 220 else 2
        panels = [dict(key="truth", label="Gold standard (gap in cyan)", img=png(panel_truth(ex), scale))]
        if "ribbon model, WM ignored" in ex and "label local geo|band" in ex:
            ex[COMB] = ex["ribbon model, WM ignored"] if r["wm_shallow"] > SWITCH else ex["label local geo|band"]
        for k, lab in SHOW:
            if k in ex:
                panels.append(dict(key=k, label=lab, img=png(panel_fill(ex, ex[k].astype(bool)), scale), dice=score(r, k)["dice"]))
        items.append(dict(id=f"{sc}-{gi}", scenario=sc, gap=int(gi), length_um=r["length_um"], width=r["width"], panels=panels,
                          switched=bool(r["wm_shallow"] > SWITCH),
                          um_per_px=4.0 / scale, h=h))
    switched = {sc: round(float(np.mean([r["wm_shallow"] > SWITCH for r in B["gaps"] if r["scenario"] == sc])), 3) for sc in ("erased", "as_wm")}
    sw_others = {sc: round(float(np.mean([r["wm_shallow"] > SWITCH for r in og if r["scenario"] == sc])), 3) if og else None
                 for sc in ("erased", "as_wm")}
    # findings, from the numbers
    allg = B["gaps"] + og
    A = lambda k, sc, G=allg: stats(G, sc, k)
    W, R, L, RI, CB = "walk (current)", "ribbon model", "label local geo|band", "ribbon model, WM ignored", COMB
    f = []
    f.append(f"<b>Both new fillers beat the walk the pipeline uses now.</b> Over all five controls, plain erased gaps are refilled with a mean Dice of "
             f"{A(L, 'erased')['mean']:.2f} by the outline filler and {A(R, 'erased')['mean']:.2f} by the ribbon gap model, against {A(W, 'erased')['mean']:.2f} for the walk, "
             f"which fills {round(100 * A(W, 'erased')['lt05'])}% of gaps badly (Dice below 0.5).")
    if og:
        f.append(f"<b>On plain gaps the outline filler and the ribbon model tie.</b> On s1_2 the ribbon model is a little better "
                 f"({A(R, 'erased', B['gaps'])['mean']:.2f} against {A(L, 'erased', B['gaps'])['mean']:.2f}); on the four 10x controls the outline filler is "
                 f"({A(L, 'erased', og)['mean']:.2f} against {A(R, 'erased', og)['mean']:.2f}).")
    f.append(f"<b>When the missing IGL carries a white-matter label, everything that trusts the labels fails.</b> The outline filler then fills nothing and "
             f"the ribbon model reaches {A(R, 'as_wm')['mean']:.2f}. The outline catches this case: real white matter almost never sits in the shallow half of "
             f"the expected IGL band, so white matter there means the label is wrong. That test fired in {round(100 * np.mean([r['wm_shallow'] > SWITCH for r in allg if r['scenario'] == 'as_wm']))}% "
             f"of mislabelled gaps and {round(100 * np.mean([r['wm_shallow'] > SWITCH for r in allg if r['scenario'] == 'erased']))}% of plain ones.")
    f.append(f"<b>Combined</b> (outline filler with the labels; when the test fires, the ribbon model with the white-matter labels ignored): "
             f"{A(CB, 'erased')['mean']:.2f} on plain gaps and {A(CB, 'as_wm')['mean']:.2f} on mislabelled ones, where the ribbon model alone gets {A(R, 'as_wm')['mean']:.2f}. "
             f"Ignoring the labels everywhere would give {A(RI, 'erased')['mean']:.2f} on both, costing accuracy on the plain gaps.")
    f.append("<b>The outline cannot replace the white-matter labels.</b> White matter predicted from the outline alone (a core along each folium's centre line) "
             "matches the labels only moderately (Dice 0.5 to 0.7 per section), and the fill built on it reaches "
             f"{A('outline local geo|band', 'erased')['mean']:.2f}. Its use is as a check on the labels, not as their substitute.")
    fits = B["whole_section_fits"]
    data = dict(findings=f, rows=rows_by, items=items, n=B["n_gaps"], fits=fits, n_others=len(others), switched=switched, switched_others=sw_others,
                others=[x["section"] for x in others])
    page = TEMPLATE.read_text().replace("__DATA__", json.dumps(data, default=float))
    dst = OUT / "gap_filler_test.html"
    dst.write_text(page)
    print(dst, round(dst.stat().st_size / 1e6, 2), "MB", len(items), "examples")


if __name__ == "__main__":
    main()
