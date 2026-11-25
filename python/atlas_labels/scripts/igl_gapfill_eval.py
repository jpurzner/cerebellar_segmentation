"""
End-to-end test of the pipeline gap filler (igl_gapfill.fill_gaps) against the walk
(igl_repair.repair_gaps): unlike igl_outline_eval.py, the gaps' ends are not given; each
method has to find and pair them itself, as on a real section.

    python scripts/igl_gapfill_eval.py

Sections: the good-IGL controls (s1_2 gold labels, s5_F, s5_G, s6_C, s6_G).
Per section and seed: 4 gaps of 50-400 um cut into the ribbon (igl_repair_eval.cut_gaps
sites); the removed IGL is relabelled as
  ML       (NeuN still IGL-bright, so the NeuN check should let it through)
  anomaly  (bleach-like, label 9)
  WM       (bleached IGL called white matter)
Scores per gap: recall (share of the removed IGL painted back), precision near the gap
(painted pixels within 150 um of it that were IGL before the damage), Dice; and painted area
far from every gap (false fills). Clean: painted area on the undamaged sections (every pixel a
false fill unless the section has a real break). Gold: s1_2 automatic labels against the gold
standard, IGL IoU before / after.
-> outputs/igl_gapfill/eval.json (env GAPFILL_PRM = filler settings as JSON, GAPFILL_EVAL_OUT = file name)
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import tifffile
from scipy import ndimage as ndi

sys.path.insert(0, str(Path(__file__).parent))
import igl_gapfill as GF  # noqa: E402
import igl_repair as IR  # noqa: E402
import igl_repair_eval as EV  # noqa: E402
import igl_ribbon as RB  # noqa: E402
import igl_ribbon_eval as RE  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/igl_gapfill"
SEEDS = (0, 1, 2)
SCEN = ("ML", "anomaly", "WM")
import os  # noqa: E402
PRM = json.loads(os.environ.get("GAPFILL_PRM", "{}"))           # filler settings to test
OUTNAME = os.environ.get("GAPFILL_EVAL_OUT", "eval.json")
METHODS = {"walk": lambda g, neun, tear: IR.repair_gaps(g, None, None, neun, tear),
           "gap filler": lambda g, neun, tear: GF.fill_gaps(g, prm=PRM, neun=neun, tear=tear)}


def damage(g, gaps, kind):
    g2 = {k: v.copy() for k, v in g.items()}
    for q in gaps:
        m = q["mask"]
        g2["igl"] &= ~m
        if kind == "ML":
            g2["ml"] |= m
        elif kind == "anomaly":
            g2["anom"] |= m
            g2["tissue"] &= ~m
        else:
            g2["wm"] |= m
    return g2


def score(painted, g, gaps):
    near_all = np.zeros(painted.shape, bool)
    res = []
    for q in gaps:
        m = q["mask"]
        near = ndi.distance_transform_edt(~m) * RB.GRID <= 150
        near_all |= near
        p = painted & near
        tp = int((p & m).sum())
        res.append(dict(length_um=q["length_um"], rec=round(tp / max(m.sum(), 1), 3),
                        prec=round(float((p & g["igl"]).sum() / max(p.sum(), 1)), 3) if p.any() else None,
                        dice=round(2 * tp / max(int(p.sum()) + int(m.sum()), 1), 3)))
    far = painted & ~near_all
    return res, round(float(far.sum() * RB.GRID ** 2)), round(float((far & ~g["igl"]).sum() * RB.GRID ** 2))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    trials, clean = [], {}
    t0 = time.time()
    for name, lab, px in RE.sections():
        g = RB.grid_labels(lab, px)
        neun, tear = EV.images(name, g["igl"].shape)
        rib0 = RB.describe(g)
        clean[name] = {}
        for meth, fn in METHODS.items():
            new, br, _ = fn(g, neun, tear)
            clean[name][meth] = dict(painted_um2=round(float((new & ~g["igl"]).sum() * RB.GRID ** 2)),
                                     joined=[{k: v for k, v in b.items() if not k.startswith("_")} for b in br if b["status"] == "joined"])
        print(f"{name}: clean painted " + ", ".join(f"{m} {v['painted_um2']} um2" for m, v in clean[name].items()), flush=True)
        for seed in SEEDS:
            rng = np.random.default_rng(100 + seed)
            _, gaps = EV.cut_gaps(g, rib0, rng)
            for kind in SCEN:
                g2 = damage(g, gaps, kind)
                for meth, fn in METHODS.items():
                    new, br, _ = fn(g2, neun, tear)
                    painted = new & ~g2["igl"]
                    res, far, far_bad = score(painted, g, gaps)
                    feats = []
                    if meth == "gap filler":
                        gm = np.any([q["mask"] for q in gaps], axis=0)
                        for b in br:
                            if "_mask" not in b or not b["status"] in ("joined", "ok") and not b["status"].startswith("superseded"):
                                continue
                            y0, y1, x0, x1 = b["_box"]
                            f = b["_mask"]
                            feats.append(dict({k: v for k, v in b.items() if not k.startswith("_") and k not in ("a", "b")},
                                              hit=round(float(gm[y0:y1, x0:x1][f].mean()), 3), was_igl=round(float(g["igl"][y0:y1, x0:x1][f].mean()), 3)))
                    trials.append(dict(section=name, seed=seed, scenario=kind, method=meth, gaps=res, far_um2=far, far_not_igl_um2=far_bad,
                                       joined=sum(b["status"] == "joined" for b in br), bridges=feats))
                print(f"  seed {seed} {kind:7s} " + " | ".join(
                    f"{t['method']}: dice " + " ".join(f"{r['dice']:.2f}" for r in t["gaps"]) + f" far {t['far_not_igl_um2']}"
                    for t in trials[-len(METHODS):]) + f"  ({time.time() - t0:.0f}s)", flush=True)
    summ = {}
    for kind in SCEN:
        for meth in METHODS:
            rr = [r for t in trials if t["scenario"] == kind and t["method"] == meth for r in t["gaps"]]
            d = np.array([r["dice"] for r in rr])
            summ[f"{kind} | {meth}"] = dict(n=len(rr), mean_dice=round(float(d.mean()), 3), restored=round(float(np.mean([r["rec"] > 0.5 for r in rr])), 2),
                                            mean_rec=round(float(np.mean([r["rec"] for r in rr])), 3),
                                            mean_prec=round(float(np.nanmean([r["prec"] if r["prec"] is not None else np.nan for r in rr])), 3),
                                            false_fill_um2_per_section=round(float(np.mean([t["far_not_igl_um2"] for t in trials if t["scenario"] == kind and t["method"] == meth]))))
    # no-harm check on the gold slide: the automatic labels of s1_2
    auto = tifffile.imread(RE.FX / "2018_05_22_s1_2_p27-0021_fixed_labels.tif")
    g = RB.grid_labels(auto, 0.5119)
    gg = RB.grid_labels(RE.gold(), 0.5119)
    neun, tear = EV.images("s1_2 auto", g["igl"].shape)
    sc = gg["tissue"]
    iou = lambda a, b: float((a & b & sc).sum() / max(((a | b) & sc).sum(), 1))
    gold = dict(igl_iou_before=round(iou(g["igl"], gg["igl"]), 4))
    for meth, fn in METHODS.items():
        new, br, _ = fn(g, neun, tear)
        add = new & ~g["igl"]
        gold[meth] = dict(igl_iou_after=round(iou(new, gg["igl"]), 4), painted_um2=round(float(add.sum() * RB.GRID ** 2)),
                          painted_that_is_gold_igl=round(float((add & gg["igl"]).sum() / max(add.sum(), 1)), 3),
                          joined=sum(b["status"] == "joined" for b in br))
    print(json.dumps(summ, indent=1))
    print("gold s1_2 auto:", gold)
    print("clean:", {k: {m: v["painted_um2"] for m, v in c.items()} for k, c in clean.items()})
    (OUT / OUTNAME).write_text(json.dumps(dict(prm=PRM, summary=summ, gold_s1_2_auto=gold, clean=clean, trials=trials), indent=1,
                                              default=lambda o: o.tolist() if hasattr(o, "tolist") else float(o)))


if __name__ == "__main__":
    main()
