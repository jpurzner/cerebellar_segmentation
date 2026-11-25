"""
Test the IGL gap repair (igl_repair.py) on synthetic gaps in the good-IGL controls.

    python scripts/igl_repair_eval.py

Per section and seed: 4 gaps of 50-400 um cut into the ribbon (sites as in
igl_ribbon_eval.py); removed IGL becomes ML (half of the gaps) or anomaly (bleach-like).
Scores:
  restored  - share of the removed IGL pixels that the repair paints back
  precision - share of all painted pixels that were IGL before the damage
  per gap   - restored if > 50 % of its removed pixels come back
  clean     - area painted when the repair runs on the undamaged section
  gold      - s1_2 automatic labels: IGL IoU and 9-class accuracy against the gold
              standard before / after repair
-> outputs/igl_ribbon/repair_eval.json
"""
import json
import sys
from pathlib import Path

import numpy as np
import tifffile
from scipy import ndimage as ndi
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).parent))
import igl_ribbon as RB  # noqa: E402
import igl_repair as IR  # noqa: E402
import igl_ribbon_eval as RE  # noqa: E402

OUT = RE.OUT
SEEDS = (0, 1, 2)


def cut_gaps(g, rib, rng, n=4):
    sites = RE.pick_sites(rib, rng, n=n)
    cl = np.vstack([b["path"] for b in rib["branches"]])
    owner = np.concatenate([[(bi, j) for j in range(len(b["path"]))] for bi, b in enumerate(rib["branches"])])
    igl = g["igl"].copy()
    yy, xx = np.nonzero(igl)
    _, nn = cKDTree(cl).query(np.stack([yy, xx], 1))
    ob, oj = owner[nn, 0], owner[nn, 1]
    g2 = {k: v.copy() for k, v in g.items()}
    gaps = []
    for k, (bi, j) in enumerate(sites):
        L = rng.uniform(50, 400)
        h = int(round(L / 2 / (RB.GRID * 1.1)))
        b = rib["branches"][bi]
        j0, j1 = max(j - h, 0), min(j + h, len(b["path"]) - 1)
        sel = (ob == bi) & (oj >= j0) & (oj <= j1)
        m = np.zeros(igl.shape, bool)
        m[yy[sel], xx[sel]] = True
        g2["igl"] &= ~m
        if k % 2 == 0:
            g2["ml"] |= m
        else:
            g2["anom"] |= m
            g2["tissue"] &= ~m
        gaps.append(dict(mask=m, length_um=round(L), filled_with="ML" if k % 2 == 0 else "anomaly"))
    return g2, gaps


SB20 = RE.SEG / "matlab_anomaly/sandbox"
SB10 = RE.SEG / "matlab_anomaly_10x/sandbox"


def images(name, shape):
    """NeuN (smoothed) and tear mask on the 4 um grid for a section."""
    if name.startswith("s1_2"):
        raw, ch, px = SB20 / "2018_05_22_s1_2_p27-0021_native.tif", 1, 0.5119
        tear = SB20 / "2018_05_22_s1_2_p27-0021_iter6b_native_anomaly_tear.tif"
    else:
        sid = name + "_p27"
        raw, ch, px = SB10 / f"{sid}_native.tif", 2, 0.65
        tear = SB10 / f"{sid}_native_anomaly_tear.tif"
    z = px / RB.GRID
    n = tifffile.imread(raw)[ch].astype(np.float32)
    n = ndi.zoom(ndi.gaussian_filter(n, 3.0 / px), z, order=1)[:shape[0], :shape[1]]
    t = (ndi.zoom((tifffile.imread(tear) > 0).astype(np.float32), z, order=1) >= 0.25)[:shape[0], :shape[1]] \
        if tear.exists() else np.zeros(shape, bool)
    def fit(a):
        out = np.zeros(shape, a.dtype)
        out[:a.shape[0], :a.shape[1]] = a
        return out
    return fit(n), fit(t)


def main(prm=None):
    rows, tot_removed, tot_restored, tot_painted, tot_painted_ok = [], 0, 0, 0, 0
    gaps_ok = gaps_n = 0
    clean = {}
    for name, lab, px in RE.sections():
        g = RB.grid_labels(lab, px)
        neun, tear = images(name, g["igl"].shape)
        rib0 = RB.describe(g)
        new0, br0, _ = IR.repair_gaps(g, prm, rib0, neun, tear)
        clean[name] = dict(painted_um2=round(float((new0 & ~g["igl"]).sum() * RB.GRID ** 2)),
                           ends=len(rib0["ends"]), statuses=[b["status"] for b in br0])
        for seed in SEEDS:
            rng = np.random.default_rng(100 + seed)
            g2, gaps = cut_gaps(g, rib0, rng)
            new, br, _ = IR.repair_gaps(g2, prm, None, neun, tear)
            painted = new & ~g2["igl"]
            removed = np.any([q["mask"] for q in gaps], axis=0)
            tot_removed += int(removed.sum())
            tot_restored += int((painted & removed).sum())
            tot_painted += int(painted.sum())
            tot_painted_ok += int((painted & g["igl"]).sum())
            res = []
            for q in gaps:
                fr = float((painted & q["mask"]).sum() / max(q["mask"].sum(), 1))
                gaps_n += 1
                gaps_ok += fr > 0.5
                res.append(f"{q['length_um']}um/{q['filled_with']} {100 * fr:.0f}%")
            rows.append(dict(section=name, seed=seed, gaps=res,
                             statuses=sorted({b["status"] for b in br}),
                             painted_um2=round(float(painted.sum() * RB.GRID ** 2)),
                             painted_outside_old_igl_um2=round(float((painted & ~g["igl"]).sum() * RB.GRID ** 2))))
            print(f"{name:12s} seed {seed}: restored per gap: {', '.join(res)} | painted outside original IGL "
                  f"{rows[-1]['painted_outside_old_igl_um2']} um2", flush=True)
    summ = dict(gaps=gaps_n, gaps_restored_over_half=gaps_ok,
                restored_share=round(tot_restored / max(tot_removed, 1), 3),
                precision=round(tot_painted_ok / max(tot_painted, 1), 3),
                clean_sections=clean)
    # no-harm check on the gold slide: repair the automatic labels of s1_2
    SEG = RE.SEG
    auto = tifffile.imread(RE.FX / "2018_05_22_s1_2_p27-0021_fixed_labels.tif")
    gold = RE.gold()
    g = RB.grid_labels(auto, 0.5119)
    gg = RB.grid_labels(gold, 0.5119)
    neun, tear = images("s1_2 auto", g["igl"].shape)
    new, br, _ = IR.repair_gaps(g, None, None, neun, tear)
    sc = gg["tissue"]
    def iou(a, b):
        return float((a & b & sc).sum() / max(((a | b) & sc).sum(), 1))
    summ["gold_s1_2"] = dict(igl_iou_before=round(iou(g["igl"], gg["igl"]), 4), igl_iou_after=round(iou(new, gg["igl"]), 4),
                             painted_um2=round(float((new & ~g["igl"]).sum() * RB.GRID ** 2)),
                             painted_that_is_gold_igl=round(float((new & ~g["igl"] & gg["igl"]).sum() / max((new & ~g["igl"]).sum(), 1)), 3),
                             bridges=[b for b in br if b["status"] != "too_long" and not b["status"].startswith("blocked")])
    print(json.dumps({k: v for k, v in summ.items() if k != "clean_sections"}, indent=1, default=float))
    print("clean sections:", {k: (v["painted_um2"], v["ends"], v["statuses"]) for k, v in clean.items()})
    (OUT / "repair_eval.json").write_text(json.dumps(dict(params={**IR.DEF, **(prm or {})}, summary=summ, trials=rows),
                                                     indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else float(o)))
    return summ


if __name__ == "__main__":
    main()
