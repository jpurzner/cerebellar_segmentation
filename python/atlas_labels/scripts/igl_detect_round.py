"""
Round 2 of the new IGL step (2026-11-20): compose MATLAB's or the best current IGL with the
automatic additions and removals, and score round 1 (trained on the 12 controls) against
round 2 (+ the user's 6 corrected cKO sections) on sections each model had not seen.

  add:    cells the model calls IGL (p >= 0.5) that the base IGL does not have, in pieces >= 10,000 um2
  remove: base IGL the model calls EGL (p >= 0.5), in pieces >= 5,000 um2 (tuned on the user's
          round-1 erasures: 63 % of them, 68 % precision; the inner EGL pulled into the cKO IGL)

    python scripts/igl_detect_round.py           # -> outputs/igl_detect/round2_eval.json
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi

sys.path.insert(0, str(Path(__file__).parent))
import igl_detect as D  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/igl_detect"
G = D.G
RULES = dict(t_add=0.5, add_min_um2=10000.0, t_egl=0.5, remove_min_um2=5000.0)


def pieces(m, min_um2):
    lab, n = ndi.label(m, structure=np.ones((3, 3), bool))
    if not n:
        return m
    size = np.bincount(lab.ravel()) * G * G
    keep = size >= min_um2
    keep[0] = False
    return keep[lab]


def compose(P, base, tissue, rules=RULES):
    """P: class probabilities (EGL, ML/PCL, IGL, WM); base: IGL to start from."""
    add = pieces((P[2] >= rules["t_add"]) & ~base & tissue, rules["add_min_um2"])
    rem = pieces(base & (P[0] >= rules["t_egl"]), rules["remove_min_um2"])
    return (base | add) & ~rem, add, rem


def probs(short, rnd):
    """held-out probabilities for a section: round 1 = loso2 (controls) / apply2 (cKO, trained on
    controls only); round 2 = loso2_r2 (the 18 training sections) / apply2_r2 (s1_5)."""
    if rnd == 1:
        f = OUT / ("loso2" if short in D.controls() else "apply2") / f"{short}.npz"
    else:
        f = OUT / ("loso2_r2" if short in D.training_sections() else "apply2_r2") / f"{short}.npz"
    return np.load(f)["P"].astype(np.float32)


def evaluate():
    rows = []
    for short in D.training_sections():
        _, y, auto, tissue = D.load(short)
        truth = y == 4                               # reviewed / corrected IGL
        mat = auto == 4
        known = (y > 0) | mat
        changed = (truth != mat) & known             # cells the review / the user changed
        r = dict(section=short, group="control" if short in D.controls() else "cKO",
                 iou_matlab=round(D.iou(mat & known, truth & known), 3), changed_um2=int(changed.sum() * G * G))
        for rnd in (1, 2):
            new, add, rem = compose(probs(short, rnd), mat, tissue)
            r[f"r{rnd}"] = dict(iou=round(D.iou(new & known, truth & known), 3),
                                edits_reproduced=round(float((new[changed] == truth[changed]).mean()), 2) if changed.any() else None,
                                extra_changes_um2=int(((new != mat) & ~changed & known).sum() * G * G),
                                added_um2=int(add.sum() * G * G), removed_um2=int(rem.sum() * G * G))
        rows.append(r)
        print(f"{short:5s} {r['group']:7s} MATLAB {r['iou_matlab']:.3f} | round 1 {r['r1']['iou']:.3f} reproduces {r['r1']['edits_reproduced']} "
              f"(+{r['r1']['extra_changes_um2']} um2 other changes) | round 2 {r['r2']['iou']:.3f} reproduces {r['r2']['edits_reproduced']} "
              f"(+{r['r2']['extra_changes_um2']} um2 other changes)", flush=True)
    summ = {}
    for grp in ("control", "cKO"):
        rr = [r for r in rows if r["group"] == grp]
        summ[grp] = dict(n=len(rr), iou_matlab=round(float(np.mean([r["iou_matlab"] for r in rr])), 3),
                         **{f"r{k}_iou": round(float(np.mean([r[f"r{k}"]["iou"] for r in rr])), 3) for k in (1, 2)},
                         **{f"r{k}_reproduced": round(float(np.sum([r[f"r{k}"]["edits_reproduced"] * r["changed_um2"] for r in rr if r["changed_um2"]]) /
                                                            max(np.sum([r["changed_um2"] for r in rr if r["changed_um2"]]), 1)), 2) for k in (1, 2)},
                         **{f"r{k}_extra_um2": int(np.mean([r[f"r{k}"]["extra_changes_um2"] for r in rr])) for k in (1, 2)})
    print(json.dumps(summ, indent=1))
    (OUT / "round2_eval.json").write_text(json.dumps(dict(rules=RULES, summary=summ, sections=rows), indent=1))
    return rows, summ


if __name__ == "__main__":
    evaluate()
