"""
Evaluate ribbon_join.py on the gold slide (s1_2) with synthetic breaks, then run it
on the eight 20x slides.

    python scripts/ribbon_join_eval.py

Input labels: depth_model/out/anatomy_fixes/<slide>_fixed_labels.tif (R1-R4 applied),
downsampled 4x (2.05 um/px), as anatomy_fixes.py decides.

Synthetic breaks (copy what the fragmented slides show between pieces):
  EGL -> ML 0.6 / IGL 0.4      ML -> IGL 0.5 / DWL 0.3 / PCL 0.2
  IGL -> ML 0.5 / DWL 0.3 / PCL 0.2
  each break removes a 30-150 um stretch of the ribbon (measured along its skeleton)
  at a spot where the gold also has that layer; 12 EGL, 12 ML, 6 IGL breaks per trial.

Scores (gold accuracy as in depth_model/tune_relabel.py: gold pixels labelled 2-7,
9 classes, gold IGL/ML swap undone):
  clean, clean+join (false-bridge cost), broken, broken+join (recovery),
  each with and without the template veto.
-> outputs/ribbon_join/eval.json, run on all slides -> outputs/ribbon_join/slides.json
"""
import json
import sys
import time
from collections import deque
from pathlib import Path

import numpy as np
import tifffile
from scipy import ndimage as ndi
from skimage import morphology

sys.path.insert(0, str(Path(__file__).parent))
import ribbon_join as RJ  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
SEG = Path("/Users/jpurzner/Dropbox/imaging_analysis/cerebellar_segmentation/python")
FIXED = SEG / "depth_model/out/anatomy_fixes"
GOLD_SLIDE = "2018_05_22_s1_2_p27-0021"
OUT = ROOT / "outputs/ribbon_join"
SLIDES = ["2018_05_22_s1_2_p27-0021", "2018_05_22_s1_4_p27-0020", "2018_05_22_s1_5_p27-0018",
          "2018_05_22_s2_2_p27", "2018_05_22_s2_4_p27-0024", "2018_05_22_s2_5_p27-0026",
          "2018_05_22_s3_2_p27-0001", "2018_05_22_s3_5_p27-0006"]
DS = 4
PX = 0.5119 * DS
SCORED = [2, 3, 4, 5, 6, 7]
BREAKS = {"EGL": (12, {5: 0.6, 4: 0.4}), "ML": (12, {4: 0.5, 6: 0.3, 7: 0.2}),
          "IGL": (6, {5: 0.5, 6: 0.3, 7: 0.2})}
LEN_UM = (30.0, 150.0)
SEEDS = range(5)


def load_gold():
    g = tifffile.imread(SEG / "labelled" / f"{GOLD_SLIDE}_labelled.tiff").astype(np.int16)[::DS, ::DS]
    gs = g.copy()
    gs[g == 4], gs[g == 5] = 5, 4
    return gs.astype(np.uint8)


def acc(lab, gold):
    m = np.isin(gold, SCORED)
    return float((lab[m] == gold[m]).mean())


def ious(lab, gold):
    m = np.isin(gold, SCORED)
    out = {}
    for name, ls in (("EGL", (2, 3)), ("ML", (5,)), ("IGL", (4,)), ("DWL", (6,)), ("PCL", (7,))):
        a, b = np.isin(lab, ls) & m, np.isin(gold, ls) & m
        out[name] = round(float((a & b).sum() / max((a | b).sum(), 1)), 3)
    return out


def continuity(lab):
    out = {}
    for name, ls in (("EGL", (2, 3)), ("ML", (5,)), ("IGL", (4,))):
        cc, n = ndi.label(np.isin(lab, ls), structure=RJ.ST)
        a = np.bincount(cc.ravel())[1:]
        out[name] = dict(largest_pct=round(100 * float(a.max() / a.sum()), 1),
                         n_ge500um2=int((a * PX * PX >= 500).sum()))
    return out


def geodesic(sk, c, maxd):
    H, W = sk.shape
    seen = {tuple(c): 0}
    q = deque([tuple(c)])
    while q:
        p = q.popleft()
        d = seen[p]
        if d >= maxd:
            continue
        for dy in (-1, 0, 1):
            for dx in (-1, 0, 1):
                r = (p[0] + dy, p[1] + dx)
                if r not in seen and 0 <= r[0] < H and 0 <= r[1] < W and sk[r]:
                    seen[r] = d + 1
                    q.append(r)
    return np.array(list(seen))


def damage(lab, gold, seed):
    rng = np.random.default_rng(seed)
    out = lab.copy()
    hit = np.zeros(lab.shape, bool)
    for layer, (n, into) in BREAKS.items():
        codes = RJ.RIBBONS[layer]["codes"]
        M = np.isin(out, codes)
        sk = morphology.skeletonize(M)
        idx = ndi.distance_transform_edt(~sk, return_distances=False, return_indices=True)
        cand = np.argwhere(sk & np.isin(gold, codes))
        rng.shuffle(cand)
        cen, made = [], 0
        for c in cand:
            if made == n:
                break
            if any(np.hypot(*(c - q)) < 300 / PX for q in cen) or hit[tuple(c)]:
                continue
            L = rng.uniform(*LEN_UM) / PX
            seg = geodesic(sk, c, L / 2)
            segm = np.zeros(lab.shape, bool)
            segm[seg[:, 0], seg[:, 1]] = True
            reg = M & segm[idx[0], idx[1]] & ~hit
            if reg.sum() * PX * PX < 100:
                continue
            k = rng.choice(list(into), p=list(into.values()))
            out[reg] = k
            hit |= ndi.binary_dilation(reg, iterations=3)
            cen.append(c)
            made += 1
    return out, out != lab


def contacts(f):
    """layer-order violations (% of the relevant edge)."""
    ext = ndi.binary_fill_holes(f > 0) == 0
    d = lambda m: ndi.binary_dilation(m, structure=RJ.ST)
    tis = (f > 0) & (f != 9)
    egl, ml, igl, wm, pcl = np.isin(f, (2, 3)), f == 5, f == 4, np.isin(f, (6, 8)), f == 7
    pia = tis & d(ext)
    eb = egl & ~ndi.binary_erosion(egl)
    mb = ml & ~ndi.binary_erosion(ml)
    # EGL gaps inside fissures: non-EGL pixels within 20 um of two large EGL pieces, > 40 um deep
    cc, n = ndi.label(egl, structure=RJ.ST)
    ar = np.bincount(cc.ravel())
    R = int(round(20 / PX))
    disk = morphology.disk(R)
    cov = np.zeros(f.shape, np.int16)
    for i, sl in enumerate(ndi.find_objects(cc), 1):
        if ar[i] * PX * PX < 500:
            continue
        sl2 = tuple(slice(max(a.start - R - 1, 0), a.stop + R + 1) for a in sl)
        cov[sl2] += ndi.binary_dilation(cc[sl2] == i, structure=disk)
    dO = ndi.distance_transform_edt(~ext) * PX
    fis_gap = (cov >= 2) & ~egl & tis & (dO > 40)
    return {"fissure EGL gap mm2": round(float(fis_gap.sum() * PX * PX / 1e6), 4),
            "ML/IGL/PCL at pia %": round(100 * float((pia & (ml | igl | pcl)).sum() / pia.sum()), 1),
            "EGL touching IGL %": round(100 * float((eb & d(igl)).sum() / eb.sum()), 1),
            "ML touching WM %": round(100 * float((mb & d(wm)).sum() / mb.sum()), 1)}


NAMES = {2: "EGL", 3: "EGL", 4: "IGL", 5: "ML", 6: "DWL", 7: "PCL"}


def evaluate(fixed, gold, rules, broken):
    prm = {"rules": rules}
    j0, _ = RJ.join(fixed, PX, prm)
    r = dict(rules=list(rules), clean=acc(fixed, gold), clean_join=acc(j0, gold),
             clean_changed_mm2=float((j0 != fixed).sum() * PX * PX / 1e6),
             clean_iou_before=ious(fixed, gold), clean_iou_after=ious(j0, gold), trials=[])
    per = {}
    for s, (br, H) in enumerate(broken):
        jb, _ = RJ.join(br, PX, prm)
        r["trials"].append(dict(seed=s, broken=acc(br, gold), joined=acc(jb, gold),
                                damaged_mm2=float(H.sum() * PX * PX / 1e6),
                                restored_frac=float((jb[H] == fixed[H]).mean()),
                                collateral_mm2=float(((jb != br) & ~H).sum() * PX * PX / 1e6)))
        for a in (2, 3, 4, 5):
            for b in (4, 5, 6, 7):
                m = H & (fixed == a) & (br == b)
                if m.any():
                    k = f"{NAMES[a]}->{NAMES[b]}"
                    t = per.setdefault(k, [0, 0])
                    t[0] += int(m.sum())
                    t[1] += int((np.isin(jb[m], (2, 3)) if a in (2, 3) else jb[m] == fixed[m]).sum())
    T = r["trials"]
    r["mean"] = {k: round(float(np.mean([t[k] for t in T])), 4)
                 for k in ("broken", "joined", "restored_frac", "damaged_mm2", "collateral_mm2")}
    gap = r["clean"] - r["mean"]["broken"]
    r["mean"]["recovered_share_of_lost_accuracy"] = round((r["mean"]["joined"] - r["mean"]["broken"]) / gap, 3)
    r["restored_by_damage"] = {k: round(v[1] / v[0], 2) for k, v in sorted(per.items())}
    m = r["mean"]
    print(f"  rules {'+'.join(rules):6s} clean {r['clean']:.4f} -> {r['clean_join']:.4f} "
          f"({r['clean_changed_mm2']:.4f} mm2 changed) | broken {m['broken']:.4f} -> {m['joined']:.4f}, "
          f"recovered {m['recovered_share_of_lost_accuracy']:.2f} of the lost accuracy", flush=True)
    print("     restored by damage type: " + "  ".join(f"{k} {v:.2f}" for k, v in r["restored_by_damage"].items()),
          flush=True)
    return r


RULESETS = (("P", "M", "W"), ("P", "M", "W", "F", "S"))


def main():
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    gold = load_gold()
    fixed = tifffile.imread(FIXED / f"{GOLD_SLIDE}_fixed_labels.tif")[::DS, ::DS]
    print(f"s1_2 anatomy-fixed labels: gold acc {acc(fixed, gold):.4f}", flush=True)
    broken = [damage(fixed, gold, s) for s in SEEDS]
    res = dict(defaults=RJ.DEFAULTS, len_um=LEN_UM,
               breaks={k: [v[0], {str(a): b for a, b in v[1].items()}] for k, v in BREAKS.items()},
               gold_contacts=contacts(gold),
               rulesets={"+".join(r): evaluate(fixed, gold, r, broken) for r in RULESETS})
    (OUT / "eval.json").write_text(json.dumps(res, indent=1, default=float))

    print(f"\nall 20x slides (gold slide contacts: {res['gold_contacts']})", flush=True)
    rows = {}
    for s in SLIDES:
        f = tifffile.imread(FIXED / f"{s}_fixed_labels.tif")[::DS, ::DS]
        rows[s] = dict(before=dict(contacts=contacts(f), continuity=continuity(f)))
        save = {"before": f}
        for rules in RULESETS:
            j, _ = RJ.join(f, PX, {"rules": rules})
            key = "+".join(rules)
            rows[s][key] = dict(contacts=contacts(j), continuity=continuity(j),
                                changed_mm2=float((j != f).sum() * PX * PX / 1e6),
                                changed_pct_tissue=round(100 * float((j != f).sum() / ((f > 0) & (f != 9)).sum()), 2))
            if s == GOLD_SLIDE:
                rows[s][key]["gold_acc"] = acc(j, gold)
            save[key.replace("+", "")] = j
        np.savez_compressed(OUT / f"{s}_joined_ds4.npz", **save)
        b = rows[s]["before"]
        line = f"  {s[-14:]:14s} before: {b['contacts']}  EGL/ML/IGL largest " + \
               "/".join(f"{b['continuity'][L]['largest_pct']:.0f}" for L in ("EGL", "ML", "IGL"))
        for rules in RULESETS:
            a = rows[s]["+".join(rules)]
            line += f"\n      {'+'.join(rules):6s} changed {a['changed_pct_tissue']:.2f}% of tissue: {a['contacts']}  largest " + \
                    "/".join(f"{a['continuity'][L]['largest_pct']:.0f}" for L in ("EGL", "ML", "IGL"))
        print(line, flush=True)
    (OUT / "slides.json").write_text(json.dumps(rows, indent=1, default=float))
    print(f"done in {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
