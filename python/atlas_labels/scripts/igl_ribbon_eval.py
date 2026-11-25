"""
Synthetic-defect test of the IGL ribbon description (igl_ribbon.py) on the good-IGL
control sections (s1_2 gold labels, s5_F, s5_G, s6_C, s6_G).

    python scripts/igl_ribbon_eval.py

Per section and seed: 2 gaps (50-400 um of ribbon removed), 2 thinnings (ribbon cut to
the central 40 % of its width over 80-300 um), 2 bulges (ribbon widened to ~2x over
80-300 um). Sites lie on branches >= 600 um long, >= 150 um from branch ends, >= 300 um
apart. Detection:
  gap    - new ribbon ends within 60 um of both gap edges (one = partial)
  thin   - a thinning flag within 60 um of the stretch
  bulge  - a bulge flag within 60 um of the stretch
False alarms: ends / flags after damage that are not within 60 um of a defect and were
not already present (within 60 um) before damage.
-> outputs/igl_ribbon/synthetic_eval.json
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

ROOT = Path(__file__).resolve().parents[1]
SEG = Path("/Users/jpurzner/Dropbox/imaging_analysis/cerebellar_segmentation/python")
FX = SEG / "depth_model/out/anatomy_fixes"
OUT = ROOT / "outputs/igl_ribbon"
SEEDS = (0, 1, 2)
NEAR = 60.0 / RB.GRID


def gold():
    g = tifffile.imread(SEG / "labelled" / "2018_05_22_s1_2_p27-0021_labelled.tiff").astype(np.int16)
    gs = g.copy()
    gs[g == 4], gs[g == 5] = 5, 4
    return gs.astype(np.uint8)


def sections():
    yield "s1_2 (gold)", gold(), 0.5119
    for s in ("s5_F_p27", "s5_G_p27", "s6_C_p27", "s6_G_p27"):
        yield s.replace("_p27", ""), tifffile.imread(FX / f"{s}_fixed_labels.tif"), 0.65


def pick_sites(rib, rng, n=6, min_branch_um=600, margin_um=150, sep_um=300):
    cands = []
    for bi, b in enumerate(rib["branches"]):
        if b["length_um"] < min_branch_um:
            continue
        m = int(margin_um / (RB.GRID * 1.1))
        for j in range(m, len(b["path"]) - m, 5):
            cands.append((bi, j))
    rng.shuffle(cands)
    sites = []
    for bi, j in cands:
        p = rib["branches"][bi]["path"][j]
        if all(np.hypot(*(p - rib["branches"][q[0]]["path"][q[1]])) * RB.GRID > sep_um for q in sites):
            sites.append((bi, j))
        if len(sites) == n:
            break
    return sites


def damage(g, rib, sites, rng):
    igl = rib["igl"].copy()
    cl = np.vstack([b["path"] for b in rib["branches"]])
    owner = np.concatenate([[(bi, j) for j in range(len(b["path"]))] for bi, b in enumerate(rib["branches"])])
    tree = cKDTree(cl)
    yy, xx = np.nonzero(igl)
    d, nn = tree.query(np.stack([yy, xx], 1))
    own_b, own_j = owner[nn, 0], owner[nn, 1]
    kinds = ["gap", "gap", "thin", "thin", "bulge", "bulge"]
    defects = []
    tissue = g["tissue"]
    for (bi, j), kind in zip(sites, kinds):
        b = rib["branches"][bi]
        L = rng.uniform(50, 400) if kind == "gap" else rng.uniform(80, 300)
        h = int(round(L / 2 / (RB.GRID * 1.1)))
        j0, j1 = max(j - h, 0), min(j + h, len(b["path"]) - 1)
        sel = (own_b == bi) & (own_j >= j0) & (own_j <= j1)
        if kind == "gap":
            igl[yy[sel], xx[sel]] = False
        elif kind == "thin":
            wloc = b["width"][own_j[sel]]
            cut = d[sel] * RB.GRID > 0.2 * wloc
            igl[yy[sel][cut], xx[sel][cut]] = False
        else:
            seg = b["path"][j0:j1 + 1]
            m = np.zeros(igl.shape, bool)
            m[seg[:, 0], seg[:, 1]] = True
            wl = float(np.median(b["width"][j0:j1 + 1]))
            igl |= (ndi.distance_transform_edt(~m) * RB.GRID <= wl) & tissue
        defects.append(dict(kind=kind, branch=bi, j0=int(j0), j1=int(j1), length_um=round(L),
                            edges=[b["path"][max(j0 - 1, 0)].tolist(), b["path"][min(j1 + 1, len(b["path"]) - 1)].tolist()],
                            stretch=b["path"][j0:j1 + 1]))
    g2 = dict(g)
    g2["igl"] = igl
    return g2, defects


def near(points, q, r=NEAR):
    return len(points) > 0 and np.min(np.hypot(*(np.asarray(points) - q).T)) <= r


def main(prm=None, quiet=False):
    OUT.mkdir(parents=True, exist_ok=True)
    rows, per_kind = [], {"gap": [0, 0, 0], "thin": [0, 0], "bulge": [0, 0]}
    fa_ends = fa_flags = 0
    clean_stats = {}
    for name, lab, px in sections():
        g = RB.grid_labels(lab, px)
        r0 = RB.describe(g, {**RB.P, **(prm or {})})
        e0 = [e["pos"] for e in r0["ends"]]
        f0 = [p for f in r0["flags"] for p in f["pts"]]
        clean_stats[name] = dict(ends=len(e0), thin=sum(f["kind"] == "thin" for f in r0["flags"]),
                                 bulge=sum(f["kind"] == "bulge" for f in r0["flags"]), length_mm=round(r0["length_mm"], 1))
        for seed in SEEDS:
            rng = np.random.default_rng(seed)
            sites = pick_sites(r0, rng)
            g2, defects = damage(g, r0, sites, rng)
            r1 = RB.describe(g2, {**RB.P, **(prm or {})})
            e1 = [e["pos"] for e in r1["ends"]]
            fl = {k: [p for f in r1["flags"] if f["kind"] == k for p in f["pts"]] for k in ("thin", "bulge")}
            for d in defects:
                if d["kind"] == "gap":
                    hits = sum(near(e1, np.array(q)) for q in d["edges"])
                    per_kind["gap"][0] += 1
                    per_kind["gap"][1] += hits == 2
                    per_kind["gap"][2] += hits == 1
                    d["found"] = {2: "both ends", 1: "one end", 0: "missed"}[hits]
                else:
                    hit = any(near(fl[d["kind"]], q) for q in d["stretch"][::3])
                    per_kind[d["kind"]][0] += 1
                    per_kind[d["kind"]][1] += hit
                    d["found"] = "yes" if hit else "missed"
            def at_defect(q):
                return any(near(d["stretch"], q, 2 * NEAR) for d in defects)
            new_e = [q for q in e1 if not near(e0, q) and not at_defect(q)]
            new_f = [q for f in r1["flags"] for q in [f["pts"][len(f["pts"]) // 2]] if not near(f0, q) and not at_defect(q)]
            fa_ends += len(new_e)
            fa_flags += len(new_f)
            rows.append(dict(section=name, seed=seed, defects=[{k: v for k, v in d.items() if k != "stretch"} for d in defects],
                             false_ends=len(new_e), false_flags=len(new_f)))
            if not quiet: print(f"{name:12s} seed {seed}: " + ", ".join(f"{d['kind']} {d['length_um']}um {d['found']}" for d in defects)
                  + f" | false ends {len(new_e)}, false flags {len(new_f)}", flush=True)
    n_trials = len(rows)
    summ = dict(gaps=dict(n=per_kind["gap"][0], both_ends=per_kind["gap"][1], one_end=per_kind["gap"][2]),
                thinning=dict(n=per_kind["thin"][0], found=per_kind["thin"][1]),
                bulges=dict(n=per_kind["bulge"][0], found=per_kind["bulge"][1]),
                false_ends_per_section=round(fa_ends / n_trials, 2), false_flags_per_section=round(fa_flags / n_trials, 2),
                clean_sections=clean_stats)
    if not quiet:
        print(json.dumps(summ, indent=1, default=int))
    (OUT / "synthetic_eval.json").write_text(json.dumps(dict(params={**RB.P, **(prm or {})}, summary=summ, trials=rows),
                                                        indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else float(o)))
    return summ


if __name__ == "__main__":
    main()
