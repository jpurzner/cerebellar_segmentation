"""
Hole-size sweep for the P7 template fit (run fit_p7_template.py first).

    python scripts/p7_template_holdout_sweep.py [slide_id]

Disks of the automatic labels are hidden (treated like anomalies), the template is
refitted without them, and each predictor is scored inside the disks against the
gold standard (pixels the gold calls 'cereb, no layer' are skipped):

  nearest_fill        label of the nearest visible automatic-label pixel
  template            similarity + 100 um B-spline template
  template+demons     + dense demons warp (8 um smoothing); hidden pixels exert no force
  hybrid              tissue vs background from template+demons; inside tissue the
                      layer of the nearest visible tissue pixel
-> outputs/p7_template_fit/<slide>/<slide>_holdout_sweep.json
"""
import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import SimpleITK as sitk
from scipy import ndimage as ndi

sys.path.insert(0, str(Path(__file__).parent))
import fit_p7_template as F  # noqa: E402

RADII = {100.0: 10, 200.0: 6, 300.0: 4, 400.0: 3}
SEEDS = (0, 1, 2)
_S = {}


def demons(fixed_np, moving_np, hidden, iters=60, sigma_px=2.0):
    f = np.where(hidden, moving_np, fixed_np)          # no force inside hidden disks
    fi, mi = F.to_img(f), F.to_img(moving_np)
    d = sitk.FastSymmetricForcesDemonsRegistrationFilter()
    d.SetNumberOfIterations(iters)
    d.SetStandardDeviations(sigma_px)
    d.SetSmoothDisplacementField(True)
    field = d.Execute(fi, mi)
    return sitk.DisplacementFieldTransform(sitk.Cast(field, sitk.sitkVectorFloat64))


def setup(slide):
    rep = json.loads((F.OUT_ROOT / slide / f"{slide}_fit.json").read_text())
    auto, gold = F.load_target(slide, rep["fit_input"][len(slide):])
    tgt = F.fractions(auto, F.BLOCK)
    _S.update(rep=rep, tgt=tgt, v=1.0 - tgt[F.CH.index(F.ANOM)], auto_g=F.grid_labels(tgt),
              gold_g=F.grid_labels(F.fractions(gold, F.BLOCK)))


def one(args):
    slide, r_um, seed = args
    if not _S:
        setup(slide)
    sitk.ProcessObject.SetGlobalDefaultNumberOfThreads(2)
    rep, tgt, v, auto_g, gold_g = _S["rep"], _S["tgt"], _S["v"], _S["auto_g"], _S["gold_g"]
    best = rep["best"]
    hid, cen = F.holdout_disks(gold_g, v > 0.5, seed, n=RADII[r_um], r_um=r_um)
    vh = v * (~hid)
    fh = F.fit_similarity(best["sec"], 0, best["mirror"], levels=(1,), tgt=tgt, v=vh, x0=best["params"])
    G = tgt.shape[1:]
    src = F.tpl_grid(best["sec"])
    wsh = F.warp(src, G, F.GT_PX, F.GT_ORG, np.array(fh["params"]), fh["mirror"],
                 np.array(fh["ct"]), np.array(fh["cs"]))
    mask = (v > 0.5) & (tgt[F.CH.index(F.UNK)] < 0.5)
    ref = F.to_img(np.zeros(G, np.float32))
    tx, _ = F.fit_bspline(tgt, wsh, mask & ~hid, 100.0)
    wb = F.apply_tx(tx, wsh, ref)
    fixed_np = F.encode(np.stack([tgt[i] for i in (1, 2, 3, 4)]) /
                        np.maximum(tgt[0:5].sum(0), 1e-6)[None])
    dtx = demons(fixed_np, F.encode(wb), hid | ~mask)
    wd = F.apply_tx(dtx, wb, ref)
    known = ~hid & (auto_g != F.ANOM)
    idx = ndi.distance_transform_edt(~known, return_distances=False, return_indices=True)
    nf = auto_g[tuple(idx)]
    tis_known = known & (auto_g > 0) & (auto_g != F.UNK)
    idx2 = ndi.distance_transform_edt(~tis_known, return_distances=False, return_indices=True)
    tdem = F.argmax_labels(wd)
    hyb = np.where(tdem > 0, auto_g[tuple(idx2)], 0)
    preds = {"nearest_fill": nf, "template": F.argmax_labels(wb), "template+demons": tdem,
             "hybrid": hyb}
    sel = hid & (gold_g != F.UNK)
    # disks reaching the pia vs fully inside the tissue
    edge = np.zeros_like(hid)
    rr = r_um / F.GT_PX
    yy, xx = np.indices(G)
    surf = []
    for cy, cx in cen:
        dk = (yy - cy) ** 2 + (xx - cx) ** 2 <= rr * rr
        surf.append(bool((gold_g[dk] == 0).any()))
        if surf[-1]:
            edge |= dk
    out = dict(radius_um=r_um, seed=seed, n_disks=len(cen), n_surface=int(sum(surf)), n_px=int(sel.sum()))
    for k, p in preds.items():
        ok = p[sel] == gold_g[sel]
        out[k] = dict(acc=round(float(ok.mean()), 3),
                      acc_surface=round(float((p[sel & edge] == gold_g[sel & edge]).mean()), 3) if (sel & edge).any() else None,
                      acc_interior=round(float((p[sel & ~edge] == gold_g[sel & ~edge]).mean()), 3) if (sel & ~edge).any() else None)
    # observed-region Dice of the demons-refined template (no hold-out influence)
    out["template+demons_dice_visible"] = F.dice_table(np.where(hid, 0, tdem), np.where(hid, 0, gold_g),
                                                       hid | (v <= 0.5))
    print(f"r {r_um:.0f} seed {seed}: " + ", ".join(f"{k} {out[k]['acc']:.3f}" for k in preds), flush=True)
    return out


def main(slide="2018_05_22_s1_2_p27-0021"):
    jobs = [(slide, r, s) for r in RADII for s in SEEDS]
    with Pool(4) as pool:
        rows = pool.map(one, jobs, chunksize=1)
    keys = ["nearest_fill", "template", "template+demons", "hybrid"]
    summary = {}
    for r in RADII:
        rs = [x for x in rows if x["radius_um"] == r]
        summary[f"{r:.0f}"] = {k: round(float(np.mean([x[k]["acc"] for x in rs])), 3) for k in keys}
    print("mean accuracy inside hidden disks, by radius (um):")
    for r, d in summary.items():
        print(f"  r {r:>3s}: " + "  ".join(f"{k} {v:.3f}" for k, v in d.items()))
    p = F.OUT_ROOT / slide / f"{slide}_holdout_sweep.json"
    p.write_text(json.dumps(dict(summary=summary, runs=rows), indent=1))
    print(p)


if __name__ == "__main__":
    main(*sys.argv[1:2])
