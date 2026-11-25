"""
Benchmark of IGL ribbon completion models on synthetic gaps in intact ribbons of the
good-IGL controls (s1_2 gold, s5_F, s5_G, s6_C, s6_G).

    python scripts/igl_model_eval.py

Per section 36 gaps of 50-500 um, each cut alone into a skeleton branch with >= 120 um of
ribbon kept on both sides. From the flanks only (smoothed, ignoring the gap) we take each
end's position, direction (last 60 um), width and depth below the EGL (median 40-200 um
behind the end). Models predict the missing centre line:
  line, hermite, elastic (igl_model.elastic, several weightings), elastic_nodepth (no
  depth term), walk (igl_repair.walk, the current method; tried from each end, IGL in
  the gap removed).
Scores per gap:
  err_um   - symmetric mean distance between predicted and true centre line
  max_um   - largest distance (Hausdorff)
  dice     - overlap of the predicted band (centre line + flank width, Hermite-interpolated)
             with the removed IGL
  width    - mean abs error of the interpolated width against the true width
Elastic weightings are chosen on the even-numbered gaps and scored on the odd ones.
-> outputs/igl_model/benchmark.json
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).parent))
import igl_model as M  # noqa: E402
import igl_repair as IR  # noqa: E402
import igl_ribbon as RB  # noqa: E402
import igl_ribbon_eval as RE  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/igl_model"
G = RB.GRID
N_GAPS = 36
MARGIN_UM = 120.0
DIR_UM = 60.0
FLANK = (40.0, 200.0)
BINS = [(50, 150), (150, 250), (250, 350), (350, 500)]
VARIANTS = {f"elastic s{sd} b{bs}": dict(sigma_depth_um=sd, bend_scale=bs)
            for sd in (8.0, 15.0, 30.0) for bs in (0.3, 1.0, 3.0)}
VARIANTS["elastic_nodepth"] = dict(sigma_depth_um=np.inf, bend_scale=1.0)
FULL = {f"model s{sd} b{bs}": dict(sigma_depth_um=sd, bend_scale=bs) for sd in (10.0, 15.0, 25.0) for bs in (0.3, 1.0)}


def arclen(P):
    return np.concatenate([[0], np.cumsum(np.hypot(*np.diff(P, axis=0).T))]) * G


def resample(C, step=0.5):
    C = np.asarray(C, float)
    if len(C) < 2:
        return C
    s = arclen(C) / G
    n = max(int(s[-1] / step), 2)
    t = np.linspace(0, s[-1], n)
    return np.stack([np.interp(t, s, C[:, k]) for k in range(2)], 1)


def curve_err(C, T):
    C, T = resample(C), resample(T)
    d1, _ = cKDTree(T).query(C)
    d2, _ = cKDTree(C).query(T)
    return float((d1.mean() + d2.mean()) / 2 * G), float(max(d1.max(), d2.max()) * G)


def band(C, wa, wb, shape, keep):
    """tube of Hermite-interpolated width around the centre line, inside keep."""
    C = resample(C)
    w = M.hermite_width(wa, wb, len(C))
    r = int(np.ceil(max(wa, wb) / 2 / G)) + 2
    y0, x0 = np.maximum(np.floor(C.min(0)).astype(int) - r, 0)
    y1, x1 = np.minimum(np.ceil(C.max(0)).astype(int) + r + 1, shape)
    yy, xx = np.mgrid[y0:y1, x0:x1]
    d, k = cKDTree(C).query(np.stack([yy.ravel(), xx.ravel()], 1))
    m = np.zeros(shape, bool)
    sel = d * G <= w[k] / 2
    m[yy.ravel()[sel], xx.ravel()[sel]] = True
    return m & keep


def dice(a, b):
    return float(2 * (a & b).sum() / max(a.sum() + b.sum(), 1))


def flank(P, from_end):
    """smoothed flank, its end point, direction pointing into the gap, and the raw points
    40-200 um behind the end."""
    Q = P if from_end else P[::-1]          # flank ordered so that the gap follows Q[-1]
    S = np.stack([ndi.gaussian_filter1d(Q[:, k].astype(float), 2.0, mode="nearest") for k in range(2)], 1)
    s = arclen(S)
    back = s[-1] - s
    k = int(np.argmin(np.abs(back - DIR_UM)))
    v = S[-1] - S[k]
    v = v / max(np.hypot(*v), 1e-9)
    sel = (back >= FLANK[0]) & (back <= FLANK[1])
    return S[-1], v, Q[sel]


def run_walk(start, v, depth, sec, igl_d):
    status, pts = IR.walk(start, v, max(depth, 8.0), sec["dEs"], sec["gy"], sec["gx"], sec["ok"], igl_d,
                          sec["g"]["egl"], sec["g"]["wm"], IR.DEF)
    if status != "joined":
        return None
    P = [tuple(start)] + [p for p in pts if p and p[0] != "v"]
    return np.array(P, float)


def section_data(name, lab, px):
    g = RB.grid_labels(lab, px)
    rib = RB.describe(g)
    dE = rib["depth_egl"]
    dEs = ndi.gaussian_filter(dE, 1.0)
    gy, gx = np.gradient(dEs)
    forbid = g["egl"] | g["wm"] | ~(g["tissue"] | g["anom"])
    sec = dict(name=name, g=g, rib=rib, dE=dE, dEs=dEs, gy=gy, gx=gx, ok=g["tissue"] | g["anom"],
               field=M.Field(dE, forbid), keep=(g["tissue"] | g["anom"]) & ~g["egl"] & ~g["wm"])
    cl = np.vstack([b["path"] for b in rib["branches"]])
    owner = np.concatenate([[(bi, j) for j in range(len(b["path"]))] for bi, b in enumerate(rib["branches"])])
    yy, xx = np.nonzero(rib["igl"])
    _, nn = cKDTree(cl).query(np.stack([yy, xx], 1))
    sec.update(yy=yy, xx=xx, ob=owner[nn, 0], oj=owner[nn, 1])
    return sec


def make_gaps(sec, rng):
    br = [(bi, b) for bi, b in enumerate(sec["rib"]["branches"])]
    gaps, tries = [], 0
    while len(gaps) < N_GAPS and tries < 5000:
        tries += 1
        L = rng.uniform(50, 500)
        cands = [(bi, b) for bi, b in br if arclen(b["path"])[-1] >= L + 2 * MARGIN_UM]
        if not cands:
            continue
        lens = np.array([arclen(b["path"])[-1] for _, b in cands])
        bi, b = cands[rng.choice(len(cands), p=lens / lens.sum())]
        s = arclen(b["path"])
        s0 = rng.uniform(MARGIN_UM, s[-1] - MARGIN_UM - L)
        j0 = int(np.searchsorted(s, s0))
        j1 = int(np.searchsorted(s, s0 + L, side="right")) - 1
        if j1 - j0 < 3:
            continue
        gaps.append(dict(branch=bi, j0=j0, j1=j1, length_um=float(s[j1] - s[j0])))
    return gaps


def evaluate_gap(sec, gp, variants):
    b = sec["rib"]["branches"][gp["branch"]]
    P = b["path"]
    j0, j1 = gp["j0"], gp["j1"]
    a, ta, fa = flank(P[:j0], True)
    bb, tb_out, fb = flank(P[j1 + 1:], False)     # direction pointing back into the gap
    tb_in = -tb_out
    W, D = sec["rib"]["width"], sec["dE"]
    wa, da = float(np.median(W[fa[:, 0], fa[:, 1]])), float(np.median(D[fa[:, 0], fa[:, 1]]))
    wb, db = float(np.median(W[fb[:, 0], fb[:, 1]])), float(np.median(D[fb[:, 0], fb[:, 1]]))
    Ps = np.stack([ndi.gaussian_filter1d(P[:, k].astype(float), 2.0, mode="nearest") for k in range(2)], 1)
    T = np.vstack([a, Ps[j0:j1 + 1], bb])
    gm = np.zeros(sec["rib"]["igl"].shape, bool)
    sel = (sec["ob"] == gp["branch"]) & (sec["oj"] >= j0) & (sec["oj"] <= j1)
    gm[sec["yy"][sel], sec["xx"][sel]] = True
    igl_d = sec["rib"]["igl"] & ~gm
    keep = sec["keep"] & ~igl_d
    tw = W[P[j0:j1 + 1, 0], P[j0:j1 + 1, 1]]
    tdep = D[np.round(T[:, 0]).astype(int), np.round(T[:, 1]).astype(int)]
    row = dict(section=sec["name"], length_um=round(gp["length_um"]), chord_um=round(float(np.hypot(*(bb - a)) * G)),
               turn_deg=round(float(np.degrees(np.arccos(np.clip(ta @ tb_in, -1, 1))))),
               width_ends=[round(wa), round(wb)], depth_ends=[round(da), round(db)],
               width_mae=round(float(np.abs(M.hermite_width(wa, wb, len(tw)) - tw).mean()), 1),
               width_mae_const=round(float(np.abs((wa + wb) / 2 - tw).mean()), 1),
               depth_dev_from_linear=round(float(np.abs(tdep - np.linspace(da, db, len(tdep))).mean()), 1),
               models={})
    preds = {"line": M.line(a, bb, 40), "hermite": M.hermite(a, ta, bb, tb_in, 40)}
    for k, v in variants.items():
        preds[k] = M.elastic(a, ta, bb, tb_in, sec["field"], da, db, v)
    avoid = igl_d & ~(ndi.binary_dilation(gm, iterations=1))
    preds["geodesic"] = M.geodesic(a, bb, sec["field"], (da + db) / 2, dict(sigma_depth_um=15.0), avoid)
    starts = {}
    for k, v in FULL.items():
        preds[k], info = M.complete(a, ta, bb, tb_in, sec["field"], da, db, v, avoid)
        starts[k] = info["start"]
    wk = run_walk(a, ta, da, sec, igl_d)
    if wk is None:
        wk2 = run_walk(bb, tb_out, db, sec, igl_d)
        wk = None if wk2 is None else wk2[::-1]
    preds["walk"] = wk
    for k, C in preds.items():
        if C is None or len(C) < 2:
            row["models"][k] = dict(err_um=None, max_um=None, dice=0.0)
            continue
        e, mx = curve_err(C, T)
        row["models"][k] = dict(err_um=round(e, 1), max_um=round(mx, 1), dice=round(dice(band(C, wa, wb, gm.shape, keep), gm), 3))
    row["start"] = starts
    ex = dict(T=T, preds={k: v for k, v in preds.items() if v is not None}, a=a, b=bb, gap=gm, ta=ta, tb_in=tb_in,
              wa=wa, wb=wb, da=da, db=db, keep=keep, avoid=avoid, igl_d=igl_d)
    return row, ex


def summarise(rows, models):
    out = {}
    for k in models:
        per = {}
        for lo, hi in BINS + [(0, 10000)]:
            rr = [r["models"][k] for r in rows if lo <= r["length_um"] < hi]
            if not rr:
                continue
            err = np.array([np.inf if r["err_um"] is None else r["err_um"] for r in rr])
            per["all" if hi == 10000 else f"{lo}-{hi}"] = dict(
                n=len(rr), median_err_um=round(float(np.median(err)), 1),
                within_10um=round(float((err <= 10).mean()), 2), within_20um=round(float((err <= 20).mean()), 2),
                failed=int(np.isinf(err).sum()), mean_dice=round(float(np.mean([r["dice"] for r in rr])), 3))
        out[k] = per
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows, examples = [], []
    t0 = time.time()
    for si, (name, lab, px) in enumerate(RE.sections()):
        sec = section_data(name, lab, px)
        gaps = make_gaps(sec, np.random.default_rng(200 + si))
        for gi, gp in enumerate(gaps):
            row, ex = evaluate_gap(sec, gp, VARIANTS)
            row["split"] = "tune" if gi % 2 == 0 else "test"
            rows.append(row)
            if gi < 6:
                examples.append(dict(section=name, **{k: v for k, v in ex.items() if k != "gap"}))
        print(f"{name}: {len(gaps)} gaps, {time.time() - t0:.0f}s", flush=True)
    models = list(rows[0]["models"])
    tune = [r for r in rows if r["split"] == "tune"]
    test = [r for r in rows if r["split"] == "test"]
    st = summarise(tune, models)
    best = min(VARIANTS, key=lambda k: st[k]["all"]["median_err_um"] if k != "elastic_nodepth" else np.inf)
    bestf = min(FULL, key=lambda k: (st[k]["all"]["median_err_um"], -st[k]["all"]["within_20um"]))
    report = ["line", "hermite", "elastic_nodepth", best, "geodesic", bestf, "walk"]
    res = dict(best_elastic=best, best_model=bestf, best_params={**M.DEF, **FULL[bestf]},
               start_used={s_: sum(r["start"][bestf] == s_ for r in test) for s_ in ("hermite", "geodesic")}, tune_all={k: st[k]["all"] for k in models},
               test=summarise(test, report),
               width=dict(hermite_mae=round(float(np.mean([r["width_mae"] for r in test])), 1),
                          const_mae=round(float(np.mean([r["width_mae_const"] for r in test])), 1)),
               depth_dev_from_linear_um=dict(median=round(float(np.median([r["depth_dev_from_linear"] for r in rows])), 1),
                                             p90=round(float(np.percentile([r["depth_dev_from_linear"] for r in rows], 90)), 1)))
    print(json.dumps(res, indent=1, default=float))
    (OUT / "benchmark.json").write_text(json.dumps(dict(summary=res, gaps=rows), indent=1,
                                                   default=lambda o: o.tolist() if hasattr(o, "tolist") else float(o)))
    np.save(OUT / "benchmark_examples.npy", np.array(examples, dtype=object), allow_pickle=True)
    return res


if __name__ == "__main__":
    main()
