"""
Which generative parameters make a sensible IGL? Approximate Bayesian computation over the
synthetic samples (igl_synth_gen.py) against the real reference sets (igl_synth_stats.py).

    python scripts/igl_synth_abc.py       # -> outputs/igl_synth/abc.json

Targets: "IF" = the 12 IF controls (5 intact + 7 others), "Allen" = Allen P7 within 600 um.
Distance of a section to a target: root-mean-square over the statistics of
  z_j = (s_j - median_j) / scale_j,   scale_j = max(1.4826 * MAD_j, floor_j)
computed on the target's own sections. A synthetic section is "sensible" for a target when
its distance is no larger than the 90th percentile of the target's own sections' distances
(each real section measured against the other real sections of its group): the synthetic
IGL is then as typical as the real ones are, by every statistic at once.
For each family and target: share of samples that are sensible, parameter ranges of the
sensible ones (posterior), sensible rate as a function of each parameter, and which
statistics fail most often among the samples that are not sensible.
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import igl_synth_stats as S  # noqa: E402

OUT = S.OUT
USE = [n for n in S.NAMES]
FLOOR = dict(igl_frac=0.02, egl_frac=0.02, ml_frac=0.02, wm_frac=0.02, solidity=0.02, igl_len=0.5, w_med=5, w_iqr=5, d_med=6,
             d_iqr=5, ml_med=5, egl_med=4, egl_contact=0.01, wm_contact=0.04, kap_med=0.5, convex=0.03, concave=0.03,
             ends_per_mm=0.15, pieces=1.0, corr_w=8)
TARGETS = {"IF": ("IF control", "IF other"), "Allen": ("Allen P7",), "cKO": ("IF EZH2 cKO",)}


def load_real():
    return json.loads((OUT / "real_stats.json").read_text())


def target_model(rows):
    X = np.array([[r[k] for k in USE] for r in rows], float)
    med = np.nanmedian(X, axis=0)
    mad = 1.4826 * np.nanmedian(np.abs(X - med), axis=0)
    sc = np.maximum(mad, [FLOOR[k] for k in USE])
    return med, sc


def dist(X, med, sc):
    Z = (X - med) / sc
    Z = np.where(np.isnan(Z), 3.0, Z)                 # a missing statistic counts as far off
    return np.sqrt((Z ** 2).mean(1)), Z


def load_samples(fam):
    p = OUT / f"samples_{fam}.jsonl"
    if not p.exists():
        return []
    return [json.loads(l) for l in p.read_text().splitlines() if l.strip()]


def main():
    real = load_real()
    res = dict(stats=USE, targets={}, families={})
    tmods = {}
    for tname, groups in TARGETS.items():
        R = [r for r in real if r["group"] in groups]
        med, sc = target_model(R)
        # each real section against the others of its target (leave one out)
        dr = []
        for i in range(len(R)):
            m_i, s_i = target_model([r for j, r in enumerate(R) if j != i])
            d, _ = dist(np.array([[R[i][k] for k in USE]], float), m_i, s_i)
            dr.append(float(d[0]))
        thr = float(np.percentile(dr, 90))
        tmods[tname] = (med, sc, thr)
        res["targets"][tname] = dict(n=len(R), median=dict(zip(USE, med.round(4).tolist())), scale=dict(zip(USE, sc.round(4).tolist())),
                                     real_distances=sorted(round(v, 2) for v in dr), threshold=round(thr, 2))
        # how far the other real groups are from this target
        for grp in ("IF control", "IF other", "IF EZH2 cKO", "Allen P7"):
            Q = [r for r in real if r["group"] == grp]
            d, _ = dist(np.array([[r[k] for k in USE] for r in Q], float), med, sc)
            res["targets"][tname][f"distance_of_{grp}"] = dict(median=round(float(np.median(d)), 2),
                                                               sensible=f"{int((d <= thr).sum())}/{len(d)}")
        print(tname, "threshold", round(thr, 2), {g: res["targets"][tname][f"distance_of_{g}"] for g in ("IF control", "IF other", "IF EZH2 cKO", "Allen P7")})
    for fam in ("fold", "fold2", "fold3", "foldma", "finger", "arbor", "fourier"):
        smp = [s for s in load_samples(fam)]
        ok = [s for s in smp if s["stats"] is not None]
        if not ok:
            continue
        pnames = [k for k in ok[0]["params"] if isinstance(ok[0]["params"][k], (int, float))]
        Pm = np.array([[s["params"][k] for k in pnames] for s in ok], float)
        X = np.array([[s["stats"][k] for k in USE] for s in ok], float)
        fres = dict(n=len(smp), valid=len(ok), params=pnames, targets={})
        for tname, (med, sc, thr) in tmods.items():
            d, Z = dist(X, med, sc)
            sens = d <= thr
            order = np.argsort(d)
            # parameter posterior (sensible samples, or the best 1 % if none are sensible)
            acc = sens if sens.sum() >= 20 else np.isin(np.arange(len(d)), order[:max(20, len(d) // 100)])
            post = {}
            for j, k in enumerate(pnames):
                lo, hi = float(Pm[:, j].min()), float(Pm[:, j].max())
                edges = np.linspace(lo, hi, 11)
                rate = []
                for a, b in zip(edges[:-1], edges[1:]):
                    inb = (Pm[:, j] >= a) & (Pm[:, j] <= b)
                    rate.append(round(float(acc[inb].mean()) if inb.sum() else 0.0, 4))
                q = np.percentile(Pm[acc, j], [5, 25, 50, 75, 95])
                post[k] = dict(prior=[round(lo, 3), round(hi, 3)], q=q.round(3).tolist(), rate=rate, edges=edges.round(3).tolist())
            # which statistics keep the rest from being sensible
            bad = ~acc
            fail = {k: round(float((np.abs(Z[bad, j]) > 2).mean()), 3) for j, k in enumerate(USE)}
            # posterior predictive: statistics of the accepted samples
            pp = {k: np.percentile(X[acc, j][~np.isnan(X[acc, j])], [10, 50, 90]).round(3).tolist() if np.isfinite(X[acc, j]).any() else None
                  for j, k in enumerate(USE)}
            # which parameters matter: random forest predicting log distance, permutation importance
            from sklearn.ensemble import RandomForestRegressor
            from sklearn.inspection import permutation_importance
            yv = np.log(d)
            idx = np.random.default_rng(0).permutation(len(yv))
            ntr = int(0.75 * len(yv))
            rf = RandomForestRegressor(n_estimators=150, min_samples_leaf=5, n_jobs=1, random_state=0).fit(Pm[idx[:ntr]], yv[idx[:ntr]])
            r2 = float(rf.score(Pm[idx[ntr:]], yv[idx[ntr:]]))
            pi = permutation_importance(rf, Pm[idx[ntr:]], yv[idx[ntr:]], n_repeats=5, random_state=0, n_jobs=1)
            imp = sorted(zip(pnames, pi.importances_mean), key=lambda t: -t[1])
            fres_imp = dict(r2=round(r2, 3), ranking=[[k, round(float(v), 4)] for k, v in imp])
            fres["targets"][tname] = dict(importance=fres_imp, sensible=int(sens.sum()), sensible_pct=round(100 * float(sens.mean()), 2),
                                          accepted=int(acc.sum()), best_distance=round(float(d.min()), 2),
                                          posterior=post, fails=dict(sorted(fail.items(), key=lambda t: -t[1])),
                                          predictive=pp, best_seeds=[int(ok[i]["seed"]) for i in order[:40]],
                                          worst_seeds=[int(ok[i]["seed"]) for i in order[-10:]],
                                          distances=dict(zip([int(s["seed"]) for s in ok], d.round(3).tolist())))
            print(f"{fam:8s} vs {tname:6s}: {sens.sum()}/{len(d)} sensible ({100 * sens.mean():.2f}%), best distance {d.min():.2f}; "
                  f"importance (R2 {r2:.2f}): {[k for k, _ in imp[:5]]}; "
                  f"most failed: {list(fres['targets'][tname]['fails'].items())[:4]}", flush=True)
        res["families"][fam] = fres
    (OUT / "abc.json").write_text(json.dumps(res, default=float))


if __name__ == "__main__":
    main()
