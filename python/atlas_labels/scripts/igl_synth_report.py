"""
Report page for the synthetic-IGL study: outputs/igl_synth/igl_synth_report.html

    python scripts/igl_synth_report.py

Collects igl_synth_stats (real statistics), igl_synth_laminar (smoothed approximations),
igl_synth_abc (which parameters make a sensible IGL) and the refinement samples, and draws
real and synthetic sections as SVG paths (sections from the user's slides are turned to the
atlas frame: anterior left, dorsal up).
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi

sys.path.insert(0, str(Path(__file__).parent))
import igl_check_page as P  # noqa: E402
import igl_synth_abc as ABC  # noqa: E402
import igl_synth_laminar as L  # noqa: E402
import igl_synth_render as R  # noqa: E402
import igl_synth_stats as S  # noqa: E402
import wm_library as WL  # noqa: E402

OUT = S.OUT
FAMS = ("fold", "fold2", "fold3", "foldma", "finger", "arbor", "fourier")
TGTS = ("IF", "Allen", "cKO")
LGROUPS = ("IF control", "IF other", "Allen P7", "IF EZH2 cKO")
KEY_STATS = ["igl_frac", "solidity", "igl_len", "w_med", "d_med", "ml_med", "egl_med", "kap_med", "ends_per_mm"]


def orient_if(name, g):
    lib = {e["id"]: e for e in json.loads((WL.OUT / "library.json").read_text())["entries"]}
    lid = "if_2018_05_22_s1_2_p27-0021" if name.startswith("s1_2") else f"if_{name}_p27"
    o = lib[lid]["orientation"]
    A = P.rot_matrix(o["rot_deg"], o["mirror"])
    H, W = g["tissue"].shape
    n = int(np.ceil(np.hypot(H, W))) + 4
    out = {k: P.rotate(g[k].astype(np.uint8), A, n, 0)[0] > 0 for k in ("tissue", "egl", "ml", "igl", "wm")}
    ys, xs = np.nonzero(out["tissue"])
    sl = (slice(ys.min() - 4, ys.max() + 5), slice(xs.min() - 4, xs.max() + 5))
    z = {k: v[sl] for k, v in out.items()}
    z["anom"] = np.zeros_like(z["tissue"])
    z["unass"] = np.zeros_like(z["tissue"])
    return z


def real_grid(name):
    if name.startswith("Allen"):
        for n, p, off in S.allen_sections():
            if n == name:
                return S.allen_grid(p)
    return orient_if(name, S.if_grid(name))


def load_jsonl(p):
    return [json.loads(l) for l in p.read_text().splitlines() if l.strip()] if p.exists() else []


def stat_line(s):
    return {k: (None if s.get(k) is None or not np.isfinite(s[k]) else round(float(s[k]), 3)) for k in KEY_STATS}


def main():
    real = json.loads((OUT / "real_stats.json").read_text())
    lam = json.loads((OUT / "laminar.json").read_text())
    abc = json.loads((OUT / "abc.json").read_text())
    data = dict(stat_names=dict(S.STATS), stats=S.NAMES, key_stats=KEY_STATS)
    # 1. real statistics by group
    groups = ["IF control", "IF other", "Allen P7", "IF EZH2 cKO"]
    data["real"] = {g: dict(n=sum(r["group"] == g for r in real),
                            q={k: np.nanpercentile([r[k] for r in real if r["group"] == g], [10, 50, 90]).round(3).tolist() for k in S.NAMES})
                    for g in groups}
    # 2. smoothed approximations
    models = ["wm_const", "band_const"] + [f"wm_smooth_{s}" for s in L.SCALES] + [f"band_smooth_{s}" for s in L.SCALES]
    data["laminar"] = dict(models=models, scales=L.SCALES,
                           dice={g: {m: np.percentile([r[m]["dice"] for r in lam if r["group"] == g], [10, 50, 90]).round(3).tolist() for m in models}
                                 for g in LGROUPS if any(r["group"] == g for r in lam)},
                           m={g: [r["wm_const"]["m"] for r in lam if r["group"] == g] for g in LGROUPS},
                           band={g: [[r["band_const"]["a"], r["band_const"]["b"]] for r in lam if r["group"] == g] for g in LGROUPS})
    ex = {}
    for name in ("s5_G", "Allen P7 +0"):
        g = real_grid(name)
        res, _ = L.fit(g)
        e = L.example_masks(g, res)
        e["dice"] = {k: v["dice"] for k, v in res.items()}
        e["m"] = res["wm_const"]["m"]
        e["band"] = [res["band_const"]["a"], res["band_const"]["b"]]
        ex[name] = e
    data["laminar"]["examples"] = ex
    print("laminar examples done", flush=True)
    al = OUT / "allometry.json"
    data["allometry"] = json.loads(al.read_text()) if al.exists() else None
    # 3. families: sensible shares, posteriors, failures, galleries
    data["targets"] = {t: dict(n=v["n"], threshold=v["threshold"], real_distances=v["real_distances"],
                               others={g: v.get(f"distance_of_{g}") for g in groups}) for t, v in abc["targets"].items()}
    data["families"] = {}
    tmods = {}
    for tname, grp in ABC.TARGETS.items():
        Rt = [r for r in real if r["group"] in grp]
        med, sc = ABC.target_model(Rt)
        tmods[tname] = (med, sc, abc["targets"][tname]["threshold"])
    for fam in FAMS:
        if fam not in abc["families"]:
            continue
        fa = abc["families"][fam]
        smp = {r["seed"]: r for r in load_jsonl(OUT / f"samples_{fam}.jsonl") if r["stats"]}
        fd = dict(n=fa["n"], valid=fa["valid"], params=fa["params"], targets={}, gallery=[])
        for t, ft in fa["targets"].items():
            fd["targets"][t] = dict(importance=ft.get("importance"), sensible=ft["sensible"], pct=ft["sensible_pct"], best=ft["best_distance"], accepted=ft["accepted"],
                                    fails=dict(list(ft["fails"].items())[:6]),
                                    posterior={k: dict(prior=v["prior"], q=v["q"], rate=v["rate"]) for k, v in ft["posterior"].items()})
        # gallery: best for IF, best for Allen, and failures with different reasons
        dIF = {int(k): v for k, v in fa["targets"]["IF"]["distances"].items()}
        dAl = {int(k): v for k, v in fa["targets"]["Allen"]["distances"].items()}
        picks = [("best for IF", s_) for s_ in sorted(dIF, key=dIF.get)[:3]] + [("best for Allen", s_) for s_ in sorted(dAl, key=dAl.get)[:2]]
        med, sc, thr = tmods["IF"]
        far = [s_ for s_ in dIF if dIF[s_] > 2 * thr]
        rng = np.random.default_rng(5)
        for s_ in rng.choice(far, size=min(4, len(far)), replace=False) if far else []:
            picks.append(("not sensible", int(s_)))
        for why, s_ in picks:
            g, p = R.regenerate(fam, s_)
            st = smp[s_]["stats"]
            z = (np.array([st[k] for k in ABC.USE], float) - med) / sc
            worst = [ABC.USE[i] for i in np.argsort(-np.abs(np.nan_to_num(z, nan=3)))[:3]]
            fd["gallery"].append(dict(why=why, seed=int(s_), d_IF=round(dIF[s_], 2), d_Allen=round(dAl[s_], 2), svg=R.svg_data(g),
                                      stats=stat_line(st), worst=[[k, round(float(z[ABC.USE.index(k)]), 1)] for k in worst],
                                      params={k: round(v, 2) if isinstance(v, float) else v for k, v in p.items()}))
        data["families"][fam] = fd
        print(fam, "gallery done", flush=True)
    # 4. refinement rounds (local search around the closest parameter sets), per family and target
    data["refine"] = []
    refined = {}                                   # (fam, target) -> (params matrix of sensible sets, keys)
    for fam in FAMS:
        for t in TGTS:
            files = sorted(OUT.glob(f"samples_{fam}_refine*_{t}.jsonl"), key=lambda f_: int(f_.stem.split("_")[-2][6:] or 1))
            if not files:
                continue
            med, sc, thr = tmods[t]
            rounds, rows = [], []
            for f_ in files:
                rr = [r for r in load_jsonl(f_) if r["stats"]]
                nall = len(load_jsonl(f_))
                if not rr:
                    continue
                Xr = np.array([[r["stats"][k] for k in ABC.USE] for r in rr], float)
                dr, _ = ABC.dist(Xr, med, sc)
                tag = f_.stem.split("_")[-2]
                rounds.append(dict(tag=tag, n=nall, valid=len(rr), best=round(float(dr.min()), 2), median=round(float(np.median(dr)), 2),
                                   sensible=int((dr <= thr).sum())))
                rows += rr
            if not rows:
                continue
            X = np.array([[r["stats"][k] for k in ABC.USE] for r in rows], float)
            d, Z = ABC.dist(X, med, sc)
            sens = d <= thr
            keys = [k for k in rows[0]["params"] if isinstance(rows[0]["params"][k], (int, float)) and not isinstance(rows[0]["params"][k], bool)]
            Pm = np.array([[r["params"][k] for k in keys] for r in rows], float)
            q = {k: np.percentile(Pm[sens, j], [5, 25, 50, 75, 95]).round(3).tolist() for j, k in enumerate(keys)} if sens.sum() >= 10 else {}
            if sens.sum() >= 10:
                refined[(fam, t)] = (Pm[sens], keys)
            if fam in data["families"] and t in data["families"][fam]["targets"]:
                data["families"][fam]["targets"][t]["refined"] = dict(n=len(rows), sensible=int(sens.sum()), q=q)
            order = np.argsort(d)
            gal = []
            for i in order[:3]:
                r = rows[i]
                g, _ = R.regenerate(fam, r["seed"], r["params"])
                gal.append(dict(seed=r["seed"], d=round(float(d[i]), 2), svg=R.svg_data(g), stats=stat_line(r["stats"])))
            still = Z[sens] if sens.sum() else Z[order[:30]]
            off = sorted(((k, float(np.median(np.abs(still[:, j])))) for j, k in enumerate(ABC.USE)), key=lambda t_: -t_[1])[:4]
            data["refine"].append(dict(target=t, family=fam, n=len(rows), sensible=int(sens.sum()), pct=round(100 * float(sens.mean()), 1), q=q,
                                       best=round(float(d.min()), 2), rounds=rounds, gallery=gal, still_off=[[k, round(v, 2)] for k, v in off],
                                       predictive={k: np.percentile(X[sens, j][np.isfinite(X[sens, j])], [10, 50, 90]).round(3).tolist()
                                                   if sens.sum() >= 10 else None for j, k in enumerate(ABC.USE)}))
            print("refine", fam, t, int(sens.sum()), "/", len(rows), flush=True)
    # parameter shifts between targets, for families refined against both
    data["contrasts"] = []
    for fam in FAMS:
        for a, b in (("IF", "cKO"), ("IF", "Allen")):
            if (fam, a) in refined and (fam, b) in refined:
                Pa, ka = refined[(fam, a)]
                Pb, kb = refined[(fam, b)]
                rows_ = []
                for j, k in enumerate(ka):
                    if k not in kb:
                        continue
                    va, vb = Pa[:, j], Pb[:, kb.index(k)]
                    ia = np.subtract(*np.percentile(va, [75, 25]))
                    ib = np.subtract(*np.percentile(vb, [75, 25]))
                    pooled = max((ia + ib) / 2, 1e-9)
                    rows_.append([k, round(float(np.median(va)), 3), round(float(ia), 3), round(float(np.median(vb)), 3), round(float(ib), 3),
                                  round(float((np.median(vb) - np.median(va)) / pooled), 2)])
                rows_.sort(key=lambda r_: -abs(r_[5]))
                data["contrasts"].append(dict(family=fam, a=a, b=b, na=int(len(Pa)), nb=int(len(Pb)), rows=rows_))
    # 5. synthetic twins: nearest synthetic section to each real one
    pool = []
    for fam in FAMS:
        pool += [r for r in load_jsonl(OUT / f"samples_{fam}.jsonl") if r["stats"]]
    for fam in FAMS:
        for t in ("IF", "Allen"):
            pool += [r for f_ in sorted(OUT.glob(f"samples_{fam}_refine*_{t}.jsonl")) for r in load_jsonl(f_) if r["stats"]]
    PX = np.array([[r["stats"][k] for k in ABC.USE] for r in pool], float)
    twins = []
    allen_mid = [n for n, p, off in S.allen_sections() if off in (0.0, -200.0, 400.0)]
    for name in S.IF_CONTROL + allen_mid:
        rr = [r for r in real if r["name"] == name][0]
        med, sc, thr = tmods["Allen" if name.startswith("Allen") else "IF"]
        x = np.array([rr[k] for k in ABC.USE], float)
        dd = np.sqrt(np.nanmean(np.where(np.isnan(PX), 9.0, ((PX - x) / sc) ** 2), axis=1))
        i = int(np.argmin(dd))
        r = pool[i]
        g, _ = R.regenerate(r["family"], r["seed"], r["params"] if r.get("parent") is not None else None)
        twins.append(dict(real=name, real_svg=R.svg_data(real_grid(name)), real_stats=stat_line(rr), family=r["family"],
                          refined=r.get("parent") is not None, d=round(float(dd[i]), 2), syn_svg=R.svg_data(g), syn_stats=stat_line(r["stats"])))
        print("twin", name, r["family"], round(float(dd[i]), 2), flush=True)
    data["twins"] = twins
    # 6. free ribbon (no surface): how often does a ribbon of real length cross itself?
    data["free_ribbon"] = free_ribbon()
    cj = OUT / "conclusions.json"
    data["conclusions"] = json.loads(cj.read_text()) if cj.exists() else []
    (OUT / "report_data.json").write_text(json.dumps(data, default=float))
    tpl = (Path(__file__).parent / "igl_synth_template.html").read_text()
    p = OUT / "igl_synth_report.html"
    p.write_text(tpl.replace("__DATA__", json.dumps(data, default=float)))
    print(p, f"{p.stat().st_size / 1e6:.2f} MB")


def _crosses(P):
    """does a polyline cross itself? (segment intersection test, skipping neighbours)"""
    from itertools import combinations
    n = len(P) - 1
    A, B = P[:-1], P[1:]
    for i in range(0, n, 1):
        a, b = A[i], B[i]
        c, d = A[i + 2:], B[i + 2:]
        if not len(c):
            continue
        def orient(p, q, r):
            return np.sign((q[..., 0] - p[..., 0]) * (r[..., 1] - p[..., 1]) - (q[..., 1] - p[..., 1]) * (r[..., 0] - p[..., 0]))
        o1 = orient(a, b, c)
        o2 = orient(a, b, d)
        o3 = orient(c, d, a[None])
        o4 = orient(c, d, b[None])
        if np.any((o1 != o2) & (o3 != o4)):
            return True
    return False


def free_ribbon(n=300, ds=8.0):
    rng = np.random.default_rng(11)
    out = {}
    for beta in (1.0, 0.3):
        sk = 0.012 * np.sqrt(beta)
        res = {}
        for L in (500, 2000, 13000):
            k = int(L / ds)
            hits = 0
            contact = 0
            ex = None
            for t in range(n if L < 13000 else 60):
                a = np.exp(-ds / 12.0)
                kap = np.empty(k)
                kap[0] = sk * rng.standard_normal()
                for i in range(1, k):
                    kap[i] = kap[i - 1] * a + sk * np.sqrt(1 - a * a) * rng.standard_normal()
                th = np.cumsum(kap * ds)
                Pp = np.cumsum(np.stack([np.cos(th), np.sin(th)], 1) * ds, 0)
                c = _crosses(Pp)
                hits += c
                # self-contact: two stretches > 150 um apart along the ribbon come closer than its width (57 um)
                from scipy.spatial import cKDTree
                pairs = cKDTree(Pp).query_pairs(57.0, output_type="ndarray")
                contact += bool(len(pairs) and (np.abs(pairs[:, 0] - pairs[:, 1]) * ds > 150).any())
                if ex is None and L == 2000:
                    ex = Pp[::2].round(1).tolist()
            res[str(L)] = round(hits / (n if L < 13000 else 60), 3)
            res["contact_" + str(L)] = round(contact / (n if L < 13000 else 60), 3)
            if ex is not None:
                res["example"] = ex
        out[str(beta)] = res
    print("free ribbon", {b: {k: v for k, v in r.items() if k != "example"} for b, r in out.items()}, flush=True)
    return out


if __name__ == "__main__":
    main()
