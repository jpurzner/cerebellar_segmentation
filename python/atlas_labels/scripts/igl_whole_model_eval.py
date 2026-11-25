"""
Tests of the whole-IGL anomaly model (igl_whole_model.py).

    python scripts/igl_whole_model_eval.py synthetic   # artificial defects in the controls
    python scripts/igl_whole_model_eval.py breaks      # break threshold: synthetic gaps vs natural ends
    python scripts/igl_whole_model_eval.py real        # flags vs the user's marked problems on s4_A / s4_C

synthetic: each control is scored by a model fitted on the other four. Per control and seed,
four sites (>= 300 um apart, on branches >= 600 um) get one defect each, over 150-250 um:
  thin   ribbon cut to its central 40 %
  bulge  ribbon widened to ~2x
  smear  ML within 100 um of the stretch relabelled IGL (the IGL reaches the EGL)
  dim    all three channels x 0.35 in a 226 x 169 um field centred on the stretch (bleach)
A defect is found when a flag lies within 60 um of the stretch. False flags: flags further
than 150 um from every defect that were not already flagged on the undamaged section.
breaks: positives = the IGL end at one side of a synthetic gap (50-500 um cut, as in
igl_model_eval), searched for its cheapest continuation without being told the partner;
negatives = every end of the intact controls. Threshold = best balanced accuracy.
"""
import json
import pickle
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).parent))
import igl_ribbon as RB  # noqa: E402
import igl_ribbon_eval as RE  # noqa: E402
import igl_whole_model as W  # noqa: E402

G = RB.GRID
OUT = W.OUT
SEEDS = (0, 1, 2)
KINDS = ("thin", "bulge", "smear", "dim")


def control_inputs():
    for name, src, px, img, ch, geno in W.sections():
        if name in W.TRAIN:
            lab = W.load_labels(src)
            g = RB.grid_labels(lab, px)
            yield name, g, W.grid_images(img, ch, px, g["igl"].shape)


def loo_models():
    meta = json.loads((OUT / "profiles.json").read_text())
    prof = {m["name"]: W.load_profile(m["name"]) for m in meta if m["name"] in W.TRAIN}
    res = json.loads((OUT / "flags.json").read_text())
    Rs = {n: W.local_features(prof[n], prof[n]["S"]) for n in W.TRAIN}
    out = {}
    for n in W.TRAIN:
        X = np.vstack([prof[m]["S"][::3] for m in W.TRAIN if m != n])
        out[n] = (W.ControlModel(X, k=res["k"]), W.LocalModel(np.vstack([Rs[m] for m in W.TRAIN if m != n])))
    return out, res["tau"], np.array(res["tau_local"])


def damage(g, imgs, rib, sites, rng):
    igl = rib["igl"].copy()
    ml = g["ml"].copy()
    imgs = {k: v.copy() for k, v in imgs.items()}
    cl = np.vstack([b["path"] for b in rib["branches"]])
    owner = np.concatenate([[(bi, j) for j in range(len(b["path"]))] for bi, b in enumerate(rib["branches"])])
    yy, xx = np.nonzero(igl)
    d, nn = cKDTree(cl).query(np.stack([yy, xx], 1))
    ob, oj = owner[nn, 0], owner[nn, 1]
    defects = []
    for (bi, j), kind in zip(sites, KINDS):
        b = rib["branches"][bi]
        L = rng.uniform(150, 250)
        h = int(round(L / 2 / (G * 1.1)))
        j0, j1 = max(j - h, 0), min(j + h, len(b["path"]) - 1)
        seg = b["path"][j0:j1 + 1]
        sel = (ob == bi) & (oj >= j0) & (oj <= j1)
        m = np.zeros(igl.shape, bool)
        m[seg[:, 0], seg[:, 1]] = True
        near = ndi.distance_transform_edt(~m) * G
        if kind == "thin":
            cut = d[sel] * G > 0.2 * b["width"][oj[sel]]
            igl[yy[sel][cut], xx[sel][cut]] = False
        elif kind == "bulge":
            wl = float(np.median(b["width"][j0:j1 + 1]))
            igl |= (near <= wl) & g["tissue"] & ~g["egl"]
        elif kind == "smear":
            add = (near <= 100) & ml
            igl |= add
            ml &= ~add
        else:
            c = b["path"][j]
            hy, hx = 169 / 2 / G, 226 / 2 / G
            box = np.zeros(igl.shape, bool)
            box[max(int(c[0] - hy), 0):int(c[0] + hy) + 1, max(int(c[1] - hx), 0):int(c[1] + hx) + 1] = True
            for k in imgs:
                imgs[k][box] *= 0.35
        defects.append(dict(kind=kind, stretch=seg, length_um=round(L)))
    g2 = dict(g)
    g2["igl"], g2["ml"] = igl, ml
    return g2, imgs, defects


def near_any(P, pts, r_um):
    if not len(pts):
        return np.zeros(len(P), bool)
    d, _ = cKDTree(np.asarray(pts, float)).query(P)
    return d * G <= r_um


def synthetic():
    models, tau, tau_l = loo_models()
    found = {k: [0, 0] for k in KINDS}
    false_flags, trials, rows = 0, 0, []
    for name, g, imgs in control_inputs():
        mo, lm = models[name]
        rib0 = RB.describe(g)
        P0 = W.profile(g, imgs, rib0)
        S0 = W.segments(P0)
        f0 = W.flags(P0, S0, mo.score(S0), tau, mo, R=W.local_features(P0, S0), lmodel=lm, tau_l=tau_l)
        base_pts = [P0["pos"][f["i0"]:f["i1"] + 1] for f in f0]
        base_pts = np.vstack(base_pts) if base_pts else np.zeros((0, 2))
        for seed in SEEDS:
            rng = np.random.default_rng(300 + seed)
            sites = RE.pick_sites(rib0, rng, n=4)
            g2, im2, defects = damage(g, imgs, rib0, sites, rng)
            P = W.profile(g2, im2)
            S = W.segments(P)
            fl = W.flags(P, S, mo.score(S), tau, mo, R=W.local_features(P, S), lmodel=lm, tau_l=tau_l)
            res = []
            for dft in defects:
                hit = [f for f in fl if near_any(P["pos"][f["i0"]:f["i1"] + 1], dft["stretch"], 60).any()]
                found[dft["kind"]][0] += 1
                found[dft["kind"]][1] += bool(hit)
                res.append(f"{dft['kind']} {'found (' + hit[0]['type'] + ')' if hit else 'missed'}")
            allst = np.vstack([d_["stretch"] for d_ in defects])
            ff = [f for f in fl if not near_any(P["pos"][f["i0"]:f["i1"] + 1], allst, 150).any()
                  and not near_any(P["pos"][f["i0"]:f["i1"] + 1], base_pts, 40).any()]
            false_flags += len(ff)
            trials += 1
            rows.append(dict(section=name, seed=seed, defects=res, false_flags=[f["type"] for f in ff]))
            print(f"{name:12s} seed {seed}: " + "; ".join(res) + f" | false flags {len(ff)}", flush=True)
    summ = dict(found={k: f"{v[1]}/{v[0]}" for k, v in found.items()}, false_flags_per_section=round(false_flags / trials, 2))
    print(json.dumps(summ, indent=1))
    (OUT / "synthetic_eval.json").write_text(json.dumps(dict(summary=summ, trials=rows), indent=1))


def breaks():
    import igl_model_eval as ME
    neg, pos = [], []
    br = json.loads((OUT / "breaks.json").read_text())
    for n in W.TRAIN:
        neg += [e["best"]["score"] if e["best"] else np.inf for e in br[n]]
    secs = list(RE.sections())
    for si, (name, lab, px) in enumerate(secs):
        sec = ME.section_data(name, lab, px)
        gaps = ME.make_gaps(sec, np.random.default_rng(200 + si))[1::3]          # 12 held-out-style gaps
        for gp in gaps:
            b = sec["rib"]["branches"][gp["branch"]]
            j0, j1 = gp["j0"], gp["j1"]
            gm = np.zeros(sec["rib"]["igl"].shape, bool)
            sel = (sec["ob"] == gp["branch"]) & (sec["oj"] >= j0) & (sec["oj"] <= j1)
            gm[sec["yy"][sel], sec["xx"][sel]] = True
            g2 = dict(sec["g"])
            g2["igl"] = sec["g"]["igl"] & ~gm
            rib2 = RB.describe(g2)
            edge = b["path"][max(j0 - 1, 0)]
            cand = [e for e in rib2["ends"] if np.hypot(*(e["pos"] - edge)) * G <= 40]
            if not cand:
                pos.append(np.nan)
                continue
            e = min(cand, key=lambda e: np.hypot(*(e["pos"] - edge)))
            best = W.end_completion(W.break_context(g2, rib2), e)
            pos.append(best["score"] if best else np.inf)
        print(name, "positives so far", len(pos), flush=True)
    pos_a = np.array([p for p in pos if np.isfinite(p) or p == np.inf])
    neg_a = np.array(neg)
    ths = np.unique(np.concatenate([pos_a[np.isfinite(pos_a)], neg_a[np.isfinite(neg_a)]]))
    best = max(ths, key=lambda t: (pos_a <= t).mean() - (neg_a <= t).mean())
    from sklearn.metrics import roc_auc_score
    y = np.r_[np.ones(len(pos_a)), np.zeros(len(neg_a))]
    s = -np.r_[np.minimum(pos_a, 1e6), np.minimum(neg_a, 1e6)]
    summ = dict(threshold=round(float(best), 2), auc=round(float(roc_auc_score(y, s)), 3),
                gap_ends_found=f"{int((pos_a <= best).sum())}/{len(pos_a)}", gap_edges_without_an_end=int(np.isnan(np.array(pos)).sum()),
                natural_ends_flagged=f"{int((neg_a <= best).sum())}/{len(neg_a)}",
                positives=sorted(round(float(v), 1) for v in pos_a), negatives=sorted(round(float(v), 1) for v in neg_a))
    print(json.dumps(summ, indent=1, default=float))
    (OUT / "break_threshold.json").write_text(json.dumps(summ, indent=1, default=float))


def real():
    import igl_check_page as P
    import igl_finalize as F
    trs = {t["sid"]: t for t in json.loads((P.OUT / "transforms.json").read_text())}
    res = json.loads((OUT / "flags.json").read_text())
    br = json.loads((OUT / "breaks.json").read_text())
    thr = json.loads((OUT / "break_threshold.json").read_text())["threshold"] if (OUT / "break_threshold.json").exists() else None
    marks = {"s4_A_p27": {"marker 1 (break, folium 9)": (0.6168, 0.4863), "folium 6 smear (marker 2)": (0.3504, 0.1428)},
             "s4_C_p27": {"C3 (folium 4)": (0.3321, 0.2319), "D4 (folium 4)": (0.4805, 0.3616), "E4 (folium 7)": (0.6192, 0.3807)}}
    out = {}
    for sid, mk in marks.items():
        n = sid.replace("_p27", "")
        P_ = W.load_profile(n)
        fl = res["sections"][n]["flags"]
        bx = np.load(W.ROOT / f"outputs/igl_ribbon/{sid}_ribbon_repair_grid.npz")["bleach"]
        bpts = np.argwhere(bx)
        items = {}
        for lab_, xy in list(mk.items()) + [("bleach box (confirmed)", None)]:
            ref = bpts if xy is None else F.disp_to_grid(trs[sid], *xy)[None]
            hits = [f"{f['type']} ({f['length_um']} um)" for f in fl if near_any(P_["pos"][f["i0"]:f["i1"] + 1], ref, 150).any()]
            bh = [f"break: continuation {b['best']['length_um']} um, score {b['best']['score']}" for b in br[n]
                  if b["best"] and near_any(np.array([b["pos"]], float), ref, 200).any() and (thr is None or b["best"]["score"] <= thr)]
            items[lab_] = hits + bh
        out[n] = items
        print(n)
        for k, v in items.items():
            print(f"   {k:28s}: {v if v else 'nothing flagged'}")
    (OUT / "real_check.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    {"synthetic": synthetic, "breaks": breaks, "real": real}[sys.argv[1]]()
