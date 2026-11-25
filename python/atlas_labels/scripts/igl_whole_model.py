"""
Whole-IGL model and anomaly detection.

    python scripts/igl_whole_model.py profiles     # describe every section's IGL ribbon (cached)
    python scripts/igl_whole_model.py fit          # fit the control model, leave-one-out threshold, flag all sections
    python scripts/igl_whole_model.py breaks       # cheapest continuation from every ribbon end

Every IGL ribbon (skeleton branch of the cleaned IGL, 4 um grid, resampled every 4 um)
is described point by point:
  w      width (2 x distance to the IGL edge)
  d      depth below the EGL (distance from the centre line to the nearest EGL)
  ml     ML thickness: closest approach of the ribbon's cross-section to the EGL
  wmgap  1 where the cross-section's inner edge is > 12 um from white matter (lifted off)
  kappa  signed curvature (per mm): + where the ribbon wraps around a folium tip,
         - where it wraps around a fissure base
  neun, dapi, p27   local IGL intensity / the section's median IGL intensity (log2)
Segments are 100 um windows (+-48 um) along the ribbon, summarised as
  [log w, d, ml, wmgap, kappa, width CV, neun, dapi, p27]
The control model is a Gaussian mixture over standardised segment features, fitted on the
five intact-IGL controls (s1_2 gold, s5_F, s5_G, s6_C, s6_G). A segment's anomaly score is
-log p(segment); the threshold is the 99th percentile of leave-one-control-out scores (each
control scored by a model fitted on the other four). Runs of anomalous segments >= 40 um
become flags; the feature furthest from the responsible mixture component names the flag.
Breaks: every ribbon end is completed by the gap model (igl_model.complete) to the best
other end or ribbon side within 500 um; an end is a likely break when that completion is
cheap (threshold calibrated on synthetic gaps vs. natural ends of the controls,
igl_whole_model_eval.py).
Outputs in outputs/igl_whole/.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import tifffile
from scipy import ndimage as ndi
from sklearn.mixture import GaussianMixture

sys.path.insert(0, str(Path(__file__).parent))
import igl_model as M  # noqa: E402
import igl_ribbon as RB  # noqa: E402
import igl_ribbon_eval as RE  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/igl_whole"
PROF = OUT / "profiles"
G = RB.GRID
FX = RE.FX
SB20 = RE.SEG / "matlab_anomaly/sandbox"
SB10 = RE.SEG / "matlab_anomaly_10x/sandbox"
TRAIN = ["s1_2 (gold)", "s5_F", "s5_G", "s6_C", "s6_G"]
CH20 = dict(p27=0, neun=1, dapi=2)
CH10 = dict(p27=0, dapi=1, neun=2)
POINT = ["w", "d", "ml", "wmgap", "kappa", "neun", "dapi", "p27"]
SEG = ["logw", "d", "ml", "wmgap", "kappa", "wcv", "neun", "dapi", "p27", "contact"]
HALF = 12                       # window half-length in points (48 um)
Q_GLOBAL, Q_LOCAL = 99.5, 99.8  # leave-one-control-out percentiles used as thresholds
MIN_RUN = 10                    # 40 um
TYPES = {("logw", 1): "IGL too thick", ("logw", -1): "IGL too thin", ("d", 1): "IGL too deep below the EGL",
         ("d", -1): "IGL too close to the EGL", ("ml", -1): "ML missing (IGL reaches the EGL)", ("ml", 1): "ML too thick",
         ("wmgap", 1): "IGL detached from white matter", ("wmgap", -1): "IGL detached from white matter",
         ("kappa", 1): "sharp bend", ("kappa", -1): "sharp bend", ("wcv", 1): "uneven width", ("wcv", -1): "uneven width",
         ("neun", -1): "dim NeuN", ("neun", 1): "bright NeuN", ("dapi", -1): "dim DAPI", ("dapi", 1): "bright DAPI",
         ("p27", -1): "dim p27", ("p27", 1): "bright p27", ("contact", 1): "IGL touches the EGL (no ML)",
         ("contact", -1): "IGL touches the EGL (no ML)"}


def sections():
    """(name, labels_path_or_'gold', px, image_path, channels, genotype)"""
    out = [("s1_2 (gold)", "gold", 0.5119, SB20 / "2018_05_22_s1_2_p27-0021_native.tif", CH20, "control")]
    for f in sorted(FX.glob("2018_05_22_*_fixed_labels.tif")):
        sid = f.name.replace("_fixed_labels.tif", "")
        p = sid.split("_")
        short = f"{p[3]}_{p[4]}"
        out.append((short + (" (auto)" if short == "s1_2" else ""), f, 0.5119, SB20 / f"{sid}_native.tif", CH20,
                    "control" if short == "s1_2" else "EZH2 cKO"))
    for f in sorted(FX.glob("s*_p27_fixed_labels.tif")):
        sid = f.name.replace("_fixed_labels.tif", "")
        out.append((sid.replace("_p27", ""), f, 0.65, SB10 / f"{sid}_native.tif", CH10, "control"))
    return out


def load_labels(src):
    return RE.gold() if src == "gold" else tifffile.imread(src)


def grid_images(path, chmap, px, shape):
    raw = tifffile.imread(path)
    f = max(1, int(round(G / px)))
    out = {}
    for k, ci in chmap.items():
        a = raw[ci].astype(np.float32)
        H2, W2 = a.shape[0] // f * f, a.shape[1] // f * f
        b = a[:H2, :W2].reshape(H2 // f, f, W2 // f, f).mean((1, 3))
        b = ndi.zoom(b, f * px / G, order=1)
        o = np.zeros(shape, np.float32)
        hh, ww = min(shape[0], b.shape[0]), min(shape[1], b.shape[1])
        o[:hh, :ww] = b[:hh, :ww]
        out[k] = o
    return out


def profile(g, imgs, rib=None):
    """point features along every ribbon branch."""
    rib = rib or RB.describe(g)
    igl = rib["igl"]
    dE = ndi.gaussian_filter(rib["depth_egl"], 1.0)
    gy, gx = np.gradient(dE)
    dW = ndi.distance_transform_edt(~g["wm"]) * G
    wmax = ndi.maximum_filter(rib["width"], 3)
    ratio = {}
    den = ndi.gaussian_filter(igl.astype(np.float32), 4.0)
    for k, im in imgs.items():
        med = float(np.median(im[igl])) if igl.any() else 1.0
        num = ndi.gaussian_filter(im * igl, 4.0)
        ratio[k] = np.log2(np.maximum(num / np.maximum(den, 1e-3), 1e-3) / max(med, 1e-3))
    rows, pos, br, idx = [], [], [], []
    for bi, b in enumerate(rib["branches"]):
        if len(b["path"]) < 6:
            continue
        C = M.smooth_path(b["path"])
        if len(C) < 6:
            continue
        d1 = np.gradient(C, axis=0)
        d2 = np.gradient(d1, axis=0)
        sp = np.maximum(np.hypot(*d1.T), 1e-9)
        T = d1 / sp[:, None]
        K = (d2 - (d2 * T).sum(1)[:, None] * T) / sp[:, None] ** 2          # curvature vector, 1/px
        K = ndi.gaussian_filter1d(K, 3.0, axis=0, mode="nearest")
        at = lambda a: ndi.map_coordinates(a, C.T, order=1, mode="nearest")
        nout = -np.stack([at(gy), at(gx)], 1)
        nout /= np.maximum(np.hypot(*nout.T), 1e-9)[:, None]
        kap = -(K * nout).sum(1) / G * 1000.0
        ii = np.clip(np.round(C).astype(int), 0, np.array(igl.shape) - 1)
        w = wmax[ii[:, 0], ii[:, 1]]
        d = at(dE)
        f = np.stack([w, d, d - w / 2, np.clip(at(dW) - w / 2, 0, 150), kap] +
                     [at(ratio[k]) if k in ratio else np.zeros(len(C)) for k in ("neun", "dapi", "p27")], 1)
        rows.append(f)
        pos.append(C)
        br.append(np.full(len(C), bi))
        idx.append(np.arange(len(C)))
    if not rows:
        return None
    F, pos = np.vstack(rows), np.vstack(pos)
    # cross-sections: every IGL pixel belongs to its nearest centre-line point; the closest
    # approach of a cross-section to the EGL is the ML thickness on that side, and the
    # distance of its inner edge to white matter says whether it lies on the white matter
    from scipy.spatial import cKDTree
    yy, xx = np.nonzero(igl)
    dd, nn = cKDTree(pos).query(np.stack([yy, xx], 1))
    ok = dd * G <= np.maximum(F[nn, 0], 8.0) / 2 + 12
    eg = np.full(len(pos), np.inf)
    np.minimum.at(eg, nn[ok], rib["depth_egl"][yy[ok], xx[ok]])
    wmin = np.full(len(pos), np.inf)
    np.minimum.at(wmin, nn[ok], dW[yy[ok], xx[ok]])
    F[:, 2] = np.where(np.isfinite(eg), eg, F[:, 1] - F[:, 0] / 2)            # ml = closest approach to the EGL
    # (segments() adds 'contact': share of the 100 um stretch where that approach is < 8 um)
    F[:, 3] = (np.where(np.isfinite(wmin), wmin, 999) > 12).astype(float)    # 1 = lifted off the white matter
    return dict(F=F, pos=pos, branch=np.concatenate(br), idx=np.concatenate(idx))


def segments(P):
    """100 um window features for every point."""
    F, br = P["F"], P["branch"]
    out = np.zeros((len(F), len(SEG)))
    for b in np.unique(br):
        sel = np.nonzero(br == b)[0]
        x = F[sel]
        n = len(sel)
        cs = lambda a: np.concatenate([[0], np.cumsum(a)])
        lw = np.log(x[:, 0] + 1)
        cols = [lw, x[:, 1], x[:, 2], x[:, 3], x[:, 4], x[:, 0], x[:, 0] ** 2, x[:, 5], x[:, 6], x[:, 7], (x[:, 2] < 8).astype(float)]
        C = [cs(c) for c in cols]
        lo = np.clip(np.arange(n) - HALF, 0, n)
        hi = np.clip(np.arange(n) + HALF + 1, 0, n)
        cnt = hi - lo
        mean = lambda k: (C[k][hi] - C[k][lo]) / cnt
        mw, mw2 = mean(5), mean(6)
        sd = np.sqrt(np.maximum(mw2 - mw ** 2, 0))
        out[sel] = np.stack([mean(0), mean(1), mean(2), mean(3), mean(4), sd / np.maximum(mw, 1), mean(7), mean(8), mean(9), mean(10)], 1)
    return out


LOCAL = ["dlogw", "dd", "dml", "dwmgap", "kappa", "wcv", "dneun", "ddapi", "dp27", "contact"]
REL = {0: 0, 1: 1, 2: 2, 3: 3, 6: 6, 7: 7, 8: 8}           # SEG columns compared with the neighbours
FLANK = (37, 112)                                           # 150-450 um, in 4 um points
LTYPES = {("dlogw", 1): "IGL thicker than its neighbours", ("dlogw", -1): "IGL thinner than its neighbours",
          ("dd", 1): "IGL deeper than its neighbours", ("dd", -1): "IGL closer to the EGL than its neighbours",
          ("dml", -1): "IGL reaches toward the EGL (ML thinner than its neighbours)", ("dml", 1): "ML thicker than its neighbours",
          ("dwmgap", 1): "IGL lifted off the white matter", ("dwmgap", -1): "IGL lying on white matter where its neighbours do not",
          ("kappa", 1): "sharp bend", ("kappa", -1): "sharp bend", ("wcv", 1): "uneven width", ("wcv", -1): "uneven width",
          ("dneun", -1): "dim NeuN", ("dneun", 1): "bright NeuN", ("ddapi", -1): "dim DAPI", ("ddapi", 1): "bright DAPI",
          ("dp27", -1): "dim p27", ("dp27", 1): "bright p27", ("contact", 1): "IGL touches the EGL (no ML)",
          ("contact", -1): "IGL touches the EGL (no ML)"}


def local_features(P, S):
    """segment features relative to the same ribbon 150-450 um away on both sides."""
    R = S.copy()
    sec_med = np.median(S, axis=0)
    for b in np.unique(P["branch"]):
        sel = np.nonzero(P["branch"] == b)[0]
        x = S[sel]
        n = len(sel)
        for i in range(n):
            fl = np.r_[np.arange(max(0, i - FLANK[1]), max(0, i - FLANK[0])), np.arange(min(n, i + FLANK[0]), min(n, i + FLANK[1]))]
            ref = np.median(x[fl], axis=0) if len(fl) >= 10 else sec_med
            for c in REL:
                R[sel[i], c] = x[i, c] - ref[c]
    return R


class LocalModel:
    """robust per-feature z-scores of the neighbour-relative features."""

    FLOOR = np.array([0.1, 5.0, 5.0, 0.1, 2.0, 0.03, 0.05, 0.05, 0.05, 0.1])   # smallest meaningful spread per feature

    def __init__(self, R):
        self.med = np.median(R, axis=0)
        self.sc = np.maximum(1.4826 * np.median(np.abs(R - self.med), axis=0), self.FLOOR)

    def z(self, R):
        return (R - self.med) / self.sc


# --------------------------------------------------------------------------- #
class ControlModel:
    def __init__(self, X, k=None, seed=0):
        self.mu, self.sd = X.mean(0), X.std(0) + 1e-6
        Z = (X - self.mu) / self.sd
        if k is None:
            bics = {}
            for kk in range(1, 9):
                gm = GaussianMixture(kk, covariance_type="full", reg_covar=1e-3, random_state=seed).fit(Z)
                bics[kk] = gm.bic(Z)
            k = min(bics, key=bics.get)
            self.bics = bics
        self.k = k
        self.gm = GaussianMixture(k, covariance_type="full", reg_covar=1e-3, random_state=seed, n_init=2).fit(Z)

    def score(self, X):
        return -self.gm.score_samples((X - self.mu) / self.sd)

    def explain(self, x):
        """per-feature z-score against the most probable component."""
        z = (x - self.mu) / self.sd
        c = int(self.gm.predict(z[None])[0])
        m, S = self.gm.means_[c], self.gm.covariances_[c]
        return (z - m) / np.sqrt(np.diag(S))


def combined(S, R, gmodel, tau_g, lmodel, tau_l):
    """anomaly level per point: > 1 means anomalous (global density or a local feature)."""
    zl = np.abs(lmodel.z(R)) / tau_l
    g = gmodel.score(S) / tau_g
    return np.maximum(g, zl.max(1)), g, zl


def flags(P, S, score, tau, model, min_run=MIN_RUN, R=None, lmodel=None, tau_l=None):
    """score/tau: global density score and its threshold; with R/lmodel/tau_l the local
    test is added and the flag is named by whichever test is furthest over its threshold."""
    out = []
    lev_g = score / tau
    if R is not None:
        zl = lmodel.z(R)
        lev_l = np.abs(zl) / tau_l
        lev = np.maximum(lev_g, lev_l.max(1))
    else:
        lev = lev_g
    an = lev > 1
    for b in np.unique(P["branch"]):
        sel = np.nonzero(P["branch"] == b)[0]
        a = an[sel]
        a = ndi.binary_closing(np.pad(a, 3), np.ones(5, bool))[3:-3]          # bridge 1-4 point dips
        lab, n = ndi.label(a)
        for r in range(1, n + 1):
            j = np.nonzero(lab == r)[0]
            if len(j) < min_run:
                continue
            pk = sel[j[np.argmax(lev[sel[j]])]]
            if R is not None and lev_l[pk].max() >= lev_g[pk]:
                order = np.argsort(-lev_l[pk])
                why = [dict(test="local", feature=LOCAL[i], z=round(float(zl[pk, i]), 1), over=round(float(lev_l[pk, i]), 2),
                            type=LTYPES[(LOCAL[i], 1 if zl[pk, i] > 0 else -1)]) for i in order[:3]]
            else:
                z = model.explain(S[pk])
                order = np.argsort(-np.abs(z))
                why = [dict(test="global", feature=SEG[i], z=round(float(z[i]), 1), over=round(float(lev_g[pk]), 2),
                            type=TYPES[(SEG[i], 1 if z[i] > 0 else -1)]) for i in order[:3]]
            out.append(dict(branch=int(b), i0=int(sel[j[0]]), i1=int(sel[j[-1]]), length_um=round(len(j) * G),
                            peak=int(pk), level=round(float(lev[pk]), 2), score=round(float(score[pk]), 1), type=why[0]["type"], why=why,
                            centre=P["pos"][sel[j[len(j) // 2]]].round(1).tolist(),
                            values=dict(zip(SEG, [round(float(v), 2) for v in S[pk]]))))
    return out


# --------------------------------------------------------------------------- #
# breaks: can the ribbon cheaply continue from an end?
R_END, R_SIDE, MAX_LEN, N_FULL = 500.0, 400.0, 600.0, 8
PRM = None


def model_params():
    global PRM
    if PRM is None:
        PRM = json.loads((ROOT / "outputs/igl_model/benchmark.json").read_text())["summary"]["best_params"]
    return PRM


WM_SCALE = 50.0                 # white matter crossed is charged (um crossed / WM_SCALE)^2


def break_context(g, rib, wm_ok=False):
    """wm_ok: white matter may be crossed (bleached IGL is often labelled white matter);
    the crossing is charged as (um of white matter crossed / 50 um)^2 instead."""
    import igl_model_apply as IA
    forbid = g["egl"] | ~(g["tissue"] | g["anom"]) | (g["wm"] if not wm_ok else False)
    return dict(g=g, rib=rib, igl=rib["igl"], field=M.Field(rib["depth_egl"], forbid), wm_ok=wm_ok,
                sk=IA.skeleton_targets(rib, rib["depth_egl"], rib["width"]))


def end_completion(ctx, e):
    """cheapest continuation of the ribbon from end e (gap model energy + length / 300 um)."""
    import igl_model_apply as IA
    prm = model_params()
    a = np.asarray(e["pos"], float)
    ta = np.asarray(e["dir"], float)
    route = IA.route_map(ctx["igl"], a, a, int((R_END + 40) / G))
    cands = []
    for e2 in ctx["rib"]["ends"]:
        if e2 is e:
            continue
        cands.append(("end", np.asarray(e2["pos"], float), -np.asarray(e2["dir"], float), e2["depth"], e2["width"]))
    for t in ctx["sk"]:
        if np.hypot(*(t["pos"] - a)) * G <= R_SIDE:
            for sgn in (1, -1):
                cands.append(("ribbon", t["pos"], sgn * t["t"], t["depth"], t["width"]))
    pre = []
    for kind, b, tb_in, db, wb in cands:
        v = b - a
        dist = np.hypot(*v)
        if dist * G < 20 or dist * G > R_END or v @ ta < -0.17 * dist or v @ tb_in < -0.17 * dist:
            continue
        iy, ix = np.clip(np.round(b).astype(int), 0, np.array(route.shape) - 1)
        if np.isfinite(route[iy, ix]) and route[iy, ix] < 3 * dist:
            continue
        Hc = M.hermite(a, ta, b, tb_in, 30)
        t = M.energy_terms(Hc, a, ta, b, tb_in, ctx["field"], e["depth"], db, prm)
        if t["length_um"] > MAX_LEN:
            continue
        pre.append((t["bend"] + t["depth"] + t["barrier"] + t["length_um"] / 300, kind, b, tb_in, db, wb))
    pre.sort(key=lambda q: q[0])
    best = None
    for _, kind, b, tb_in, db, wb in pre[:N_FULL]:
        avoid = ctx["igl"].copy()
        yy, xx = np.ogrid[:avoid.shape[0], :avoid.shape[1]]
        for p in (a, b):
            avoid &= (yy - p[0]) ** 2 + (xx - p[1]) ** 2 > 16
        C, info = M.complete(a, ta, b, tb_in, ctx["field"], e["depth"], db, prm, avoid)
        L = float(np.hypot(*np.diff(C, axis=0).T).sum() * G)
        if L > MAX_LEN:
            continue
        ii = np.clip(np.round(C).astype(int), 0, np.array(ctx["igl"].shape) - 1)
        wm_um = float(ctx["g"]["wm"][ii[:, 0], ii[:, 1]].mean() * L)
        sc = info["energy"] + L / 300 + ((wm_um / WM_SCALE) ** 2 if ctx["wm_ok"] else 0.0)
        if best is None or sc < best["score"]:
            best = dict(score=round(sc, 2), energy=round(info["energy"], 2), length_um=round(L), target=kind,
                        wm_um=round(wm_um), curve=C.round(1).tolist(), width=[round(e["width"]), round(wb)])
    return best


def section_breaks(g, rib):
    ctx = break_context(g, rib)
    ctx_wm = break_context(g, rib, wm_ok=True)
    out = []
    for e in rib["ends"]:
        out.append(dict(pos=np.asarray(e["pos"]).tolist(), dir=np.round(e["dir"], 3).tolist(),
                        width=round(e["width"]), depth=round(e["depth"]), best=end_completion(ctx, e),
                        best_wm=end_completion(ctx_wm, e)))
    return out


def build_breaks():
    meta = json.loads((OUT / "profiles.json").read_text())
    res = {}
    for name, src, px, img, ch, geno in sections():
        t0 = time.time()
        g = RB.grid_labels(load_labels(src), px)
        rib = RB.describe(g)
        res[name] = section_breaks(g, rib)
        sc = sorted(b["best"]["score"] for b in res[name] if b["best"])
        sw = sorted(b["best_wm"]["score"] for b in res[name] if b["best_wm"])
        print(f"{name:14s} {len(res[name])} ends; cheapest: {sc[:4]}; crossing white matter: {sw[:4]}  {time.time() - t0:.0f}s", flush=True)
    (OUT / "breaks.json").write_text(json.dumps(res))


# --------------------------------------------------------------------------- #
def build_profiles():
    PROF.mkdir(parents=True, exist_ok=True)
    meta = []
    for name, src, px, img, ch, geno in sections():
        t0 = time.time()
        lab = load_labels(src)
        g = RB.grid_labels(lab, px)
        rib = RB.describe(g)
        imgs = grid_images(img, ch, px, g["igl"].shape) if img.exists() else {}
        P = profile(g, imgs, rib)
        if P is None:
            print(name, "no ribbon")
            continue
        S = segments(P)
        ends = np.array([e["pos"] for e in rib["ends"]]) if rib["ends"] else np.zeros((0, 2))
        np.savez_compressed(PROF / f"{name.replace(' ', '_').replace('(', '').replace(')', '')}.npz", S=S, ends=ends,
                            **P)
        meta.append(dict(name=name, px=px, genotype=geno, labels=str(src), image=str(img), points=len(S),
                         ribbon_mm=round(len(S) * G / 1000, 1), ends=len(ends), shape=list(g["igl"].shape)))
        print(f"{name:14s} {len(S):6d} points, {len(ends)} ends, {time.time() - t0:.0f}s", flush=True)
    (OUT / "profiles.json").write_text(json.dumps(meta, indent=1))


def load_profile(name):
    z = np.load(PROF / f"{name.replace(' ', '_').replace('(', '').replace(')', '')}.npz")
    return {k: z[k] for k in z.files}


def fit():
    meta = json.loads((OUT / "profiles.json").read_text())
    prof = {m["name"]: load_profile(m["name"]) for m in meta}
    Xtr = np.vstack([prof[n]["S"][::3] for n in TRAIN])
    full = ControlModel(Xtr)
    print("mixture components (BIC):", full.k, {k: round(v) for k, v in full.bics.items()})

    # leave-one-control-out scores -> threshold
    loo = {}
    for n in TRAIN:
        Xo = np.vstack([prof[m]["S"][::3] for m in TRAIN if m != n])
        mo = ControlModel(Xo, k=full.k)
        loo[n] = (mo, mo.score(prof[n]["S"]))
    pooled = np.concatenate([s for _, s in loo.values()])
    tau = float(np.percentile(pooled, Q_GLOBAL))
    Rs = {n: local_features(prof[n], prof[n]["S"]) for n in prof}
    lfull = LocalModel(np.vstack([Rs[n] for n in TRAIN]))
    lloo = {n: LocalModel(np.vstack([Rs[m] for m in TRAIN if m != n])) for n in TRAIN}
    zp = np.abs(np.vstack([lloo[n].z(Rs[n]) for n in TRAIN]))
    tau_l = np.percentile(zp, Q_LOCAL, axis=0)
    res = dict(k=full.k, tau=round(tau, 2), features=SEG, local_features=LOCAL, tau_local=tau_l.round(2).tolist(),
               local_med=lfull.med.round(3).tolist(), local_scale=lfull.sc.round(3).tolist(), q_global=Q_GLOBAL, q_local=Q_LOCAL,
               sections={})
    for m in meta:
        n = m["name"]
        P = prof[n]
        mo, sc = loo[n] if n in loo else (full, full.score(P["S"]))
        lm = lloo.get(n, lfull)
        fl = flags(P, P["S"], sc, tau, mo, R=Rs[n], lmodel=lm, tau_l=tau_l)
        res["sections"][n] = dict(genotype=m["genotype"], ribbon_mm=m["ribbon_mm"], training=n in TRAIN,
                                  flagged_um=int(sum(f["length_um"] for f in fl)), n_flags=len(fl),
                                  flagged_per_mm=round(sum(f["length_um"] for f in fl) / 1000 / max(m["ribbon_mm"], 1e-6) * 100, 1),
                                  flags=fl)
        np.save(OUT / f"score_{n.replace(' ', '_').replace('(', '').replace(')', '')}.npy", sc)
        print(f"{n:14s} {m['genotype']:9s} {'train' if n in TRAIN else '     '} ribbon {m['ribbon_mm']:5.1f} mm  "
              f"flagged {res['sections'][n]['flagged_per_mm']:5.1f}%  ({len(fl)} flags)  types: "
              + ", ".join(sorted({f['type'] for f in fl})), flush=True)
    (OUT / "flags.json").write_text(json.dumps(res, indent=1, default=float))
    import pickle
    (OUT / "control_model.pkl").write_bytes(pickle.dumps(dict(model=full, tau=tau)))
    return res


if __name__ == "__main__":
    {"profiles": build_profiles, "fit": fit, "breaks": build_breaks}[sys.argv[1]]()
