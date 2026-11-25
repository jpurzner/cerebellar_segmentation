"""
Error regions in the flattened cerebellum, and their correction (2026-11-21).

In the flattened strips (flatten.py --repaired: columns at right angles to the reference line, 2 um apart,
sampled every 2 um in depth) the layers should come out as horizontal bands in the order EGL -> ML -> IGL ->
white matter, at depths that change slowly along the line. Where the labels break that, something is wrong.

Per section:
1. Expected depths. EGL bottom, IGL top and IGL bottom per column from the labels (flatten.bounds), smoothed
   along the line with flatten.fit_line (rolling median, Gaussian, gaps <= 300 um bridged): the course of
   each boundary from the neighbouring columns.
2. The section's own stain model. Trusted columns = their boundaries within a few um of the expected ones,
   with ML between EGL and IGL. On their samples, a quadratic discriminant (QDA) learns EGL / ML / IGL /
   white matter from nuclei, NeuN and p27, each as log intensity against its running median along the line
   (bleaching cancels), at two smoothings. Each section is learned on itself, so a cKO section is not held
   to the controls' appearance. White matter and ML look alike in the stains; only the order separates them.
3. Column segmentation (Viterbi on a left-to-right model EGL -> ML -> IGL -> WM). Each sample scores the
   stain log-likelihood plus a bonus w_label where it agrees with its label; a layer boundary at depth d
   costs ((d - expected) / sigma)^2 / 2 (flat_cost where no expected depth); the EGL bottom's sigma is loose
   (30 um) because a thick, uneven EGL is the cKO phenotype. A search band around the labels' course (dev_*_um, as in
   OCT segmentation) is available but off: tried at IGL top +-24 um, it blocked the large boundary moves the user's
   own corrections make (gold s1_2: 84 -> 71 % of changes matching). Skipping a layer (no ML) costs skip_cost. The
   model covers each column from the first sample past the leading meninges / unassigned cells to the
   expected IGL bottom + below_bot_um (at most max_model_um). Two passes: the second takes its expected
   depths from the first segmentation.
4. Error regions = runs of columns where the segmentation and the labels disagree by >= min_run_um of
   contiguous depth (shorter disagreements are boundary jitter), closed over close_um, at least min_err_um
   long. Each is typed by its main change (IGL missing, IGL in the ML, inner EGL called IGL, ...).
   Geometry errors are found separately and only flagged: runs of columns with EGL under the line that end
   before they reach the IGL (shallower than the expected IGL top, or 60 um where none is expected) and have
   no IGL. There the reference line or the columns are wrong (a sliver between two cuts, a cut inside the
   EGL), and no relabelling can fix it.
   Not cortex (v2, the user's flatten v4 notes): stretches whose columns have no IGL label within 160 um,
   white matter from within 80 um making up 30 % of the column, no EGL band in the nuclei stain at the top and
   no EGL stain in that white matter, facing open background (a furrow base faces a narrow channel), are the
   deep cerebellar nuclei boundary or a cut edge: their EGL labels become unassigned, nothing else there is
   corrected, and they leave the EGL-covered line.
5. Correction: inside error regions the changed samples take the segmentation's layer; unassigned,
   anomaly and background samples are never changed. Back on the 4 um label grid by vote: a cell changes
   when at least half the column samples that fall in it (bilinear) agree. The gold s1_2 is detected
   but not changed.

    python scripts/flat_fix.py eval            # from the MATLAB labels: how much of the user's corrections it makes
    python scripts/flat_fix.py                 # all sections, on the final EGL repair -> outputs/flat_fix/
    python scripts/flat_fix.py s3_5 s4_A       # some sections
    python scripts/flat_fix.py export          # corrected label grids -> outputs/flat_fix/final/ (corrected(short))
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis

sys.path.insert(0, str(Path(__file__).parent))
import flatten as F  # noqa: E402
import egl_repair as ER  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/flat_fix"
STEP, G = F.STEP, F.G
GOLD = "s1_2"
LAYERS = (2, 3, 4, 5)                       # EGL, ML, IGL, WM label codes = model states 0..3
NAMES = {2: "EGL", 3: "ML", 4: "IGL", 5: "WM"}
PRM = dict(ref_win_um=400.0, sm_lo=1.0, sm_hi=3.0, lead_um=40.0, below_bot_um=24.0, below_top_um=100.0, max_model_um=160.0,
           trust_eb_um=12.0, trust_top_um=12.0, trust_bot_um=20.0, trust_ml_um=8.0, max_train=200000, qda_reg=0.05,
           ll_floor=-8.0, w_label=2.0, sig_eb_um=30.0, sig_top_um=12.0, sig_bot_um=20.0, flat_cost=4.0,
           skip_cost=6.0, start_ml_cost=3.0, dev_eb_um=1e9, dev_top_um=1e9, dev_bot_um=1e9, min_run_um=16.0, close_um=20.0, min_err_um=100.0, vote=0.5, geo_min_um=60.0,
           ncx_igl_um=160.0, ncx_wm_um=160.0, ncx_wm0_um=80.0, ncx_wm_frac=0.3, ncx_top_um=30.0, ncx_nuc_z=0.7, ncx_egl_frac=0.3,
           ncx_min_len_um=60.0, ncx_len_um=60.0, ncx_open_um=40.0,
           band=False, band_lam=1.5, band_free=1, stain_top_um=20.0, stain_top_frac=0.5, min_egl_run_um=8.0, egl_gap_um=8.0,
           ref_min_len_um=0.0)   # v10 work in progress (2026-11-24), off: defaults reproduce v9 / page v12


def strip_features(Z, L, prm=PRM):
    """(N, Dn, 6): nuclei, NeuN, p27 as log intensity minus the running median (ref_win_um along the line) of
    the column medians over tissue samples; smoothed sm_lo and sm_hi samples. v10: only columns at least
    ref_min_len_um long give the reference (short columns are mostly EGL and raised it: s3_2 at 8.2 mm, an outer
    EGL with the section's EGL nuclei brightness came out between ML and EGL)."""
    valid = Z["valid"]
    tis = valid & (L > 0)
    w = int(prm["ref_win_um"] / STEP) | 1
    long_ = valid.sum(1) * STEP >= prm.get("ref_min_len_um", 0.0)
    out = []
    for k in ("S_nuclei", "S_neun", "S_p27"):
        S = np.log(np.maximum(Z[k].astype(np.float32), 1.0))
        med = np.array([np.median(S[i][tis[i]]) if tis[i].any() and long_[i] else np.nan for i in range(len(S))])
        ref = pd.Series(med).rolling(w, center=True, min_periods=20).median().bfill().ffill().to_numpy()
        A = np.where(valid, S - ref[:, None], 0.0).astype(np.float32)
        for s in (prm["sm_lo"], prm["sm_hi"]):
            out.append(ndi.gaussian_filter(A, s))
    return np.stack(out, -1)


def ranges(L, valid, ok, exp_top, exp_bot, prm=PRM):
    """per column the modelled samples [lead, end): past the leading unassigned / anomaly cells (at most
    lead_um), to the column end, the expected IGL bottom + below_bot_um, or (no bottom expected) the expected
    IGL top + below_top_um, at most max_model_um; only columns with EGL under the line. Deeper, a column can
    cross white matter into the IGL of the next sub-folium, which the one-way model cannot follow."""
    N, Dn = L.shape
    skip = np.isin(L, (1, 6)) & valid
    lead = np.minimum(np.argmin(skip, axis=1), int(prm["lead_um"] / STEP))
    end = valid.sum(1)
    lim_um = np.where(np.isfinite(exp_bot), exp_bot + prm["below_bot_um"],
                      np.where(np.isfinite(exp_top), exp_top + prm["below_top_um"], prm["max_model_um"]))
    lim = np.ceil(np.minimum(lim_um, prm["max_model_um"]) / STEP).astype(int)
    end = np.minimum(end, lim)
    end = np.where(ok & (end > lead), end, lead)
    return lead, end


def expected(L, valid, ok, prm=PRM):
    """expected EGL bottom, IGL top, IGL bottom (um) per column from the neighbours, and the raw ones."""
    eb, top, bot = F.bounds(L, valid, F.PRM)
    t_ok = ok & valid[:, 3]
    fe, _ = F.fit_line(eb, t_ok, F.PRM)
    ft, _ = F.fit_line(top, t_ok, F.PRM)
    fb, _ = F.fit_line(bot, t_ok, F.PRM)
    return dict(eb=eb, top=top, bot=bot, feb=fe, ftop=ft, fbot=fb)


def stain_model(Fx, L, valid, ok, ex, lead, end, prm=PRM):
    """per-sample log-likelihood (N, Dn, 4) of EGL / ML / IGL / WM from the section's trusted columns."""
    tr = (ok & np.isfinite(ex["top"]) & np.isfinite(ex["eb"]) & (np.abs(ex["top"] - ex["ftop"]) <= prm["trust_top_um"])
          & (np.abs(ex["bot"] - ex["fbot"]) <= prm["trust_bot_um"]) & (np.abs(ex["eb"] - ex["feb"]) <= prm["trust_eb_um"])
          & (ex["top"] - ex["eb"] >= prm["trust_ml_um"]))
    N, Dn = L.shape
    k = np.arange(Dn)[None, :]
    inr = (k >= lead[:, None]) & (k < end[:, None])
    sel = tr[:, None] & inr & np.isin(L, LAYERS)
    idx = np.flatnonzero(sel.ravel())
    rng = np.random.default_rng(0)
    if len(idx) > prm["max_train"]:
        idx = rng.choice(idx, prm["max_train"], replace=False)
    Xf = Fx.reshape(-1, Fx.shape[-1])
    y = L.ravel()[idx]
    q = QuadraticDiscriminantAnalysis(reg_param=prm["qda_reg"]).fit(Xf[idx], y)
    ll = np.full((N * Dn, 4), prm["ll_floor"], np.float32)
    flat_in = np.flatnonzero(inr.ravel())
    lp = q.predict_log_proba(Xf[flat_in]) - np.log(q.priors_)[None, :]
    order = [list(q.classes_).index(c) for c in LAYERS]
    lp = lp[:, order]
    lp -= lp.max(1, keepdims=True)
    ll[flat_in] = np.maximum(lp, prm["ll_floor"])
    return ll.reshape(N, Dn, 4), dict(trusted_frac=round(float(tr[ok].mean()), 3) if ok.any() else 0.0, train_n=int(len(idx)))


def viterbi(ll, L, lead, end, exp3, prm=PRM, band3=None):
    """best layer sequence per column; states 0..3 (EGL, ML, IGL, WM), -1 outside the modelled range."""
    N, Dn, _ = ll.shape
    sig = (prm["sig_eb_um"], prm["sig_top_um"], prm["sig_bot_um"])
    R = end - lead
    K = int(R.max()) if len(R) else 0
    st = np.full((N, Dn), -1, np.int8)
    if K <= 0:
        return st
    rows = np.arange(N)

    dev = (prm["dev_eb_um"], prm["dev_top_um"], prm["dev_bot_um"])
    band3 = exp3 if band3 is None else band3

    def bcost(j, d):                                   # boundary j (0: EGL/ML, 1: ML/IGL, 2: IGL/WM) at sample d
        e, b = exp3[j], band3[j]
        c = np.where(np.isfinite(e), ((d * STEP - e) / sig[j]) ** 2 / 2, prm["flat_cost"])
        # searched only within dev of the course the labels give (the second pass may not drift further)
        return np.where(np.isfinite(b) & (np.abs(d * STEP - b) > dev[j]), 1e3, c)

    lab_bonus = np.stack([(L == c) for c in LAYERS], -1) * prm["w_label"]
    em_all = ll + lab_bonus
    bp = np.zeros((N, K, 4), np.int8)
    d0 = lead
    em0 = em_all[rows, np.minimum(d0, Dn - 1)]
    dp = np.stack([np.zeros(N), np.full(N, -prm["start_ml_cost"]) - bcost(0, d0), np.full(N, -np.inf), np.full(N, -np.inf)], 1) + em0
    for k in range(1, K):
        d = lead + k
        live = k < R
        em = np.where(live[:, None], em_all[rows, np.minimum(d, Dn - 1)], 0.0)
        bc = [bcost(j, d) for j in range(3)]
        new = np.empty_like(dp)
        for b in range(4):
            best, arg = dp[:, b].copy(), np.full(N, b, np.int8)
            for a in range(b):
                c = dp[:, a] - sum(bc[j] for j in range(a, b)) - prm["skip_cost"] * (b - a - 1)
                c = np.where(live, c, -np.inf)
                better = c > best
                best = np.where(better, c, best)
                arg = np.where(better, a, arg).astype(np.int8)
            new[:, b] = best + em[:, b]
            bp[:, k, b] = arg
        dp = new
    s = np.argmax(dp, 1).astype(np.int8)
    for k in range(K - 1, -1, -1):
        live = k < R
        d = np.minimum(lead + k, Dn - 1)
        st[rows[live], d[live]] = s[live]
        if k > 0:
            s = np.where(live, bp[rows, k, s], s).astype(np.int8)
    return st


def _trans(dp_from, bc, live, prm):
    """best predecessor per state (a <= b) and its score; bc[j] = boundary j cost at this depth."""
    N = dp_from.shape[0]
    best = np.full((N, 4), -np.inf)
    arg = np.zeros((N, 4), np.int8)
    for b in range(4):
        for a in range(b + 1):
            c = dp_from[:, a] - (sum(bc[j] for j in range(a, b)) + prm["skip_cost"] * (b - a - 1) if b > a else 0.0)
            if b > a:
                c = np.where(live, c, -np.inf)
            better = c > best[:, b]
            best[:, b] = np.where(better, c, best[:, b])
            arg[:, b] = np.where(better, a, arg[:, b])
    return best, arg


def egl_band(ll, L, lead, end, exp3, prm=PRM):
    """third pass (v10, the user's flatten v10 / v11 notes: wrong EGL boundaries under a smooth line). The same
    column model as viterbi, but the EGL bottom is chosen for all columns together: for every column and every
    depth the best score of the column with its EGL ending there (forward-backward max-marginals), then a
    dynamic programme along the line that charges band_lam per sample of change in the EGL bottom between
    neighbouring columns beyond band_free (the EGL bottom, like the surface, does not make sharp turns; a speckle
    of a few columns costs more than it gains, a change held over several columns is followed). Each column is
    then segmented with its EGL bottom fixed there. Returns the states (as viterbi) and the EGL bottom (samples)."""
    N, Dn, _ = ll.shape
    sig = (prm["sig_eb_um"], prm["sig_top_um"], prm["sig_bot_um"])
    R = end - lead
    K = int(R.max()) if len(R) else 0
    st = np.full((N, Dn), -1, np.int8)
    eb = lead.copy()
    if K <= 0:
        return st, eb
    rows = np.arange(N)

    def bcost(j, d):
        e = exp3[j]
        return np.where(np.isfinite(e), ((d * STEP - e) / sig[j]) ** 2 / 2, prm["flat_cost"])

    em_all = ll + np.stack([(L == c) for c in LAYERS], -1) * prm["w_label"]
    dep = lambda k: np.minimum(lead + k, Dn - 1)
    em = np.stack([em_all[rows, dep(k)] for k in range(K)], 1)                 # (N, K, 4)
    BC = np.stack([np.stack([bcost(j, lead + k) for j in range(3)], -1) for k in range(K)], 1)   # (N, K, 3)
    live = np.arange(K)[None, :] < R[:, None]
    # forward
    A = np.full((N, K, 4), -np.inf)
    A[:, 0] = np.stack([np.zeros(N), -prm["start_ml_cost"] - BC[:, 0, 0], np.full(N, -np.inf), np.full(N, -np.inf)], 1) + em[:, 0]
    for k in range(1, K):
        best, _ = _trans(A[:, k - 1], [BC[:, k, j] for j in range(3)], live[:, k], prm)
        A[:, k] = np.where(live[:, k, None], best + em[:, k], A[:, k - 1])
    # backward (score of depths k+1.. given state at k), with pointers
    Bt = np.zeros((N, K, 4))
    bp = np.zeros((N, K, 4), np.int8)
    for k in range(K - 2, -1, -1):
        nxt = Bt[:, k + 1] + em[:, k + 1]                                       # (N, 4) score from state b at k+1 on
        bc = [BC[:, k + 1, j] for j in range(3)]
        for a in range(4):
            best, arg = np.full(N, -np.inf), np.full(N, a, np.int8)
            for b in range(a, 4):
                c = nxt[:, b] - (sum(bc[j] for j in range(a, b)) + prm["skip_cost"] * (b - a - 1) if b > a else 0.0)
                better = c > best
                best, arg = np.where(better, c, best), np.where(better, b, arg).astype(np.int8)
            ok_k = live[:, k + 1]
            Bt[:, k, a] = np.where(ok_k, best, 0.0)
            bp[:, k, a] = np.where(ok_k, arg, a)
    # S[i, b]: EGL over k < b, a later state at k = b (b = R: EGL to the end of the range)
    S = np.full((N, K + 1), -np.inf)
    S[:, 0] = A[:, 0, 1] + Bt[:, 0, 1]
    first = np.zeros((N, K + 1), np.int8)
    first[:, 0] = 1
    for b in range(1, K):
        tb = np.stack([A[:, b - 1, 0] - sum(BC[:, b, j] for j in range(0, s)) - prm["skip_cost"] * (s - 1) + em[:, b, s] + Bt[:, b, s]
                       for s in (1, 2, 3)], 1)
        S[:, b] = np.where(live[:, b], tb.max(1), -np.inf)
        first[:, b] = (tb.argmax(1) + 1).astype(np.int8)
    S[rows, R] = A[rows, np.maximum(R - 1, 0), 0]
    # along the line, in absolute depth D = lead + b
    Dm = int((lead + R).max()) + 1
    G = np.full((N, Dm), -np.inf)
    for b in range(K + 1):
        m = b <= R
        G[rows[m], (lead + b)[m]] = S[m, b]
    G -= np.nanmax(np.where(np.isfinite(G), G, np.nan), 1, keepdims=True) if np.isfinite(G).any() else 0.0
    G[R <= 0] = 0.0                                                           # unmodelled columns: free
    G = np.where(np.isfinite(G), G, -1e6)
    t = np.arange(Dm)
    C = prm["band_lam"] * np.maximum(np.abs(t[:, None] - t[None, :]) - prm["band_free"], 0)
    dp = G[0].copy()
    ptr = np.zeros((N, Dm), np.int32)
    for i in range(1, N):
        M = dp[:, None] - C
        ptr[i] = np.argmax(M, 0)
        dp = M[ptr[i], t] + G[i]
    D = np.zeros(N, int)
    D[-1] = int(np.argmax(dp))
    for i in range(N - 1, 0, -1):
        D[i - 1] = ptr[i, D[i]]
    # segment each modelled column with its EGL bottom at D
    b = np.clip(D - lead, 0, np.maximum(R, 0))
    for i in np.flatnonzero(R > 0):
        bi, Ri = int(b[i]), int(R[i])
        st[i, lead[i]:lead[i] + bi] = 0
        if bi < Ri:
            s = int(first[i, bi])
            for k in range(bi, Ri):
                st[i, lead[i] + k] = s
                if k < Ri - 1:
                    s = int(bp[i, k, s])
    eb = np.where(R > 0, lead + b, lead)
    return st, eb


def seg_bounds(st, ok, valid):
    """expected depths from a segmentation (first sample of each layer), smoothed along the line."""
    N, Dn = st.shape
    out = []
    for j in (1, 2, 3):
        has = (st >= j).any(1)
        v = np.where(has, np.argmax(st >= j, 1) * STEP, np.nan)
        f, _ = F.fit_line(v, ok & valid[:, 3], F.PRM)
        out.append(f)
    return out


def _runs_close(m, close, min_len):
    m = ndi.binary_closing(np.pad(m, close), np.ones(2 * close + 1, bool))[close:-close] if close else m
    out = np.zeros_like(m)
    for a, b in F._runs(m):
        if b - a >= min_len:
            out[a:b] = True
    return out


def stain_egl(Z, short):
    """the EGL repair's stain EGL (bright nuclei, low NeuN) sampled on the strip."""
    E = np.load(ER.OUT / f"{short}.npz")
    valid = Z["valid"]
    N, Dn = valid.shape
    return ndi.map_coordinates(E["egl_stain"].astype(np.uint8), Z["X"].reshape(-1, 2).T, order=0, mode="constant").reshape(N, Dn).astype(bool) & valid


def stain_top(L, valid, egl_s, prm=PRM):
    """columns whose first stain_top_um (past the leading unassigned / anomaly cells) are at least stain_top_frac
    stain EGL: modelled even without an EGL label at the top (s3_2 at 8.2 mm: ML labelled over the EGL)."""
    skip = np.isin(L, (1, 6)) & valid
    lead = np.minimum(np.argmin(skip, axis=1), int(prm["lead_um"] / STEP))
    k = int(prm["stain_top_um"] / STEP)
    idx = np.minimum(lead[:, None] + np.arange(k)[None, :], L.shape[1] - 1)
    r = np.arange(L.shape[0])[:, None]
    return (egl_s[r, idx] & valid[r, idx]).mean(1) >= prm["stain_top_frac"]


def segment(Z, L, prm=PRM, egl_s=None):
    valid, ok = Z["valid"], Z["ok"]
    ex = expected(L, valid, ok, prm)
    mod = ok | (stain_top(L, valid, egl_s, prm) if egl_s is not None and prm.get("band") else False)
    lead, end = ranges(L, valid, mod, ex["ftop"], ex["fbot"], prm)
    Fx = strip_features(Z, L, prm)
    ll, info = stain_model(Fx, L, valid, ok, ex, lead, end, prm)
    st = viterbi(ll, L, lead, end, (ex["feb"], ex["ftop"], ex["fbot"]), prm)
    e2 = seg_bounds(st, ok, valid)
    exp2 = tuple(np.where(np.isfinite(a), a, b) for a, b in zip(e2, (ex["feb"], ex["ftop"], ex["fbot"])))
    if prm.get("band"):
        st, _ = egl_band(ll, L, lead, end, exp2, prm)
    else:
        st = viterbi(ll, L, lead, end, exp2, prm, band3=(ex["feb"], ex["ftop"], ex["fbot"]))
    info["modelled_extra_um"] = int((mod & ~ok).sum() * STEP)
    return st, ex, info, (lead, end)


def deep_egl(L, st, prm=PRM):
    """labelled EGL below the column's top EGL zone: past the model's EGL bottom and separated from the top run of
    EGL labels (gaps up to egl_gap_um allowed) by more than egl_gap_um. A fissure wall the column runs along, or the
    next folium's EGL: not the top band's business, left as labelled."""
    N, Dn = L.shape
    out = np.zeros((N, Dn), bool)
    gap = int(prm["egl_gap_um"] / STEP)
    mod = st >= 0
    for i in np.flatnonzero(((L == 2) & mod).any(1)):
        d = np.flatnonzero(mod[i])
        a, b = d[0], d[-1] + 1
        eb = a + int(np.argmax(np.r_[st[i, a:b] != 0, True]))          # model EGL bottom
        e = np.flatnonzero(L[i, a:b] == 2) + a
        z = a - 1
        for x in e:
            if x - z - 1 > gap:
                break
            z = x
        z = max(z + 1, eb)
        out[i, z:b] = L[i, z:b] == 2
    return out


def errors(L, st, prm=PRM):
    """changed samples (contiguous >= min_run_um in depth), error-region columns, the corrected strip. v10: a
    column's whole EGL band removed or added counts even when thinner than min_run_um (s2_4 at 14.1 mm: a 10 um
    EGL label over ML-like nuclei on a cut wall), and an EGL boundary moved by at least min_egl_run_um."""
    model = np.where(st >= 0, np.array(LAYERS)[np.maximum(st, 0)], 0).astype(np.uint8)
    ch0 = (st >= 0) & np.isin(L, LAYERS) & (model != L)
    k = max(1, int(prm["min_run_um"] / STEP))
    ch = ndi.binary_opening(ch0, np.ones((1, k), bool))
    if prm.get("band"):
        lab_e = ((L == 2) & (st >= 0)).any(1)
        mod_e = (st == 0).any(1)
        whole = ((lab_e & ~mod_e)[:, None] & ch0 & (L == 2)) | ((~lab_e & mod_e)[:, None] & ch0 & (model == 2))
        ke = max(1, int(prm["min_egl_run_um"] / STEP))
        egl_ch = ndi.binary_opening(ch0 & ((L == 2) | (model == 2)), np.ones((1, ke), bool))
        ch |= whole | egl_ch
        ch &= ~deep_egl(L, st, prm)
    col = ch.any(1)
    reg = _runs_close(col, int(prm["close_um"] / STEP / 2), int(prm["min_err_um"] / STEP))
    fix = ch & reg[:, None]
    Lf = L.copy()
    Lf[fix] = model[fix]
    return Lf, fix, reg


TYPES = {(3, 4): "IGL missing (ML -> IGL)", (5, 4): "IGL missing (WM -> IGL)", (4, 3): "IGL in the ML (IGL -> ML)",
         (4, 2): "inner EGL called IGL (IGL -> EGL)", (2, 4): "IGL called EGL (EGL -> IGL)", (3, 2): "EGL missing (ML -> EGL)",
         (2, 3): "EGL too deep (EGL -> ML)", (4, 5): "WM in the IGL band (IGL -> WM)", (5, 3): "WM above the IGL (WM -> ML)",
         (3, 5): "ML below the IGL (ML -> WM)", (2, 5): "EGL -> WM", (5, 2): "WM -> EGL"}


def regions(L, Lf, fix, reg):
    out = []
    for a, b in F._runs(reg):
        f = fix[a:b]
        pairs = pd.Series(list(zip(L[a:b][f].tolist(), Lf[a:b][f].tolist()))).value_counts()
        main = pairs.index[0] if len(pairs) else None
        out.append(dict(a=int(a), b=int(b), um=int((b - a) * STEP), changed_um2=int(f.sum() * STEP * STEP),
                        type=TYPES.get(main, f"{NAMES.get(main[0])} -> {NAMES.get(main[1])}") if main else "",
                        pairs={f"{NAMES[p[0]]}->{NAMES[p[1]]}": int(n * STEP * STEP) for p, n in pairs.items()}))
    return out


def geometry(L, valid, ok, ex, prm=PRM):
    """columns with EGL under the line that end before reaching the IGL, in runs >= min_err_um."""
    ext = valid.sum(1) * STEP
    need = np.where(np.isfinite(ex["ftop"]), np.maximum(ex["ftop"], prm["geo_min_um"]), prm["geo_min_um"])
    m = ok & ~np.isfinite(ex["top"]) & (ext < need)
    return _runs_close(m, int(prm["close_um"] / STEP / 2), int(prm["min_err_um"] / STEP))


def noncortex(L, valid, ok, nuc, egl_s, prm=PRM):
    """stretches of line that are not a cortical surface (the deep cerebellar nuclei boundary, a cut edge; the
    user's flatten v4 notes on s5_A-C and s5_E). Per column: no IGL label within ncx_igl_um; white matter starts
    within ncx_wm0_um and makes up at least ncx_wm_frac of the first ncx_wm_um; no EGL band at the top (the
    nuclei over the first ncx_top_um below ncx_nuc_z of the way from the section's ML to its EGL level, on the
    strip; the thin EGL labels there have ML-like nuclei); and the white matter does not look like EGL (a deep
    furrow's fused EGL labelled white matter, s5_F). Runs of at least ncx_len_um."""
    N, Dn = L.shape
    d = (np.arange(Dn) * STEP)[None, :]
    m_ = valid & ok[:, None]
    egl_lvl, ml_lvl = float(np.median(nuc[m_ & (L == 2)])), float(np.median(nuc[m_ & (L == 3)]))
    top = valid & (d < prm["ncx_top_um"]) & (L > 0)
    z = np.where(top.any(1), (np.where(top, nuc, 0).sum(1) / np.maximum(top.sum(1), 1) - ml_lvl) / max(egl_lvl - ml_lvl, 1e-6), 9.0)
    has_igl = ((L == 4) & valid & (d < prm["ncx_igl_um"])).any(1)
    W = (L == 5) & valid & (d < prm["ncx_wm_um"])
    w0 = np.where(W.any(1), W.argmax(1) * STEP, np.inf)
    wfrac = W.sum(1) / np.maximum((valid & (d < prm["ncx_wm_um"])).sum(1), 1)
    wegl = np.where(W.any(1), (egl_s & W).sum(1) / np.maximum(W.sum(1), 1), 0.0)
    m = (ok & ~has_igl & (w0 < prm["ncx_wm0_um"]) & (wfrac >= prm["ncx_wm_frac"]) & (z < prm["ncx_nuc_z"])
         & (wegl < prm["ncx_egl_frac"]) & (valid.sum(1) * STEP >= prm["ncx_min_len_um"]))
    return _runs_close(m, int(prm["close_um"] / STEP / 2), int(prm["ncx_len_um"] / STEP))


def open_front(Q, n, tissue, max_um=320.0):
    """per line point, how much background lies in front of it: the largest distance to tissue met walking
    outward (against the column direction) until tissue is reached again; inf if the walk leaves the image.
    A furrow base faces a narrow channel, the deep nuclei boundary or a cut edge open background."""
    bg = ~tissue
    dt = ndi.distance_transform_edt(bg) * G
    H, W = tissue.shape
    v = -n / np.maximum(np.hypot(*n.T), 1e-6)[:, None]
    best = np.zeros(len(Q))
    alive = np.ones(len(Q), bool)
    for k in range(1, int(max_um / G)):
        q = np.round(Q + v * k).astype(int)
        out = (q[:, 0] < 0) | (q[:, 0] >= H) | (q[:, 1] < 0) | (q[:, 1] >= W)
        best[alive & out] = np.inf
        alive &= ~out
        qc = np.clip(q, 0, [H - 1, W - 1])
        inbg = bg[qc[:, 0], qc[:, 1]]
        alive &= inbg | (k <= 3)
        best = np.where(alive & inbg, np.maximum(best, dt[qc[:, 0], qc[:, 1]]), best)
    return best


def paint(lab, tissue, X, valid, fix, Lf, prm=PRM):
    """the strip's changes back on the label grid: bilinear votes of the changed samples against all samples."""
    H, W = lab.shape
    tot = np.zeros((H, W))
    TARGETS = (1,) + LAYERS
    votes = {c: np.zeros((H, W)) for c in TARGETS}

    def splat(acc, pts, w=None):
        y, x = pts[:, 0], pts[:, 1]
        y0, x0 = np.floor(y).astype(int), np.floor(x).astype(int)
        fy, fx = y - y0, x - x0
        for dy, dx, ww in ((0, 0, (1 - fy) * (1 - fx)), (1, 0, fy * (1 - fx)), (0, 1, (1 - fy) * fx), (1, 1, fy * fx)):
            yy, xx = np.clip(y0 + dy, 0, H - 1), np.clip(x0 + dx, 0, W - 1)
            np.add.at(acc, (yy, xx), ww if w is None else ww * w)

    splat(tot, X[valid])
    for c in TARGETS:
        m = fix & (Lf == c)
        if m.any():
            splat(votes[c], X[m])
    new = lab.copy()
    best_c, best_v = np.zeros((H, W), np.uint8), np.zeros((H, W))
    for c in TARGETS:
        better = votes[c] > best_v
        best_c[better], best_v[better] = c, votes[c][better]
    chg = tissue & np.isin(lab, LAYERS) & (best_v >= prm["vote"] * np.maximum(tot, 1e-9)) & (tot > 0) & (best_c != lab)
    new[chg] = best_c[chg]
    return new, chg


def grids(short, Z=None):
    """the label grid and tissue mask the flattening used: the EGL repair's, with the shape check's changes (flatten
    v9: meninges trimmed -> unassigned, cracks filled -> tissue)."""
    lab, tissue = ER.repaired(short)
    Z = Z if Z is not None else np.load(F.OUT2 / f"{short}.npz")
    if "lab_grid" in Z.files:
        lab, tissue = Z["lab_grid"].astype(np.uint8), tissue | Z["shape_fill"]
    return lab, tissue


def run(short, prm=PRM, save=True):
    t0 = time.time()
    Z = np.load(F.OUT2 / f"{short}.npz")
    L, valid, X = Z["L_repaired"], Z["valid"], Z["X"]
    flat = X.reshape(-1, 2).T
    N, Dn = valid.shape
    egl_s = stain_egl(Z, short)
    st, ex, info, (lead, end) = segment(Z, L, prm, egl_s)
    Lf, fix, reg = errors(L, st, prm)
    ncx = noncortex(L, valid, Z["ok"], strip_features(Z, L, prm)[..., 1], egl_s, prm)
    if ncx.any():                                           # a furrow base faces a narrow channel: that is cortex (s5_F)
        front = open_front(Z["Q"], Z["n"], Z["pia"] if "pia" in Z.files else ER.pia_tissue(short))
        for a, b in F._runs(ncx):
            if np.median(front[a:b]) < prm["ncx_open_um"]:
                ncx[a:b] = False
    # not cortex: no layer corrections there; its EGL labels become unassigned
    fix &= ~ncx[:, None]
    reg &= ~ncx
    Lf = np.where(ncx[:, None], L, Lf)
    rm = ncx[:, None] & (L == 2) & valid
    Lf[rm] = 1
    fix |= rm
    regs = regions(L, Lf, fix & ~ncx[:, None], reg)
    regs += [dict(a=int(a), b=int(b), um=int((b - a) * STEP), changed_um2=int(rm[a:b].sum() * STEP * STEP),
                  type="not cortex (deep nuclei boundary or cut edge): EGL label removed", pairs={"EGL->unassigned": int(rm[a:b].sum() * STEP * STEP)},
                  noncortex=True) for a, b in F._runs(ncx)]
    geo = geometry(L, valid, Z["ok"], ex, prm) & ~ncx
    regs += [dict(a=int(a), b=int(b), um=int((b - a) * STEP), changed_um2=0, type="columns end before the IGL (geometry)", pairs={}, geometry=True)
             for a, b in F._runs(geo)]
    lab, tissue = grids(short, Z)
    if short == GOLD:
        new, chg = lab.copy(), np.zeros_like(tissue)
    else:
        new, chg = paint(lab, tissue, X, valid, fix, Lf, prm)
    Lnew = ndi.map_coordinates(new, flat, order=0, mode="constant").reshape(N, Dn).astype(np.uint8)
    Lnew[~valid] = 0
    ok_new = Z["ok"] & ~ncx                                  # not-cortex stretches leave the EGL-covered line
    before = F.igl_line(L, valid, Z["ok"], F.PRM)["stats"]
    after = F.igl_line(Lnew, valid, ok_new, F.PRM)["stats"]
    a2 = lambda m: int(m.sum() * G * G)
    trans = {f"{NAMES[a]}->{NAMES[b]}": a2(chg & (lab == a) & (new == b)) for a in LAYERS for b in LAYERS if a != b}
    res = dict(section=short, group="control" if short in F.D.controls() else "cKO", **info,
               error_um=int(reg.sum() * STEP), error_frac=round(float(reg[Z["ok"]].mean()), 3) if Z["ok"].any() else 0.0,
               regions=int(sum(not r.get("geometry") and not r.get("noncortex") for r in regs)), geometry_um=int(geo.sum() * STEP),
               noncortex_um=int(ncx.sum() * STEP),
               geometry_regions=int(sum(bool(r.get("geometry")) for r in regs)), changed_um2=a2(chg), changes_um2={k: v for k, v in trans.items() if v},
               igl_present=dict(before=round(before["igl_present_um"] / max(before["covered_um"], 1), 3),
                                after=round(after["igl_present_um"] / max(after["covered_um"], 1), 3)),
               flags_um=dict(before={k: before[f"{k}_um"] for k in ("missing", "off", "no_ml")} | dict(short=before["short_columns_um"]),
                             after={k: after[f"{k}_um"] for k in ("missing", "off", "no_ml")} | dict(short=after["short_columns_um"])),
               region_types=pd.Series([r["type"] for r in regs if not r.get("geometry") and not r.get("noncortex")]).value_counts().to_dict() if regs else {})
    if save:
        OUT.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(OUT / f"{short}.npz", st=st, fix=fix, reg=reg, geo=geo, ncx=ncx, ok_new=ok_new, L=L, L_fix=Lf, L_new=Lnew, lab_new=new, changed=chg,
                            lead=lead, end=end, **{f"exp_{k}": v for k, v in ex.items()})
        (OUT / f"{short}_regions.json").write_text(json.dumps(regs))
    res["seconds"] = round(time.time() - t0, 1)
    print(json.dumps(res), flush=True)
    return res


FINAL = OUT / "final"


def export(shorts=None):
    """corrected label grids (lab_new of each <short>.npz) -> final/<short>_labels.tif (same grid and codes as
    egl_repair.repaired; the gold s1_2 unchanged)."""
    import tifffile
    FINAL.mkdir(parents=True, exist_ok=True)
    for short in shorts or F.all_shorts():
        new = np.load(OUT / f"{short}.npz")["lab_new"].astype(np.uint8)
        lab, tissue = grids(short)
        assert new.shape == lab.shape and (new[~tissue] == lab[~tissue]).all()
        tifffile.imwrite(FINAL / f"{short}_labels.tif", new, compression="zlib")
        tifffile.imwrite(FINAL / f"{short}_tissue.tif", tissue.astype(np.uint8), compression="zlib")
    (FINAL / "manifest.json").write_text(json.dumps(dict(version="flat_fix v9 (2026-11-24), on flatten v11 (the walk: one loop clockwise from the cut surface, corners cut by arcs, no jumps or hops over cuts, no closed slits, cuts continued to the furrow base, thin ML caps relabelled) / EGL repair v10", grid_um=G,
                                                          codes=ER.CODES, script="scripts/flat_fix.py", gold_unchanged=GOLD), indent=1))


def corrected(short):
    import tifffile
    return tifffile.imread(FINAL / f"{short}_labels.tif"), tifffile.imread(FINAL / f"{short}_tissue.tif") > 0


egl_depth = F.egl_depth                     # moved to flatten (v7: the peaks also drive flatten.peak_cuts)


def egl_peaks(L, valid, ok, thr_um=40.0, win_um=200.0, sep_um=60.0):
    """flatten.egl_peaks: [(index, excess um)], excess per column."""
    pk, exc, _, _ = F.egl_peaks(L, valid, ok, thr_um, win_um, sep_um)
    return pk, exc


def evaluate(prm=PRM, shorts=None):
    """start from the MATLAB labels on the same columns: how much of the user's corrections (reviewed controls,
    corrected cKO sections, the gold s1_2 tracing) the correction makes, and how much else it changes."""
    rows = []
    for short in shorts or [s for s in F.all_shorts() if s != "s1_5"]:
        Z = np.load(F.OUT2 / f"{short}.npz")
        auto, best, tissue = F.labels(short)
        valid, X = Z["valid"], Z["X"]
        N, Dn = valid.shape
        flat = X.reshape(-1, 2).T
        samp = lambda g: np.where(valid, ndi.map_coordinates(g, flat, order=0, mode="constant").reshape(N, Dn), 0).astype(np.uint8)
        La, Lb = samp(auto), samp(best)
        st, ex, info, (lead, end) = segment(Z, La, prm, stain_egl(Z, short))
        Lf, fix, reg = errors(La, st, prm)
        k = np.arange(Dn)[None, :]
        inr = (k >= lead[:, None]) & (k < end[:, None]) & np.isin(La, LAYERS) & np.isin(Lb, LAYERS)
        user = inr & (La != Lb)
        mine = fix & inr
        good = mine & (Lf == Lb)
        a2 = lambda m: int(m.sum() * STEP * STEP)
        r = dict(section=short, group="control" if short in F.D.controls() else "cKO", user_changed_um2=a2(user), fix_changed_um2=a2(mine),
                 agree_um2=a2(good), recall=round(a2(good) / max(a2(user), 1), 3), precision=round(a2(good) / max(a2(mine), 1), 3),
                 against_user_um2=a2(mine & ~user), error_frac=round(float(reg[Z["ok"]].mean()), 3))
        rows.append(r)
        print(json.dumps(r), flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "eval.json").write_text(json.dumps(dict(prm=prm, sections=rows), indent=1))
    return rows


if __name__ == "__main__":
    if sys.argv[1:2] == ["eval"]:
        evaluate(shorts=sys.argv[2:] or None)
        sys.exit()
    if sys.argv[1:2] == ["export"]:
        export(sys.argv[2:] or None)
        sys.exit()
    shorts = sys.argv[1:] or F.all_shorts()
    rows = [run(s) for s in shorts]
    if not sys.argv[1:]:
        (OUT / "summary.json").write_text(json.dumps(dict(prm=PRM, sections=rows), indent=1))
