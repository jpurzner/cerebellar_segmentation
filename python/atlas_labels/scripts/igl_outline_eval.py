"""
Benchmark of the outline-driven IGL gap filler (igl_outline_fill.py) against the ribbon
gap model (igl_model.complete, the best setting of igl_model_eval.py) and the current
walk (igl_repair.walk), on synthetic gaps cut into the intact IGL of the gold standard.

    python scripts/igl_outline_eval.py [n_gaps] [section ...]     # default: 80 gaps, s1_2 (gold)

Gaps are made exactly as in igl_model_eval.py (50-500 um, cut into one skeleton branch,
>= 120 um of ribbon kept on both sides; the removed IGL = IGL pixels whose nearest
centre-line point lies in the cut stretch). One gap at a time; the rest stays intact.
Two scenarios:
  erased  the gap's IGL is simply missing (unlabelled tissue)
  as_wm   the gap's IGL is labelled white matter (as when bleached IGL is called WM)
The outline rule's parameters (m; rho, h0, margin) are fitted per gap on the section's
own IGL with the gap removed, so the fit never sees the answer. m is used either as the
section's fitted value ("global") or as the ML thickness at the two flanks ("local":
median of depth - width/2 over the flank, 40-200 um behind each end).
Fill: "geo" = the predicted IGL in the geodesic corridor between the two ends;
"band" = the predicted IGL inside a tube around the ribbon model's centre line.
Scores per gap: Dice, precision and recall of the fill against the removed IGL.
-> outputs/igl_outline/benchmark_<section>.json, examples_<section>.npz
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi

sys.path.insert(0, str(Path(__file__).parent))
import igl_model as M  # noqa: E402
import igl_model_eval as ME  # noqa: E402
import igl_outline_fill as OF  # noqa: E402
import igl_ribbon_eval as RE  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/igl_outline"
MODEL = dict(sigma_depth_um=15.0, bend_scale=0.3)     # best ribbon model (igl_model_eval)
MODES = ("label", "gated", "outline", "union")
N_EX = 14
WM_SHARE = 0.35      # corridor counted as covered by white-matter labels above this share


def prf(fill, gm):
    tp = np.count_nonzero(fill & gm)
    p = tp / max(np.count_nonzero(fill), 1)
    r = tp / max(np.count_nonzero(gm), 1)
    return dict(dice=round(2 * tp / max(np.count_nonzero(fill) + np.count_nonzero(gm), 1), 3), prec=round(p, 3), rec=round(r, 3))


def gap_setup(sec, gp):
    b = sec["rib"]["branches"][gp["branch"]]
    P = b["path"]
    j0, j1 = gp["j0"], gp["j1"]
    a, ta, fa = ME.flank(P[:j0], True)
    bb, tb_out, fb = ME.flank(P[j1 + 1:], False)
    W, D = sec["rib"]["width"], sec["dE"]
    wa, da = float(np.median(W[fa[:, 0], fa[:, 1]])), float(np.median(D[fa[:, 0], fa[:, 1]]))
    wb, db = float(np.median(W[fb[:, 0], fb[:, 1]])), float(np.median(D[fb[:, 0], fb[:, 1]]))
    ma = float(np.median(D[fa[:, 0], fa[:, 1]] - W[fa[:, 0], fa[:, 1]] / 2))
    mb = float(np.median(D[fb[:, 0], fb[:, 1]] - W[fb[:, 0], fb[:, 1]] / 2))
    gm = np.zeros(sec["rib"]["igl"].shape, bool)
    sel = (sec["ob"] == gp["branch"]) & (sec["oj"] >= j0) & (sec["oj"] <= j1)
    gm[sec["yy"][sel], sec["xx"][sel]] = True
    igl_d = sec["rib"]["igl"] & ~gm
    return dict(a=a, ta=ta, b=bb, tb_out=tb_out, tb_in=-tb_out, wa=wa, wb=wb, da=da, db=db,
                m_local=max((ma + mb) / 2, 8.0), gm=gm, igl_d=igl_d)


def run_gap(sec, geo, gp, scenario, want_ex):
    s = gap_setup(sec, gp)
    g = sec["g"]
    gm, igl_d = s["gm"], s["igl_d"]
    wm_lab = g["wm"] | gm if scenario == "as_wm" else g["wm"]
    ok = sec["ok"]
    forbid = g["egl"] | wm_lab | ~ok
    field = M.Field(sec["dE"], forbid) if scenario == "as_wm" else sec["field"]
    avoid = igl_d & ~ndi.binary_dilation(gm, iterations=1)
    C, info = M.complete(s["a"], s["ta"], s["b"], s["tb_in"], field, s["da"], s["db"], MODEL, avoid)
    # the ribbon model's own band (as benchmarked before); in as_wm the white-matter label is
    # not allowed to cut its band, which favours it
    keep_m = (ok & ~g["egl"] & ~igl_d) if scenario == "as_wm" else (sec["keep"] & ~igl_d)
    fills = {"ribbon model": ME.band(C, s["wa"], s["wb"], gm.shape, keep_m)}
    if scenario == "erased":
        wk = ME.run_walk(s["a"], s["ta"], s["da"], sec, igl_d)
        if wk is None:
            wk2 = ME.run_walk(s["b"], s["tb_out"], s["db"], sec, igl_d)
            wk = None if wk2 is None else wk2[::-1]
        fills["walk (current)"] = ME.band(wk, s["wa"], s["wb"], gm.shape, keep_m) if wk is not None else np.zeros_like(gm)
    allowed = geo["inner"]
    fits, status = {}, {}
    for mode in MODES:
        prm = OF.fit(geo, igl_d, wm_lab, mode)
        fits[mode] = prm
        for mname, mval in (("global", prm["m"]), ("local", s["m_local"])):
            P = OF.predict(geo, prm, wm_lab, mode, m=mval)
            F, st = OF.fill_between(P, igl_d, allowed, s["a"], s["b"], s["wa"], s["wb"])
            fills[f"{mode} {mname} geo"] = F
            status[f"{mode} {mname} geo"] = st
            Fb = OF.band(P, igl_d, allowed, C, s["wa"], s["wb"])
            fills[f"{mode} {mname} band"] = Fb
            # practical: corridor when the ends connect, otherwise the band
            fills[f"{mode} {mname} geo|band"] = F if st == "ok" else Fb
    # no white matter at all: depth band from the flanks, alone and minus the outline core
    w_gap = (s["wa"] + s["wb"]) / 2
    Pd = OF.predict_depth(geo, s["m_local"], w_gap)
    for nm, PP in (("depth local", Pd), ("depth+core local", Pd & ~OF.wm_outline(geo, fits["outline"]))):
        F, st = OF.fill_between(PP, igl_d, allowed, s["a"], s["b"], s["wa"], s["wb"])
        fills[f"{nm} geo"] = F
        status[f"{nm} geo"] = st
        fills[f"{nm} geo|band"] = F if st == "ok" else OF.band(PP, igl_d, allowed, C, s["wa"], s["wb"])
    # combined: is the corridor blocked when the white-matter labels are trusted? (no ground truth used)
    chord_px = float(np.hypot(*(np.asarray(s["b"]) - s["a"])))
    exp_px = chord_px * w_gap / OF.G
    lab = fills["label local geo"]
    blocked = status["label local geo"] != "ok" or np.count_nonzero(lab) < 0.5 * exp_px
    # second test: share of the depth-band corridor (drawn without any white matter) that carries
    # a white-matter label; real white matter lies deeper than the IGL band
    Fdc = fills["depth+core local geo|band"]
    wm_share = float(np.count_nonzero(Fdc & wm_lab) / max(np.count_nonzero(Fdc), 1))
    shallow = Fdc & (geo["dS"] <= s["m_local"] + w_gap / 2)
    wm_shallow = float(np.count_nonzero(shallow & wm_lab) / max(np.count_nonzero(shallow), 1))
    covered = wm_shallow > WM_SHARE
    # the ribbon model with the white-matter labels ignored (barrier = EGL and background only)
    C2, _ = M.complete(s["a"], s["ta"], s["b"], s["tb_in"], M.Field(sec["dE"], g["egl"] | ~ok), s["da"], s["db"], MODEL, avoid)
    fills["ribbon model, WM ignored"] = ME.band(C2, s["wa"], s["wb"], gm.shape, ok & ~g["egl"] & ~igl_d)
    fills["ribbon model, WM ignored, depth+core"] = OF.band(Pd & ~OF.wm_outline(geo, fits["outline"]), igl_d, allowed, C2, s["wa"], s["wb"])
    for fb in ("outline local geo|band", "depth local geo|band", "depth+core local geo|band", "ribbon model, WM ignored",
               "ribbon model, WM ignored, depth+core"):
        fills[f"auto: labels | {fb}"] = fills[fb] if blocked else fills["label local geo|band"]
        fills[f"auto: ribbon | {fb}"] = fills[fb] if blocked else fills["ribbon model"]
        fills[f"auto2: ribbon | {fb}"] = fills[fb] if covered else fills["ribbon model"]
        fills[f"auto2: labels | {fb}"] = fills[fb] if covered else fills["label local geo|band"]
    row = dict(section=sec["name"], scenario=scenario, blocked=bool(blocked), wm_share=round(wm_share, 3), wm_shallow=round(wm_shallow, 3), length_um=round(gp["length_um"]), width=[round(s["wa"]), round(s["wb"])],
               m_local=round(s["m_local"], 1), start=info["start"], fits=fits, status=status,
               scores={k: prf(v, gm) for k, v in fills.items()})
    ex = None
    if want_ex:
        ys, xs = np.nonzero(ndi.binary_dilation(gm, iterations=1))
        pad = 45
        y0, y1 = max(ys.min() - pad, 0), min(ys.max() + pad + 1, gm.shape[0])
        x0, x1 = max(xs.min() - pad, 0), min(xs.max() + pad + 1, gm.shape[1])
        crop = (slice(y0, y1), slice(x0, x1))
        keys = ["ribbon model", "walk (current)", "label local geo|band", "outline local geo|band", "depth local geo|band",
                "depth+core local geo|band", "auto: ribbon | outline local geo|band", "auto: ribbon | depth local geo|band",
                "auto: ribbon | depth+core local geo|band", "auto: labels | depth local geo|band", "auto2: ribbon | depth+core local geo|band",
                "ribbon model, WM ignored", "auto2: ribbon | ribbon model, WM ignored", "auto2: ribbon | ribbon model, WM ignored, depth+core"]
        ex = dict(box=[y0, y1, x0, x1], tissue=g["tissue"][crop], egl=g["egl"][crop], wm=wm_lab[crop], igl=igl_d[crop], gap=gm[crop],
                  **{k: fills[k][crop] for k in keys if k in fills})
    return row, ex


def summarise(rows):
    out = {}
    for sc in ("erased", "as_wm"):
        rr = [r for r in rows if r["scenario"] == sc]
        if not rr:
            continue
        methods = list(rr[0]["scores"])
        res = {}
        for k in methods:
            d = np.array([r["scores"][k]["dice"] for r in rr])
            p = np.array([r["scores"][k]["prec"] for r in rr])
            rc = np.array([r["scores"][k]["rec"] for r in rr])
            res[k] = dict(n=len(rr), mean_dice=round(float(d.mean()), 3), median_dice=round(float(np.median(d)), 3),
                          dice_ge_08=round(float((d >= 0.8).mean()), 2), dice_lt_05=round(float((d < 0.5).mean()), 2),
                          prec=round(float(p.mean()), 3), rec=round(float(rc.mean()), 3),
                          by_len={f"{lo}-{hi}": round(float(np.mean([r["scores"][k]["dice"] for r in rr if lo <= r["length_um"] < hi])), 3)
                                  for lo, hi in ME.BINS if any(lo <= r["length_um"] < hi for r in rr)})
        out[sc] = dict(sorted(res.items(), key=lambda t: -t[1]["mean_dice"]))
        out[sc + "_blocked_share"] = round(float(np.mean([r["blocked"] for r in rr])), 3)
    return out


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 80
    want = sys.argv[2:] or ["s1_2 (gold)"]
    OUT.mkdir(parents=True, exist_ok=True)
    ME.N_GAPS = n
    for si, (name, lab, px) in enumerate(RE.sections()):
        if name not in want:
            continue
        t0 = time.time()
        sec = ME.section_data(name, lab, px)
        geo = OF.geometry(sec["g"])
        full = {m: OF.fit(geo, sec["rib"]["igl"], sec["g"]["wm"], m) for m in MODES}
        print(name, "whole-section fits", full, flush=True)
        gaps = ME.make_gaps(sec, np.random.default_rng(900 + si))
        rows, exs = [], {}
        for gi, gp in enumerate(gaps):
            for sc in ("erased", "as_wm"):
                row, ex = run_gap(sec, geo, gp, sc, gi < N_EX)
                row["gap"] = gi
                rows.append(row)
                if ex is not None:
                    for k, v in ex.items():
                        exs[f"{sc}|{gi}|{k}"] = np.asarray(v)
            if gi % 10 == 9:
                print(f"  {gi + 1}/{len(gaps)} gaps, {time.time() - t0:.0f}s", flush=True)
        summ = summarise(rows)
        tag = name.split(" ")[0]
        (OUT / f"benchmark_{tag}.json").write_text(json.dumps(dict(section=name, n_gaps=len(gaps), whole_section_fits=full,
                                                                   summary=summ, gaps=rows), indent=1, default=float))
        np.savez_compressed(OUT / f"examples_{tag}.npz", **exs)
        for sc, res in summ.items():
            if not isinstance(res, dict):
                print(f"{name} {sc}: {res}")
                continue
            print(f"{name} [{sc}]")
            for k, v in list(res.items())[:30]:
                print(f"   {k:28s} dice {v['mean_dice']:.3f} (median {v['median_dice']:.3f}, >=0.8 {v['dice_ge_08']:.2f}, <0.5 {v['dice_lt_05']:.2f})"
                      f"  prec {v['prec']:.2f} rec {v['rec']:.2f}  by length {v['by_len']}")
        print(f"{name}: done in {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
