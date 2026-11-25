"""
Data for the IGL ribbon model explainer page (outputs/igl_model/igl_model_explainer.html).

    python scripts/igl_model_explainer.py

  stats     - per-control ribbon statistics (ribbon_statistics.json)
  ribbon    - one real ribbon stretch (s5_G): label map as SVG paths, centre line with
              curvature, width and depth along it, cross-sections, iso-depth line
  accuracy  - median centre-line error by gap length for each model (held-out gaps)
  examples  - synthetic gaps from the benchmark: what was cut, the model's centre line
              and painted band, the old repair, the smooth curve and the straight line
  explorer  - one long folium-tip gap completed with a grid of weights, with each
              result's energy terms
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi
from skimage import measure

sys.path.insert(0, str(Path(__file__).parent))
import igl_model as M  # noqa: E402
import igl_model_eval as ME  # noqa: E402
import igl_ribbon as RB  # noqa: E402
import igl_ribbon_eval as RE  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/igl_model"
G = RB.GRID
STEP = 2.0                                   # um per crop pixel
BEST = "model s15.0 b0.3"
ELASTIC = "elastic s15.0 b0.3"


def paths(mask, tol=0.8):
    """mask -> SVG path string (even-odd fill) in crop pixels."""
    m = np.pad(mask.astype(float), 1)
    out = []
    for c in measure.find_contours(m, 0.5):
        c = measure.approximate_polygon(c - 1, tol)
        if len(c) < 3:
            continue
        out.append("M" + "L".join(f"{x:.1f},{y:.1f}" for y, x in c) + "Z")
    return "".join(out)


def crop_maps(lab, px, y0, x0, n, gap_grid=None):
    """native labels sampled at STEP um on an n x n crop whose corner is (y0, x0) um."""
    yy, xx = np.mgrid[0:n, 0:n].astype(float)
    nat = np.stack([(y0 + yy * STEP) / px, (x0 + xx * STEP) / px])
    L = ndi.map_coordinates(lab, nat, order=0, mode="constant", cval=0)
    maps = dict(tissue=L > 0, egl=np.isin(L, (2, 3)), igl=L == 4, wm=np.isin(L, (6, 8)))
    if gap_grid is not None:
        gq = np.stack([(y0 + yy * STEP) / G, (x0 + xx * STEP) / G])
        gap = ndi.map_coordinates(gap_grid.astype(np.uint8), gq, order=0, mode="constant", cval=0) > 0
        maps["gap"] = gap & maps["igl"]
        maps["igl"] = maps["igl"] & ~maps["gap"]
    return {k: paths(v) for k, v in maps.items()}


def to_crop(C, y0, x0):
    C = np.asarray(C, float) * G
    return [[round(float((p[1] - x0) / STEP), 1), round(float((p[0] - y0) / STEP), 1)] for p in C]


def band_path(C, wa, wb, shape, keep, y0, x0, n):
    b = ME.band(C, wa, wb, shape, keep)
    yy, xx = np.mgrid[0:n, 0:n].astype(float)
    gq = np.stack([(y0 + yy * STEP) / G, (x0 + xx * STEP) / G])
    m = ndi.map_coordinates(b.astype(np.uint8), gq, order=0, mode="constant", cval=0) > 0
    return paths(m)


def ribbon_figure():
    secs = {name: (lab, px) for name, lab, px in RE.sections()}
    name = "s5_G"
    lab, px = secs[name]
    g = RB.grid_labels(lab, px)
    rib = RB.describe(g)
    dE = rib["depth_egl"]
    best = None
    for b in rib["branches"]:
        if b["length_um"] < 700:
            continue
        S = M.smooth_path(b["path"])
        th, k, w, d = M.features(S, rib["width"], dE)
        ks = ndi.gaussian_filter1d(k, 3)
        L = len(S) * G
        sc = float(np.percentile(np.abs(ks), 95)) - abs(L - 1100) / 5000
        if best is None or sc > best[0]:
            best = (sc, S, ks, w, d)
    _, S, ks, w, d = best
    if len(S) * G > 1300:                    # keep the most curved 1.3 mm
        m = int(1300 / G)
        j = int(np.argmax(np.convolve(np.abs(ks), np.ones(m), "valid")))
        S, ks, w, d = S[j:j + m], ks[j:j + m], w[j:j + m], d[j:j + m]
    lo, hi = S.min(0) * G - 120, S.max(0) * G + 120
    size = float(max(hi - lo))
    c = (lo + hi) / 2
    y0, x0 = c - size / 2
    n = int(size / STEP)
    maps = crop_maps(lab, px, y0, x0, n)
    # iso-depth line at the stretch's median depth, inside tissue
    yy, xx = np.mgrid[0:n, 0:n].astype(float)
    gq = np.stack([(y0 + yy * STEP) / G, (x0 + xx * STEP) / G])
    dmap = ndi.map_coordinates(ndi.gaussian_filter(dE, 1.0), gq, order=1, mode="nearest")
    tis = ndi.map_coordinates(g["tissue"].astype(np.uint8), gq, order=0) > 0
    dmed = float(np.median(d))
    iso = []
    for cc in measure.find_contours(np.where(tis, dmap, np.nan), dmed):
        if len(cc) > 20:
            cc = measure.approximate_polygon(cc, 0.8)
            iso.append("M" + "L".join(f"{x:.1f},{y:.1f}" for y, x in cc))
    T = np.gradient(S, axis=0)
    T /= np.maximum(np.hypot(*T.T), 1e-9)[:, None]
    Nrm = np.stack([-T[:, 1], T[:, 0]], 1)
    secs_x = []
    for j in range(5, len(S) - 5, 10):
        half = w[j] / 2 / G
        secs_x.append(to_crop([S[j] - Nrm[j] * half, S[j] + Nrm[j] * half], y0, x0))
    s_um = np.arange(len(S)) * G
    return dict(section=name, n=n, maps=maps, iso="".join(iso), depth_med=round(dmed),
                line=to_crop(S, y0, x0), s=[round(float(v)) for v in s_um],
                kappa=[round(float(v) * 1000, 2) for v in ks], width=[round(float(v), 1) for v in w],
                depth=[round(float(v), 1) for v in d], cross=secs_x)


def accuracy(rows):
    bins = [(50, 100), (100, 200), (200, 300), (300, 400), (400, 501)]
    test = [r for r in rows if r["split"] == "test"]
    models = [("model", BEST), ("old", "walk"), ("hermite", "hermite"), ("line", "line"), ("elastic", ELASTIC)]
    out = dict(bins=[f"{a}–{min(b, 500)}" for a, b in bins], mid=[(a + min(b, 500)) / 2 for a, b in bins], series={}, n=[], within20={})
    for key, m in models:
        med, w20 = [], []
        for a, b in bins:
            e = np.array([np.inf if r["models"][m]["err_um"] is None else r["models"][m]["err_um"]
                          for r in test if a <= r["length_um"] < b])
            med.append(round(float(np.median(e)), 1))
            w20.append(round(float((e <= 20).mean() * 100)))
        out["series"][key] = med
        out["within20"][key] = w20
    out["n"] = [sum(a <= r["length_um"] < b for r in test) for a, b in bins]
    out["failed_old"] = [sum(a <= r["length_um"] < b and r["models"]["walk"]["err_um"] is None for r in test) for a, b in bins]
    return out


def pick_examples(rows):
    test = [(i, r) for i, r in enumerate(rows) if r["split"] == "test"]
    err = lambda r, m: np.inf if r["models"][m]["err_um"] is None else r["models"][m]["err_um"]
    chosen, used = [], set()

    def take(label, cond, key):
        c = [(i, r) for i, r in test if cond(r) and i not in used]
        if c:
            i, r = min(c, key=lambda t: key(t[1]))
            used.add(i)
            chosen.append((label, i))
    take("Short gap", lambda r: 90 <= r["length_um"] <= 160, lambda r: abs(err(r, BEST) - 2))
    take("Gap on a bend", lambda r: 200 <= r["length_um"] <= 320 and r["turn_deg"] >= 60, lambda r: err(r, BEST))
    take("Folium tip removed", lambda r: r["length_um"] >= 350 and r["turn_deg"] >= 140, lambda r: err(r, BEST))
    take("Long gap that bows out", lambda r: r["length_um"] >= 380 and r["turn_deg"] <= 45, lambda r: err(r, BEST))
    take("Where the model fails", lambda r: True, lambda r: -err(r, BEST))
    return chosen


def example(sec, lab, px, gp, label, rows_r):
    row, ex = ME.evaluate_gap(sec, gp, {})          # no elastic sweep needed here
    pts = np.vstack([ex["T"]] + [v for v in ex["preds"].values()])
    lo, hi = pts.min(0) * G - 140, pts.max(0) * G + 140
    size = float(max(hi - lo))
    c = (lo + hi) / 2
    y0, x0 = c - size / 2
    n = int(size / STEP)
    maps = crop_maps(lab, px, y0, x0, n, ex["gap"])
    curves = {"truth": to_crop(ex["T"], y0, x0)}
    for key, m in (("model", BEST), ("old", "walk"), ("hermite", "hermite"), ("line", "line")):
        if m in ex["preds"]:
            curves[key] = to_crop(ex["preds"][m], y0, x0)
    shape = ex["gap"].shape
    bands = {"model": band_path(ex["preds"][BEST], ex["wa"], ex["wb"], shape, ex["keep"], y0, x0, n)}
    if "walk" in ex["preds"]:
        bands["old"] = band_path(ex["preds"]["walk"], ex["wa"], ex["wb"], shape, ex["keep"], y0, x0, n)
    metr = {k: row["models"][m] for k, m in (("model", BEST), ("old", "walk"), ("hermite", "hermite"), ("line", "line"))}
    return dict(label=label, section=sec["name"], length_um=row["length_um"], turn_deg=row["turn_deg"],
                chord_um=row["chord_um"], width_ends=row["width_ends"], depth_ends=row["depth_ends"], n=n,
                maps=maps, curves=curves, bands=bands, metrics=metr,
                ends=[to_crop([ex["a"]], y0, x0)[0], to_crop([ex["b"]], y0, x0)[0]]), ex, (y0, x0, n)


def explorer(sec, ex, crop):
    y0, x0, n = crop
    out = []
    for sd in (5.0, 15.0, 40.0, np.inf):
        for bs in (0.1, 0.3, 1.0, 3.0):
            prm = dict(ME.FULL[BEST], sigma_depth_um=sd, bend_scale=bs)
            C, info = M.complete(ex["a"], ex["ta"], ex["b"], ex["tb_in"], sec["field"], ex["da"], ex["db"], prm, ex["avoid"])
            t = M.energy_terms(C, ex["a"], ex["ta"], ex["b"], ex["tb_in"], sec["field"], ex["da"], ex["db"], prm)
            e, _ = ME.curve_err(C, ex["T"])
            out.append(dict(sigma=None if not np.isfinite(sd) else sd, bend=bs, start=info["start"],
                            curve=to_crop(C, y0, x0), err_um=round(e, 1),
                            terms={k: round(v, 2) for k, v in t.items()}))
    h = M.hermite(ex["a"], ex["ta"], ex["b"], ex["tb_in"], 60)
    gd = M.geodesic(ex["a"], ex["b"], sec["field"], (ex["da"] + ex["db"]) / 2, dict(sigma_depth_um=15.0), ex["avoid"])
    return dict(grid=out, starts=dict(hermite=to_crop(h, y0, x0), geodesic=to_crop(gd, y0, x0) if gd is not None else None))


def main():
    bench = json.loads((OUT / "benchmark.json").read_text())
    rows = bench["gaps"]
    data = dict(stats=json.loads((ROOT / "outputs/igl_ribbon/ribbon_statistics.json").read_text()),
                params={k: (None if isinstance(v, float) and not np.isfinite(v) else v) for k, v in bench["summary"]["best_params"].items()},
                width=bench["summary"]["width"], depth_dev=bench["summary"]["depth_dev_from_linear_um"])
    data["ribbon"] = ribbon_figure()
    print("ribbon figure done", flush=True)
    data["accuracy"] = accuracy(rows)
    picks = pick_examples(rows)
    # rebuild each section's gaps exactly as the benchmark did (same seeds, same order)
    secs = list(RE.sections())
    exs, explore = [], None
    for label, idx in picks:
        r = rows[idx]
        si = [s[0] for s in secs].index(r["section"])
        name, lab, px = secs[si]
        sec = ME.section_data(name, lab, px)
        gaps = ME.make_gaps(sec, np.random.default_rng(200 + si))
        gi = idx - next(i for i, q in enumerate(rows) if q["section"] == name)
        e, ex, crop = example(sec, lab, px, gaps[gi], label, r)
        exs.append(e)
        print(label, name, e["length_um"], {k: v["err_um"] for k, v in e["metrics"].items()}, flush=True)
        if label == "Folium tip removed":
            explore = dict(index=len(exs) - 1, **explorer(sec, ex, crop))
    data["examples"] = exs
    data["explorer"] = explore
    (OUT / "explainer_data.json").write_text(json.dumps(data, default=float))
    tpl = (Path(__file__).parent / "igl_model_explainer_template.html").read_text()
    p = OUT / "igl_model_explainer.html"
    p.write_text(tpl.replace("__DATA__", json.dumps(data, default=float)))
    print(p, f"{p.stat().st_size / 1e6:.2f} MB")


if __name__ == "__main__":
    main()
