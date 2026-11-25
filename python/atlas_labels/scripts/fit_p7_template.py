"""
Fit the Allen P7 Nissl cerebellum template (outputs/allen_p7, made by
allen_nissl_template.py) to one of the user's P7 immunofluorescence sections.

    python scripts/fit_p7_template.py [slide_id] [labels_suffix]

  slide_id       default 2018_05_22_s1_2_p27-0021
  labels_suffix  default _iter6b_native_labels.tif (MATLAB automatic labels)

Fit input:   <sandbox>/<slide><suffix>        automatic MATLAB labels
Scored on:   labelled/<slide>_labelled.tiff   hand-painted gold standard (if present)
Output:      outputs/p7_template_fit/<slide>/

Stages
  1. Both label maps -> template codes (0 BG, 3 EGL, 5 ML, 4 IGL, 6 WM; 1 = tissue
     of unknown layer; 9 = anomaly, tear/bleach, left out of every fit).
     PCL -> IGL, DCN -> WM, iEGL+oEGL -> EGL. The gold standard has IGL/ML swapped
     (PROJECT_HISTORY section 4) and is swapped back here.
     IF labels are block-averaged 8x to class-fraction maps at 4.10 um/px.
  2. Orientation + mirror: midline template, every 5 deg x mirror x 3 scales.
  3. Every template section: similarity fit (rotation, scale, shift) by Powell on
     soft Dice (tissue + 4 layers) -> best section ~ medial-lateral position.
  4. Best section: B-spline warp (SimpleITK) of the depth-coded label image
     (BG 0, EGL 1, ML 2, IGL 3, WM 4); anomaly pixels masked out of the metric.
     Three control-point spacings are compared.
  5. Hold-out test: hide disks of the fit input, refit, and score the template's
     prediction inside the disks against the gold standard, next to a
     nearest-label fill from the surrounding automatic labels.
"""
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import SimpleITK as sitk
import tifffile
from scipy import ndimage as ndi
from scipy import optimize

ROOT = Path(__file__).resolve().parents[1]
SEG = Path("/Users/jpurzner/Dropbox/imaging_analysis/cerebellar_segmentation/python")
SANDBOX = SEG / "matlab_anomaly/sandbox"
TPL_DIR = ROOT / "outputs/allen_p7"
OUT_ROOT = ROOT / "outputs/p7_template_fit"

IF_PX = 0.5119            # 20x, um/px
BLOCK = 8
GT_PX = IF_PX * BLOCK     # 4.095 um/px working grid
GT_ORG = (BLOCK - 1) / 2 * IF_PX
TPL_PX = 3.96
TPL_MID = 218

UNK, ANOM = 1, 9
LAYERS = (3, 5, 4, 6)
NAMES = {0: "BG", 3: "EGL", 5: "ML", 4: "IGL", 6: "WM"}
DEPTH = np.array([0.0, 1.0, 2.0, 3.0, 4.0], np.float32)   # BG, EGL, ML, IGL, WM
CH = (0, 3, 5, 4, 6, UNK, ANOM)                            # fraction-stack channels
# MATLAB raw labels: 1 cereb, 2 iEGL, 3 oEGL, 4 IGL, 5 ML, 6 DWL, 7 PCL, 8 DCN, 9 anomaly
AUTO_LUT = {0: 0, 1: UNK, 2: 3, 3: 3, 4: 4, 5: 5, 6: 6, 7: 4, 8: 6, 9: ANOM}
# gold standard: same scheme but 4 = ML, 5 = IGL (napari swap baked in)
GOLD_LUT = {0: 0, 1: UNK, 2: 3, 3: 3, 4: 5, 5: 4, 6: 6, 7: 4, 8: 6, 9: ANOM}

BSPLINE_UM = (400.0, 200.0, 100.0)    # final control-point spacing variants
HOLD_R_UM, HOLD_N, HOLD_SEEDS = 100.0, 10, (0, 1, 2)


# --------------------------------------------------------------------------- #
# grids: (stack (C,H,W), px um, origin um).  physical = index * px + origin
# --------------------------------------------------------------------------- #
def remap(a, lut):
    t = np.zeros(max(int(a.max()) + 1, max(lut) + 1), np.uint8)
    for k, v in lut.items():
        t[k] = v
    return t[a]


def fractions(lab, block, chans=CH):
    H, W = lab.shape
    h, w = H // block, W // block
    lab = lab[:h * block, :w * block]
    out = np.empty((len(chans), h, w), np.float32)
    for i, c in enumerate(chans):
        out[i] = (lab == c).reshape(h, block, w, block).mean((1, 3))
    return out


def coarsen(grid, k):
    st, px, org = grid
    if k == 1:
        return grid
    C, H, W = st.shape
    h, w = H // k, W // k
    st = st[:, :h * k, :w * k].reshape(C, h, k, w, k).mean((2, 4))
    return st, px * k, org + (k - 1) / 2 * px


def centroid_um(tissue, px, org):
    r, c = np.indices(tissue.shape)
    m = tissue.sum()
    return np.array([(r * tissue).sum() / m, (c * tissue).sum() / m]) * px + org


def lin(params, mirror):
    th, logs = params[0], params[1]
    c, s = np.cos(th), np.sin(th)
    M = np.diag([1.0, -1.0]) if mirror else np.eye(2)
    return np.exp(logs) * np.array([[c, -s], [s, c]]) @ M


def warp(src, tgt_shape, tgt_px, tgt_org, params, mirror, ct, cs, order=1):
    """Pull template layer channels into target space.
    template_um = L (target_um - ct) + cs + d"""
    st, spx, sorg = src
    L = lin(params, mirror)
    d = np.array(params[2:4])
    A = L * (tgt_px / spx)
    o = (L @ (np.array([tgt_org, tgt_org]) - ct) + cs + d - sorg) / spx
    return np.stack([ndi.affine_transform(ch, A, o, output_shape=tgt_shape,
                                          order=order, cval=0.0) for ch in st])


def soft_dice(p, q, v):
    a = (v * p).sum() + (v * q).sum()
    return 2 * (v * p * q).sum() / a if a > 0 else 1.0


def score_parts(tgt, w, v):
    """tgt: (7,H,W) target fractions; w: (4,H,W) warped template layers; v: weights."""
    pt = tgt[1:6].sum(0)
    qt = w.sum(0)
    parts = [soft_dice(pt, qt, v)] + [soft_dice(tgt[1 + i], w[i], v) for i in range(4)]
    return np.array(parts)


# --------------------------------------------------------------------------- #
# data
# --------------------------------------------------------------------------- #
def load_target(slide, suffix):
    auto = remap(tifffile.imread(SANDBOX / f"{slide}{suffix}"), AUTO_LUT)
    gp = SEG / "labelled" / f"{slide}_labelled.tiff"
    gold = remap(tifffile.imread(gp).astype(np.int64), GOLD_LUT) if gp.exists() else None
    return auto, gold


def tpl_sections():
    st = json.loads((TPL_DIR / "p7_nissl_stats.json").read_text())
    return [(r["sec"], r["offset_um"]) for r in st]


def tpl_grid(sec):
    lab = tifffile.imread(TPL_DIR / f"sections/sec{sec:03d}_labels.tif")
    st = np.stack([(lab == c).astype(np.float32) for c in LAYERS])
    return st, TPL_PX, 0.0


# --------------------------------------------------------------------------- #
# similarity fit
# --------------------------------------------------------------------------- #
_W = {}


def _init(tgt, v):
    _W["tgt"], _W["v"] = tgt, v


def fit_similarity(sec, th0, mirror, levels=(4, 2), th_span=30, tgt=None, v=None, x0=None):
    tgt = _W["tgt"] if tgt is None else tgt
    v = _W["v"] if v is None else v
    src0 = tpl_grid(sec)
    T0 = (tgt, GT_PX, GT_ORG)
    ct = centroid_um(tgt[1:6].sum(0) * v, GT_PX, GT_ORG)
    cs = centroid_um(src0[0].sum(0), TPL_PX, 0.0)
    s0 = np.log(np.sqrt(src0[0].sum() * TPL_PX ** 2 / ((tgt[1:6].sum(0)).sum() * GT_PX ** 2)))
    best = None
    if x0 is None:
        tg, tpx, torg = coarsen(T0, levels[0])
        vg = coarsen((v[None], GT_PX, GT_ORG), levels[0])[0][0]
        src = coarsen(src0, levels[0])
        for th in np.deg2rad(th0 + np.arange(-th_span, th_span + 1, 5)):
            for ls in s0 + np.array([-0.1, 0.0, 0.1]):
                p = np.array([th, ls, 0.0, 0.0])
                sc = score_parts(tg, warp(src, tg.shape[1:], tpx, torg, p, mirror, ct, cs), vg).mean()
                if best is None or sc > best[0]:
                    best = (sc, p)
        x = best[1]
    else:
        x = np.array(x0, float)
    scl = np.array([0.1, 0.1, 100.0, 100.0])
    for k in levels[1:] if x0 is None else levels:
        tg, tpx, torg = coarsen(T0, k)
        vg = coarsen((v[None], GT_PX, GT_ORG), k)[0][0]
        src = coarsen(src0, k)

        def f(z):
            return -score_parts(tg, warp(src, tg.shape[1:], tpx, torg, z * scl, mirror, ct, cs), vg).mean()

        r = optimize.minimize(f, x / scl, method="Powell",
                              options=dict(xtol=1e-3, ftol=1e-5, maxfev=600))
        x = r.x * scl
    tg = coarsen(T0, k)[0]
    vg = coarsen((v[None], GT_PX, GT_ORG), k)[0][0]
    parts = score_parts(tg, warp(src, tg.shape[1:], GT_PX * k, coarsen(T0, k)[2], x, mirror, ct, cs), vg)
    return dict(sec=sec, params=x.tolist(), mirror=bool(mirror), ct=ct.tolist(), cs=cs.tolist(),
                score=float(parts.mean()), parts=parts.round(4).tolist(), level=k)


def _fit_task(a):
    return fit_similarity(*a)


def orientation_sweep(tgt, v, sec=TPL_MID):
    src0 = tpl_grid(sec)
    ct = centroid_um(tgt[1:6].sum(0) * v, GT_PX, GT_ORG)
    cs = centroid_um(src0[0].sum(0), TPL_PX, 0.0)
    s0 = np.log(np.sqrt(src0[0].sum() * TPL_PX ** 2 / (tgt[1:6].sum(0).sum() * GT_PX ** 2)))
    tg, tpx, torg = coarsen((tgt, GT_PX, GT_ORG), 4)
    vg = coarsen((v[None], GT_PX, GT_ORG), 4)[0][0]
    src = coarsen(src0, 4)
    res = []
    for mirror in (False, True):
        for th in range(0, 360, 5):
            for ls in s0 + np.array([-0.1, 0.0, 0.1]):
                p = np.array([np.deg2rad(th), ls, 0.0, 0.0])
                sc = score_parts(tg, warp(src, tg.shape[1:], tpx, torg, p, mirror, ct, cs), vg).mean()
                res.append((sc, th, mirror, ls))
    res.sort(key=lambda r: -r[0])
    return res


# --------------------------------------------------------------------------- #
# B-spline
# --------------------------------------------------------------------------- #
def encode(fr_layers):
    """depth code from (4,H,W) layer fractions (BG = remainder)."""
    bg = np.clip(1 - fr_layers.sum(0), 0, 1)
    tot = bg + fr_layers.sum(0)
    return (DEPTH[1:, None, None] * fr_layers).sum(0) / np.maximum(tot, 1e-6)


def to_img(a):
    im = sitk.GetImageFromArray(np.ascontiguousarray(a, np.float32))
    im.SetSpacing((GT_PX, GT_PX))
    im.SetOrigin((GT_ORG, GT_ORG))
    return im


def fit_bspline(tgt, wsim, mask, final_um):
    fixed_np = encode(np.stack([tgt[i] for i in (1, 2, 3, 4)]) /
                      np.maximum(tgt[0:5].sum(0), 1e-6)[None])
    fixed, moving = to_img(fixed_np), to_img(encode(wsim))
    m = sitk.GetImageFromArray(mask.astype(np.uint8))
    m.CopyInformation(fixed)
    size_um = np.array(fixed_np.shape) * GT_PX
    mesh0 = [max(1, int(round(s / (final_um * 4)))) for s in size_um[::-1]]
    tx = sitk.BSplineTransformInitializer(fixed, mesh0, 3)
    R = sitk.ImageRegistrationMethod()
    R.SetMetricAsMeanSquares()
    R.SetMetricFixedMask(m)
    R.SetMetricSamplingStrategy(R.NONE)
    R.SetInterpolator(sitk.sitkLinear)
    # (LBFGSB cannot change its parameter count between B-spline levels)
    R.SetOptimizerAsGradientDescentLineSearch(5.0, 100, convergenceMinimumValue=1e-5,
                                              convergenceWindowSize=5)
    R.SetInitialTransformAsBSpline(tx, inPlace=True, scaleFactors=[1, 2, 4])
    R.SetShrinkFactorsPerLevel([4, 2, 1])
    R.SetSmoothingSigmasPerLevel([4, 2, 1])
    R.SmoothingSigmasAreSpecifiedInPhysicalUnitsOff()
    out = R.Execute(fixed, moving)
    return out, R.GetMetricValue()


def apply_tx(tx, layers, ref):
    return np.stack([sitk.GetArrayFromImage(sitk.Resample(to_img(ch), ref, tx, sitk.sitkLinear, 0.0))
                     for ch in layers])


def argmax_labels(layers):
    full = np.concatenate([np.clip(1 - layers.sum(0), 0, 1)[None], layers])
    return np.array([0, 3, 5, 4, 6], np.uint8)[full.argmax(0)]


# --------------------------------------------------------------------------- #
# scoring
# --------------------------------------------------------------------------- #
def dice_table(pred, ref, ignore):
    out = {}
    ok = ~ignore
    x, y = (pred > 0) & ok, (ref > 0) & ok
    out["tissue"] = 2 * (x & y).sum() / max(x.sum() + y.sum(), 1)
    ok2 = ok & (ref != UNK)
    for c in LAYERS:
        x, y = (pred == c) & ok2, (ref == c) & ok2
        out[NAMES[c]] = 2 * (x & y).sum() / max(x.sum() + y.sum(), 1)
    return {k: round(float(v), 3) for k, v in out.items()}


def grid_labels(fr):
    """argmax label on the working grid from a 7-channel fraction stack."""
    return np.array(CH, np.uint8)[fr.argmax(0)]


def holdout_disks(gold_lab, valid, seed, n=HOLD_N, r_um=HOLD_R_UM):
    rng = np.random.default_rng(seed)
    r = r_um / GT_PX
    H, W = gold_lab.shape
    cand = np.argwhere((gold_lab > 0) & (gold_lab != UNK) & valid)
    cand = cand[(cand[:, 0] > r) & (cand[:, 0] < H - r) & (cand[:, 1] > r) & (cand[:, 1] < W - r)]
    cen = []
    for i in rng.permutation(len(cand)):
        p = cand[i]
        if all(np.hypot(*(p - q)) > 2.5 * r for q in cen):
            cen.append(p)
        if len(cen) == n:
            break
    yy, xx = np.indices((H, W))
    hid = np.zeros((H, W), bool)
    for cy, cx in cen:
        hid |= (yy - cy) ** 2 + (xx - cx) ** 2 <= r * r
    return hid, np.array(cen)


# --------------------------------------------------------------------------- #
def main(slide="2018_05_22_s1_2_p27-0021", suffix="_iter6b_native_labels.tif"):
    t0 = time.time()
    out = OUT_ROOT / slide
    out.mkdir(parents=True, exist_ok=True)
    auto, gold = load_target(slide, suffix)
    native_shape = auto.shape
    tgt = fractions(auto, BLOCK)
    v = 1.0 - tgt[CH.index(ANOM)]
    auto_g = grid_labels(tgt)
    gold_fr = fractions(gold, BLOCK) if gold is not None else None
    gold_g = grid_labels(gold_fr) if gold is not None else None
    tissue_mm2 = float(tgt[1:6].sum() * GT_PX ** 2 / 1e6)
    print(f"{slide}: native {native_shape}, grid {tgt.shape[1:]} at {GT_PX:.3f} um/px, "
          f"tissue {tissue_mm2:.2f} mm2, anomaly {tgt[6].sum() * GT_PX ** 2 / 1e6:.3f} mm2", flush=True)
    for c in LAYERS:
        print(f"  {NAMES[c]:4s} auto {tgt[CH.index(c)].sum() * GT_PX ** 2 / 1e6:.2f} mm2"
              + (f"   gold {gold_fr[CH.index(c)].sum() * GT_PX ** 2 / 1e6:.2f} mm2" if gold is not None else ""))

    # 2. orientation
    sw = orientation_sweep(tgt, v)
    th0, mirror = sw[0][1], sw[0][2]
    bm = {m: max(r[0] for r in sw if r[2] == m) for m in (False, True)}
    print(f"orientation (midline template): rot {th0} deg, mirror {mirror}; "
          f"best score mirror=False {bm[False]:.3f}, True {bm[True]:.3f}; "
          f"runner-up angles {[r[1] for r in sw[1:6]]}", flush=True)

    # 3. all sections, similarity fit
    secs = tpl_sections()
    cache = out / "similarity_fits_cache.json"
    if cache.exists() and json.loads(cache.read_text()).get("key") == [th0, bool(mirror), suffix]:
        fits = json.loads(cache.read_text())["fits"]
        print("similarity fits loaded from cache", flush=True)
    else:
        with Pool(8, initializer=_init, initargs=(tgt, v)) as pool:
            fits = pool.map(_fit_task, [(s, th0, mirror) for s, _ in secs], chunksize=4)
        cache.write_text(json.dumps(dict(key=[th0, bool(mirror), suffix], fits=fits)))
    off = dict(secs)
    for f in fits:
        f["offset_um"] = off[f["sec"]]
    fits.sort(key=lambda f: f["sec"])
    best8 = sorted(fits, key=lambda f: -f["score"])[:5]
    print("similarity fit, 8 um grid, top 5 template sections:", flush=True)
    for f in best8:
        print(f"  sec {f['sec']} offset {f['offset_um']:+.0f} um  score {f['score']:.3f}  "
              f"rot {np.rad2deg(f['params'][0]):.1f} deg  scale {np.exp(f['params'][1]):.3f}  "
              f"parts(tissue,EGL,ML,IGL,WM) {f['parts']}", flush=True)
    # refine the top 3 on the 4 um grid
    fin = [fit_similarity(f["sec"], 0, f["mirror"], levels=(1,), tgt=tgt, v=v, x0=f["params"])
           for f in best8[:3]]
    for f in fin:
        f["offset_um"] = off[f["sec"]]
    fin.sort(key=lambda f: -f["score"])
    best = fin[0]
    print("4 um refine:", [(f["sec"], f["offset_um"], round(f["score"], 4)) for f in fin], flush=True)

    # warped template (similarity) on the working grid
    G = tgt.shape[1:]
    src = tpl_grid(best["sec"])
    wsim = warp(src, G, GT_PX, GT_ORG, np.array(best["params"]), best["mirror"],
                np.array(best["ct"]), np.array(best["cs"]))
    lab_sim = argmax_labels(wsim)

    # 4. B-spline variants
    mask = (v > 0.5) & (tgt[CH.index(UNK)] < 0.5)
    ref = to_img(np.zeros(G, np.float32))
    variants = {}
    for um in BSPLINE_UM:
        t1 = time.time()
        tx, mv = fit_bspline(tgt, wsim, mask, um)
        wb = apply_tx(tx, wsim, ref)
        variants[um] = dict(tx=tx, layers=wb, labels=argmax_labels(wb), metric=mv)
        print(f"B-spline {um:.0f} um: metric {mv:.4f}  ({time.time() - t1:.0f} s)", flush=True)

    # 5a. scores vs gold
    report = dict(slide=slide, fit_input=f"{slide}{suffix}", grid_um=GT_PX, tissue_mm2=tissue_mm2,
                  orientation=dict(rot_deg=th0, mirror=bool(mirror), score_by_mirror={str(k): v_ for k, v_ in bm.items()}),
                  sections=fits, best=best, top_refined=fin)
    if gold is not None:
        ign = tgt[CH.index(ANOM)] > 0.5
        report["dice_vs_gold"] = {
            "auto_segmentation": dice_table(auto_g, gold_g, ign),
            "template_similarity": dice_table(lab_sim, gold_g, ign),
            **{f"template_bspline_{um:.0f}um": dice_table(variants[um]["labels"], gold_g, ign)
               for um in BSPLINE_UM}}
        print("Dice vs gold (4 um grid, anomaly px ignored):", flush=True)
        for k, d in report["dice_vs_gold"].items():
            print(f"  {k:28s} {d}", flush=True)

        # 5b. hold-out
        valid = v > 0.5
        rows = []
        for seed in HOLD_SEEDS:
            hid, cen = holdout_disks(gold_g, valid, seed)
            v_h = v * (~hid)
            fh = fit_similarity(best["sec"], 0, best["mirror"], levels=(1,), tgt=tgt, v=v_h,
                                x0=best["params"])
            wsh = warp(src, G, GT_PX, GT_ORG, np.array(fh["params"]), fh["mirror"],
                       np.array(fh["ct"]), np.array(fh["cs"]))
            preds = {"template_similarity": argmax_labels(wsh)}
            for um in BSPLINE_UM:
                txh, _ = fit_bspline(tgt, wsh, mask & ~hid, um)
                preds[f"template_bspline_{um:.0f}um"] = argmax_labels(apply_tx(txh, wsh, ref))
            known = ~hid & (auto_g != ANOM)
            idx = ndi.distance_transform_edt(~known, return_distances=False, return_indices=True)
            preds["nearest_label_fill"] = auto_g[tuple(idx)]
            preds["auto_seen (upper ref)"] = auto_g
            sel = hid & (gold_g != UNK)
            row = {"seed": seed, "n_px": int(sel.sum())}
            for k, p in preds.items():
                row[k] = dict(acc=round(float((p[sel] == gold_g[sel]).mean()), 3),
                              **dice_table(np.where(hid, p, 0), np.where(hid, gold_g, 0), ~hid))
            rows.append(row)
            print(f"hold-out seed {seed}: " + ", ".join(f"{k} acc {row[k]['acc']:.3f}"
                                                        for k in preds), flush=True)
        keys = [k for k in rows[0] if k not in ("seed", "n_px")]
        report["holdout"] = dict(radius_um=HOLD_R_UM, n_disks=HOLD_N, seeds=rows,
                                 mean_acc={k: round(float(np.mean([r[k]["acc"] for r in rows])), 3) for k in keys})
        print("hold-out mean accuracy inside hidden disks:", report["holdout"]["mean_acc"], flush=True)

    # pick the variant used for the saved native map: best hold-out accuracy among template variants
    if gold is not None:
        cand = {um: report["holdout"]["mean_acc"][f"template_bspline_{um:.0f}um"] for um in BSPLINE_UM}
        use_um = max(cand, key=cand.get)
    else:
        use_um = 200.0
    report["saved_variant_um"] = use_um
    print(f"saving native-resolution template labels with the {use_um:.0f} um B-spline", flush=True)

    # 6. native-resolution outputs
    nat = sitk.Image(int(native_shape[1]), int(native_shape[0]), sitk.sitkFloat32)
    nat.SetSpacing((IF_PX, IF_PX))
    best_v, best_lab = None, None
    codes = np.array([3, 5, 4, 6], np.uint8)
    tx = variants[use_um]["tx"]
    acc = None
    for i, ch in enumerate(wsim):
        a = sitk.GetArrayFromImage(sitk.Resample(to_img(ch), nat, tx, sitk.sitkLinear, 0.0))
        acc = a if acc is None else acc + a
        if best_v is None:
            best_v, best_lab = a, np.full(a.shape, codes[i], np.uint8)
        else:
            upd = a > best_v
            best_v[upd] = a[upd]
            best_lab[upd] = codes[i]
    native = np.where(acc > 0.5, best_lab, 0).astype(np.uint8)
    tifffile.imwrite(out / f"{slide}_p7template_labels.tif", native, compression="zlib")
    for um in BSPLINE_UM:
        variants[um]["layers"] = variants[um]["layers"].astype(np.float16)
    np.savez_compressed(out / f"{slide}_p7template_grid.npz",
                        grid_um=GT_PX, grid_origin_um=GT_ORG, auto=auto_g,
                        gold=gold_g if gold is not None else np.zeros(0),
                        similarity=lab_sim,
                        **{f"bspline_{um:.0f}um": variants[um]["labels"] for um in BSPLINE_UM})
    for um in BSPLINE_UM:
        sitk.WriteTransform(variants[um]["tx"], str(out / f"{slide}_bspline_{um:.0f}um.tfm"))
    (out / f"{slide}_fit.json").write_text(json.dumps(report, indent=1, default=float))
    print(f"native labels {native.shape}: " + ", ".join(
        f"{NAMES[c]} {(native == c).sum() * IF_PX ** 2 / 1e6:.2f} mm2" for c in LAYERS), flush=True)
    print(f"done in {time.time() - t0:.0f} s -> {out}", flush=True)
    return report


if __name__ == "__main__":
    main(*sys.argv[1:3])
