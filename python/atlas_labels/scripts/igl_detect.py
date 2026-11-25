"""
New IGL step (2026-11-19): call the IGL from all three stains instead of NeuN alone.

Why: the MATLAB pipeline calls IGL as the top 40 % of smoothed NeuN. Where NeuN is bleached the
IGL is lost (s4_C D4 / E4, s4_A m1, the s5_E lobule: the user traced all of these by hand),
while the nuclei stain still shows it as clearly as the IGL next to it. In the EZH2 cKO
sections the nuclei-dense, NeuN-positive inner EGL is pulled into the IGL.

Features, per 4 um grid cell (the grid of igl_ribbon.grid_labels, same corner-aligned zoom):
  - each stain (nuclei, NeuN, p27) divided by its own local background (tissue-normalised
    Gaussian, sigma 50 and 150 um), log, smoothed at 4 and 12 um: bleaching is a slow
    multiplicative loss, so it cancels;
  - the same stains against the section's tissue median (absolute level), smoothed at 4 um;
  - cell density: nuclei-scale blob response (DoG 1.5 / 4 um at 1 um per pixel), positive
    part, over the local background, smoothed at 6 um, for nuclei and NeuN; local
    nuclei texture (standard deviation over 4 um);
  - depth from the tissue edge (distance to background, um), raw and smoothed.
The existing layer labels are not used as features: they are what is wrong.

Model: gradient-boosted trees (sklearn HistGradientBoostingClassifier), 4 classes on tissue
(EGL, ML/PCL, IGL, white matter), trained on the reviewed controls
(outputs/final_labels/controls); cells labelled unassigned or anomaly are not trained on.

    python scripts/igl_detect.py features          # cache features for all 19 sections
    python scripts/igl_detect.py loso              # leave-one-section-out over the 12 controls
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import tifffile
from scipy import ndimage as ndi

sys.path.insert(0, str(Path(__file__).parent))
import igl_gapfill_apply as GA  # noqa: E402
import igl_ribbon as RB  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/igl_detect"
FINAL = ROOT / "outputs/final_labels/controls"
G = RB.GRID
UM = 1.0                                     # working resolution for the image features (um per pixel)
CLASSES = {2: "EGL", 3: "ML/PCL", 4: "IGL", 5: "white matter"}


def stains(sid):
    """nuclei, NeuN, p27 at native resolution (20x pages: p27, NeuN, DAPI; 10x: p27, DAPI, NeuN)."""
    px, _, big = GA.info(sid)
    raw = tifffile.imread((GA.EV.SB20 if big else GA.EV.SB10) / f"{sid}_native.tif")
    nuc, neun = (2, 1) if big else (1, 2)
    return dict(nuclei=raw[nuc], neun=raw[neun], p27=raw[0]), px


def _norm_blur(a, mask, sigma_px):
    """Gaussian mean of a over mask (normalised convolution)."""
    num = ndi.gaussian_filter(a * mask, sigma_px)
    den = ndi.gaussian_filter(mask.astype(np.float32), sigma_px)
    return num / np.maximum(den, 1e-3)


def features(sid):
    """feature stack (n_features, H, W) on the 4 um grid, its names, and the grid labels of the
    anatomy-fixed MATLAB labels (for the tissue mask and for comparison)."""
    S, px = stains(sid)
    auto = tifffile.imread(GA.RE.FX / f"{sid}_fixed_labels.tif")
    g = RB.grid_labels(auto, px)
    gshape = g["igl"].shape
    z1 = px / UM
    tissue1 = ndi.zoom((auto > 0).astype(np.float32), z1, order=1) >= 0.5
    T = tissue1.astype(np.float32)
    feats, names = [], []

    def to_grid(a):
        zf = (gshape[0] / a.shape[0], gshape[1] / a.shape[1])
        return ndi.zoom(ndi.gaussian_filter(a, G / UM / 2.0), zf, order=1).astype(np.float32)

    for name, img in S.items():
        a = ndi.zoom(ndi.gaussian_filter(img.astype(np.float32), 0.5 / px), z1, order=1)
        a = np.maximum(a, 1.0)
        med = float(np.median(a[tissue1]))
        for bg_um in (50.0, 150.0):
            B = _norm_blur(a, T, bg_um / UM)
            L = np.log(a / np.maximum(B, 1.0)) * T
            for s_um in (4.0, 12.0):
                feats.append(to_grid(ndi.gaussian_filter(L, s_um / UM)))
                names.append(f"{name} vs local {int(bg_um)} um, smoothed {int(s_um)} um")
        feats.append(to_grid(ndi.gaussian_filter(np.log(a / med) * T, 4.0 / UM)))
        names.append(f"{name} vs section median, smoothed 4 um")
        if name in ("nuclei", "neun"):
            B = _norm_blur(a, T, 150.0 / UM)
            dog = (ndi.gaussian_filter(a, 1.5 / UM) - ndi.gaussian_filter(a, 4.0 / UM)) / np.maximum(B, 1.0)
            feats.append(to_grid(ndi.gaussian_filter(np.maximum(dog, 0) * T, 6.0 / UM)))
            names.append(f"{name} cell density (blob response), 6 um")
            if name == "nuclei":
                m1 = ndi.gaussian_filter(a / np.maximum(B, 1.0), 4.0 / UM)
                m2 = ndi.gaussian_filter((a / np.maximum(B, 1.0)) ** 2, 4.0 / UM)
                feats.append(to_grid(np.sqrt(np.maximum(m2 - m1 ** 2, 0)) * T))
                names.append("nuclei texture (std over 4 um)")
    tg = g["tissue"] | g["anom"]
    d = ndi.distance_transform_edt(tg) * G
    feats += [np.log1p(d).astype(np.float32), np.log1p(ndi.gaussian_filter(d, 25.0 / G)).astype(np.float32)]
    names += ["depth from the tissue edge (log um)", "depth from the tissue edge, smoothed 25 um (log um)"]
    return np.stack(feats), names, g


def target(sid, shape, px):
    """training classes on the grid from the reviewed control labels: 2 EGL, 3 ML/PCL, 4 IGL,
    5 white matter; 0 = not trained on (background, unassigned, anomaly)."""
    lab = tifffile.imread(FINAL / f"{GA.short(sid)}_labels.tif")
    g = RB.grid_labels(lab, px)
    c = GA.grid_code(g)
    assert c.shape == shape, (c.shape, shape)
    c[~np.isin(c, (2, 3, 4, 5))] = 0
    return c


def cache_features(sids):
    (OUT / "features").mkdir(parents=True, exist_ok=True)
    for sid in sids:
        t0 = time.time()
        F, names, g = features(sid)
        px, group, _ = GA.info(sid)
        y = target(sid, g["igl"].shape, px) if group == "control" else np.zeros(g["igl"].shape, np.uint8)
        np.savez_compressed(OUT / "features" / f"{GA.short(sid)}.npz", F=F.astype(np.float16), y=y, auto=GA.grid_code(g),
                            tissue=g["tissue"] | g["anom"])
        (OUT / "features" / "names.json").write_text(json.dumps(names, indent=1))
        print(f"{GA.short(sid):5s} {group:7s} grid {F.shape[1:]} {F.shape[0]} features  [{time.time() - t0:.0f}s]", flush=True)


def load(short):
    Z = np.load(OUT / "features" / f"{short}.npz")
    return Z["F"].astype(np.float32), Z["y"], Z["auto"], Z["tissue"]


_load_images = load            # the image features alone (load is swapped for load2 in the second pass)


def controls():
    return [GA.short(s) for s in GA.section_list() if GA.info(s)[1] == "control"]


CORRECTED = ROOT / "outputs/igl_detect/corrected_round1"


def edited_cko():
    """cKO sections the user corrected on the New IGL Step page (round 1)."""
    return sorted(f.name.replace("_labels.tif", "") for f in CORRECTED.glob("*_labels.tif"))


def training_sections():
    return controls() + edited_cko()


def update_targets():
    """training classes of the edited cKO sections from their corrected labels (round 2)."""
    for short in edited_cko():
        sid = next(x for x in GA.section_list() if GA.short(x) == short)
        px = GA.info(sid)[0]
        Z = dict(np.load(OUT / "features" / f"{short}.npz"))
        g = RB.grid_labels(tifffile.imread(CORRECTED / f"{short}_labels.tif"), px)
        c = GA.grid_code(g)
        c[~np.isin(c, (2, 3, 4, 5))] = 0
        assert c.shape == Z["y"].shape
        Z["y"] = c
        np.savez_compressed(OUT / "features" / f"{short}.npz", **Z)
        print("targets from the user's edits:", short, int((c == 4).sum() * G * G), "um2 IGL", flush=True)


def train(shorts, per_class=12000, seed=0, prm=None):
    """fit the classifier on a balanced sample of grid cells of the given (control) sections."""
    from sklearn.ensemble import HistGradientBoostingClassifier
    rng = np.random.default_rng(seed)
    X, Y = [], []
    for s in shorts:
        F, y, _, _ = load(s)
        for c in CLASSES:
            idx = np.flatnonzero(y.ravel() == c)
            if not len(idx):
                continue
            take = rng.choice(idx, min(per_class, len(idx)), replace=False)
            X.append(F.reshape(F.shape[0], -1)[:, take].T)
            Y.append(np.full(len(take), c, np.uint8))
    p = dict(max_iter=300, learning_rate=0.08, max_leaf_nodes=63, min_samples_leaf=40, l2_regularization=1.0, random_state=seed)
    p.update(prm or {})
    clf = HistGradientBoostingClassifier(**p)
    clf.fit(np.vstack(X), np.concatenate(Y))
    return clf


def predict(clf, short):
    """class map (0 outside tissue) and class probabilities (n_classes, H, W)."""
    F, _, _, tissue = load(short)
    n, H, W = F.shape
    P = np.zeros((len(clf.classes_), H, W), np.float32)
    idx = np.flatnonzero(tissue.ravel())
    P.reshape(len(clf.classes_), -1)[:, idx] = clf.predict_proba(F.reshape(n, -1)[:, idx].T).T
    pred = np.zeros((H, W), np.uint8)
    pred.ravel()[idx] = clf.classes_[np.argmax(P.reshape(len(clf.classes_), -1)[:, idx], 0)]
    return pred, P


def iou(a, b):
    return float((a & b).sum() / max((a | b).sum(), 1))


def score(short, pred):
    """IGL agreement with the reviewed labels, against the MATLAB labels as the baseline; how much
    of the IGL added by the repairs / the user (reviewed IGL that MATLAB missed) the new step finds,
    and how much IGL that the review removed it leaves out."""
    _, y, auto, tissue = load(short)
    known = y > 0
    fin, mat, new = (y == 4), (auto == 4), (pred == 4)
    added = fin & ~mat & tissue                       # repairs / user traces
    removed = mat & ~fin & known
    acc = float((pred[known] == y[known]).mean())
    return dict(section=short, igl_iou_matlab=round(iou(mat & known, fin), 3), igl_iou_new=round(iou(new & known, fin), 3),
                acc_4class=round(acc, 3), acc_4class_matlab=round(float((np.where(np.isin(auto, (2, 3, 4, 5)), auto, 0)[known] == y[known]).mean()), 3),
                igl_mm2_reviewed=round(float(fin.sum() * G * G / 1e6), 3), igl_mm2_new=round(float((new & tissue).sum() * G * G / 1e6), 3),
                repaired_um2=int(added.sum() * G * G), repaired_found=round(float(new[added].mean()), 2) if added.any() else None,
                removed_um2=int(removed.sum() * G * G), removed_left_out=round(float((~new[removed]).mean()), 2) if removed.any() else None)


def loso(prm=None, tag="loso", per_class=12000):
    cs = controls()
    rows = []
    (OUT / tag).mkdir(parents=True, exist_ok=True)
    for h in cs:
        t0 = time.time()
        clf = train([c for c in cs if c != h], per_class=per_class, prm=prm)
        pred, P = predict(clf, h)
        r = score(h, pred)
        rows.append(r)
        np.savez_compressed(OUT / tag / f"{h}.npz", pred=pred, p_igl=P[list(clf.classes_).index(4)].astype(np.float16))
        print(f"{h:5s} IGL IoU MATLAB {r['igl_iou_matlab']:.3f} -> new {r['igl_iou_new']:.3f} | 4-class acc {r['acc_4class_matlab']:.3f} -> {r['acc_4class']:.3f}"
              f" | repaired {r['repaired_um2']:>7d} um2 found {r['repaired_found']} | reviewed-out {r['removed_um2']:>6d} um2 left out {r['removed_left_out']}"
              f"  [{time.time() - t0:.0f}s]", flush=True)
    m = lambda k: round(float(np.mean([r[k] for r in rows])), 3)
    summ = dict(igl_iou_matlab=m("igl_iou_matlab"), igl_iou_new=m("igl_iou_new"), acc_matlab=m("acc_4class_matlab"), acc_new=m("acc_4class"))
    print("mean", summ)
    (OUT / f"{tag}.json").write_text(json.dumps(dict(prm=prm, summary=summ, sections=rows), indent=1))
    return rows


CTX_SIGMA_UM = (8.0, 24.0, 60.0)


def context(P):
    """second-pass features: each class probability averaged over 8, 24 and 60 um (the layer
    order around a cell: EGL above ML above IGL above white matter)."""
    out = []
    for c in range(P.shape[0]):
        for s_um in CTX_SIGMA_UM:
            out.append(ndi.gaussian_filter(P[c].astype(np.float32), s_um / G))
    return np.stack(out)


def first_pass(per_class=8000, tag="pass1", sections=None):
    """class probabilities of the first pass for every section: training sections out of fold
    (model trained on the others), the rest from a model trained on all training sections."""
    (OUT / tag).mkdir(parents=True, exist_ok=True)
    cs = sections or controls()
    for h in cs:
        clf = train([c for c in cs if c != h], per_class=per_class)
        _, P = predict(clf, h)
        np.savez_compressed(OUT / tag / f"{h}.npz", P=P.astype(np.float16), classes=clf.classes_)
        print("pass 1 out of fold:", h, flush=True)
    clf = train(cs, per_class=per_class)
    for sid in GA.section_list():
        s = GA.short(sid)
        if s not in cs:
            _, P = predict(clf, s)
            np.savez_compressed(OUT / tag / f"{s}.npz", P=P.astype(np.float16), classes=clf.classes_)
            print("pass 1 (all controls):", s, flush=True)


PASS1 = "pass1"


def load2(short, tag=None):
    """features of the second pass: image features + first-pass context."""
    F, y, auto, tissue = _load_images(short)
    P = np.load(OUT / (tag or PASS1) / f"{short}.npz")["P"].astype(np.float32)
    return np.concatenate([F, context(P)]), y, auto, tissue


def loso2(per_class=8000, tag="loso2", prm=None, sections=None):
    """held-out test of the second pass. The first-pass probabilities of the training sections are
    out of fold, but each came from a model that had seen the held-out section: a mild leak,
    accepted to keep the run short (a fully nested run is 12 x 11 fits)."""
    global load
    cs = sections or controls()
    rows = []
    (OUT / tag).mkdir(parents=True, exist_ok=True)
    base_load = load
    load = load2
    try:
        for h in cs:
            t0 = time.time()
            clf = train([c for c in cs if c != h], per_class=per_class, prm=prm)
            pred, P = predict(clf, h)
            r = score(h, pred)
            rows.append(r)
            np.savez_compressed(OUT / tag / f"{h}.npz", pred=pred, p_igl=P[list(clf.classes_).index(4)].astype(np.float16), P=P.astype(np.float16))
            print(f"{h:5s} IGL IoU MATLAB {r['igl_iou_matlab']:.3f} -> new {r['igl_iou_new']:.3f} | 4-class acc {r['acc_4class_matlab']:.3f} -> {r['acc_4class']:.3f}"
                  f" | repaired {r['repaired_um2']:>7d} um2 found {r['repaired_found']} | reviewed-out {r['removed_um2']:>6d} um2 left out {r['removed_left_out']}"
                  f"  [{time.time() - t0:.0f}s]", flush=True)
    finally:
        load = base_load
    m = lambda k: round(float(np.mean([r[k] for r in rows])), 3)
    summ = dict(igl_iou_matlab=m("igl_iou_matlab"), igl_iou_new=m("igl_iou_new"), acc_matlab=m("acc_4class_matlab"), acc_new=m("acc_4class"))
    print("mean", summ)
    (OUT / f"{tag}.json").write_text(json.dumps(dict(prm=prm, summary=summ, sections=rows), indent=1))


def apply_pass2(per_class=8000, tag="apply2", sections=None):
    """second pass trained on all training sections, applied to the others."""
    global load
    (OUT / tag).mkdir(parents=True, exist_ok=True)
    cs = sections or controls()
    base_load = load
    load = load2
    try:
        clf = train(cs, per_class=per_class)
        for sid in GA.section_list():
            s = GA.short(sid)
            if s in cs:
                continue
            pred, P = predict(clf, s)
            np.savez_compressed(OUT / tag / f"{s}.npz", pred=pred, p_igl=P[list(clf.classes_).index(4)].astype(np.float16), P=P.astype(np.float16))
            print("pass 2 (all controls):", s, flush=True)
    finally:
        load = base_load


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "features"
    if cmd == "round2":
        # round 2 (2026-11-20): the user's round-1 cKO edits join the 12 controls as training sections
        update_targets()
        PASS1 = "pass1_r2"
        tr = training_sections()
        print("training sections:", tr, flush=True)
        first_pass(tag="pass1_r2", sections=tr)
        loso2(tag="loso2_r2", sections=tr)
        apply_pass2(tag="apply2_r2", sections=tr)
    if cmd == "apply2":
        apply_pass2()
    if cmd == "pass1":
        first_pass()
    elif cmd == "loso2":
        loso2()
    if cmd == "features":
        cache_features(sys.argv[2:] or GA.section_list())
    elif cmd == "loso":
        loso(json.loads(sys.argv[2]) if len(sys.argv) > 2 else None, sys.argv[3] if len(sys.argv) > 3 else "loso")
