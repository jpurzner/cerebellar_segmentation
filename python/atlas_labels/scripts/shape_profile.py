"""
Generalized cerebellum shape profile (outline + white-matter tree), built from the
curated good shapes and grown as more sections are fixed and accepted.

    python scripts/shape_profile.py stage0

Stage 0 - prepare the shape set
  * IF sections (curation.json grades): labels -> 4 um grid; tissue, white matter
    (DWL + DCN), EGL. Base check: tissue outside the convex hull of the main EGL
    (beyond the lobule I-X line) is reported and cut.
    20x 0.5119 um/px; 10x 0.65 um/px (not the 1.0239 the MATLAB port assumes).
  * Allen P56 within +-0.6 mm (white matter checked against the atlas there).
  * All rotated into atlas orientation (anterior left, dorsal up) by the library's
    orientation search, then aligned to the s1_2 gold standard by a similarity
    transform (rotation, scale, shift) maximizing the overlap of tissue and white
    matter. Scale is recorded separately (shape vs size).
Output: outputs/shape_profile/stage0_shape_set.npz, stage0.json
"""
import json
import sys
from pathlib import Path

import numpy as np
import tifffile
from scipy import ndimage as ndi
from scipy import optimize
from skimage import morphology

sys.path.insert(0, str(Path(__file__).parent))
import fit_p7_template as F  # noqa: E402
import wm_library as WL  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/shape_profile"
LIB = ROOT / "outputs/wm_library"
GRID = 4.0
ST = np.ones((3, 3), bool)


def slide_px(sid):
    return 0.5119 if sid.startswith("2018_") else WL.PX_10X


def load_ds(sid):
    """labels -> tissue / WM / EGL at the 4 um grid (majority of each block)."""
    lab = tifffile.imread(WL.FIXED / f"{sid}_fixed_labels.tif")
    px = slide_px(sid)
    out = {}
    for name, codes in (("tissue", (1, 2, 3, 4, 5, 6, 7, 8)), ("wm", (6, 8)), ("egl", (2, 3))):
        m = np.isin(lab, codes).astype(np.float32)
        out[name] = ndi.zoom(m, px / GRID, order=1) >= 0.5
    return out


def base_cut(tissue, egl, max_frac=0.4):
    l, n = ndi.label(egl, structure=ST)
    if n == 0:
        return tissue, 0.0
    sz = np.bincount(l.ravel())[1:]
    main = np.isin(l, 1 + np.nonzero((sz >= 0.1 * sz.max()) & (sz >= 100))[0])
    hull = morphology.convex_hull_image(main)
    removed = tissue & ~hull
    if removed.sum() > max_frac * tissue.sum():
        return tissue, -1.0
    return tissue & hull, float(removed.sum() * GRID ** 2 / 1e6)


def align(src, ref, x0=None):
    """similarity transform pulling src (2,H,W) onto ref (2,Hr,Wr); returns params, warped, score."""
    sgrid = (src.astype(np.float32), GRID, 0.0)
    ct = F.centroid_um(ref[0].astype(float), GRID, 0.0)
    cs = F.centroid_um(src[0].astype(float), GRID, 0.0)
    s0 = np.log(np.sqrt(src[0].sum() / ref[0].sum()))
    v = np.ones(ref.shape[1:], np.float32)
    refs = ref.astype(np.float32)

    def score(p):
        w = F.warp(sgrid, ref.shape[1:], GRID, 0.0, p, False, ct, cs)
        return np.mean([F.soft_dice(refs[i], w[i], v) for i in range(2)])

    best = None
    for th in np.deg2rad([-20, -10, 0, 10, 20]):
        p = np.array([th, s0, 0.0, 0.0])
        sc = score(p)
        if best is None or sc > best[0]:
            best = (sc, p)
    scl = np.array([0.05, 0.05, 50.0, 50.0])
    r = optimize.minimize(lambda z: -score(z * scl), best[1] / scl, method="Powell",
                          options=dict(xtol=1e-3, ftol=1e-5, maxfev=800))
    p = r.x * scl
    w = F.warp(sgrid, ref.shape[1:], GRID, 0.0, p, False, ct, cs)
    d = [float(F.soft_dice(refs[i], w[i], v)) for i in range(2)]
    return p, w >= 0.5, d


def pad_to(m, shape):
    H, W = shape
    h, w = m.shape[-2:]
    out = np.zeros(m.shape[:-2] + (H, W), m.dtype)
    y0, x0 = (H - h) // 2, (W - w) // 2
    out[..., y0:y0 + h, x0:x0 + w] = m
    return out


def dice(a, b):
    return float(2 * (a & b).sum() / max(a.sum() + b.sum(), 1))


def align_tissue(src3, ref_t):
    """similarity transform (rotation, scale, shift) maximizing tissue overlap with ref_t,
    optimized on an 8 um grid; returns params and all channels of src3 warped at 4 um."""
    S = (src3.astype(np.float32), GRID, 0.0)
    ct = F.centroid_um(ref_t.astype(float), GRID, 0.0)
    cs = F.centroid_um(src3[0].astype(float), GRID, 0.0)
    s0 = np.log(np.sqrt(src3[0].sum() / ref_t.sum()))
    Sc = F.coarsen(S, 2)
    rt, rpx, rorg = F.coarsen((ref_t[None].astype(np.float32), GRID, 0.0), 2)
    v = np.ones(rt.shape[1:], np.float32)
    st = (Sc[0][:1], Sc[1], Sc[2])

    def score(p):
        return F.soft_dice(rt[0], F.warp(st, rt.shape[1:], rpx, rorg, p, False, ct, cs)[0], v)

    best = max(((score(np.array([th, s0, 0.0, 0.0])), np.array([th, s0, 0.0, 0.0]))
                for th in np.deg2rad([-20, -10, 0, 10, 20])), key=lambda r: r[0])[1]
    scl = np.array([0.05, 0.05, 50.0, 50.0])
    r = optimize.minimize(lambda z: -score(z * scl), best / scl, method="Powell",
                          options=dict(xtol=1e-3, ftol=1e-5, maxfev=800))
    p = r.x * scl
    W = F.warp(S, ref_t.shape, GRID, 0.0, p, False, ct, cs) >= 0.5
    return p, W


def stage0():
    """Control shape set (10x IF controls + Allen P56 within +-0.6 mm), tissue / WM / IGL,
    aligned on the tissue outline into a mean-shape frame (two rounds, started from s5_G).
    Per-feature membership follows the user's review (curation.json)."""
    OUT.mkdir(parents=True, exist_ok=True)
    cur = json.loads((LIB / "curation.json").read_text())
    wz = np.load(LIB / "library.npz")
    iz = np.load(ROOT / "outputs/igl_library/library.npz")
    lib = {e["id"]: e for e in json.loads((LIB / "library.json").read_text())["entries"]}
    shapes = []
    for sid in WL.SLIDES_10X:                       # 10x controls
        eid = f"if_{sid}"
        shapes.append(dict(id=eid, source="IF 10x control",
                           in_igl=cur["igl_grade"][eid] in ("good", "usable"),
                           in_wm=cur["wm_grade"][eid] in ("good", "usable"),
                           igl_grade=cur["igl_grade"][eid], wm_grade=cur["wm_grade"][eid]))
    # s1_2 is a control (the other 20x slides are EZH2 cKO); use its hand-painted gold labels
    g = "if_2018_05_22_s1_2_p27-0021"
    shapes.append(dict(id="gold_s1_2", source="IF 20x control (gold labels)",
                       in_igl=cur["igl_grade"][g] in ("good", "usable"), in_wm=cur["wm_grade"][g] in ("good", "usable"),
                       igl_grade=cur["igl_grade"][g], wm_grade=cur["wm_grade"][g]))
    for e in sorted(lib.values(), key=lambda e: e.get("offset_um", 0)):
        if e["source"] == "Allen P56" and abs(e["offset_um"]) <= cur["allen_shape_max_offset_um"]:
            shapes.append(dict(id=e["id"], source="Allen P56", offset_um=e["offset_um"],
                               in_igl=e["id"] in cur["p56_igl_good"], in_wm=e["id"] in cur["p56_wm_good"]))
    for s in shapes:
        s["stack"] = np.stack([wz[f"tissue_{s['id']}"], wz[f"wm_{s['id']}"], iz[f"igl_{s['id']}"]])
    seed = [s for s in shapes if s["id"] == "if_s5_G_p27"][0]["stack"][0]
    m = 80
    ref = np.pad(seed, m)
    is_if = np.array([s["source"].startswith("IF") for s in shapes])
    w = np.where(is_if, 0.5 / is_if.sum(), 0.5 / (~is_if).sum())
    for rnd in (1, 2):
        out = [align_tissue(s["stack"], ref) for s in shapes]
        A = np.stack([o[1] for o in out])
        mean_t = (w[:, None, None] * A[:, 0]).sum(0)
        print(f"round {rnd}: tissue Dice to frame, IF median "
              f"{np.median([dice(A[i, 0], ref) for i in np.nonzero(is_if)[0]]):.3f}, P56 median "
              f"{np.median([dice(A[i, 0], ref) for i in np.nonzero(~is_if)[0]]):.3f}", flush=True)
        ref = mean_t >= 0.5
    meta = []
    for s, (p, W) in zip(shapes, out):
        mm = {k: v for k, v in s.items() if k != "stack"}
        mm.update(scale_to_frame=round(float(np.exp(-p[1])), 3), rot_deg=round(float(np.rad2deg(p[0])), 1),
                  dice_tissue_to_mean=round(dice(W[0], ref), 3),
                  tissue_mm2=round(float(s["stack"][0].sum() * GRID ** 2 / 1e6), 3))
        meta.append(mm)
    inI = np.array([m_["in_igl"] for m_ in meta])
    inW = np.array([m_["in_wm"] for m_ in meta])

    def pmean(ch, sel):
        ww = np.where(is_if[sel], 0.5 / max(is_if[sel].sum(), 1), 0.5 / max((~is_if[sel]).sum(), 1))
        return (ww[:, None, None] * A[sel, ch]).sum(0).astype(np.float32)

    np.savez_compressed(OUT / "stage0_shape_set.npz", aligned=A, ids=np.array([m_["id"] for m_ in meta]),
                        ref_tissue=ref, mean_tissue=mean_t.astype(np.float32),
                        p_wm=pmean(1, inW), p_igl=pmean(2, inI),
                        p_wm_if=A[inW & is_if, 1].mean(0).astype(np.float32),
                        p_igl_if=A[inI & is_if, 2].mean(0).astype(np.float32))
    (OUT / "stage0.json").write_text(json.dumps(dict(
        grid_um=GRID, frame="mean control outline (2 rounds of alignment on tissue, seeded with s5_G)",
        excluded="20x slides except s1_2 (EZH2 cKO)", weights="IF and P56 groups 50/50",
        n=dict(outline=len(meta), outline_if=int(is_if.sum()), igl_if=int((inI & is_if).sum()),
               igl_p56=int((inI & ~is_if).sum()), wm_if=int((inW & is_if).sum()), wm_p56=int((inW & ~is_if).sum())),
        shapes=meta), indent=1, default=float))
    print(f"control shape set: {len(meta)} outlines ({is_if.sum()} IF, {(~is_if).sum()} P56); "
          f"IGL {int((inI & is_if).sum())} IF + {int((inI & ~is_if).sum())} P56; "
          f"WM {int((inW & is_if).sum())} IF + {int((inW & ~is_if).sum())} P56", flush=True)
    for mm in meta:
        if mm["source"].startswith("IF"):
            print(f"  {mm['id'].replace('if_', ''):10s} scale {mm['scale_to_frame']:.2f} rot {mm['rot_deg']:+5.1f} "
                  f"Dice to mean {mm['dice_tissue_to_mean']:.3f}  IGL {mm['igl_grade']:7s} WM {mm['wm_grade']}", flush=True)


# --------------------------------------------------------------------------- #
# Stage 1a - outline model (option A)
# --------------------------------------------------------------------------- #
ENV_CLOSE_UM = 40.0     # envelope: tissue closed over fissures, holes filled
SDF_CLIP_UM = 300.0


def envelope(tissue):
    r = int(round(ENV_CLOSE_UM / GRID))
    e = ndi.binary_closing(np.pad(tissue, r + 2), structure=morphology.disk(r))[r + 2:-r - 2, r + 2:-r - 2]
    e = ndi.binary_fill_holes(e)
    l, n = ndi.label(e, structure=ST)
    if n > 1:
        e = l == 1 + np.argmax(np.bincount(l.ravel())[1:])
    return e


def sdf(mask):
    """signed distance in um: negative inside, positive outside, clipped."""
    d = (ndi.distance_transform_edt(~mask) - ndi.distance_transform_edt(mask)) * GRID
    return np.clip(d, -SDF_CLIP_UM, SDF_CLIP_UM).astype(np.float32)


def pca(X, w):
    """weighted PCA of rows of X; returns mean, components (k, d), variances."""
    w = w / w.sum()
    mu = (w[:, None] * X).sum(0)
    Xc = (X - mu) * np.sqrt(w)[:, None]
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    return mu, Vt, S ** 2


def dice(a, b):
    return float(2 * (a & b).sum() / max(a.sum() + b.sum(), 1))


def stage1_outline():
    z = np.load(OUT / "stage0_shape_set.npz")
    meta = json.loads((OUT / "stage0.json").read_text())["shapes"]
    A = z["aligned"]
    ids = list(z["ids"])
    is_if = np.array([m["source"].startswith("IF") for m in meta])
    env = np.stack([envelope(a[0]) for a in A])
    D = np.stack([sdf(e) for e in env])
    shp = D.shape[1:]
    X = D.reshape(len(D), -1)
    w = np.where(is_if, 0.5 / is_if.sum(), 0.5 / (~is_if).sum())
    mu, V, var = pca(X, w)
    k_keep = int(np.searchsorted(np.cumsum(var) / var.sum(), 0.95) + 1)
    print(f"outline model: {len(D)} shapes ({is_if.sum()} IF, {(~is_if).sum()} P56); "
          f"variance explained by modes 1-5: {np.round(np.cumsum(var)[:5] / var.sum(), 3).tolist()}; "
          f"{k_keep} modes reach 95 %", flush=True)
    mean_env = mu.reshape(shp) < 0
    print("Dice of each shape's envelope with the mean envelope: IF " +
          ", ".join(f"{ids[i][3:].replace('2018_05_22_', '')} {dice(env[i], mean_env):.3f}" for i in np.nonzero(is_if)[0]) +
          f"; P56 median {np.median([dice(env[i], mean_env) for i in np.nonzero(~is_if)[0]]):.3f}", flush=True)
    # leave-one-out on the IF shapes: model from the rest, reconstruct the held-out envelope
    loo = {}
    for i in np.nonzero(is_if)[0]:
        keep = np.arange(len(D)) != i
        wi = np.where(is_if[keep], 0.5 / is_if[keep].sum(), 0.5 / (~is_if[keep]).sum())
        m_i, V_i, _ = pca(X[keep], wi)
        row = {"mean only": dice(env[i], m_i.reshape(shp) < 0)}
        for k in (2, 4, 8):
            c = V_i[:k] @ (X[i] - m_i)
            row[f"{k} modes"] = dice(env[i], (m_i + c @ V_i[:k]).reshape(shp) < 0)
        loo[ids[i]] = {a: round(b, 3) for a, b in row.items()}
        print(f"  leave-one-out {ids[i][3:].replace('2018_05_22_', ''):18s} " +
              "  ".join(f"{a}: {b:.3f}" for a, b in loo[ids[i]].items()), flush=True)
    # does mode 1-3 track medial-lateral position in P56?
    off = np.array([m.get("offset_um", np.nan) for m in meta])
    C = (X - mu) @ V[:3].T
    corr = {f"mode {j + 1}": round(float(np.corrcoef(np.abs(off[~is_if]), C[~is_if, j])[0, 1]), 2) for j in range(3)}
    print("P56: correlation of mode scores with |offset|:", corr, flush=True)
    np.savez_compressed(OUT / "stage1_outline.npz", mean_sdf=mu.reshape(shp), modes=V[:12].reshape((-1,) + shp),
                        variances=var[:12], scores=(X - mu) @ V[:12].T, ids=np.array(ids), envelopes=env)
    (OUT / "stage1_outline.json").write_text(json.dumps(dict(
        grid_um=GRID, envelope_close_um=ENV_CLOSE_UM, sdf_clip_um=SDF_CLIP_UM, weights="IF and P56 groups 50/50",
        variance_explained=np.round(np.cumsum(var)[:12] / var.sum(), 4).tolist(), modes_for_95pct=k_keep,
        leave_one_out_dice=loo, p56_mode_offset_corr=corr), indent=1))


def stage2_outline_rank():
    """Fit the outline model to every IF section (graded ones too, for reference):
    rotate to atlas orientation, align the tissue to the s1_2 frame (similarity),
    then report envelope Dice with the mean and with the best 8-mode fit, and the
    area the section has beyond / is missing against its own 8-mode fit."""
    lib = json.loads((LIB / "library.json").read_text())
    ent = {e["id"]: e for e in lib["entries"]}
    cur = json.loads((LIB / "curation.json").read_text())
    st0 = np.load(OUT / "stage0_shape_set.npz")
    om = np.load(OUT / "stage1_outline.npz")
    mu, V = om["mean_sdf"], om["modes"][:8]
    shp = mu.shape
    ref_t = st0["ref_tissue"]
    mean_env = mu < 0
    rows = []
    for sid in WL.SLIDES_20X + WL.SLIDES_10X:
        d = load_ds(sid)
        tis, _ = base_cut(d["tissue"], d["egl"])
        o = cur.get("orientation_overrides", {}).get(f"if_{sid}") or ent[f"if_{sid}"]["orientation"]
        t = WL.crop(WL.rotate_to_atlas(tis, o["rot_deg"], o["mirror"]))[0]
        p, w = align_tissue(np.stack([t]), ref_t)
        env = envelope(w[0])
        x = sdf(env).ravel() - mu.ravel()
        c = V.reshape(len(V), -1) @ x
        fit = (mu.ravel() + c @ V.reshape(len(V), -1)).reshape(shp) < 0
        r = dict(slide=sid.replace("2018_05_22_", ""),
                 genotype="EZH2 cKO" if sid.startswith("2018_") and "s1_2" not in sid else "control",
                 igl_grade=cur["igl_grade"][f"if_{sid}"], wm_grade=cur["wm_grade"][f"if_{sid}"],
                 orientation_margin=o.get("mirror_margin", -1), scale_to_ref=round(float(np.exp(-p[1])), 2),
                 dice_mean=round(dice(env, mean_env), 3), dice_fit8=round(dice(env, fit), 3),
                 extra_mm2=round(float((env & ~fit).sum() * GRID ** 2 / 1e6), 3),
                 missing_mm2=round(float((fit & ~env).sum() * GRID ** 2 / 1e6), 3),
                 mode_z=np.round(c[:3] / np.sqrt(om["variances"][:3] / om["variances"][:3].sum()
                                                 * 1.0), 2).tolist())
        rows.append(r)
        print(f"{r['slide']:16s} {r['genotype']:9s} IGL {r['igl_grade']:7s} WM {r['wm_grade']:7s} Dice mean {r['dice_mean']:.3f} fit {r['dice_fit8']:.3f} "
              f"extra {r['extra_mm2']:.3f} missing {r['missing_mm2']:.3f} mm2  orient margin {r['orientation_margin']:.3f}",
              flush=True)
    (OUT / "stage2_outline_rank.json").write_text(json.dumps(rows, indent=1))


# --------------------------------------------------------------------------- #
# Stage 1b - IGL and white-matter shape models from the viable (graded) shapes
# --------------------------------------------------------------------------- #
FEATURES = {"igl": dict(ch=2, flag="in_igl", clip=100.0), "wm": dict(ch=1, flag="in_wm", clip=150.0)}


def sdf_clip(mask, clip):
    d = (ndi.distance_transform_edt(~mask) - ndi.distance_transform_edt(mask)) * GRID
    return np.clip(d, -clip, clip).astype(np.float32)


def stage1_features():
    z = np.load(OUT / "stage0_shape_set.npz")
    meta = json.loads((OUT / "stage0.json").read_text())["shapes"]
    A = z["aligned"]
    ids = [m["id"] for m in meta]
    is_if = np.array([m["source"].startswith("IF") for m in meta])
    res, save = {}, {}
    for name, f in FEATURES.items():
        sel = np.array([m[f["flag"]] for m in meta])
        idx = np.nonzero(sel)[0]
        M = A[idx, f["ch"]] & A[idx, 0]
        w = np.where(is_if[idx], 0.5 / is_if[idx].sum(), 0.5 / (~is_if[idx]).sum())
        prob = (w[:, None, None] * M).sum(0)
        prob_if = M[is_if[idx]].mean(0)
        prob_p56 = M[~is_if[idx]].mean(0)
        X = np.stack([sdf_clip(m, f["clip"]) for m in M]).reshape(len(M), -1)
        mu, V, var = pca(X, w.copy())
        shp = M.shape[1:]
        mean_mask = mu.reshape(shp) < 0
        modes = []
        for j in range(3):
            sd = np.sqrt(var[j] / w.sum()) if var[j] > 0 else 0
            sd = np.sqrt(var[j])                      # weighted PCA: components scaled by sqrt(w)
            modes.append([(mu + k * sd * V[j]).reshape(shp) < 0 for k in (-2, 2)])
        # leave-one-out on the IF shapes; reference = Dice with each other IF shape
        loo = {}
        if_i = [i for i in range(len(idx)) if is_if[idx][i]]
        for i in if_i:
            keep = np.arange(len(idx)) != i
            wi = np.where(is_if[idx][keep], 0.5 / is_if[idx][keep].sum(), 0.5 / (~is_if[idx][keep]).sum())
            m_i, V_i, _ = pca(X[keep], wi)
            row = dict(mean=dice(M[i], m_i.reshape(shp) < 0))
            c = V_i[:8] @ (X[i] - m_i)
            row["modes8"] = dice(M[i], (m_i + c @ V_i[:8]).reshape(shp) < 0)
            pr = (wi[:, None, None] * M[keep]).sum(0)
            row["prob50"] = dice(M[i], pr >= 0.5)
            others = [dice(M[i], M[j]) for j in if_i if j != i]
            row["other_median"] = float(np.median(others))
            row["other_best"] = float(np.max(others))
            loo[ids[idx[i]]] = {k: round(v, 3) for k, v in row.items()}
        res[name] = dict(n_if=int(is_if[idx].sum()), n_p56=int((~is_if[idx]).sum()), ids=[ids[i] for i in idx],
                         variance_explained=np.round(np.cumsum(var)[:10] / var.sum(), 3).tolist(),
                         area_mean_mm2=round(float(mean_mask.sum() * GRID ** 2 / 1e6), 3),
                         area_shapes_mm2=[round(float(m.sum() * GRID ** 2 / 1e6), 3) for m in M], loo=loo)
        save.update({f"{name}_prob": prob.astype(np.float32), f"{name}_prob_if": prob_if.astype(np.float32),
                     f"{name}_prob_p56": prob_p56.astype(np.float32), f"{name}_mean": mean_mask,
                     f"{name}_modes": np.array(modes), f"{name}_shapes": M, f"{name}_idx": idx})
        L = res[name]["loo"]
        print(f"{name}: {res[name]['n_if']} IF + {res[name]['n_p56']} P56 viable shapes; modes 1-5 explain "
              f"{res[name]['variance_explained'][:5]}; mean shape {res[name]['area_mean_mm2']} mm2", flush=True)
        for k, v in L.items():
            print(f"   held out {k:12s} Dice: mean shape {v['mean']:.3f}, 8-mode fit {v['modes8']:.3f}, "
                  f"probability>=0.5 {v['prob50']:.3f} | vs another section: median {v['other_median']:.3f}, best {v['other_best']:.3f}",
                  flush=True)
    np.savez_compressed(OUT / "stage1_features.npz", tissue_mean=z["mean_tissue"], ref_tissue=z["ref_tissue"], **save)
    (OUT / "stage1_features.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    {"stage0": stage0, "outline": stage1_outline, "rank": stage2_outline_rank, "features": stage1_features}[sys.argv[1] if len(sys.argv) > 1 else "stage0"]()
