"""
IGL problem check page: IF image, IGL / white-matter outlines and labels of s4_A and s4_C
(plus s5_G as the good-IGL reference), with click-to-mark problem spots saved to the
page's db (collection "igl_check", one document per section).

    python scripts/igl_check_page.py  -> outputs/igl_check/igl_check.html
                                         outputs/igl_check/transforms.json (display -> native)

Images: 10x native pages p27, DAPI, NeuN; drawn at 3x downsampling (1.95 um/px) in atlas
orientation (anterior left, dorsal up) using outputs/wm_library/library.json rotations.
Marker coordinates are stored as fractions of the displayed image; transforms.json maps
them back to native pixels.
"""
import base64
import io
import json
import sys
from pathlib import Path

import numpy as np
import tifffile
from PIL import Image
from scipy import ndimage as ndi

sys.path.insert(0, str(Path(__file__).parent))
import wm_library as WL  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/igl_check"
SB10 = WL.SEG / "matlab_anomaly_10x/sandbox"
DS = 3
PALETTE = {0: (0, 0, 0, 0), 1: (150, 152, 160, 255), 2: (110, 175, 235, 255), 3: (0, 196, 214, 255),
           4: (232, 186, 28, 255), 5: (146, 212, 140, 255), 6: (226, 122, 50, 255), 7: (214, 64, 64, 255),
           8: (196, 84, 196, 255), 9: (80, 80, 86, 255)}
SECTIONS = [("s4_A_p27", True), ("s4_C_p27", True), ("s5_G_p27", False)]


def rot_matrix(rot_deg, mirror):
    th = np.deg2rad(rot_deg)
    c, s = np.cos(th), np.sin(th)
    return (np.diag([1.0, -1.0]) if mirror else np.eye(2)) @ np.array([[c, -s], [s, c]]).T


def rotate(arr, A, n, order):
    H, W = arr.shape
    o = np.array([H, W]) / 2.0 - A @ (np.array([n, n]) / 2.0)
    return ndi.affine_transform(arr, A, o, output_shape=(n, n), order=order), o


def uri(img, fmt, **kw):
    buf = io.BytesIO()
    img.save(buf, fmt, **kw)
    return f"data:image/{fmt.lower()};base64," + base64.b64encode(buf.getvalue()).decode()


def build(sid, lib):
    o = lib[f"if_{sid}"]["orientation"]
    lab = tifffile.imread(WL.FIXED / f"{sid}_fixed_labels.tif")[::DS, ::DS]
    raw = tifffile.imread(SB10 / f"{sid}_native.tif")[:, ::DS, ::DS].astype(np.float32)
    H, W = lab.shape
    A = rot_matrix(o["rot_deg"], o["mirror"])
    n = int(np.ceil(np.hypot(H, W))) + 4
    labr, off = rotate(lab, A, n, 0)
    tis = labr > 0
    ys, xs = np.nonzero(tis)
    m = 12
    y0, y1, x0, x1 = max(ys.min() - m, 0), min(ys.max() + m + 1, n), max(xs.min() - m, 0), min(xs.max() + m + 1, n)
    labr = labr[y0:y1, x0:x1]
    tis_raw = lab > 0
    chans = []
    for ci in (0, 2, 1):                               # R = p27, G = NeuN, B = DAPI
        ch = raw[ci]
        lo, hi = np.percentile(ch[tis_raw], [1, 99.7])
        r, _ = rotate(np.clip((ch - lo) / max(hi - lo, 1e-6), 0, 1), A, n, 1)
        chans.append(r[y0:y1, x0:x1])
    rgb = (np.stack(chans, -1) * 255).astype(np.uint8)
    rgb[labr == 0] = (rgb[labr == 0] * 0.35).astype(np.uint8)
    # outlines
    igl, wm = labr == 4, np.isin(labr, (6, 8))
    ov = np.zeros(labr.shape + (4,), np.uint8)
    ov[wm & ~ndi.binary_erosion(wm)] = (255, 140, 40, 255)
    ov[igl & ~ndi.binary_erosion(igl)] = (255, 236, 60, 255)
    lr = np.zeros(labr.shape + (4,), np.uint8)
    for k, col in PALETTE.items():
        lr[labr == k] = col
    fol, base = folia(labr, 0.65 * DS)
    segs = {f["n"]: f.pop("seg") for f in fol}
    np.savez_compressed(OUT / f"{sid}_folium_segments.npz", **{f"f{k}": np.array(v) for k, v in segs.items()})
    tr = dict(sid=sid, ds=DS, rot_deg=o["rot_deg"], mirror=o["mirror"], A=A.tolist(), offset=off.tolist(),
              folia=fol, base=base,
              crop=[int(y0), int(x0)], shape=list(labr.shape), canvas=n,
              note="display (row, col) = crop + (r, c); native_ds = A @ (row, col) + offset; native = native_ds * ds")
    return dict(id=sid, w=int(labr.shape[1]), h=int(labr.shape[0]),
                folia=[dict(n=f["n"], x=f["col"] / labr.shape[1], y=f["row"] / labr.shape[0]) for f in fol],
                cell_px=250.0 / (0.65 * DS),
                img=uri(Image.fromarray(rgb, "RGB"), "JPEG", quality=86),
                outline=uri(Image.fromarray(ov, "RGBA"), "PNG", optimize=True),
                labels=uri(Image.fromarray(lr, "RGBA"), "PNG", optimize=True)), tr, (A, off, y0, x0, labr.shape)


def _endpoints(sk):
    nb = ndi.convolve(sk.astype(np.uint8), np.ones((3, 3), np.uint8), mode="constant")
    return sk & (nb == 2)


def _prune(sk, min_len_px, rounds=4):
    """remove end branches shorter than min_len_px (walk from each end to the first junction)."""
    sk = sk.copy()
    H, W = sk.shape

    def nbrs(p):
        return [(p[0] + dy, p[1] + dx) for dy in (-1, 0, 1) for dx in (-1, 0, 1)
                if (dy or dx) and 0 <= p[0] + dy < H and 0 <= p[1] + dx < W and sk[p[0] + dy, p[1] + dx]]
    for _ in range(rounds):
        changed = False
        for e in map(tuple, np.argwhere(_endpoints(sk))):
            path, prev, cur = [e], None, e
            while True:
                nb = [q for q in nbrs(cur) if q != prev and q not in path]
                if len(nb) != 1:
                    break
                prev, cur = cur, nb[0]
                if len(nbrs(cur)) > 2:
                    if len(path) < min_len_px:
                        for q in path:
                            sk[q] = False
                        changed = True
                    break
                path.append(cur)
        if not changed:
            break
    return sk


def folia(labr, px, prm=None):
    """Number the folia of a display-frame label map (anterior left, dorsal up).
    Fissure openings along the envelope (tissue closed over 60 um and filled):
      closed fissures - EGL skeleton deeper than 40 um below the tissue surface, pieces
                        >= 100 um long; opening = its point nearest the surface;
      open fissures   - stretches of the tissue outline more than 60 um inside the envelope.
    Base = the longest stretch of envelope with white matter within 25 um (peduncle cut).
    Folia = envelope stretches between consecutive openings on the dorsal way from the
    anterior to the posterior end of the base; label placed 40 um inside the midpoint."""
    from skimage import measure, morphology as mo
    q = dict(env_um=60.0, deep_um=40.0, fis_min_um=100.0, notch_um=60.0, wm_um=25.0, merge_um=100.0, inset_um=40.0,
             min_folium_um=80.0, **(prm or {}))
    tissue = labr > 0
    filled = ndi.binary_fill_holes(tissue)
    r = int(round(q["env_um"] / px))
    env = ndi.binary_fill_holes(ndi.binary_closing(np.pad(filled, r + 2), structure=mo.disk(r))[r + 2:-r - 2, r + 2:-r - 2])
    dT = ndi.distance_transform_edt(filled) * px           # depth below the tissue surface
    dE = ndi.distance_transform_edt(env) * px              # depth below the envelope
    cont = max(measure.find_contours(env.astype(float), 0.5), key=len)
    ci = np.clip(np.round(cont).astype(int), 0, np.array(env.shape) - 1)
    mouths = []
    egl = np.isin(labr, (2, 3))
    sk = mo.skeletonize(egl) & (dT > q["deep_um"])
    cc, n = ndi.label(sk, structure=np.ones((3, 3), bool))
    for i in range(1, n + 1):
        m = np.argwhere(cc == i)
        if len(m) * px < q["fis_min_um"]:
            continue
        top = m[np.argmin(dT[tuple(m.T)])]
        mouths.append(int(np.argmin(np.hypot(*(cont - top).T))))
    # open notches: tissue outline points deep inside the envelope
    tc = max(measure.find_contours(filled.astype(float), 0.5), key=len)
    tci = np.clip(np.round(tc).astype(int), 0, np.array(env.shape) - 1)
    deep = dE[tci[:, 0], tci[:, 1]] > q["notch_um"]
    lab_runs, nr = ndi.label(deep)
    for i in range(1, nr + 1):
        pts = tc[lab_runs == i]
        bottom = pts[np.argmax(dE[tci[lab_runs == i, 0], tci[lab_runs == i, 1]])]
        mouths.append(int(np.argmin(np.hypot(*(cont - bottom).T))))
    # fissures predicted from the white matter (midline between neighbouring branches)
    import wm_scaffold as WS
    sc = WS.fissure_lines(labr, px)
    lc, nl = ndi.label(sc["lines"], structure=np.ones((3, 3), bool))
    for i in range(1, nl + 1):
        m = np.argwhere(lc == i)
        if len(m) * px < q["fis_min_um"]:
            continue
        top = m[np.argmin(dT[tuple(m.T)])]
        if dT[tuple(top)] <= 150:
            mouths.append(int(np.argmin(np.hypot(*(cont - top).T))))
    # base: the longest envelope stretch with no EGL nearby that has white matter close to it
    wm = np.isin(labr, (6, 8))
    dW = ndi.distance_transform_edt(~wm) * px
    dEg = ndi.distance_transform_edt(~egl) * px
    noegl = dEg[ci[:, 0], ci[:, 1]] > q["wm_um"]
    wmnear = dW[ci[:, 0], ci[:, 1]] <= 80.0
    N = len(cont)
    runs, k = [], 0
    lab1, nrun = ndi.label(np.concatenate([noegl, noegl]))
    seen = set()
    for i in range(1, nrun + 1):
        idx = np.nonzero(lab1 == i)[0]
        if len(idx) > N or idx[0] >= N:
            continue
        key = tuple(sorted(set((idx % N).tolist())))[:3]
        if key in seen:
            continue
        seen.add(key)
        runs.append((len(idx), float(wmnear[idx % N].mean()), idx[0] % N, idx[-1] % N))
    # orientation is verified (dorsal up): the base is the most ventral EGL-free stretch of
    # at least 150 um (the peduncle cut); fall back to the longest stretch
    def mean_row(r_):
        idx = (np.arange(r_[2], r_[2] + r_[0])) % N
        return float(cont[idx, 0].mean())
    cand = [r_ for r_ in runs if r_[0] * px >= 150.0] or runs
    if not cand:
        return [], dict(error="no base found")
    best_run = max(cand, key=mean_row)
    best = best_run[0]
    b0, b1 = best_run[2], best_run[3]
    ant, post = (b0, b1) if cont[b0, 1] < cont[b1, 1] else (b1, b0)
    fwd = np.arange(ant, ant + ((post - ant) % N) + 1) % N
    bwd = np.arange(ant, ant - ((ant - post) % N) - 1, -1) % N
    path = fwd if cont[fwd, 0].mean() < cont[bwd, 0].mean() else bwd
    where = {int(v): k for k, v in enumerate(path)}
    pos = sorted({where[m] for m in mouths if m in where})
    step = np.hypot(*np.diff(cont[path], axis=0).T)
    arc = np.concatenate([[0.0], np.cumsum(step)]) * px
    merged = []
    for pp in pos:
        if not merged or arc[pp] - arc[merged[-1]] > q["merge_um"]:
            merged.append(pp)
    bounds = [0] + [pp for pp in merged if q["min_folium_um"] < arc[pp] < arc[-1] - q["min_folium_um"]] + [len(path) - 1]
    out = []
    inner = ndi.distance_transform_edt(env) * px
    for k in range(len(bounds) - 1):
        a_, b_ = bounds[k], bounds[k + 1]
        if arc[b_] - arc[a_] < q["min_folium_um"]:
            continue
        mid = path[int(np.argmin(np.abs(arc[a_:b_ + 1] - (arc[a_] + arc[b_]) / 2))) + a_]
        p0 = cont[mid]
        # move inward: step down the envelope depth gradient
        gy, gx = np.gradient(inner)
        pnt = p0.copy()
        for _ in range(int(q["inset_um"] / px)):
            iy, ix = np.clip(np.round(pnt).astype(int), 0, np.array(env.shape) - 1)
            g = np.array([gy[iy, ix], gx[iy, ix]])
            if np.hypot(*g) < 1e-6:
                break
            pnt = pnt + g / np.hypot(*g)
        out.append(dict(n=len(out) + 1, row=float(pnt[0]), col=float(pnt[1]),
                        length_um=round(float(arc[b_] - arc[a_])),
                        seg=cont[path[a_:b_ + 1]].tolist()))
    return out, dict(base_ant=cont[ant].tolist(), base_post=cont[post].tolist(),
                     base_um=round(float(best * px)), openings=len(merged))


def s4c_gap(lib, geom):
    """the ML gap between IGL pieces found in s4_C, in display fractions."""
    lab = tifffile.imread(WL.FIXED / "s4_C_p27_fixed_labels.tif")[::DS, ::DS]
    ST = np.ones((3, 3), bool)
    cc, _ = ndi.label(lab == 4, structure=ST)
    a = np.bincount(cc.ravel())
    a[0] = 0
    i, j = np.argsort(-a)[:2]
    da, (iy, ix) = ndi.distance_transform_edt(cc != i, return_indices=True)
    k = np.argmin(np.where(cc == j, da, np.inf))
    pb = np.array(np.unravel_index(k, cc.shape))
    pa = np.array([iy[tuple(pb)], ix[tuple(pb)]])
    mid = (pa + pb) / 2.0
    A, off, y0, x0, shp = geom
    out = np.linalg.solve(A, mid - off) - [y0, x0]
    return dict(x=float(out[1] / shp[1]), y=float(out[0] / shp[0]), gap_um=float(da[tuple(pb)] * 0.65 * DS))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    lib = {e["id"]: e for e in json.loads((WL.OUT / "library.json").read_text())["entries"]}
    items, trs, geoms = [], [], {}
    for sid, ask in SECTIONS:
        it, tr, geom = build(sid, lib)
        it["ask"] = ask
        items.append(it)
        trs.append(tr)
        geoms[sid] = geom
    import bleach_boxes as BB
    for it in items:
        sid = it["id"]
        lab = tifffile.imread(WL.FIXED / f"{sid}_fixed_labels.tif")
        raw = tifffile.imread(SB10 / f"{sid}_native.tif")
        bx = BB.find_boxes(raw, (lab > 0) & (lab != 9), WL.PX_10X)[:5]
        A, off, y0, x0, shp = geoms[sid]
        polys = []
        for k, b in enumerate(bx):
            corners = [(b["row0"], b["col0"]), (b["row0"], b["col0"] + b["cols"]),
                       (b["row0"] + b["rows"], b["col0"] + b["cols"]), (b["row0"] + b["rows"], b["col0"])]
            pts = []
            for (r_, c_) in corners:
                d = np.linalg.solve(A, np.array([r_ / DS, c_ / DS]) - off) - [y0, x0]
                pts.append([float(d[1] / shp[1]), float(d[0] / shp[0])])
            polys.append(dict(id=f"B{k + 1}", pts=pts, score=round(b["score"], 2), native=b))
        it["boxes"] = [dict(id=q["id"], pts=q["pts"], score=q["score"]) for q in polys]
        trs[[t["sid"] for t in trs].index(sid)]["bleach_candidates"] = [dict(id=q["id"], **q["native"]) for q in polys]
    # proposed repairs (folium_repair.py v2: smear depth fix + shape-based bleach fix B)
    FIXES = {"s4_A_p27": [("smear6", "Folium 6: smear fix (layers rebuilt from neighbouring folia)"),
                          ("bleachB5", "Box B5: bleach fix (IGL shell around white matter)")],
             "s4_C_p27": [("bleachB2", "Box B2: bleach fix (IGL shell around white matter)")]}
    for it in items:
        it["fixes"] = []
        fz = ROOT / f"outputs/folium_repair/{it['id']}_repair_display.npz"
        if it["id"] in FIXES and fz.exists():
            z = np.load(fz)
            new = z["fix_B"]
            old_l = z["before"]
            if new.shape != (it["h"], it["w"]):
                continue
            igl, wm = new == 4, np.isin(new, (6, 8))
            ov = np.zeros(new.shape + (4,), np.uint8)
            ov[wm & ~ndi.binary_erosion(wm)] = (255, 140, 40, 255)
            ov[igl & ~ndi.binary_erosion(igl)] = (255, 236, 60, 255)
            ch = ndi.binary_dilation(new != old_l, iterations=2)
            ov[ch & ~ndi.binary_erosion(ch, iterations=2)] = (120, 220, 255, 255)
            lr = np.zeros(new.shape + (4,), np.uint8)
            for k, col in PALETTE.items():
                lr[new == k] = col
            lr[ch & ~ndi.binary_erosion(ch, iterations=2)] = (120, 220, 255, 255)
            it["fix_outline"] = uri(Image.fromarray(ov, "RGBA"), "PNG", optimize=True)
            it["fix_labels"] = uri(Image.fromarray(lr, "RGBA"), "PNG", optimize=True)
            it["fixes"] = [dict(id=k, text=t) for k, t in FIXES[it["id"]]]
    # ribbon repair (igl_repair.py): bridged IGL gaps, as numbered repairs R1, R2, ...
    rj = ROOT / "outputs/igl_ribbon/ribbon_repair_items.json"
    ritems = json.loads(rj.read_text()) if rj.exists() else {}
    for it in items:
        it["ribbon"] = []
        fz = ROOT / f"outputs/igl_ribbon/{it['id']}_ribbon_repair_display.npz"
        if it["id"] in ritems and fz.exists():
            z = np.load(fz)
            new, old_l = z["after"], z["before"]
            if new.shape != (it["h"], it["w"]):
                continue
            igl, wm = new == 4, np.isin(new, (6, 8))
            ch = ndi.binary_dilation(new != old_l, iterations=2)
            ov = np.zeros(new.shape + (4,), np.uint8)
            ov[wm & ~ndi.binary_erosion(wm)] = (255, 140, 40, 255)
            ov[igl & ~ndi.binary_erosion(igl)] = (255, 236, 60, 255)
            ov[ch & ~ndi.binary_erosion(ch, iterations=2)] = (90, 255, 150, 255)
            lr = np.zeros(new.shape + (4,), np.uint8)
            for k, col in PALETTE.items():
                lr[new == k] = col
            lr[ch & ~ndi.binary_erosion(ch, iterations=2)] = (90, 255, 150, 255)
            it["ribbon_outline"] = uri(Image.fromarray(ov, "RGBA"), "PNG", optimize=True)
            it["ribbon_labels"] = uri(Image.fromarray(lr, "RGBA"), "PNG", optimize=True)
            it["ribbon"] = ritems[it["id"]]["items"]
    # final repaired labels (igl_finalize.py)
    for it in items:
        fz = ROOT / f"outputs/igl_fixed/{it['id']}_final_display.npz"
        if not fz.exists():
            continue
        z = np.load(fz)
        new, old_l = z["after"], z["before"]
        if new.shape != (it["h"], it["w"]):
            continue
        igl, wm = new == 4, np.isin(new, (6, 8))
        ch = ndi.binary_dilation(new != old_l, iterations=2)
        ov = np.zeros(new.shape + (4,), np.uint8)
        ov[wm & ~ndi.binary_erosion(wm)] = (255, 140, 40, 255)
        ov[igl & ~ndi.binary_erosion(igl)] = (255, 236, 60, 255)
        ov[ch & ~ndi.binary_erosion(ch, iterations=2)] = (255, 255, 255, 255)
        lr = np.zeros(new.shape + (4,), np.uint8)
        for k, col in PALETTE.items():
            lr[new == k] = col
        lr[ch & ~ndi.binary_erosion(ch, iterations=2)] = (255, 255, 255, 255)
        it["final_outline"] = uri(Image.fromarray(ov, "RGBA"), "PNG", optimize=True)
        it["final_labels"] = uri(Image.fromarray(lr, "RGBA"), "PNG", optimize=True)
    # IGL ribbon model: candidate completions at the user's markers (igl_model_apply.py)
    import igl_model_page as MP
    trmap = {t["sid"]: t for t in trs}
    for it in items:
        it["model"] = MP.build(it["id"], trmap[it["id"]]) if it["ask"] else []
    gap = s4c_gap(lib, geoms["s4_C_p27"])
    for it in items:
        it["hint"] = ([dict(x=gap["x"], y=gap["y"], text=f"my guess: {gap['gap_um']:.0f} µm ML gap between IGL pieces")]
                      if it["id"] == "s4_C_p27" else [])
    (OUT / "transforms.json").write_text(json.dumps(trs, indent=1))
    page = (Path(__file__).parent / "igl_check_template.html").read_text()
    page = page.replace("__DATA__", json.dumps(items))
    p = OUT / "igl_check.html"
    p.write_text(page)
    print(p, f"{p.stat().st_size / 1e6:.2f} MB", [(i["id"], i["w"], i["h"], len(i["folia"])) for i in items], gap)


if __name__ == "__main__":
    main()
