"""
Repair layer labels inside user-flagged folia, from the section's own layer thicknesses.

    python scripts/folium_repair.py

Folium numbers and territories come from igl_check_page.py (display frame: anterior
left, dorsal up; territory = tissue whose nearest surface point lies on that folium's
stretch of surface). Thicknesses are measured in the section's other folia.

  "smear" (IGL called over the ML and EGL): IGL / PCL / unassigned / anomaly pixels
      shallower than the local EGL depth -> EGL (oEGL above the local oEGL depth, iEGL
      below); shallower than the local EGL + ML depth -> ML; deeper anomaly pixels
      within the local IGL thickness of white matter -> IGL.
  "bleach" (IGL breached): within 60 um of the bleached (anomaly) patch, anomaly / ML /
      unassigned / PCL pixels within the local IGL thickness of white matter and deeper
      than the local EGL + ML depth -> IGL.
  ML thickness = median EGL-to-IGL gap in the neighbouring folia (fissure-safe).

Changes are decided at 3x downsampling in the display frame, mapped back to the
section's own frame, and stamped onto native pixels carrying the source label.
Output: outputs/folium_repair/<slide>_repaired_labels.tif, <slide>_repair.json,
        <slide>_repair_display.npz (before / after display-frame labels for the page)
"""
import json
import sys
from pathlib import Path

import numpy as np
import tifffile
from scipy import ndimage as ndi
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).parent))
import igl_check_page as P  # noqa: E402
import wm_library as WL  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/folium_repair"
REQUESTS = {"s4_A_p27": [(6, "smear"), (8, "bleach")]}       # from the user, 2026-11-18
PX = 0.65 * P.DS
BLEACH_REACH_UM = 60.0


def display_labels(sid, tr):
    lab = tifffile.imread(WL.FIXED / f"{sid}_fixed_labels.tif")
    small = lab[::P.DS, ::P.DS]
    A = np.array(tr["A"])
    labr, _ = P.rotate(small, A, tr["canvas"], 0)
    y0, x0 = tr["crop"]
    h, w = tr["shape"]
    return lab, small, labr[y0:y0 + h, x0:x0 + w].copy()


def territories(labr, sid):
    segs = np.load(P.OUT / f"{sid}_folium_segments.npz")
    pts, ids = [], []
    for k in segs.files:
        s = segs[k]
        pts.append(s)
        ids += [int(k[1:])] * len(s)
    tree = cKDTree(np.vstack(pts))
    tis = ndi.binary_fill_holes(labr > 0)
    yy, xx = np.nonzero(tis)
    _, nn = tree.query(np.stack([yy, xx], 1))
    terr = np.zeros(labr.shape, np.int16)
    terr[yy, xx] = np.array(ids)[nn]
    return terr


def thickness(labr, terr, use, dO, dW):
    """EGL depth below the surface, oEGL depth, ML thickness (EGL-to-IGL gap, so ML in
    closed fissures does not inflate it) and IGL thickness from white matter, in the
    neighbouring folia."""
    m = np.isin(terr, use)
    egl_m = np.isin(labr, (2, 3))
    dE = ndi.distance_transform_edt(~egl_m) * PX
    dI = ndi.distance_transform_edt(labr != 4) * PX
    egl = m & egl_m & (dO <= 100)
    oegl = m & (labr == 3) & (dO <= 100)
    ml = m & (labr == 5) & (dE + dI <= 150)
    igl = m & (labr == 4) & (dW <= 200)
    # oEGL / iEGL boundary: the depth where oEGL stops being the majority of the EGL
    bins = np.arange(0, 101, 2.0)
    d = dO[egl]
    is_o = labr[egl] == 3
    t_o = 0.0
    for b0 in bins[:-1]:
        sel = (d >= b0) & (d < b0 + 2)
        if sel.sum() >= 20 and is_o[sel].mean() < 0.5:
            t_o = float(b0)
            break
    t = dict(egl=float(np.percentile(dO[egl], 90)), oegl=t_o,
             ml=float(np.median(dE[ml] + dI[ml])), igl=float(np.percentile(dW[igl], 90)))
    t["outer"] = t["egl"] + t["ml"]
    return t


def repair(labr, terr, requests):
    tis = labr > 0
    filled = ndi.binary_fill_holes(tis)
    dO = ndi.distance_transform_edt(filled) * PX
    dW = ndi.distance_transform_edt(~np.isin(labr, (6, 8))) * PX
    out = labr.copy()
    report = []
    flagged = [f for f, _ in requests]
    for fol, kind in requests:
        nbr = [k for k in (fol - 2, fol - 1, fol + 1, fol + 2) if k > 0 and k not in flagged]
        t = thickness(labr, terr, nbr, dO, dW)
        R = (terr == fol) & tis
        before = np.bincount(labr[R], minlength=10)
        if kind == "smear":
            src = R & np.isin(labr, (1, 4, 7, 9))
            e = src & (dO <= t["egl"])
            out[e] = np.where(dO[e] <= t["oegl"], 3, 2)
            m_ = src & ~e & (dO <= t["outer"])
            out[m_] = 5
            deep_anom = R & (labr == 9) & (dO > t["outer"]) & (dW <= t["igl"])
            out[deep_anom] = 4
        else:
            near_bleach = ndi.distance_transform_edt(~(R & (labr == 9))) * PX <= BLEACH_REACH_UM
            src = R & near_bleach & np.isin(labr, (1, 5, 7, 9)) & (dW <= t["igl"]) & (dO > t["outer"])
            out[src] = 4
        after = np.bincount(out[R], minlength=10)
        surf = R & (dO <= 60)
        report.append(dict(folium=fol, kind=kind, neighbours=nbr, thickness_um={k: round(v, 1) for k, v in t.items()},
                           territory_mm2=round(float(R.sum() * PX * PX / 1e6), 4),
                           changed_mm2=round(float(((out != labr) & R).sum() * PX * PX / 1e6), 4),
                           before={P_NAMES[k]: round(float(before[k] / R.sum()), 3) for k in range(10) if before[k]},
                           after={P_NAMES[k]: round(float(after[k] / R.sum()), 3) for k in range(10) if after[k]},
                           igl_in_outer_60um=[round(float((labr[surf] == 4).mean()), 3), round(float((out[surf] == 4).mean()), 3)]))
    return out, report


P_NAMES = {0: "bg", 1: "unassigned", 2: "iEGL", 3: "oEGL", 4: "IGL", 5: "ML", 6: "WM", 7: "PCL", 8: "DCN", 9: "anomaly"}


def back_to_section(new_disp, old_disp, small, tr):
    """map display-frame changes onto the section's own 3x grid (nearest pixel)."""
    A = np.array(tr["A"])
    off = np.array(tr["offset"])
    y0, x0 = tr["crop"]
    H, W = small.shape
    yy, xx = np.indices((H, W))
    q = np.stack([yy.ravel(), xx.ravel()], 0).astype(float)
    disp = np.linalg.solve(A, q - off[:, None]) - np.array([[y0], [x0]])
    r = np.round(disp[0]).astype(int)
    c = np.round(disp[1]).astype(int)
    ok = (r >= 0) & (r < new_disp.shape[0]) & (c >= 0) & (c < new_disp.shape[1])
    out = small.copy().ravel()
    idx = np.nonzero(ok)[0]
    nv, ov = new_disp[r[idx], c[idx]], old_disp[r[idx], c[idx]]
    ch = (nv != ov) & (small.ravel()[idx] == ov)
    out[idx[ch]] = nv[ch]
    return out.reshape(H, W)


def apply_native(native, ds, decided):
    H, W = native.shape
    ys = np.minimum(np.arange(H) // ds, decided.shape[0] - 1)
    xs = np.minimum(np.arange(W) // ds, decided.shape[1] - 1)
    orig_up = native[::ds, ::ds][np.ix_(ys, xs)]
    new_up = decided[np.ix_(ys, xs)]
    change = (new_up != orig_up) & (native == orig_up)
    out = native.copy()
    out[change] = new_up[change]
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    trs = {t["sid"]: t for t in json.loads((P.OUT / "transforms.json").read_text())}
    for sid, req in REQUESTS.items():
        tr = trs[sid]
        native, small, labr = display_labels(sid, tr)
        terr = territories(labr, sid)
        new, rep = repair(labr, terr, req)
        small_new = back_to_section(new, labr, small, tr)
        nat_new = apply_native(native, P.DS, small_new)
        tifffile.imwrite(OUT / f"{sid}_repaired_labels.tif", nat_new, compression="zlib")
        np.savez_compressed(OUT / f"{sid}_repair_display.npz", before=labr, after=new, territory=terr)
        (OUT / f"{sid}_repair.json").write_text(json.dumps(dict(slide=sid, requests=req, folia=rep), indent=1))
        for r in rep:
            print(json.dumps(r))
        px_n = WL.PX_10X
        ch = (nat_new != native).sum() * px_n * px_n / 1e6
        print(f"{sid}: native pixels changed {ch:.4f} mm2 -> {OUT / (sid + '_repaired_labels.tif')}")


# --------------------------------------------------------------------------- #
# v2: user-confirmed bleach boxes, two candidate repairs (A image, B shape)
# --------------------------------------------------------------------------- #
CONFIRMED_BLEACH = {"s4_A_p27": ["B5"], "s4_C_p27": ["B2"]}      # user, igl_check page, 2026-11-18


def box_mask_display(box, tr, shape):
    from skimage.draw import polygon
    A = np.array(tr["A"])
    off = np.array(tr["offset"])
    y0, x0 = tr["crop"]
    corners = [(box["row0"], box["col0"]), (box["row0"], box["col0"] + box["cols"]),
               (box["row0"] + box["rows"], box["col0"] + box["cols"]), (box["row0"] + box["rows"], box["col0"])]
    pts = np.array([np.linalg.solve(A, np.array([r / P.DS, c / P.DS]) - off) - [y0, x0] for r, c in corners])
    rr, cc = polygon(pts[:, 0], pts[:, 1], shape)
    m = np.zeros(shape, bool)
    m[rr, cc] = True
    return m


def display_image(sid, tr, ch):
    raw = tifffile.imread(P.SB10 / f"{sid}_native.tif")[ch, ::P.DS, ::P.DS].astype(np.float32)
    A = np.array(tr["A"])
    img, _ = P.rotate(raw, A, tr["canvas"], 1)
    y0, x0 = tr["crop"]
    h, w = tr["shape"]
    return img[y0:y0 + h, x0:x0 + w]


def bleach_fix(labr, R, method, img3, t):
    """img3 = display-frame (p27, DAPI, NeuN). Returns a new label map (only R changes)."""
    out = labr.copy()
    tis = labr > 0
    filled = ndi.binary_fill_holes(tis)
    dO = ndi.distance_transform_edt(filled) * PX
    dW = ndi.distance_transform_edt(~np.isin(labr, (6, 8))) * PX
    ring = (ndi.distance_transform_edt(~R) * PX <= 100) & ~R & tis & (labr != 9)
    inside = R & tis
    if method == "A":
        g = 1.0 / np.median([np.median(img3[c][inside]) / np.median(img3[c][ring]) for c in range(3)])
        n = ndi.gaussian_filter(img3[2], 3 / PX)
        out_ok = tis & ~R & (labr != 9)
        t_n = np.sqrt(np.median(n[out_ok & (labr == 4)]) * np.median(n[out_ok & (labr == 5)]))
        bright = inside & (n * g > t_n)
        bright = ndi.binary_opening(bright, structure=np.ones((3, 3), bool))
        zone = inside & (dO > t["egl"]) & ~np.isin(labr, (2, 3))
        out[bright & zone] = 4
        out[inside & (labr == 4) & ~bright] = 5
        return out, dict(gain=round(float(g), 2), neun_threshold=float(t_n))
    zone = inside & (dW <= t["igl"]) & (dO > t["outer"]) & np.isin(labr, (1, 5, 7, 9))
    out[zone] = 4
    eg = inside & np.isin(labr, (1, 4, 7, 9)) & (dO <= t["egl"])
    out[eg] = np.where(dO[eg] <= t["oegl"], 3, 2)
    ml = inside & np.isin(labr, (1, 4, 7, 9)) & (dO > t["egl"]) & (dO <= t["outer"])
    out[ml] = 5
    return out, {}


def main_v2():
    OUT.mkdir(parents=True, exist_ok=True)
    trs = {t["sid"]: t for t in json.loads((P.OUT / "transforms.json").read_text())}
    summary = {}
    for sid in CONFIRMED_BLEACH:
        tr = trs[sid]
        native, small, labr = display_labels(sid, tr)
        terr = territories(labr, sid)
        base, rep_smear = repair(labr, terr, [r for r in REQUESTS.get(sid, []) if r[1] == "smear"])
        img3 = np.stack([display_image(sid, tr, c) for c in range(3)])
        filled = ndi.binary_fill_holes(labr > 0)
        dO = ndi.distance_transform_edt(filled) * PX
        dW = ndi.distance_transform_edt(~np.isin(labr, (6, 8))) * PX
        t = thickness(labr, terr, [k for k in range(1, 40) if k not in [f for f, _ in REQUESTS.get(sid, [])]], dO, dW)
        regions, fixes = [], {}
        for meth in ("A", "B"):
            cur = base.copy()
            for bid in CONFIRMED_BLEACH[sid]:
                box = [b for b in tr["bleach_candidates"] if b["id"] == bid][0]
                R = box_mask_display(box, tr, labr.shape)
                new, info = bleach_fix(cur, R, meth, img3, t)
                chg = (new != cur) & R
                regions.append(dict(box=bid, method=meth, changed_mm2=round(float(chg.sum() * PX * PX / 1e6), 4),
                                    before={P_NAMES[k]: round(float(v), 3) for k, v in enumerate(np.bincount(cur[R], minlength=10) / max(R.sum(), 1)) if v >= 0.01},
                                    after={P_NAMES[k]: round(float(v), 3) for k, v in enumerate(np.bincount(new[R], minlength=10) / max(R.sum(), 1)) if v >= 0.01},
                                    **info))
                cur = new
            fixes[meth] = cur
        np.savez_compressed(OUT / f"{sid}_repair_display.npz", before=labr, fix_A=fixes["A"], fix_B=fixes["B"], territory=terr)
        for meth in ("A", "B"):
            nat = apply_native(native, P.DS, back_to_section(fixes[meth], labr, small, tr))
            tifffile.imwrite(OUT / f"{sid}_repaired_{meth}.tif", nat, compression="zlib")
        summary[sid] = dict(thickness_um={k: round(v, 1) for k, v in t.items()}, smear=rep_smear, bleach=regions)
        for r in rep_smear + regions:
            print(sid, json.dumps(r))
    (OUT / "repair_v2.json").write_text(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main_v2() if (len(sys.argv) > 1 and sys.argv[1] == "v2") else main()
