"""
Run the pipeline gap filler (igl_gapfill.fill_gaps) on every real section.

    python scripts/igl_gapfill_apply.py [section ...]      # default: all 19 anatomy-fixed sections

Input: depth_model/out/anatomy_fixes/<sid>_fixed_labels.tif (the same input as the walk in
igl_finalize.py), NeuN and the MATLAB tear mask from the anomaly sandboxes.
Pixel size: 20x (2018_05_22_*) 0.5119 um/px, 10x 0.65 um/px.
Controls (all 10x + s1_2): accepted bridges are stamped onto the native labels (grid
decisions written only onto pixels labelled unassigned / ML / PCL / anomaly, and onto white
matter where the bridge's outline test found the white-matter labels wrong).
EZH2 cKO (20x except s1_2): bridges are computed and shown for review but not stamped; a
break in a cKO IGL may be phenotype.
No user-confirmed bleach regions are used here; igl_finalize.py adds them for s4_A / s4_C.
-> outputs/igl_gapfill/<sid>_gapfilled_labels.tif   (controls)
   outputs/igl_gapfill/<sid>_gapfill.npz            (grid: labels, NeuN, bridges' fills and curves)
   outputs/igl_gapfill/apply.json
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import tifffile
from scipy import ndimage as ndi

sys.path.insert(0, str(Path(__file__).parent))
import igl_gapfill as GF  # noqa: E402
import igl_repair as IR  # noqa: E402
import igl_repair_eval as EV  # noqa: E402
import igl_ribbon as RB  # noqa: E402
import igl_ribbon_eval as RE  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/igl_gapfill"
G = RB.GRID


def section_list():
    return sorted(f.name.replace("_fixed_labels.tif", "") for f in RE.FX.glob("*_fixed_labels.tif"))


def short(sid):
    return sid.replace("2018_05_22_", "").split("_p27")[0]


def info(sid):
    big = sid.startswith("2018_05_22")
    px = 0.5119 if big else 0.65
    group = "control" if (not big or "_s1_2_" in sid) else "cKO"
    return px, group, big


def images(sid, shape):
    """NeuN (3 um blur) and tear mask on the 4 um grid."""
    px, _, big = info(sid)
    sb = EV.SB20 if big else EV.SB10
    raw, tear = sb / f"{sid}_native.tif", sb / f"{sid}_native_anomaly_tear.tif"
    z = px / G
    n = tifffile.imread(raw)[1 if big else 2].astype(np.float32)
    n = ndi.zoom(ndi.gaussian_filter(n, 3.0 / px), z, order=1)
    t = ndi.zoom((tifffile.imread(tear) > 0).astype(np.float32), z, order=1) >= 0.25 if tear.exists() else np.zeros(shape, bool)

    def fit(a):
        out = np.zeros(shape, a.dtype)
        h, w = min(a.shape[0], shape[0]), min(a.shape[1], shape[1])
        out[:h, :w] = a[:h, :w]
        return out
    return fit(n), fit(t)


def stamp(native, px, grid_paint, grid_allow_wm):
    f = G / px
    H, W = native.shape
    ys = np.minimum((np.arange(H) / f).astype(int), grid_paint.shape[0] - 1)
    xs = np.minimum((np.arange(W) / f).astype(int), grid_paint.shape[1] - 1)
    p_up = grid_paint[np.ix_(ys, xs)]
    w_up = grid_allow_wm[np.ix_(ys, xs)]
    allowed = np.isin(native, (1, 5, 7, 9)) | (w_up & np.isin(native, (6, 8)))
    out = native.copy()
    out[p_up & allowed] = 4
    return out


def grid_code(g):
    """one code per grid pixel for display: 0 bg, 1 other tissue, 2 EGL, 3 ML/PCL, 4 IGL, 5 WM, 6 anomaly."""
    c = np.zeros(g["igl"].shape, np.uint8)
    c[g["tissue"]] = 1
    c[g["ml"]] = 3
    c[g["wm"]] = 5
    c[g["anom"]] = 6
    c[g["egl"]] = 2
    c[g["igl"]] = 4
    return c


def run(sid):
    px, group, _ = info(sid)
    native = tifffile.imread(RE.FX / f"{sid}_fixed_labels.tif")
    g = RB.grid_labels(native, px)
    neun, tear = images(sid, g["igl"].shape)
    t0 = time.time()
    new, br, rib = GF.fill_gaps(g, neun=neun, tear=tear, verbose=True)
    t_fill = time.time() - t0
    # the walk on the same input, for comparison
    new_w, br_w, _ = IR.repair_gaps(g, None, None, neun, tear)
    joined = [b for b in br if b["status"] == "joined"]
    paint = new & ~g["igl"]
    allow_wm = np.zeros_like(paint)
    arrays = dict(code=grid_code(g), neun=neun.astype(np.float16), tear=tear, igl_clean=rib["igl"], walk_paint=new_w & ~g["igl"])
    for k, b in enumerate(br):
        b["id"] = k
        if "_mask" in b:
            arrays[f"fill_{k}"] = b["_mask"]
            arrays[f"box_{k}"] = np.array(b["_box"])
        arrays[f"curve_{k}"] = np.asarray(b["_curve"], np.float32)
        if b["status"] == "joined" and b["method"].startswith("ribbon model, white"):
            y0, y1, x0, x1 = b["_box"]
            allow_wm[y0:y1, x0:x1] |= b["_mask"]
    np.savez_compressed(OUT / f"{sid}_gapfill.npz", **arrays)
    entry = dict(sid=sid, short=short(sid), group=group, px_um=px, grid_shape=list(g["igl"].shape),
                 fits=rib["fill_setup"], fill_seconds=round(t_fill, 1),
                 candidates=[{k: v for k, v in b.items() if not k.startswith("_")} for b in br],
                 joined=len(joined), painted_mm2=round(float(paint.sum() * G * G / 1e6), 4),
                 walk=dict(joined=sum(b["status"] == "joined" for b in br_w), painted_mm2=round(float((new_w & ~g["igl"]).sum() * G * G / 1e6), 4)))
    if group == "control":
        out = stamp(native, px, paint, allow_wm)
        tifffile.imwrite(OUT / f"{sid}_gapfilled_labels.tif", out, compression="zlib")
        entry.update(stamped=True, igl_added_mm2=round(float(((out == 4) & (native != 4)).sum() * px * px / 1e6), 4),
                     from_labels={str(c): round(float(((native == c) & (out == 4)).sum() * px * px / 1e6), 4) for c in (1, 5, 6, 7, 8, 9)})
    else:
        entry["stamped"] = False
    print(f"{short(sid):6s} {group:7s} ends {rib['fill_setup']['n_ends']:3d}  pairs {len(br):3d}  joined {len(joined)} "
          f"({entry['painted_mm2']} mm2; walk {entry['walk']['joined']}, {entry['walk']['painted_mm2']} mm2)  "
          + "; ".join(f"{b['chord_um']:.0f}um {b['method'].split(',')[0]} wm_sh {b['wm_shallow']}" for b in joined)
          + f"  [{t_fill:.0f}s]", flush=True)
    return entry


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    want = sys.argv[1:]
    secs = [s for s in section_list() if not want or short(s) in want or s in want]
    res = {}
    p = OUT / "apply.json"
    if want and p.exists():
        res = json.loads(p.read_text())
    for sid in secs:
        res[sid] = run(sid)
    p.write_text(json.dumps(res, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else float(o)))


if __name__ == "__main__":
    main()
