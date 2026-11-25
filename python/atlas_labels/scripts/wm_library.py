"""
Library of cerebellar white-matter shapes from all our data.

    python scripts/wm_library.py

Sources
  user IF, 20x (8 slides) and 10x (11 slides): depth_model anatomy-fixed labels,
      white matter = DWL (6) + DCN (8); rotated into the atlas orientation
      (anterior left, dorsal up) by the orientation search of fit_p7_template.py
  user gold s1_2: hand-painted DWL + DCN
  Allen P4, P7, P56 reference Nissl: outputs/allen_<age>/sections (label 6).
      P14 has no separable white matter (merged into IGL) and is not included.
      P56 white matter is reliable only within +-0.6 mm of the midline.

Output: outputs/wm_library/
  library.npz   per entry: wm_<id>, tissue_<id> (bool, 4 um/px, cropped)
  library.json  per entry: source, age, offset, areas, branch tips/junctions, orientation
  thumbs/<id>.png  RGBA thumbnails at 15 um/px, common canvas (same scale for all)
"""
import json
import sys
from pathlib import Path

import numpy as np
import tifffile
from PIL import Image
from scipy import ndimage as ndi
from skimage import morphology

sys.path.insert(0, str(Path(__file__).parent))
import fit_p7_template as F  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
SEG = Path("/Users/jpurzner/Dropbox/imaging_analysis/cerebellar_segmentation/python")
FIXED = SEG / "depth_model/out/anatomy_fixes"
OUT = ROOT / "outputs/wm_library"
SLIDES_20X = ["2018_05_22_s1_2_p27-0021", "2018_05_22_s1_4_p27-0020", "2018_05_22_s1_5_p27-0018",
              "2018_05_22_s2_2_p27", "2018_05_22_s2_4_p27-0024", "2018_05_22_s2_5_p27-0026",
              "2018_05_22_s3_2_p27-0001", "2018_05_22_s3_5_p27-0006"]
SLIDES_10X = ["s4_A_p27", "s4_C_p27", "s4_F_p27", "s5_A_p27", "s5_B_p27", "s5_C_p27",
              "s5_E_p27", "s5_F_p27", "s5_G_p27", "s6_C_p27", "s6_G_p27"]
PX_10X = 0.65       # um/px of the 10x panoramas
GRID = 4.0          # library grid, um/px
THUMB = 15.0        # thumbnail scale, um/px
ALLEN = {"P4": (3.96, None), "P7": (3.96, None), "P56": (4.188, 600.0)}
GALLERY_STEP = {"P4": 200, "P7": 200, "P56": 150}
GALLERY_MAX = {"P4": 2000, "P7": 2000, "P56": 1200}


def resample(mask, px_from, px_to):
    z = px_from / px_to
    return ndi.zoom(mask.astype(np.float32), z, order=1) >= 0.5


def crop(*ms, pad=4):
    t = ms[0]
    ys, xs = np.nonzero(t)
    sl = (slice(max(ys.min() - pad, 0), ys.max() + pad + 1), slice(max(xs.min() - pad, 0), xs.max() + pad + 1))
    return [m[sl] for m in ms]


def branch_stats(wm, px):
    sk = morphology.skeletonize(wm)
    nb = ndi.convolve(sk.astype(np.uint8), np.ones((3, 3), np.uint8), mode="constant") - 1
    s = sk.copy()
    for _ in range(int(round(20 / px))):        # prune spurs < 20 um
        e = s & (ndi.convolve(s.astype(np.uint8), np.ones((3, 3), np.uint8), mode="constant") == 2)
        if not e.any():
            break
        s &= ~e
    nb = ndi.convolve(s.astype(np.uint8), np.ones((3, 3), np.uint8), mode="constant") - 1
    tips = int((s & (nb == 1)).sum())
    junc = int(ndi.label(s & (nb >= 3), structure=np.ones((3, 3), bool))[1])
    return dict(tips=tips, junctions=junc, skeleton_mm=round(float(s.sum() * px / 1000), 2))


def orient_user(lab, block):
    tgt = F.fractions(F.remap(lab, F.AUTO_LUT), block)
    v = 1.0 - tgt[F.CH.index(F.ANOM)]
    sw = F.orientation_sweep(tgt, v)
    best = sw[0]
    other = max(r[0] for r in sw if r[2] != best[2])
    return dict(rot_deg=int(best[1]), mirror=bool(best[2]), score=round(float(best[0]), 3),
                mirror_margin=round(float(best[0] - other), 3))


def rotate_to_atlas(mask, rot_deg, mirror):
    th = np.deg2rad(rot_deg)
    c, s = np.cos(th), np.sin(th)
    R = np.array([[c, -s], [s, c]])
    M = np.diag([1.0, -1.0]) if mirror else np.eye(2)
    A = M @ R.T
    H, W = mask.shape
    n = int(np.ceil(np.hypot(H, W))) + 4
    cin = np.array([H, W]) / 2.0
    cout = np.array([n, n]) / 2.0
    o = cin - A @ cout
    return ndi.affine_transform(mask.astype(np.float32), A, o, output_shape=(n, n), order=1) >= 0.5


def user_entries():
    out = []
    # 10x panoramas are 0.65 um/px (camera 6.5 um / 10x; confirmed by the user), not the
    # 1.0239 um/px the MATLAB 10x port assumes. Block 6 -> 3.9 um grid for the orientation search.
    for mag, slides, px, block in (("20x", SLIDES_20X, 0.5119, 8), ("10x", SLIDES_10X, PX_10X, 6)):
        for s in slides:
            p = FIXED / f"{s}_fixed_labels.tif"
            if not p.exists():
                continue
            lab = tifffile.imread(p)
            ori = orient_user(lab, block)
            yield dict(id=f"if_{s}", source=f"IF {mag}", age="P7", name=s.replace("2018_05_22_", ""),
                       px=px, lab=lab, orientation=ori, notes=[])
    g = tifffile.imread(SEG / "labelled" / "2018_05_22_s1_2_p27-0021_labelled.tiff").astype(np.int16)
    gl = g.copy()
    gl[g == 4], gl[g == 5] = 5, 4
    lab20 = tifffile.imread(FIXED / "2018_05_22_s1_2_p27-0021_fixed_labels.tif")
    ori = orient_user(lab20, 8)
    yield dict(id="gold_s1_2", source="IF gold", age="P7", name="s1_2 gold standard", px=0.5119,
               lab=gl.astype(np.uint8), orientation=ori, notes=["hand-painted"])


def allen_entries():
    for age, (px, reliable) in ALLEN.items():
        d = ROOT / f"outputs/allen_{age.lower()}"
        st = json.loads((d / f"{age.lower()}_nissl_stats.json").read_text())
        for r in st:
            lab = tifffile.imread(d / f"sections/sec{r['sec']:03d}_labels.tif")
            notes = []
            if reliable and abs(r["offset_um"]) > reliable:
                notes.append("WM unreliable here")
            yield dict(id=f"allen_{age}_sec{r['sec']:03d}", source=f"Allen {age}", age=age,
                       name=f"sec {r['sec']}", offset_um=float(r["offset_um"]), px=px, lab=lab,
                       orientation=None, notes=notes)


def thumb(tissue, wm, px, canvas):
    t = resample(tissue, px, THUMB)
    w = resample(wm, px, THUMB) & t
    sk = morphology.skeletonize(w)
    H, W = canvas
    img = np.zeros((H, W, 4), np.uint8)
    ys, xs = np.nonzero(t)
    cy, cx = (ys.min() + ys.max()) // 2, (xs.min() + xs.max()) // 2
    oy, ox = H // 2 - cy, W // 2 - cx

    def put(mask, rgba):
        yy, xx = np.nonzero(mask)
        yy, xx = yy + oy, xx + ox
        k = (yy >= 0) & (yy < H) & (xx >= 0) & (xx < W)
        img[yy[k], xx[k]] = rgba

    edge = t & ~ndi.binary_erosion(t)
    put(t, (128, 132, 140, 60))
    put(edge, (128, 132, 140, 200))
    put(w, (226, 122, 50, 255))
    put(sk, (90, 38, 6, 255))
    return Image.fromarray(img, "RGBA")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "thumbs").mkdir(exist_ok=True)
    meta, arrays, keep = [], {}, []
    for e in list(user_entries()) + list(allen_entries()):
        lab = e.pop("lab")
        tissue = (lab > 0) & (lab != 9)
        wm = np.isin(lab, (6, 8)) & tissue
        wm = morphology.remove_small_objects(wm, max_size=int(500 / e["px"] ** 2))
        tissue_g = resample(tissue, e["px"], GRID)
        wm_g = resample(wm, e["px"], GRID) & tissue_g
        if e["orientation"]:
            o = e["orientation"]
            tissue_g, wm_g = rotate_to_atlas(tissue_g, o["rot_deg"], o["mirror"]), \
                rotate_to_atlas(wm_g, o["rot_deg"], o["mirror"])
            if o["mirror_margin"] < 0.02:
                e["notes"].append("orientation uncertain")
        tissue_g, wm_g = crop(tissue_g, wm_g)
        e.update(wm_mm2=round(float(wm_g.sum() * GRID ** 2 / 1e6), 3),
                 tissue_mm2=round(float(tissue_g.sum() * GRID ** 2 / 1e6), 3),
                 extent_um=[int(tissue_g.shape[0] * GRID), int(tissue_g.shape[1] * GRID)],
                 **(branch_stats(wm_g, GRID) if wm_g.any() else dict(tips=0, junctions=0, skeleton_mm=0)))
        arrays[f"wm_{e['id']}"] = wm_g
        arrays[f"tissue_{e['id']}"] = tissue_g
        meta.append(e)
        print(f"{e['id']:42s} WM {e['wm_mm2']:.2f} mm2, tissue {e['tissue_mm2']:.2f}, tips {e['tips']:3d}, "
              f"junctions {e['junctions']:3d} {e.get('orientation') or ''} {e['notes']}", flush=True)
    np.savez_compressed(OUT / "library.npz", **arrays)
    # common thumbnail canvas
    ext = np.array([m["extent_um"] for m in meta])
    canvas = tuple(int(np.ceil(v / THUMB)) + 6 for v in ext.max(0))
    for m in meta:
        thumb(arrays[f"tissue_{m['id']}"], arrays[f"wm_{m['id']}"], GRID, canvas).save(OUT / "thumbs" / f"{m['id']}.png")
    (OUT / "library.json").write_text(json.dumps(dict(grid_um=GRID, thumb_um=THUMB, canvas_px=canvas,
                                                      entries=meta), indent=1))
    print(f"{len(meta)} entries; thumbnail canvas {canvas} px at {THUMB} um/px")


if __name__ == "__main__":
    main()
