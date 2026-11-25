"""
Library of cerebellar internal granule layer (IGL) shapes from all our data.
Companion to wm_library.py (same grid, crop, thumbnail canvas approach).

    python scripts/igl_library.py

Sources
  user IF, 20x (8 slides) and 10x (11 slides): depth_model anatomy-fixed labels,
      IGL = code 4 only (PCL, code 7, is recorded separately as pcl_mm2, not merged);
      rotated into the atlas orientation (anterior left, dorsal up) with the per-section
      rotation already found for the white-matter library (outputs/wm_library/library.json);
      the orientation search is not rerun.
  user gold s1_2: hand-painted; its 4 and 5 are swapped (4 = ML, 5 = IGL) and are swapped back.
  Allen P4, P7, P14, P56 reference Nissl: outputs/allen_<age>/sections (label 4).
      P14 label 4 is IGL + white matter (not separable in P14 Nissl): every P14 entry is flagged.
      P56 IGL was checked against 13 atlas plates (outputs/allen_p56/p56_validation_vs_atlas.json).
      P4, P7 IGL comes from Nissl grey-level classes and is not checked against any annotation.

Output: outputs/igl_library/
  library.npz   per entry: igl_<id>, tissue_<id> (bool, 4 um/px, cropped to tissue)
  library.json  per entry: source, age, offset/orientation, areas, pieces, thickness, skeleton
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
from wm_library import (FIXED, GRID, PX_10X, SEG, SLIDES_10X, SLIDES_20X, THUMB,  # noqa: E402
                        crop, resample, rotate_to_atlas)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/igl_library"
WM_LIB = ROOT / "outputs/wm_library/library.json"
P56_CHECK = ROOT / "outputs/allen_p56/p56_validation_vs_atlas.json"
GOLD = SEG / "labelled" / "2018_05_22_s1_2_p27-0021_labelled.tiff"
PX_20X = 0.5119
ALLEN = {"P4": 3.96, "P7": 3.96, "P14": 4.188, "P56": 4.188}
ALLEN_NOTES = {"P4": ["IGL not validated"], "P7": ["IGL not validated"],
               "P14": ["IGL + white matter merged"], "P56": []}
IGL, PCL, ANOM = 4, 7, 9
UNCERTAIN = 0.02        # mirror_margin below this -> "orientation uncertain"
SPECK_UM2 = 500.0       # specks removed at native resolution (as wm_library)
PIECE_UM2 = 1000.0      # connected components counted as pieces (8-connectivity)
SPUR_UM = 20.0          # skeleton end pruning (as wm_library.branch_stats)
HOLE_UM2 = 500.0        # holes in the IGL below this are filled for thickness_filled_um only
IGL_RGBA = (226, 180, 30, 255)
S1_2 = "if_2018_05_22_s1_2_p27-0021"


def orientations():
    lib = json.loads(WM_LIB.read_text())
    return {e["id"]: {k: e["orientation"][k] for k in ("rot_deg", "mirror", "mirror_margin")}
            for e in lib["entries"] if e.get("orientation")}


def user_entries(ori):
    for mag, slides, px in (("20x", SLIDES_20X, PX_20X), ("10x", SLIDES_10X, PX_10X)):
        for s in slides:
            p = FIXED / f"{s}_fixed_labels.tif"
            if not p.exists():
                print("missing", p, flush=True)
                continue
            yield dict(id=f"if_{s}", source=f"IF {mag}", age="P7", name=s.replace("2018_05_22_", ""),
                       px=px, lab=tifffile.imread(p), orientation=ori[f"if_{s}"], notes=[])
    g = tifffile.imread(GOLD).astype(np.int16)
    gl = g.copy()
    gl[g == 4], gl[g == 5] = 5, 4
    yield dict(id="gold_s1_2", source="IF gold", age="P7", name="s1_2 gold standard", px=PX_20X,
               lab=gl.astype(np.uint8), orientation=ori["gold_s1_2"], notes=["hand-painted"])


def allen_entries():
    for age, px in ALLEN.items():
        d = ROOT / f"outputs/allen_{age.lower()}"
        st = json.loads((d / f"{age.lower()}_nissl_stats.json").read_text())
        for r in sorted(st, key=lambda r: r["offset_um"]):
            yield dict(id=f"allen_{age}_sec{r['sec']:03d}", source=f"Allen {age}", age=age,
                       name=f"sec {r['sec']}", offset_um=float(r["offset_um"]), px=px,
                       lab=tifffile.imread(d / f"sections/sec{r['sec']:03d}_labels.tif"),
                       notes=list(ALLEN_NOTES[age]))


def skeleton(mask):
    s = morphology.skeletonize(mask)
    raw = s.copy()
    k = np.ones((3, 3), np.uint8)
    for _ in range(int(round(SPUR_UM / GRID))):
        e = s & (ndi.convolve(s.astype(np.uint8), k, mode="constant") == 2)
        if not e.any():
            break
        s &= ~e
    return s if s.any() else raw


def metrics(igl, tissue):
    a2 = GRID ** 2
    igl_mm2, tissue_mm2 = igl.sum() * a2 / 1e6, tissue.sum() * a2 / 1e6
    m = dict(igl_mm2=round(float(igl_mm2), 3), tissue_mm2=round(float(tissue_mm2), 3),
             igl_pct=round(float(100 * igl_mm2 / tissue_mm2), 1) if tissue_mm2 else 0.0)
    if not igl.any():
        m.update(pieces=0, largest_piece_pct=0.0, thickness_um=0.0, skeleton_mm=0.0, holes=0, thickness_filled_um=0.0)
        return m
    cc, _ = ndi.label(igl, structure=np.ones((3, 3), bool))
    a = np.bincount(cc.ravel())[1:] * a2
    sk = skeleton(igl)
    dt = ndi.distance_transform_edt(igl)
    # holes enclosed by the IGL (4-connected background); small ones filled for a hole-insensitive thickness
    hc, nh = ndi.label(ndi.binary_fill_holes(igl) & ~igl)
    small = np.bincount(hc.ravel()) * a2 < HOLE_UM2
    small[0] = False
    filled = igl | small[hc]
    m.update(pieces=int((a >= PIECE_UM2).sum()), largest_piece_pct=round(float(100 * a.max() / a.sum()), 1),
             thickness_um=round(float(2 * np.median(dt[sk]) * GRID), 1),
             skeleton_mm=round(float(sk.sum() * GRID / 1000), 2), holes=int(nh),
             thickness_filled_um=round(float(2 * np.median(ndi.distance_transform_edt(filled)[skeleton(filled)]) * GRID), 1))
    return m


def iou(a, b):
    return float((a & b).sum() / max((a | b).sum(), 1))


def thumb(tissue, igl, px, canvas):
    t = resample(tissue, px, THUMB)
    g = resample(igl, px, THUMB) & t
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
        return int((~k).sum())

    edge = t & ~ndi.binary_erosion(t)
    lost = put(t, (128, 132, 140, 60))
    put(edge, (128, 132, 140, 200))
    put(g, IGL_RGBA)
    return Image.fromarray(img, "RGBA"), int(g.sum()), lost


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "thumbs").mkdir(exist_ok=True)
    ori = orientations()
    meta, arrays, s12 = [], {}, {}
    for e in list(user_entries(ori)) + list(allen_entries()):
        lab = e.pop("lab")
        tissue = (lab > 0) & (lab != ANOM)
        raw = (lab == IGL) & tissue
        igl = morphology.remove_small_objects(raw, max_size=int(SPECK_UM2 / e["px"] ** 2))
        if "orientation" in e:
            e["pcl_mm2"] = round(float(((lab == PCL) & tissue).sum() * e["px"] ** 2 / 1e6), 3)
        tissue_g = resample(tissue, e["px"], GRID)
        igl_g = resample(igl, e["px"], GRID) & tissue_g
        if e["id"] in (S1_2, "gold_s1_2"):
            s12[e["id"]] = (raw, igl_g, lab != 1)       # 1 = tissue left unassigned
        if "orientation" in e:
            o = e["orientation"]
            tissue_g = rotate_to_atlas(tissue_g, o["rot_deg"], o["mirror"])
            igl_g = rotate_to_atlas(igl_g, o["rot_deg"], o["mirror"]) & tissue_g
            if o["mirror_margin"] < UNCERTAIN:
                e["notes"].append("orientation uncertain")
        tissue_g, igl_g = crop(tissue_g, igl_g)
        e.update(metrics(igl_g, tissue_g), extent_um=[int(tissue_g.shape[0] * GRID), int(tissue_g.shape[1] * GRID)])
        arrays[f"igl_{e['id']}"] = igl_g
        arrays[f"tissue_{e['id']}"] = tissue_g
        meta.append(e)
        print(f"{e['id']:40s} IGL {e['igl_mm2']:.3f} mm2 ({e['igl_pct']:4.1f}% of {e['tissue_mm2']:.2f}), "
              f"pieces {e['pieces']:3d}, largest {e['largest_piece_pct']:5.1f}%, thick {e['thickness_um']:5.1f} um, "
              f"skel {e['skeleton_mm']:5.2f} mm {e['notes']}", flush=True)
        del lab, tissue, raw, igl

    # s1_2: automatic IGL against the gold standard
    (ra, ga, _), (rb, gb, assigned) = s12[S1_2], s12["gold_s1_2"]
    checks = dict(s1_2_igl_iou_native=round(iou(ra, rb), 3),
                  s1_2_igl_iou_native_gold_assigned=round(iou(ra & assigned, rb & assigned), 3),
                  s1_2_igl_iou_grid=round(iou(ga, gb), 3))
    v = json.loads(P56_CHECK.read_text())
    d = [r["dice_igl"] for r in v]
    checks.update(p56_igl_dice_min=round(min(d), 2), p56_igl_dice_max=round(max(d), 2), p56_plates=len(v),
                  p56_plate_offset_um=[min(r["offset_um"] for r in v), max(r["offset_um"] for r in v)])
    print("checks", checks, flush=True)

    np.savez_compressed(OUT / "library.npz", **arrays)
    ext = np.array([m["extent_um"] for m in meta])
    canvas = tuple(int(np.ceil(x / THUMB)) + 6 for x in ext.max(0))
    ratios = []
    for m in meta:
        im, n_igl, lost = thumb(arrays[f"tissue_{m['id']}"], arrays[f"igl_{m['id']}"], GRID, canvas)
        im.save(OUT / "thumbs" / f"{m['id']}.png")
        if m["igl_mm2"] > 0:
            ratios.append(n_igl * THUMB ** 2 / 1e6 / m["igl_mm2"])
        if lost:
            print("thumbnail clipped", m["id"], lost, "px", flush=True)
    r = np.array(ratios)
    print(f"thumb IGL area / library IGL area: median {np.median(r):.3f}, range {r.min():.3f}-{r.max():.3f}")
    (OUT / "library.json").write_text(json.dumps(dict(
        grid_um=GRID, thumb_um=THUMB, canvas_px=canvas,
        definitions=dict(
            igl="IF: label 4 only (PCL, label 7, kept apart as pcl_mm2); gold: label 5 (swapped back to 4); "
                "Allen: label 4 (P14: IGL + white matter)",
            speck_um2=SPECK_UM2, piece_min_um2=PIECE_UM2, uncertain_mirror_margin=UNCERTAIN, spur_um=SPUR_UM,
            pieces="8-connected IGL components >= piece_min_um2",
            largest_piece_pct="share of IGL area in the largest component",
            thickness_um="2 x median Euclidean distance transform of the IGL sampled on its skeleton "
                         "(skeleton ends pruned by spur_um), on the grid",
            skeleton_mm="pruned skeleton pixel count x grid_um",
            hole_fill_um2=HOLE_UM2,
            holes="background components (4-connected) enclosed by the IGL",
            thickness_filled_um="thickness_um after filling IGL holes < hole_fill_um2 (masks are not filled); "
                                "thickness_um drops wherever the IGL is perforated",
            s1_2_igl_iou="automatic vs gold IGL on s1_2 before resampling: _native over all pixels, "
                         "_native_gold_assigned excluding tissue the gold standard leaves unassigned (label 1), "
                         "_grid on the grid before rotation"),
        checks=checks, entries=meta), indent=1))
    print(f"{len(meta)} entries; thumbnail canvas {canvas} px at {THUMB} um/px")


if __name__ == "__main__":
    main()
