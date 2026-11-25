"""
Cerebellum outline + internal layers from the Allen Developing Mouse Brain
reference sagittal Nissl series, one age at a time:

    python scripts/allen_nissl_template.py P7 [max_offset_um]

  P7  SectionDataSet 100073790, 20 um sections, allen_reference/P7_nissl
  P4  SectionDataSet 100033295, 20 um sections, allen_reference/P4_nissl
  (d2 images = 4x downsampled = 3.96 um/px)
Seeds: P7 from the hand-picked cerebellum crop; ages with a reference atlas
(P4, P14, P56) from the atlas plate's r1R (cerebellum) polygon, which is drawn
on the same Nissl image.

Per section:
  1. track the cerebellum outward from the midline seed (template match of the
     previous section, keep the dense-cell component overlapping the prior)
  2. envelope = dense component, closed + filled, restricted to tissue
  3. 4-class multi-Otsu on envelope pixels (adaptive per section):
        darkest -> EGL | next -> IGL | next -> WM | palest -> ML / gaps
  4. fissures = medial line of EGL that lies deeper than FISSURE_MIN_DEPTH_UM
     below the envelope (fused EGL pairs); cut out of the envelope -> pia
  5. layers re-assigned against depth from the carved pia

Outputs (outputs/allen_<age>/):
  sections/secNNN_crop.tif        RGB crop of the Nissl image
  sections/secNNN_labels.tif      uint8 labels 0 BG, 3 EGL, 5 ML, 4 IGL, 6 WM
  sections/secNNN_boundaries.tif  uint8 1-px  1 pia (incl. fissures),
                                  3 EGL inner edge, 4 IGL outer edge, 6 WM edge
  <age>_nissl_stats.json, <age>_nissl_contours.npz, <age>_nissl_overview.png,
  <age>_nissl_midline.png
Label codes follow ribbon_unfold_laplace.LBL (EGL not split into o/i: code 3).
Orientation as in the Allen image: dorsal up, anterior LEFT.
"""
import json
import os
import re
import sys
import warnings
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage as ndi
from skimage import filters, measure, morphology
from skimage.feature import match_template
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

ROOT = Path("/Users/jpurzner/Dropbox/Q_research/py_projects/cerebellar_atlas_labels")
# ADMBA prosomeric annotation: cerebellum = isthmic roof (isR, vermis) + r1 roof (r1R)
CB_SEED_IDS = {"16811", "16916"}
CCF_CB_ROOT = 512                 # P56 plates use the adult CCF ontology (CB subtree)

AGES = {
    "P7": dict(src="allen_reference/P7_nissl", px_um=3.96, sec_um=20.0, mid_sec=218,
               seed=dict(kind="crop", sec=217, img_id=101229768, offset=(98, 1993))),
    "P4": dict(src="allen_reference/P4_nissl", px_um=3.96, sec_um=20.0, mid_sec=None,
               seed=dict(kind="atlas", sec=185, atlas_image_id=100711175)),
    "P14": dict(src="allen_reference/P14_nissl_series", px_um=4.188, sec_um=25.0, mid_sec=None,
                seed=dict(kind="atlas", sec=184, atlas_image_id=100519041),
                mode="adult", egl=True, ml_max_um=100.0, ml_typ_um=70.0, wm="none"),
    "P56": dict(src="allen_reference/P56_nissl_series", px_um=4.188, sec_um=25.0, mid_sec=None,
                seed=dict(kind="atlas", sec=201, atlas_image_id=100883921, ontology="ccf"),
                mode="adult", ml_max_um=220.0, ml_typ_um=150.0, wm="pale"),
}
MODE = "egl"          # "egl" (EGL present: P4-P14) or "adult" (no EGL: P56)
ML_MAX_UM = 250.0     # adult: ML taken as tissue within this distance of the IGL/WM core
ML_TYP_UM = 120.0     # adult: typical ML thickness (base-cut hull expansion)
EGL_FLAG = False      # P14: thin EGL on the outer ML surface is labelled too
WM_KIND = "pale"      # adult: WM pale in Nissl; P14 "none": WM as dark as IGL, not separable
OPEN_PX = 3           # opening radius that strips thin dark lines (EGL) from the IGL core
FISS_MARGIN_UM = 10.0 # adult/P14 fissure axis: envelope edge must be this much farther than core
WM_CLOSE_UM = 34.0    # adult: IGL closing radius before filling the WM inside it
EGL_MAX_UM = 40.0     # P14: EGL = dark pixels within this distance of the pia

# set by configure(age)
AGE = None
SRC = OUT = None
PX_UM = SEC_UM = None
MID_SEC = None
DS = 4
WIN = 300                         # px around prior when searching
MARGIN = 40                       # px around cerebellum in saved crops
FISSURE_MIN_DEPTH_UM = 50.0       # a fissure must reach at least this deep
FISSURE_EXT_DEPTH_UM = 26.0       # skeleton deeper than a single EGL band
EGL_ADJ_UM = 60.0                 # EGL+ML zone below pia (partial-volume fix)
# brainstem cut = convex hull of the main EGL network (lobule I -> X line)
HULL_MIN_COMP_FRAC = 0.10         # EGL pieces >= 10% of the largest define the hull
BASE_MAX_REMOVE_FRAC = 0.4        # never remove more than this fraction of the envelope
LBL = dict(BG=0, EGL=3, IGL=4, ML=5, WM=6)


def configure(age):
    global AGE, SRC, OUT, PX_UM, SEC_UM, MID_SEC, MODE, EGL_FLAG, ML_MAX_UM, ML_TYP_UM, WM_KIND
    cfg = AGES[age]
    AGE, SRC, OUT = age, ROOT / cfg["src"], ROOT / f"outputs/allen_{age.lower()}"
    PX_UM, SEC_UM = cfg["px_um"], cfg["sec_um"]
    MODE = cfg.get("mode", "egl")
    EGL_FLAG = cfg.get("egl", False)
    ML_MAX_UM = cfg.get("ml_max_um", ML_MAX_UM)
    ML_TYP_UM = cfg.get("ml_typ_um", ML_TYP_UM)
    WM_KIND = cfg.get("wm", "pale")
    (OUT / "sections").mkdir(parents=True, exist_ok=True)
    if cfg["mid_sec"] is not None:
        MID_SEC = cfg["mid_sec"]
    elif (SRC / "midline.json").exists():      # written by download_allen_nissl.py
        MID_SEC = int(round(json.loads((SRC / "midline.json").read_text())["mid_section"]))
        print(f"{AGE}: midline section {MID_SEC} (from full-series thumbnails)", flush=True)
    else:
        MID_SEC = find_midline(section_files())
    return cfg


def find_midline(files):
    """Symmetry centre of the per-section tissue-area profile."""
    secs, area = [], []
    for sec in sorted(files):
        im = Image.open(files[sec])
        W, H = im.size
        im.draft("L", (W // 8, H // 8))
        a = np.asarray(im.convert("L")).astype(float)
        secs.append(sec)
        area.append((a < 200).sum() * (W / a.shape[1]) ** 2)
    secs, area = np.array(secs), np.array(area)
    best = (-1, None)
    lo, hi = np.percentile(secs, 25), np.percentile(secs, 75)
    for c in np.arange(lo, hi, 0.5):
        m = np.array([(a1, np.interp(2 * c - s0, secs, area)) for s0, a1 in zip(secs, area)
                      if secs.min() <= 2 * c - s0 <= secs.max()])
        if len(m) > 0.6 * len(secs):
            r = np.corrcoef(m[:, 0], m[:, 1])[0, 1]
            if r > best[0]:
                best = (r, c)
    print(f"{AGE}: midline section {best[1]} (area-symmetry r={best[0]:.3f})", flush=True)
    return int(round(best[1]))


def ccf_cb_ids():
    g = json.loads((ROOT / "allen_reference/CCF_structure_graph_1.json").read_text())
    ids = set()

    def walk(n, inside):
        inside = inside or n["id"] == CCF_CB_ROOT
        if inside:
            ids.add(str(n["id"]))
        for c in n.get("children", []) or []:
            walk(c, inside)
    for t in g["msg"]:
        walk(t, False)
    return ids


def atlas_cb_mask(atlas_image_id, sec, shape, ontology="admba"):
    """Rasterise the atlas plate's cerebellum polygons (isR + r1R) in d2 pixels."""
    import urllib.request
    import xml.etree.ElementTree as ET
    from svgpathtools import parse_path
    from PIL import ImageDraw
    cache = ROOT / "allen_reference/atlas_svgs_d2" / f"{AGE}_sec{sec}_{atlas_image_id}.svg"
    cache.parent.mkdir(exist_ok=True)
    if not cache.exists():
        url = f"http://api.brain-map.org/api/v2/svg_download/{atlas_image_id}?downsample=2"
        cache.write_bytes(urllib.request.urlopen(url, timeout=60).read())
    root = ET.parse(cache).getroot()
    seed_ids = ccf_cb_ids() if ontology == "ccf" else CB_SEED_IDS
    img = Image.new("L", (shape[1], shape[0]), 0)
    draw = ImageDraw.Draw(img)

    def walk(el, scale):
        tr = el.attrib.get("transform", "")
        for m in re.finditer(r"scale\(([-\d.e]+)(?:[ ,]+([-\d.e]+))?\)", tr):
            scale = scale * float(m.group(1))
        if el.attrib.get("structure_id") in seed_ids and "d" in el.attrib:
            path = parse_path(el.attrib["d"])
            pts = [path.point(t) for t in np.linspace(0, 1, 2000)]
            draw.polygon([(p.real * scale, p.imag * scale) for p in pts], fill=1)
        for c in el:
            walk(c, scale)

    walk(root, 1.0)
    return np.asarray(img).astype(bool)


def section_files():
    out = {}
    for f in os.listdir(SRC):
        m = re.match(r"^sec(\d{4})_id(\d+)_d2\.jpg$", f)
        if m:
            out[int(m.group(1))] = SRC / f
    return out


def bbox_of(mask, margin, shape):
    ys, xs = np.nonzero(mask)
    return (max(0, ys.min() - margin), min(shape[0], ys.max() + margin + 1),
            max(0, xs.min() - margin), min(shape[1], xs.max() + margin + 1))


def sl(b):
    return (slice(b[0], b[1]), slice(b[2], b[3]))


# ----------------------------------------------------------------------------
def envelope(gs, prior=None):
    """gs: smoothed gray window. Returns envelope mask, thresholds."""
    t_dense, t_mid, t_bg = filters.threshold_multiotsu(gs, 4)
    tissue = morphology.remove_small_objects(gs < t_bg, max_size=2000)
    dense = morphology.remove_small_objects(gs < t_dense, max_size=200)
    lab = measure.label(dense)
    if lab.max() == 0:
        return None
    sizes = np.bincount(lab.ravel())[1:]
    if prior is None:
        keep = [int(np.argmax(sizes)) + 1]
    else:
        ov = np.bincount(lab[prior], minlength=lab.max() + 1)[1:]
        if ov.max() == 0:
            return None
        keep = [int(np.argmax(ov)) + 1]
        keep += [i + 1 for i in np.nonzero((ov > 0.5 * sizes) & (sizes > 2000))[0]
                 if i + 1 not in keep]
    core = np.isin(lab, keep)
    small = ndi.binary_fill_holes(ndi.binary_closing(core, morphology.disk(6)))
    big = ndi.binary_closing(core, morphology.disk(25)) & tissue
    cb = ndi.binary_fill_holes(small | big) & (tissue | small)
    cb = ndi.binary_opening(cb, morphology.disk(2))
    lab2 = measure.label(cb)
    if lab2.max() == 0:
        return None
    sizes = np.bincount(lab2.ravel())[1:]
    return np.isin(lab2, [i + 1 for i in np.nonzero(sizes > 0.05 * sizes.max())[0]])


def segment_layers(g1, env):
    """g1: lightly smoothed gray window; env: envelope. Returns pia mask,
    label image, fissure skeleton, class thresholds."""
    t = filters.threshold_multiotsu(g1[env], 4)
    cls = np.digitize(g1, t)                       # 0 EGL,1 IGL,2 WM,3 ML/gap
    cls = ndi.median_filter(cls, size=3)
    depth_env = ndi.distance_transform_edt(env) * PX_UM

    # EGL network attached to the envelope surface
    egl0 = morphology.remove_small_objects(env & (cls == 0), max_size=40)
    l0 = measure.label(egl0)
    ids = np.unique(l0[egl0 & (depth_env <= 2 * PX_UM)])
    egl_s = np.isin(l0, ids[ids > 0])

    # fissures: medial line of the (gap-closed) EGL network, deep parts only
    eglc = ndi.binary_closing(egl_s, morphology.disk(2)) & env
    sk = morphology.skeletonize(eglc)
    # skeleton deeper than a single surface EGL band (mid-band sits ~16 um deep)
    cand = sk & (depth_env > FISSURE_EXT_DEPTH_UM)
    lc = measure.label(cand, connectivity=2)
    fiss_full = np.zeros_like(env)
    surf_r, surf_c = np.nonzero(env & ~ndi.binary_erosion(env))
    for p in measure.regionprops(lc):
        rr, cc = p.coords[:, 0], p.coords[:, 1]
        if depth_env[rr, cc].max() < FISSURE_MIN_DEPTH_UM or p.area < 8:
            continue
        fiss_full[rr, cc] = True
        # join the shallowest point to the nearest surface pixel (fissure mouth)
        k = int(np.argmin(depth_env[rr, cc]))
        r0, c0 = rr[k], cc[k]
        j = int(np.argmin((surf_r - r0) ** 2 + (surf_c - c0) ** 2))
        n = int(max(abs(surf_r[j] - r0), abs(surf_c[j] - c0))) + 1
        fiss_full[np.linspace(r0, surf_r[j], n).round().astype(int),
                  np.linspace(c0, surf_c[j], n).round().astype(int)] = True
    cut = ndi.binary_dilation(fiss_full, structure=np.ones((3, 3)))
    pia = env & ~cut
    pia = ndi.binary_opening(pia, np.ones((2, 2)))
    # cuts that only shave off slivers are wrong: undo them
    lp = measure.label(pia)
    if lp.max() > 1:
        sz = np.bincount(lp.ravel())[1:]
        small = np.isin(lp, np.nonzero(sz < 0.02 * sz.max())[0] + 1)
        restore = ndi.binary_dilation(small, iterations=2) & env & cut
        pia |= restore
        fiss_full &= ~ndi.binary_dilation(restore, iterations=1)

    depth = ndi.distance_transform_edt(pia) * PX_UM
    lab = np.zeros(env.shape, np.uint8)
    lab[pia] = LBL["ML"]
    l1 = measure.label(pia & (cls == 0))
    ids = np.unique(l1[pia & (cls == 0) & (depth <= 3 * PX_UM)])
    egl = np.isin(l1, ids[ids > 0])
    igl = pia & ((cls == 1) | ((cls == 0) & ~egl))
    wm = pia & (cls == 2)
    # partial-volume rim along EGL is not IGL/WM
    near_egl = ndi.binary_dilation(egl, iterations=3) & (depth <= EGL_ADJ_UM)
    igl &= ~(near_egl & (cls >= 1))
    wm &= ~near_egl
    igl = morphology.remove_small_objects(igl, max_size=60)
    wm = morphology.remove_small_objects(wm, max_size=150)
    lab[igl] = LBL["IGL"]
    lab[wm] = LBL["WM"]
    lab[egl] = LBL["EGL"]
    return pia, lab, fiss_full, t


def envelope_adult(gs, prior):
    """IGL-based envelope (P14, adult).
    dense (dark) pixels are opened to strip thin lines (EGL at P14); the thick
    dark component overlapping the prior is the IGL(+WM at P14) core; envelope
    = tissue within ML_MAX_UM of it.  Pale regions not touching background are
    WM (adult arbor vitae).  Returns envelope, core, IGL, thin-dark, thresholds."""
    t_dense, t_mid, t_bg = filters.threshold_multiotsu(gs, 4)
    tissue = morphology.remove_small_objects(gs < t_bg, max_size=2000)
    dense = gs < t_dense
    dense_o = morphology.remove_small_objects(
        ndi.binary_opening(dense, morphology.disk(OPEN_PX)), max_size=200)
    thin = dense & ~ndi.binary_dilation(dense_o, iterations=1)
    lab = measure.label(dense_o)
    if lab.max() == 0:
        return None
    sizes = np.bincount(lab.ravel())[1:]
    ov = np.bincount(lab[prior], minlength=lab.max() + 1)[1:]
    if ov.max() == 0:
        return None
    keep = [int(np.argmax(ov)) + 1]
    keep += [i + 1 for i in np.nonzero((ov > 0.5 * sizes) & (sizes > 2000))[0] if i + 1 not in keep]
    igl = ndi.binary_closing(np.isin(lab, keep), morphology.disk(2))
    near = ndi.distance_transform_edt(~igl) * PX_UM <= ML_MAX_UM
    env = tissue & (near | igl)
    env = ndi.binary_opening(env, morphology.disk(2))
    le = measure.label(env)
    env = le == (np.argmax(np.bincount(le.ravel())[1:]) + 1)
    env = ndi.binary_fill_holes(env)
    igl &= env
    core = igl.copy()
    if WM_KIND == "pale":
        # WM (arbor vitae) is as pale as ML but sits inside the IGL ribbons:
        # close small IGL gaps (r = WM_CLOSE_UM) and fill.  Checked on P56 plate
        # 201 against the atlas: Dice 0.87 with atlas WM, ~0 leakage into ML.
        r = max(1, int(round(WM_CLOSE_UM / PX_UM)))
        wm = ndi.binary_fill_holes(ndi.binary_closing(igl, morphology.disk(r))) & env & ~igl
        core |= morphology.remove_small_objects(wm, max_size=100)
    return env, core, igl, thin & env, (t_dense, t_mid, t_bg)


def segment_adult(g1, env, core, igl, t):
    """IGL-based layers + fissures (P14, adult).
    Fissures = medial axis of the ML where the nearest boundary on both sides is
    IGL core (not the outer surface): only exists where two folia face each
    other.  Labels: IGL 4, WM 6, ML 5, EGL 3 (P14 only)."""
    D = ndi.distance_transform_edt(~core) * PX_UM          # distance to IGL core
    Dout = ndi.distance_transform_edt(env) * PX_UM         # distance to envelope edge
    ml = env & ~core
    ax = morphology.skeletonize(ml)
    fa = ax & (Dout > D + FISS_MARGIN_UM) & (D >= 15.0)
    surf_r, surf_c = np.nonzero(env & ~ndi.binary_erosion(env))
    fiss = np.zeros_like(env)
    for p in measure.regionprops(measure.label(fa, connectivity=2)):
        if p.area * PX_UM < 100:
            continue
        rr, cc = p.coords[:, 0], p.coords[:, 1]
        fiss[rr, cc] = True
        k = int(np.argmin(Dout[rr, cc]))                    # mouth end
        r0, c0 = rr[k], cc[k]
        j = int(np.argmin((surf_r - r0) ** 2 + (surf_c - c0) ** 2))
        n = int(max(abs(surf_r[j] - r0), abs(surf_c[j] - c0))) + 1
        fiss[np.linspace(r0, surf_r[j], n).round().astype(int),
             np.linspace(c0, surf_c[j], n).round().astype(int)] = True
    pia = env & ~ndi.binary_dilation(fiss, structure=np.ones((3, 3)))
    pia = ndi.binary_opening(pia, np.ones((2, 2)))
    lab = np.zeros(env.shape, np.uint8)
    lab[pia] = LBL["ML"]
    lab[pia & core] = LBL["IGL"]
    if WM_KIND == "pale":
        wm = pia & core & ~igl
        lab[morphology.remove_small_objects(wm, max_size=150)] = LBL["WM"]
    # WM_KIND "none" (P14): glia-rich WM is as dark as the IGL in Nissl (grey level
    # flat with depth inside the core), so the core is labelled IGL incl. WM
    if EGL_FLAG and ml.any():
        t_egl = 0.5 * (t[0] + np.median(g1[ml]))           # between IGL and ML grey
        egl = morphology.remove_small_objects(pia & ml & (g1 < t_egl), max_size=15)
        lab[egl] = LBL["EGL"]
    return pia, lab, fiss


def base_cut(env, egl):
    """Brainstem / peduncle cut = straight line from lobule I to lobule X.
    Implemented as envelope & convex hull of the main EGL network: the EGL
    wraps the cerebellum everywhere except the base, so the hull edge across the
    base joins the lingula tip to the nodulus tip.  Returns keep-mask, info."""
    from scipy.spatial import ConvexHull
    l = measure.label(egl, connectivity=2)
    if l.max() == 0:
        return np.ones_like(env), dict(skipped="no EGL")
    sz = np.bincount(l.ravel())[1:]
    main = np.isin(l, np.nonzero((sz >= HULL_MIN_COMP_FRAC * sz.max()) & (sz >= 100))[0] + 1)
    hull = morphology.convex_hull_image(main)
    removed = env & ~hull
    if removed.sum() > BASE_MAX_REMOVE_FRAC * env.sum():
        return np.ones_like(env), dict(skipped="cut would remove too much",
                                       removed_mm2=float(removed.sum() * (PX_UM / 1000) ** 2))
    pts = np.argwhere(main).astype(float)
    hv = pts[ConvexHull(pts).vertices]
    e = np.roll(hv, -1, axis=0) - hv
    k = int(np.argmax(np.hypot(e[:, 0], e[:, 1])))
    A, B = hv[k], hv[(k + 1) % len(hv)]
    return hull, dict(A_px=[float(A[0]), float(A[1])], B_px=[float(B[0]), float(B[1])],
                      chord_um=float(np.hypot(*(B - A)) * PX_UM),
                      removed_mm2=float(removed.sum() * (PX_UM / 1000) ** 2),
                      removed_egl_frac=float((removed & egl).sum() / max(egl.sum(), 1)))


def boundaries(pia, lab):
    b = np.zeros(pia.shape, np.uint8)
    b[pia & ~ndi.binary_erosion(pia)] = 1
    for code in (LBL["EGL"], LBL["IGL"], LBL["WM"]):
        m = lab == code
        e = m & ~ndi.binary_erosion(m) & (b == 0)
        if code == LBL["EGL"]:
            e &= ~ndi.binary_dilation(b == 1)   # inner EGL edge only
        b[e] = code
    return b


def track_shift(prev_g, prev_mask, cur_g):
    r0, r1, c0, c1 = bbox_of(prev_mask, 40, prev_g.shape)
    tpl = prev_g[r0:r1:DS, c0:c1:DS]
    img = cur_g[::DS, ::DS]
    pr, pc = max(0, tpl.shape[0] - img.shape[0]), max(0, tpl.shape[1] - img.shape[1])
    if pr or pc:
        img = np.pad(img, ((0, pr), (0, pc)), constant_values=255)
    res = match_template(img, tpl)
    k = np.unravel_index(res.argmax(), res.shape)
    return (k[0] * DS - r0, k[1] * DS - c0), float(res.max())


def shift_mask(mask, dr, dc, shape):
    out = np.zeros(shape, bool)
    ys, xs = np.nonzero(mask)
    ys, xs = ys + dr, xs + dc
    ok = (ys >= 0) & (ys < shape[0]) & (xs >= 0) & (xs < shape[1])
    out[ys[ok], xs[ok]] = True
    return out


# ----------------------------------------------------------------------------
def main(age="P7", max_offset_um=2000.0):
    cfg = configure(age)
    max_offset = int(round(max_offset_um / SEC_UM))
    files = section_files()
    results = {}

    def run_section(sec, rgb, g_full, prior_full):
        """prior_full: bool mask in this image's coords (or seed window)."""
        w = bbox_of(prior_full, WIN, g_full.shape)
        gw = g_full[sl(w)]
        gs = ndi.gaussian_filter(gw, 1.5)
        g1 = ndi.gaussian_filter(gw, 1.0)
        pw = ndi.binary_dilation(prior_full[sl(w)], iterations=10)
        if MODE == "adult":
            res = envelope_adult(gs, pw)
            if res is None:
                return None
            env, core, igl, thin, t = res
        else:
            env = envelope(gs, pw)
        if env is None or env.sum() < 2000:
            return None
        sol = env.sum() / max(morphology.convex_hull_image(env).sum(), 1)
        if prev_area[0] and (env.sum() > 1.15 * prev_area[0] or sol < 0.85 * prev_sol[0]):
            # sudden growth / appendage (choroid plexus, brainstem): clip to prior
            env &= ndi.binary_dilation(prior_full[sl(w)], iterations=15)
            env = ndi.binary_opening(env, morphology.disk(3))
            lab_e = measure.label(env)
            env = lab_e == (np.argmax(np.bincount(lab_e.ravel())[1:]) + 1)
            sol = env.sum() / max(morphology.convex_hull_image(env).sum(), 1)
        prev_sol[0] = sol
        cb_full = np.zeros(g_full.shape, bool)
        cb_full[sl(w)] = env                 # tracking uses the uncut envelope
        if MODE == "adult":
            core &= env
            igl &= env
            hull_src = ndi.binary_dilation(core, iterations=int(ML_TYP_UM / PX_UM)) & env
            keep, cut_info = base_cut(env, hull_src)
            # note: WM is reliable only within ~0.6 mm of the midline at P56
            # (Dice 0.82-0.88 vs atlas); laterally the arbor vitae/nuclei open
            # toward the peduncles and are missed - use the atlas there.
            pia, lab, fiss = segment_adult(g1, env, core, igl, t)
        else:
            pia, lab, fiss, t = segment_layers(g1, env)
            keep, cut_info = base_cut(env, lab == LBL["EGL"])
        if cut_info is not None and "skipped" not in cut_info:
            env &= keep
            lab_e = measure.label(env)
            szs = np.bincount(lab_e.ravel())[1:]
            env = np.isin(lab_e, np.nonzero(szs > 0.05 * szs.max())[0] + 1)
            pia &= env
            fiss &= env
            lab[~pia] = 0
        b = bbox_of(env, MARGIN, env.shape)
        crop_rgb = rgb[sl(w)][sl(b)]
        lab_c, pia_c, fiss_c = lab[sl(b)], pia[sl(b)], fiss[sl(b)]
        bnd = boundaries(pia_c, lab_c)
        Image.fromarray(crop_rgb).save(OUT / f"sections/sec{sec:03d}_crop.tif", compression="tiff_lzw")
        Image.fromarray(lab_c).save(OUT / f"sections/sec{sec:03d}_labels.tif", compression="tiff_lzw")
        Image.fromarray(bnd).save(OUT / f"sections/sec{sec:03d}_boundaries.tif", compression="tiff_lzw")
        a = lambda m: float(m.sum() * (PX_UM / 1000) ** 2)
        hull = morphology.convex_hull_image(pia_c)
        cont = measure.find_contours(pia_c.astype(float), 0.5)
        results[sec] = dict(
            sec=sec, offset_um=(sec - MID_SEC) * SEC_UM,
            crop_origin_px=[int(w[0] + b[0]), int(w[2] + b[2])],
            crop_shape=list(lab_c.shape), class_thr=[round(float(v), 1) for v in t],
            cb_mm2=a(pia_c), egl_mm2=a(lab_c == 3), ml_mm2=a(lab_c == 5),
            igl_mm2=a(lab_c == 4), wm_mm2=a(lab_c == 6),
            solidity_env=float(env.sum() / max(morphology.convex_hull_image(env).sum(), 1)),
            solidity_pia=float(pia_c.sum() / max(hull.sum(), 1)),
            fissure_len_mm=float(fiss_c.sum() * PX_UM / 1000),
            n_fissures=int(sum(1 for p in measure.regionprops(measure.label(fiss_c, connectivity=2)) if p.area * PX_UM >= 100)),
            n_parts=int(measure.label(pia_c).max()),
            base_cut=(cut_info if (cut_info is None or "skipped" in cut_info) else dict(
                cut_info, A_px=[cut_info["A_px"][0] - b[0], cut_info["A_px"][1] - b[2]],
                B_px=[cut_info["B_px"][0] - b[0], cut_info["B_px"][1] - b[2]])),
            _pia=(max(cont, key=len) * PX_UM) if cont else np.zeros((0, 2)),
        )
        return cb_full

    prev_area = [0]
    prev_sol = [1.0]
    # seed
    sd = cfg["seed"]
    SEED_SEC = sd["sec"]
    rgb = np.asarray(Image.open(files[SEED_SEC]).convert("RGB"))
    g = rgb.astype(float).mean(2)
    if sd["kind"] == "crop":
        crop_shape = Image.open(SRC / f"cerebellum_crop_{sd['img_id']}.jpg").size[::-1]
        seed = np.zeros(g.shape, bool)
        r0, c0 = sd["offset"]
        seed[r0 + 40:r0 + crop_shape[0] - 40, c0 + 40:c0 + crop_shape[1] - 40] = True
        w = bbox_of(seed, 0, g.shape)
        env0 = envelope(ndi.gaussian_filter(g[sl(w)], 1.5), None)
    else:
        atlas = atlas_cb_mask(sd["atlas_image_id"], SEED_SEC, g.shape, sd.get("ontology", "admba"))
        print(f"atlas isR+r1R seed: {atlas.sum() * (PX_UM / 1000) ** 2:.2f} mm2", flush=True)
        w = bbox_of(atlas, WIN, g.shape)
        env0 = envelope(ndi.gaussian_filter(g[sl(w)], 1.5),
                        ndi.binary_erosion(atlas[sl(w)], iterations=10))
    seed_mask = np.zeros(g.shape, bool)
    seed_mask[sl(w)] = env0
    cb = run_section(SEED_SEC, rgb, g, seed_mask)
    prev_area[0] = int(cb.sum())
    seed_state = (g, cb, prev_area[0], prev_sol[0])
    print(f"seed sec {SEED_SEC}: cb {results[SEED_SEC]['cb_mm2']:.2f} mm2", flush=True)

    for direction in (+1, -1):
        g_prev, cb_prev, prev_area[0], prev_sol[0] = seed_state
        s = SEED_SEC + direction
        while abs(s - MID_SEC) <= max_offset:
            if s not in files:
                s += direction
                if s < min(files) or s > max(files):
                    break
                continue
            rgb = np.asarray(Image.open(files[s]).convert("RGB"))
            g = rgb.astype(float).mean(2)
            (dr, dc), score = track_shift(g_prev, cb_prev, g)
            prior = shift_mask(cb_prev, dr, dc, g.shape)
            if not prior.any():
                print(f"sec {s}: lost track"); break
            cb = run_section(s, rgb, g, prior)
            if cb is None:
                print(f"sec {s}: cerebellum not found (score {score:.2f}); stop"); break
            results[s]["track_score"] = round(score, 3)
            g_prev, cb_prev, prev_area[0] = g, cb, int(cb.sum())
            s += direction

    done = sorted(results)
    stats = [{k: v for k, v in results[s].items() if not k.startswith("_")} for s in done]
    (OUT / f"{AGE.lower()}_nissl_stats.json").write_text(json.dumps(stats, indent=1))
    np.savez_compressed(OUT / f"{AGE.lower()}_nissl_contours.npz", secs=np.array(done),
                        offset_um=np.array([results[s]["offset_um"] for s in done]),
                        **{f"pia_{s}": results[s]["_pia"] for s in done})
    make_figures(results, done)

    print("\n sec  off_um  thr(EGL,IGL,WM)      CB   EGL   ML   IGL   WM  solEnv solPia nFiss fissmm parts | base: chord_um cut_mm2 EGLlost")
    for st in stats:
        if st["sec"] % 5 == 0 or st["sec"] == SEED_SEC:
            print(f" {st['sec']:3d} {st['offset_um']:+6.0f}  {str(st['class_thr']):20s}"
                  f"{st['cb_mm2']:5.2f} {st['egl_mm2']:5.2f} {st['ml_mm2']:5.2f} {st['igl_mm2']:5.2f} "
                  f"{st['wm_mm2']:5.2f}  {st['solidity_env']:.2f}   {st['solidity_pia']:.2f}  "
                  f"{st['n_fissures']:4d}  {st['fissure_len_mm']:5.2f}  {st['n_parts']:4d} | "
                  + ("none" if (st['base_cut'] is None or 'skipped' in st['base_cut']) else
                     f"{st['base_cut']['chord_um']:6.0f}  {st['base_cut']['removed_mm2']:5.2f}  "
                     f"{st['base_cut']['removed_egl_frac']:.3f}"))
    print(f"processed {len(stats)} sections: {done[0]}..{done[-1]}")


CMAP = matplotlib.colors.ListedColormap(
    [(1, 1, 1), (0.1, 0.1, 0.4), (0.4, 0.7, 0.95), (0, 0.75, 0.9), (1, 0.85, 0.1),
     (0.55, 0.9, 0.55), (1, 0.5, 0.1)])


def make_figures(results, done):
    def show(ax_img, ax_lab, s):
        crop = np.asarray(Image.open(OUT / f"sections/sec{s:03d}_crop.tif"))
        lab = np.asarray(Image.open(OUT / f"sections/sec{s:03d}_labels.tif"))
        bnd = np.asarray(Image.open(OUT / f"sections/sec{s:03d}_boundaries.tif"))
        ax_img.imshow(crop)
        yy, xx = np.nonzero(bnd == 1)
        ax_img.plot(xx, yy, ",", color="red")
        ax_lab.imshow(lab, cmap=CMAP, vmin=0, vmax=6, interpolation="nearest")
        for a in (ax_img, ax_lab):
            a.set_xticks([]); a.set_yticks([])
        ax_img.set_title(f"sec {s} ({results[s]['offset_um']:+.0f} um)", fontsize=10)

    pick = []
    for o in (0, 5, 10, 20, 35, 50):
        q = min(done, key=lambda s: abs(s - (MID_SEC + o)))
        if q not in pick:
            pick.append(q)
    fig, axes = plt.subplots(2, len(pick), figsize=(18, 9), dpi=100, squeeze=False)
    for j, s in enumerate(pick):
        show(axes[0, j], axes[1, j], s)
    axes[0, 0].set_ylabel("Nissl + pia (red)", fontsize=9)
    axes[1, 0].set_ylabel("cyan EGL / green ML\nyellow IGL / orange WM", fontsize=9)
    fig.suptitle(f"Allen {AGE} Nissl, sagittal: cerebellum outline + layers (anterior left, dorsal up)")
    fig.tight_layout()
    fig.savefig(OUT / f"{AGE.lower()}_nissl_overview.png")
    plt.close(fig)

    s = min(done, key=lambda q: abs(q - MID_SEC))
    fig, axes = plt.subplots(1, 2, figsize=(18, 9), dpi=100)
    show(axes[0], axes[1], s)
    fig.tight_layout()
    fig.savefig(OUT / f"{AGE.lower()}_nissl_midline.png")
    plt.close(fig)
    for f in (f"{AGE.lower()}_nissl_overview.png", f"{AGE.lower()}_nissl_midline.png"):
        print(f, Image.open(OUT / f).size)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "P7",
         float(sys.argv[2]) if len(sys.argv) > 2 else 2000.0)
