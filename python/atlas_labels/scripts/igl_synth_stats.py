"""
Summary statistics of an IGL (real or synthetic), computed identically for every source,
and the real reference sets.

    python scripts/igl_synth_stats.py      # -> outputs/igl_synth/real_stats.json

A section is a dict of boolean masks on the 4 um grid: tissue, egl, ml, igl, wm (+ anom,
unass, empty for synthetic and Allen sections).

Reference sets
  IF control       the five intact-IGL controls (s1_2 gold, s5_F, s5_G, s6_C, s6_G)
  IF other         the other 10x controls (s4_A, s4_C, s4_F, s5_A, s5_B, s5_C, s5_E)
  Allen P7         Allen Developing Mouse Brain P7 Nissl, sagittal sections within
                   600 um of the midline (layers from Nissl density classes)
"""
import json
import sys
from pathlib import Path

import numpy as np
import tifffile
from scipy import ndimage as ndi
from skimage import measure, morphology

sys.path.insert(0, str(Path(__file__).parent))
import igl_model as M  # noqa: E402
import igl_ribbon as RB  # noqa: E402
import igl_whole_model as W  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/igl_synth"
G = RB.GRID
ALLEN = ROOT / "outputs/allen_p7"
ALLEN_MAX_OFFSET = 600.0

# statistics used to compare synthetic with real IGLs (name, description)
STATS = [
    ("igl_frac", "IGL area / cerebellum area"),
    ("egl_frac", "EGL area / cerebellum area"),
    ("ml_frac", "ML area / cerebellum area"),
    ("wm_frac", "white matter area / cerebellum area"),
    ("solidity", "cerebellum area / convex hull area (folding)"),
    ("igl_len", "IGL ribbon length / sqrt(cerebellum area)"),
    ("w_med", "IGL width, median (um)"),
    ("w_iqr", "IGL width, interquartile range (um)"),
    ("d_med", "IGL centre depth below the EGL, median (um)"),
    ("d_iqr", "IGL centre depth below the EGL, IQR (um)"),
    ("ml_med", "ML thickness (closest IGL-EGL approach), median (um)"),
    ("egl_med", "EGL thickness, median (um)"),
    ("egl_contact", "share of ribbon where the IGL touches the EGL"),
    ("wm_contact", "share of ribbon lying on white matter"),
    ("kap_med", "|curvature| of the ribbon, median (per mm)"),
    ("convex", "share of ribbon wrapping a folium tip (curvature > 5/mm)"),
    ("concave", "share of ribbon wrapping a fissure base (curvature < -5/mm)"),
    ("ends_per_mm", "ribbon free ends per mm of ribbon"),
    ("pieces", "IGL pieces >= 2000 um2"),
    ("corr_w", "width correlation length along the ribbon (um)"),
]
NAMES = [s[0] for s in STATS]


def _thickness_med(mask):
    if mask.sum() < 20:
        return np.nan
    sk = morphology.skeletonize(mask)
    e = ndi.distance_transform_edt(mask)
    v = e[sk]
    return float(np.median(2 * v) * G) if len(v) else np.nan


def summarize(g):
    T = int(g["tissue"].sum())
    if T < 1000 or g["igl"].sum() < 200:
        return None
    rib = RB.describe(g)
    if not rib["branches"]:
        return None
    P = W.profile(g, {}, rib)
    if P is None or len(P["F"]) < 50:
        return None
    F = P["F"]
    w, d, ml, lift, kap = F[:, 0], F[:, 1], F[:, 2], F[:, 3], F[:, 4]
    hull = morphology.convex_hull_image(g["tissue"])
    lab, n = ndi.label(rib["igl"])
    sizes = np.bincount(lab.ravel())[1:] * G ** 2
    cw = []
    for b in np.unique(P["branch"]):
        sel = P["branch"] == b
        if sel.sum() >= 150:
            cw.append(M.corr_length(w[sel], 100))
    s = dict(tissue_mm2=T * G ** 2 / 1e6,
             igl_frac=g["igl"].sum() / T, egl_frac=g["egl"].sum() / T, ml_frac=g["ml"].sum() / T, wm_frac=g["wm"].sum() / T,
             solidity=T / max(hull.sum(), 1), igl_len=rib["length_mm"] * 1000 / np.sqrt(T * G ** 2),
             w_med=np.median(w), w_iqr=np.subtract(*np.percentile(w, [75, 25])),
             d_med=np.median(d), d_iqr=np.subtract(*np.percentile(d, [75, 25])),
             ml_med=np.median(ml), egl_med=_thickness_med(g["egl"]),
             egl_contact=float((ml < 8).mean()), wm_contact=float(1 - lift.mean()),
             kap_med=np.median(np.abs(kap)), convex=float((kap > 5).mean()), concave=float((kap < -5).mean()),
             ends_per_mm=len(rib["ends"]) / max(rib["length_mm"], 1e-3), pieces=int((sizes >= 2000).sum()),
             corr_w=float(np.nanmedian(cw)) if cw else np.nan)
    return {k: (float(v) if v is not None else np.nan) for k, v in s.items()}


def allen_grid(path):
    L = ndi.zoom(tifffile.imread(path), 3.96 / G, order=0)
    z = np.zeros(L.shape, bool)
    return dict(tissue=L > 0, egl=L == 3, igl=L == 4, ml=L == 5, wm=L == 6, anom=z, unass=z.copy())


def if_grid(name):
    for n, src, px, img, ch, geno in W.sections():
        if n == name:
            return RB.grid_labels(W.load_labels(src), px)
    raise KeyError(name)


IF_CONTROL = ["s1_2 (gold)", "s5_F", "s5_G", "s6_C", "s6_G"]
IF_OTHER = ["s4_A", "s4_C", "s4_F", "s5_A", "s5_B", "s5_C", "s5_E"]
IF_CKO = ["s1_4", "s1_5", "s2_2", "s2_4", "s2_5", "s3_2", "s3_5"]


def allen_sections():
    st = json.loads((ALLEN / "p7_nissl_stats.json").read_text())
    return [(f"Allen P7 {s['offset_um']:+.0f}", ALLEN / "sections" / f"sec{s['sec']}_labels.tif", s["offset_um"])
            for s in st if abs(s["offset_um"]) <= ALLEN_MAX_OFFSET]


def real_sets():
    for n in IF_CONTROL:
        yield "IF control", n, lambda n=n: if_grid(n)
    for n in IF_OTHER:
        yield "IF other", n, lambda n=n: if_grid(n)
    for n in IF_CKO:
        yield "IF EZH2 cKO", n, lambda n=n: if_grid(n)
    for n, p, off in allen_sections():
        yield "Allen P7", n, lambda p=p: allen_grid(p)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for group, name, load in real_sets():
        g = load()
        s = summarize(g)
        if s is None:
            print(name, "skipped")
            continue
        rows.append(dict(group=group, name=name, **s))
        print(f"{group:12s} {name:18s} " + " ".join(f"{k}={s[k]:.3g}" for k in ("igl_frac", "w_med", "d_med", "ml_med", "egl_med", "igl_len", "solidity")), flush=True)
    (OUT / "real_stats.json").write_text(json.dumps(rows, indent=1))
    # group summaries
    for grp in ("IF control", "IF other", "IF EZH2 cKO", "Allen P7"):
        R = [r for r in rows if r["group"] == grp]
        print(grp, len(R), {k: f"{np.nanmedian([r[k] for r in R]):.3g} [{np.nanpercentile([r[k] for r in R], 10):.3g}-{np.nanpercentile([r[k] for r in R], 90):.3g}]" for k in NAMES})


if __name__ == "__main__":
    main()
