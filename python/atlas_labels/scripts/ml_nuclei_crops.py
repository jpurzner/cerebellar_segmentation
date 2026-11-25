"""
ML nuclei on the old repo's 20x folium crops (2026-11-24): the same segmentation and IGL-perpendicular angle as
ml_nuclei.py, with the layers classified from the stains, because these crops (and the 20x wild type, slide
position 3: s1_3, s2_3, s3_3) have no labels.

Layers per crop, on a 4 um grid: log DAPI and log NeuN smoothed over 4 um, a 3-component Gaussian mixture; the
component with the highest NeuN is IGL, the denser of the other two EGL, the sparser ML (the crops have no empty
background: a 4th component took the ML). IGL and EGL pieces under 2000 um^2 become ML candidates; ML is kept only
where the distances to the EGL and to the IGL add up to <= ml_reach_um (130 um): between the two, which fissure
lumen, meninges and white matter are not.

    python scripts/ml_nuclei_crops.py           # -> outputs/quant/ml_nuclei_crops.csv.gz, ml_nuclei_crops_sections.csv
    python scripts/ml_nuclei_crops.py --map F   # ASCII map of one crop's classes (check)
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi
from sklearn.mixture import GaussianMixture

sys.path.insert(0, str(Path(__file__).parent))
import ml_nuclei as M  # noqa: E402
import quant_layers as QL  # noqa: E402

PIECES = Path.home() / "Dropbox/images/edu_repeat/p27/pieces"
OUT = QL.OUT
PX = 0.5119
G = M.G
CPRM = dict(sm_um=4.0, min_layer_um2=2000.0, ml_reach_um=130.0, pcl_um=16.0)


def crop_list():
    return sorted(PIECES.glob("2018_05_22_*_fused_crop_*.tif"))


def short_of(p):
    return re.search(r"2018_05_22_(s\d_\d)_", p.name).group(1)


def load(p):
    a = tifffile.imread(p).astype(np.float32)
    return dict(p27=a[0], neun=a[1], nuclei=a[2])


def stain_layers(St, px=PX, prm=CPRM):
    z = px / G
    sm = prm["sm_um"] / px
    Dg = ndi.zoom(ndi.gaussian_filter(np.log(np.maximum(St["nuclei"], 1)), sm), z, order=1)
    Ng = ndi.zoom(ndi.gaussian_filter(np.log(np.maximum(St["neun"], 1)), sm), z, order=1)
    X = np.stack([Dg.ravel(), Ng.ravel()], 1)
    gm = GaussianMixture(3, random_state=0, n_init=3).fit(X[::3])
    c = gm.predict(X).reshape(Dg.shape)
    mu = gm.means_
    igl = int(np.argmax(mu[:, 1]))
    rest = [k for k in range(3) if k != igl]
    egl = max(rest, key=lambda k: mu[k, 0])
    ml = min(rest, key=lambda k: mu[k, 0])
    bg = -1
    lab = np.full(Dg.shape, 3, np.uint8)
    lab[c == egl] = 2
    lab[c == igl] = 4
    min_cells = prm["min_layer_um2"] / G ** 2
    for code in (2, 4):
        lb, n = ndi.label(lab == code)
        if n:
            small = np.bincount(lb.ravel(), minlength=n + 1) < min_cells
            lab[(lb > 0) & small[lb]] = 3
    # dense, NeuN-low cells within pcl_um of the IGL are the Purkinje cell layer or white matter, not EGL (the EGL
    # faces the IGL across the ML)
    near_igl = ndi.distance_transform_edt(lab != 4) * G <= prm["pcl_um"]
    lab[(lab == 2) & near_igl] = 4
    lb, n = ndi.label(lab == 2)
    if n:
        small = np.bincount(lb.ravel(), minlength=n + 1) < min_cells
        lab[(lb > 0) & small[lb]] = 3
    # ML only between an EGL and an IGL: the two within ml_reach_um together, and in opposite directions (fissure
    # lumen, meninges and white matter between two IGLs are not)
    d_igl = ndi.distance_transform_edt(lab != 4) * G
    d_egl = ndi.distance_transform_edt(lab != 2) * G
    gi = np.gradient(ndi.gaussian_filter(d_igl, 2.0))
    ge = np.gradient(ndi.gaussian_filter(d_egl, 2.0))
    dot = (gi[0] * ge[0] + gi[1] * ge[1]) / np.maximum(np.hypot(*gi) * np.hypot(*ge), 1e-6)
    far = (lab == 3) & ((d_igl + d_egl > prm["ml_reach_um"]) | (dot > -0.3))
    lab[far] = 1
    return lab, mu, (egl, ml, igl, bg)


def run():
    rows, infos = [], []
    for p in crop_list():
        s = short_of(p)
        St = load(p)
        lab, mu, _ = stain_layers(St)
        df, keep, extra = M.segment_arrays(St, PX, lab)
        animal, g, mag = QL.genotype(s)
        d = df[keep].drop(columns=["area", "axis_major_length", "axis_minor_length", "in_ml"]).assign(
            crop=p.stem.split("fused_crop_")[1], section=s, animal=animal, genotype=g, imaging="20x crop")
        rows.append(d)
        frac = {k: round(float((lab == c).mean()), 3) for c, k in ((2, "egl"), (3, "ml"), (4, "igl"))}
        infos.append(dict(crop=p.name, section=s, animal=animal, genotype=g, kept=int(keep.sum()), **frac))
        print(p.name, g, int(keep.sum()), frac, flush=True)
    N = pd.concat(rows, ignore_index=True)
    N.to_csv(OUT / "ml_nuclei_crops.csv.gz", index=False)
    pd.DataFrame(infos).to_csv(OUT / "ml_nuclei_crops_info.csv", index=False)
    return N


def amap(name):
    p = next(x for x in crop_list() if name in x.name)
    St = load(p)
    lab, mu, ids = stain_layers(St)
    print(p.name, "components (logDAPI, logNeuN):", np.round(mu, 2).tolist(), "egl/ml/igl/bg =", ids)
    H, W = lab.shape
    sy, sx = max(1, H // 30), max(1, W // 90)
    ch = {0: " ", 1: ".", 2: "e", 3: "-", 4: "I"}
    for r in range(0, H, sy):
        print("".join(ch[int(v)] for v in lab[r, ::sx]))


if __name__ == "__main__":
    if sys.argv[1:2] == ["--map"]:
        amap(sys.argv[2])
    else:
        run()
