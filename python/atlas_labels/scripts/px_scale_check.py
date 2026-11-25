"""
Is the 10x pixel size right? (2026-11-25)  Three rulers, all read in the ASSUMED pixel sizes (10x 0.65, 20x 0.5119 um/px):

  1. packing   spacing of the packed nuclei in the IGL and the EGL: the first side peak of the radially averaged
               autocorrelation of DAPI inside the layer (high-passed at 12 um, contrast-normalised). A nucleus-to-
               nucleus distance does not depend on how big the section is or where it was cut.
  2. allen     similarity fit (rotation, scale, shift) of the Allen P7 Nissl template (fit_p7_template) to the final
               labels, every template section within +-1 mm of the midline; scale = true px / assumed px if the
               section matches the template in size.
  3. size      labelled cerebellum area against the Allen P7 sections near the midline.

    python scripts/px_scale_check.py [packing|allen|all]
Writes outputs/px_scale/packing.csv, allen_fit.csv, summary.txt.
"""
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi

sys.path.insert(0, str(Path(__file__).parent))
import ml_nuclei as MN  # noqa: E402
import quant_layers as QL  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/px_scale"
SECTIONS = ["s1_2", "s1_4", "s1_5", "s2_2", "s2_4", "s2_5", "s3_2", "s3_5",
            "s4_A", "s4_C", "s4_F", "s5_A", "s5_B", "s5_C", "s5_E", "s5_F", "s5_G", "s6_C", "s6_G"]
LAYER = {"EGL": 2, "IGL": 4}


# --------------------------------------------------------------------------- #
# 1. packing
# --------------------------------------------------------------------------- #
def radial_ac(img, mask, px, rmax_um=20.0):
    """radially averaged autocorrelation of img inside mask (masked, normalised); r in assumed um."""
    m = mask.astype(np.float32)
    f = np.where(mask, img, 0).astype(np.float32)
    sh = [int(2 ** np.ceil(np.log2(s * 1.0 + rmax_um / px + 2))) for s in img.shape]
    F = np.fft.rfft2(f, sh)
    M = np.fft.rfft2(m, sh)
    num = np.fft.irfft2(F * np.conj(F), sh)
    den = np.fft.irfft2(M * np.conj(M), sh)
    R = int(np.ceil(rmax_um / px)) + 1
    idx = np.r_[0:R, -R + 1:0]
    num = num[np.ix_(idx, idx)]
    den = den[np.ix_(idx, idx)]
    ac = num / np.maximum(den, 1.0)
    ac /= ac[0, 0]
    yy, xx = np.meshgrid(np.r_[0:R, -R + 1:0], np.r_[0:R, -R + 1:0], indexing="ij")
    r = np.hypot(yy, xx) * px
    # radial profile: mean in rings one pixel wide, then on a common 0.25 um axis by interpolation
    # (rings narrower than a pixel leave empty bins at small r)
    k = np.round(r.ravel() / px).astype(int)
    ok = k <= R - 1
    n = np.bincount(k[ok], minlength=R)
    rp = np.bincount(k[ok], r.ravel()[ok], R) / np.maximum(n, 1)
    pp = np.bincount(k[ok], ac.ravel()[ok], R) / np.maximum(n, 1)
    keep = n > 0
    rr = np.arange(0.125, rmax_um, 0.25)
    return rr, np.interp(rr, rp[keep], pp[keep])


def side_peak(r, prof, lo=3.0, hi=16.0):
    """first local maximum of the profile beyond its first minimum (parabolic interpolation)."""
    sel = np.where((r >= lo) & (r <= hi))[0]
    p = prof[sel]
    i0 = int(np.argmin(p[: len(p) // 2 + 1]))
    for i in range(max(i0, 1), len(p) - 1):
        if p[i] >= p[i - 1] and p[i] >= p[i + 1]:
            a, b, c = p[i - 1], p[i], p[i + 1]
            d = 0.5 * (a - c) / (a - 2 * b + c) if (a - 2 * b + c) != 0 else 0.0
            j = sel[i]
            return float(r[j] + d * (r[1] - r[0])), float(b - p[i0]), float(r[sel[i0]])
    return np.nan, np.nan, np.nan


def packing(short):
    sid = MN.sid_of(short)
    St, px = MN.D.stains(sid)
    dapi = St["nuclei"].astype(np.float32)
    lab = MN.native_labels(tifffile.imread(MN.FINAL / f"{short}_labels.tif"), dapi.shape, px)
    hp = dapi - ndi.gaussian_filter(dapi, 12.0 / px)                       # high-pass at 12 um
    sd = np.sqrt(ndi.gaussian_filter(hp ** 2, 12.0 / px))
    z = hp / np.maximum(sd, 1e-3)
    rows, profs = [], {}
    for name, code in LAYER.items():
        m = ndi.binary_erosion(lab == code, iterations=max(1, int(round(4.0 / px))))   # 4 um away from the layer edge
        r, prof = radial_ac(z, m, px)
        pk, depth, rmin = side_peak(r, prof)
        animal, g, mag = QL.genotype(short)
        _, g_proj, _ = QL.genotype(short, "project")
        rows.append(dict(section=short, animal=animal, genotype=g, genotype_project=g_proj, imaging=mag, px_assumed=px, layer=name,
                         area_mm2=float(m.sum() * px * px / 1e6), spacing_um=pk, first_min_um=rmin, peak_height=depth))
        profs[f"{short}_{name}"] = prof
        profs["r"] = r
    return rows, profs


# --------------------------------------------------------------------------- #
# 2. Allen template fit
# --------------------------------------------------------------------------- #
FINAL_TO_TPL = {0: 0, 1: 1, 2: 3, 3: 5, 4: 4, 5: 6, 6: 9}     # final codes -> fit_p7_template codes (UNK 1, ANOM 9)


def allen_fit(short, max_off_um=1000.0, step=2):
    import fit_p7_template as FT
    lab = tifffile.imread(MN.FINAL / f"{short}_labels.tif")
    t = np.zeros(10, np.uint8)
    for k, v in FINAL_TO_TPL.items():
        t[k] = v
    lab = t[lab]
    FT.GT_PX, FT.GT_ORG = 4.0, 0.0
    tgt = FT.fractions(lab, 1)
    v = 1.0 - tgt[FT.CH.index(FT.ANOM)]
    sw = FT.orientation_sweep(tgt, v)
    th0, mirror = sw[0][1], sw[0][2]
    secs = [(s, o) for s, o in FT.tpl_sections() if abs(o) <= max_off_um][::step]
    with Pool(8, initializer=_init_ft, initargs=(tgt, v)) as pool:
        fits = pool.map(_fit_task, [(s, th0, mirror) for s, _ in secs], chunksize=2)
    off = dict(secs)
    rows = []
    for f in fits:
        rows.append(dict(section=short, sec=f["sec"], offset_um=off[f["sec"]], score=f["score"], scale=float(np.exp(f["params"][1])),
                         rot_deg=float(np.rad2deg(f["params"][0])), mirror=f["mirror"], tissue_dice=f["parts"][0]))
    return rows


def _init_ft(tgt, v):
    import fit_p7_template as FT
    FT.GT_PX, FT.GT_ORG = 4.0, 0.0
    FT._init(tgt, v)


def _fit_task(a):
    import fit_p7_template as FT
    return FT.fit_similarity(*a)


# --------------------------------------------------------------------------- #
def main(what="all"):
    OUT.mkdir(parents=True, exist_ok=True)
    if what in ("packing", "all"):
        rows, P = [], {}
        for s in SECTIONS:
            t0 = time.time()
            r, p = packing(s)
            rows += r
            P.update(p)
            print(s, [(x["layer"], round(x["spacing_um"], 2), round(x["peak_height"], 3)) for x in r], f"{time.time() - t0:.0f} s", flush=True)
        df = pd.DataFrame(rows)
        df.to_csv(OUT / "packing.csv", index=False)
        np.savez_compressed(OUT / "packing_profiles.npz", **P)
    if what in ("allen", "all"):
        rows = []
        for s in SECTIONS:
            t0 = time.time()
            r = allen_fit(s)
            rows += r
            b = sorted(r, key=lambda x: -x["score"])[:5]
            print(s, "top5 (offset, score, scale):", [(x["offset_um"], round(x["score"], 3), round(x["scale"], 3)) for x in b],
                  f"{time.time() - t0:.0f} s", flush=True)
            pd.DataFrame(rows).to_csv(OUT / "allen_fit.csv", index=False)


if __name__ == "__main__":
    main(*(sys.argv[1:] or ["all"]))
