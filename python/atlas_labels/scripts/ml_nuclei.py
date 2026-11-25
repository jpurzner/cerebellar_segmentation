"""
Shape and orientation of the nuclei in the molecular layer (2026-11-24), the new-label version of the old R
pipeline's segment_ml_nuclei (R_projects/cerebellar_layer_quant/R/segment_ml_nuclei.R).

The old function worked on hand-cut folium crops (rows = depth, pia horizontal), found the ML as the rows between
the p27 minimum and the IGL, and measured each nucleus's angle against the crop's row axis. Here the ML is the
final label (outputs/flat_fix/final), whole sections are used, and each nucleus's major axis is measured against
the local line perpendicular to the IGL: the gradient of the distance to the IGL label (smoothed over igl_sigma_um),
at the nucleus. 0 deg = radial (along that line, the migration direction), 90 deg = tangential (parallel to the
IGL). The angle to the perpendicular of the EGL (the old method's pial reference) is kept as well.

Per section (native resolution; 20x 0.5119 um/px, 10x 1.0239 um/px; the 10x labels and their 4 um grid were built
at 0.65 um/px, so grid units are scaled by quant_layers.um_scale and every um parameter is in true um):
1. DAPI minus its Gaussian background (bg_sigma_um), lightly smoothed (sm_um).
2. Foreground: above the local mean over adaptive_um (as the old adaptive threshold, window 21 px at 20x =
   10.7 um) by adaptive_k robust SDs of the ML's DAPI, and above the ML's own Otsu level x otsu_frac; holes
   filled, opened with a 0.5 um disc.
3. Touching nuclei split by a watershed on the distance transform (peaks >= min_sep_um apart).
4. Kept: area area_um2 range, solidity >= min_solidity, >= in_ml_frac of the pixels in the ML label, the centroid
   at least edge_um from the EGL and the IGL labels (the dense layers either side merge with ML nuclei).
5. NeuN+ as in the old code: >= neun_frac of the nucleus above the Otsu level of NeuN over the ML.
6. Per nucleus: area, eccentricity, major / minor axis (um), the angle to the IGL perpendicular and to the
   EGL perpendicular (0-90 deg), depth in the ML (0 = EGL side, 1 = IGL side), mean NeuN / p27 / DAPI.

    python scripts/ml_nuclei.py                  # all sections -> outputs/quant/ml_nuclei_<s>.csv.gz, ml_nuclei.csv
    python scripts/ml_nuclei.py --down s1_2 ...  # 20x sections resampled to the 10x 1.0239 um/px first (resolution test)
    python scripts/ml_nuclei.py --qc s4_C ...    # outline crops for review -> figures/quant/qc/
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi
from skimage.feature import peak_local_max
from skimage.filters import threshold_otsu
from skimage.measure import regionprops_table
from skimage.morphology import disk
from skimage.segmentation import watershed

sys.path.insert(0, str(Path(__file__).parent))
import igl_detect as D  # noqa: E402
import quant_layers as QL  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/quant"
FINAL = ROOT / "outputs/flat_fix/final"
G = 4.0
GA = D.GA
PRM = dict(bg_sigma_um=10.0, sm_um=0.4, adaptive_um=10.7, adaptive_k=0.25, otsu_frac=0.6, min_sep_um=2.6,
           area_um2=(14.0, 120.0), min_solidity=0.85, in_ml_frac=0.8, edge_um=3.0, neun_frac=0.3, igl_sigma_um=10.0,
           down_px=1.0239)


def sid_of(short):
    return next(s for s in GA.section_list() if GA.short(s) == short)


def grid_fields(short):
    """on the 4 um grid: labels, tissue, the distances to the IGL and EGL labels and their unit gradients (the
    perpendiculars to the IGL and to the EGL, smoothed over igl_sigma_um = 10 um: the 4 um label edges are
    jagged, and at 6 um the two disagreed by a median 13-21 deg in the ML; the distance to the outside of the tissue is no use: across a closed
    fissure the nearest outside can be far away)."""
    lab = tifffile.imread(FINAL / f"{short}_labels.tif")
    tis = tifffile.imread(FINAL / f"{short}_tissue.tif") > 0
    return (lab, tis) + fields_from_labels(lab, ks=QL.um_scale(short))


def fields_from_labels(lab, prm=PRM, ks=1.0):
    """distances to the IGL / EGL labels (true um) on the label grid (G built um = G * ks true um per cell) and
    their smoothed unit gradients."""
    d_igl = ndi.distance_transform_edt(lab != 4) * G * ks
    d_egl = ndi.distance_transform_edt(lab != 2) * G * ks
    s = prm["igl_sigma_um"] / (G * ks)

    def unit(d):
        gr, gc = np.gradient(ndi.gaussian_filter(d, s))
        n = np.hypot(gr, gc)
        return gr / np.maximum(n, 1e-6), gc / np.maximum(n, 1e-6)
    return d_igl, d_egl, unit(d_igl), unit(d_egl)


def native_labels(lab, shape, px):
    gy = np.minimum(np.round(np.arange(shape[0]) * px / G).astype(int), lab.shape[0] - 1)
    gx = np.minimum(np.round(np.arange(shape[1]) * px / G).astype(int), lab.shape[1] - 1)
    return lab[gy][:, gx]


def _odd(x):
    x = int(round(x))
    return x + 1 - x % 2


def segment(short, down=False, prm=PRM, keep_labels=False):
    t0 = time.time()
    sid = sid_of(short)
    St, px = D.stains(sid)                       # px as the labels were built (10x 0.65)
    ks = QL.um_scale(short)
    px = px * ks                                  # true um/px (10x 1.0239)
    St = {k: v.astype(np.float32) for k, v in St.items()}
    if down and px < prm["down_px"]:
        z = px / prm["down_px"]
        St = {k: ndi.zoom(ndi.gaussian_filter(v, 0.5 / z), z, order=1) for k, v in St.items()}
        px = prm["down_px"]
    lab, tis = grid_fields(short)[:2]
    df, keep, extra = segment_arrays(St, px, lab, prm, ks)
    animal, g, mag = QL.genotype(short)
    info = dict(section=short, animal=animal, genotype=g, imaging=mag, px_um=px, down=bool(down), ml_area_mm2=extra["ml_area_mm2"],
                objects=int(len(df)), kept=int(keep.sum()), neun_thr=extra["neun_thr"], seconds=round(time.time() - t0, 1))
    df = df[keep].drop(columns=["area", "axis_major_length", "axis_minor_length", "in_ml"]).assign(section=short, animal=animal, genotype=g, imaging=mag)
    if keep_labels:
        return df, info, dict(lbl=extra["lbl"], keep=set(df.label), dapi=St["nuclei"], neun=St["neun"], Ln=extra["Ln"], px=px)
    return df, info


def segment_arrays(St, px, lab, prm=PRM, ks=1.0):
    """nuclei of the ML (lab == 3 on the label grid, aligned to the native image's top-left corner): all objects
    with their measures, the keep mask, and extras. px = true um per image pixel; ks = true / built um of the
    label grid (one grid cell = G * ks true um)."""
    d_igl, d_egl, n_igl, n_egl = fields_from_labels(lab, prm, ks)
    H, W = St["nuclei"].shape
    pg = px / ks                                  # image pixel in grid units
    Ln = native_labels(lab, (H, W), pg)
    ml = Ln == 3
    um = lambda v: v / px
    dapi = St["nuclei"]
    bs = dapi - ndi.gaussian_filter(dapi, um(prm["bg_sigma_um"]))
    bs = ndi.gaussian_filter(bs, um(prm["sm_um"]))
    v = bs[ml]
    med = np.median(v)
    mad = np.median(np.abs(v - med)) * 1.4826
    w = _odd(um(prm["adaptive_um"]))
    loc = ndi.uniform_filter(bs, w)
    otsu = threshold_otsu(v[np.isfinite(v)])
    fg = (bs > loc + prm["adaptive_k"] * mad) & (bs > med + prm["otsu_frac"] * (otsu - med))
    # restrict the work to the ML and a margin around it
    near = ndi.binary_dilation(ml, iterations=int(um(6)) + 1)
    fg &= near
    fg = ndi.binary_fill_holes(fg)
    fg = ndi.binary_opening(fg, structure=disk(max(1, int(round(um(0.5))))))
    lb_, n_ = ndi.label(fg)
    fg &= (np.bincount(lb_.ravel(), minlength=n_ + 1) >= prm["area_um2"][0] / px ** 2 / 2)[lb_] & (lb_ > 0)
    dist = ndi.distance_transform_edt(fg)
    lbl0, _ = ndi.label(fg)
    pk = peak_local_max(ndi.gaussian_filter(dist, 0.7), min_distance=max(1, int(round(um(prm["min_sep_um"])))), labels=lbl0, exclude_border=False)
    markers = np.zeros(fg.shape, np.int32)
    markers[tuple(pk.T)] = np.arange(1, len(pk) + 1)
    lbl = watershed(-dist, markers, mask=fg)
    neun_thr = threshold_otsu(St["neun"][ml])
    props = regionprops_table(lbl, intensity_image=np.stack([ml.astype(np.float32), (St["neun"] > neun_thr).astype(np.float32),
                                                             St["neun"], St["p27"], dapi], -1),
                              properties=("label", "area", "centroid", "eccentricity", "orientation", "axis_major_length",
                                          "axis_minor_length", "solidity", "intensity_mean"))
    df = pd.DataFrame(props)
    df = df.rename(columns={"intensity_mean-0": "in_ml", "intensity_mean-1": "neun_pos_frac", "intensity_mean-2": "neun",
                            "intensity_mean-3": "p27", "intensity_mean-4": "dapi", "centroid-0": "r", "centroid-1": "c"})
    df["area_um2"] = df.area * px ** 2
    df["major_um"] = df.axis_major_length * px
    df["minor_um"] = df.axis_minor_length * px
    gr, gc = df.r * pg / G, df.c * pg / G
    co = np.stack([gr, gc])
    samp = lambda a: ndi.map_coordinates(a, co, order=1, mode="nearest")
    df["d_igl_um"], df["d_egl_um"] = samp(d_igl), samp(d_egl)
    df["ml_depth_frac"] = df.d_egl_um / np.maximum(df.d_egl_um + df.d_igl_um, 1e-6)
    ur, uc = np.cos(df.orientation), np.sin(df.orientation)
    nv = {}
    for nm, (nr, nc) in (("igl", n_igl), ("egl", n_egl)):
        a, b = samp(nr), samp(nc)
        nn = np.maximum(np.hypot(a, b), 1e-6)
        nv[nm] = (a / nn, b / nn)
        df[f"angle_{nm}_deg"] = np.rad2deg(np.arccos(np.clip(np.abs(ur * a + uc * b) / nn, 0, 1)))
    # how far the IGL and EGL perpendiculars disagree at the nucleus (0 = the ML runs straight between them)
    df["normal_dev_deg"] = np.rad2deg(np.arccos(np.clip(-(nv["igl"][0] * nv["egl"][0] + nv["igl"][1] * nv["egl"][1]), -1, 1)))
    lo, hi = prm["area_um2"]
    keep = ((df.area_um2 >= lo) & (df.area_um2 <= hi) & (df.solidity >= prm["min_solidity"]) & (df.in_ml >= prm["in_ml_frac"])
            & (df.d_igl_um >= prm["edge_um"]) & (df.d_egl_um >= prm["edge_um"]))
    df["neun_positive"] = df.neun_pos_frac >= prm["neun_frac"]
    return df, keep, dict(ml_area_mm2=round(float(ml.sum() * px * px / 1e6), 4), neun_thr=float(neun_thr), lbl=lbl, Ln=Ln)


def qc(short, n=6, size_um=110.0):
    """outline crops of ML stretches for the user to review: DAPI grey, kept nuclei outlined coloured by the angle
    to the IGL perpendicular (blue radial -> orange tangential), the IGL / EGL label edges, a short tick showing
    the IGL perpendicular. PNGs <= 1900 px."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from skimage.segmentation import find_boundaries
    df, info, A = segment(short, keep_labels=True)
    px = A["px"]
    out = ROOT / "figures/quant/qc"
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(0)
    h = int(size_um / px)
    # crops centred on random kept nuclei far apart
    pick = df.sample(min(len(df), 400), random_state=1)
    centres = []
    for _, r_ in pick.iterrows():
        if all(np.hypot(r_.r - a, r_.c - b) > 2 * h for a, b in centres):
            centres.append((r_.r, r_.c))
        if len(centres) >= n:
            break
    fig, axs = plt.subplots(2, 3, figsize=(7.2, 4.9))
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("ang", ["#2a78d6", "#c3c2b7", "#eb6834"])
    for ax, (cr, cc) in zip(axs.ravel(), centres):
        r0, c0 = int(max(cr - h / 2, 0)), int(max(cc - h / 2, 0))
        sl = (slice(r0, r0 + h), slice(c0, c0 + h))
        im = A["dapi"][sl]
        lo, hi = np.percentile(im, [1, 99.7])
        ax.imshow(np.clip((im - lo) / (hi - lo), 0, 1), cmap="gray", interpolation="nearest")
        L = A["Ln"][sl]
        for code, col in ((2, "#ffd24d"), (4, "#7fe08a")):
            ax.contour(L == code, levels=[0.5], colors=[col], linewidths=0.5)
        sub = df[(df.r >= r0) & (df.r < r0 + h) & (df.c >= c0) & (df.c < c0 + h)]
        lb = A["lbl"][sl]
        for _, q in sub.iterrows():
            m = lb == q.label
            b = find_boundaries(m, mode="inner")
            yy, xx = np.nonzero(b)
            ax.scatter(xx, yy, s=0.15, c=[cmap(q.angle_igl_deg / 90)], marker="s", linewidths=0)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(f"{short} at ({int(cr * px)}, {int(cc * px)}) µm: {len(sub)} nuclei", fontsize=6.5, loc="left")
        ax.plot([5, 5 + 20 / px], [h - 6, h - 6], color="white", lw=1.2)
        ax.text(5, h - 9, "20 µm", color="white", fontsize=5.5)
    fig.suptitle(f"{short} ML nuclei: outline colour = angle to the IGL perpendicular (blue 0° radial, grey 45°, orange 90° tangential); "
                 "yellow EGL edge, green IGL edge", fontsize=6.5)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.9, bottom=0.01, wspace=0.04, hspace=0.12)
    p = out / f"qc_{short}.png"
    fig.savefig(p, dpi=min(250, int(1900 / 7.2)))
    plt.close(fig)
    print("qc", p)


def main(argv):
    down = "--down" in argv
    if "--qc" in argv:
        for s in [a for a in argv if not a.startswith("--")] or ["s4_C", "s5_A"]:
            qc(s)
        return
    shorts = [a for a in argv if not a.startswith("--")] or [GA.short(s) for s in GA.section_list()]
    OUT.mkdir(parents=True, exist_ok=True)
    infos = []
    for s in shorts:
        df, info = segment(s, down=down)
        tag = "_down" if down else ""
        df.to_csv(OUT / f"ml_nuclei_{s}{tag}.csv.gz", index=False)
        pos = df[df.neun_positive]
        info.update(neun_pos_frac=round(float(df.neun_positive.mean()), 3), density_per_mm2=round(len(df) / max(info["ml_area_mm2"], 1e-9), 1),
                    area_med=round(float(df.area_um2.median()), 1), ecc_pos=round(float(pos.eccentricity.mean()), 3),
                    angle_igl_pos=round(float(pos.angle_igl_deg.mean()), 1), angle_egl_pos=round(float(pos.angle_egl_deg.mean()), 1), r_igl_egl=round(float(np.corrcoef(df.angle_igl_deg, df.angle_egl_deg)[0, 1]), 3),
                    radial30_pos=round(float((pos.angle_igl_deg < 30).mean()), 3))
        infos.append(info)
        print(json.dumps(info), flush=True)
    if len(shorts) > 1:
        pd.DataFrame(infos).to_csv(OUT / f"ml_nuclei_sections{'_down' if down else ''}.csv", index=False)


if __name__ == "__main__":
    main(sys.argv[1:])
