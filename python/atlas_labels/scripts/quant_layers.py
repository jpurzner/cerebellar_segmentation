"""
Layer size and stain intensity per section, Ezh2 cKO vs wild type (2026-11-24).

The figures this feeds are the new-segmentation versions of the paper's Figure 3 H and I (per-layer area,
per-layer p27 intensity) and of the old R pipeline's depth profiles (R_projects/cerebellar_layer_quant,
align_folia_images_connect -> plot_layer_boundary2). Everything comes from the final labels and the flattened
strips; nothing is re-segmented here.

Genotype (GENOTYPE below, with its sources):
  10x panoramas, 2017-05-17/18, experiment IF_ABW. The bench plan (amos imaging/plans/notes/IF_ABW.txt) lists
  "Cre + A, Cre + B, Cre + E" and "Cre - C, Cre - F, Cre - G"; the old analysis notebook maps the same letters
  (A, B, E -> Ezh2_cKO; C, D, F, G -> WT). The section letter is the animal: s4_A and s5_A are one mouse.
  20x panoramas, 2018-05-22. The old notebook and segment_ml_nuclei.R map the section's position on the slide
  (the digit after s1_ / s2_ / s3_) to the animal: c(cKO, cKO, WT, cKO, cKO), so positions 2, 4 and 5 (all the
  20x sections here) are cKO. Together that is the paper's n = 7 cKO / 4 WT mice.
  The labelling project (igl_gapfill_apply.info, the 2026-11-18 note) instead treats every 10x section and the
  gold s1_2 as controls. scheme="project" reproduces that for comparison.

Pixel sizes: 10x 1.0239 um/px, 20x 0.5119 (PX_TRUE below; the 10x labels were built at 0.65 and are rescaled).

Per section, from outputs/flatten_v11/<s>.npz (columns at right angles to the surface, 2 um apart, sampled every
2 um in depth, in the units the labels were built with: 3.15 true um on the 10x; S_nuclei / S_neun / S_p27 = the native stains block-averaged to ~2 um and sampled on the columns)
and outputs/flat_fix/<s>.npz (L_new = the final labels on the columns, ok_new = EGL-covered cortex, geo =
geometry errors):
  * a column counts when it is EGL-covered cortex, not a geometry error, its labels run EGL -> ML -> IGL from
    the surface (leading unassigned / anomaly <= 40 um skipped, gaps <= 4 um in the EGL and <= 8 um in the IGL
    closed) and the IGL ends inside the column;
  * thickness of EGL, ML, IGL per column; the section value is the median over its counted columns;
  * mean stain per layer per column (EGL split into the outer and the inner half of its depth: the inner half
    is where the paper's iEGL lies), the section value the median over columns; raw intensities (comparable only
    within one imaging set) and as a ratio to the section's own ML;
  * the p27+ inner EGL: EGL samples whose p27 (smoothed 2 um in depth) is past half way from the section's ML p27
    level to its IGL level (0.4 and 0.6 as sensitivity), its thickness per column and share of the EGL; the rest is
    the p27- outer EGL (the paper's iEGL / oEGL);
  * depth profiles: the median over counted columns of each stain at each depth below the surface, and in layer
    coordinates (each column's EGL, ML and IGL resampled to 20 bins each);
  * areas of EGL, ML, IGL and white matter on the final 4 um label grid (outputs/flat_fix/final).

    python scripts/quant_layers.py            # -> outputs/quant/sections.csv, columns.parquet-like npz, profiles.npz
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi
from scipy.stats import trim_mean

sys.path.insert(0, str(Path(__file__).parent))
import flatten as F  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/quant"
FIX = ROOT / "outputs/flat_fix"
STEP, G = F.STEP, F.G
STAINS = {"p27": "S_p27", "NeuN": "S_neun", "DAPI": "S_nuclei"}
LAYERS = ("oEGL", "iEGL", "EGL", "ML", "IGL")
NBIN = 20                                   # bins per layer in layer coordinates
PIEGL = (0.5, 0.4, 0.6)                     # p27+ inner EGL: past this fraction of the way from ML to IGL p27 (first = main)
DMAX_UM = 250.0
# pixel sizes. The labels, the flattened strips and the 4 um label grid of the 10x sections were built with 0.65 um/px
# (2026-11-18 to 11-24); the microscope's is 1.0239 (AxioCam MR3 6.45 um pixels / (10x x 0.63 camera adapter), as the
# 20x is 6.45 / (20 x 0.63) = 0.5119; outputs/px_scale/summary.txt). Every length is multiplied by PX_TRUE / PX_BUILT
# and every um threshold is converted with the true size, so the 10x labels need not be rebuilt.
PX_BUILT = {"10x": 0.65, "20x": 0.5119}
PX_TRUE = {"10x": 1.0239, "20x": 0.5119}


def um_scale(short):
    """true um per built um for the section's imaging set (10x 1.575, 20x 1)."""
    m = genotype(short)[2]
    return PX_TRUE[m] / PX_BUILT[m]

GENOTYPE = dict(
    sources={
        "10x": "IF_ABW bench plan: Cre+ A, B, E; Cre- C, F, G (old notebook: A, B, E -> Ezh2_cKO; C, D, F, G -> WT)",
        "20x": "old notebook / segment_ml_nuclei.R: slide position 1-5 -> cKO, cKO, WT, cKO, cKO (positions 2, 4, 5 here)",
        "project": "labelling project note 2026-11-18: all 10x and the gold s1_2 are controls, other 20x cKO",
    })


def genotype(short, scheme="bench"):
    """(animal, genotype 'WT' / 'cKO', imaging set '10x' / '20x')."""
    slide, pos = short.split("_")
    if slide in ("s1", "s2", "s3"):
        g = {1: "cKO", 2: "cKO", 3: "WT", 4: "cKO", 5: "cKO"}[int(pos)]
        if scheme == "project":
            g = "WT" if short == "s1_2" else "cKO"
        return f"2018-{pos}", g, "20x"
    g = "cKO" if pos in ("A", "B", "E") else "WT"
    if scheme == "project":
        g = "WT"
    return f"2017-{pos}", g, "10x"


def _run_from(m, a, gap):
    """end of the run of m starting at a, bridging gaps <= gap samples (m[a] must be True)."""
    n, d, last = len(m), a, a
    while d < n:
        if m[d]:
            last = d
        elif d - last > gap:
            break
        d += 1
    return last + 1


def column_layers(L, valid, step=STEP):
    """per column sample indices: lead, EGL bottom, IGL top, IGL bottom, column end; complete flag.
    step = true um per sample: lead <= 40 um skipped, gaps <= 4 um in the EGL and <= 8 um in the IGL bridged,
    IGL runs >= 8 um."""
    N, Dn = L.shape
    n_lead, g_egl, g_igl = int(40 / step), max(1, round(4 / step)), max(1, round(8 / step))
    lead = np.full(N, -1)
    eb, it, ib = (np.full(N, -1) for _ in range(3))
    end = valid.sum(1)
    for i in range(N):
        e = end[i]
        if e < 10:
            continue
        col = L[i, :e]
        a = 0
        while a < min(e, n_lead) and col[a] in (0, 1, 6):
            a += 1
        if a >= e or col[a] != 2:
            continue
        lead[i] = a
        eb[i] = _run_from(col == 2, a, g_egl)
        igl = np.flatnonzero(col[eb[i]:] == 4)
        if not len(igl):
            continue
        t = eb[i] + igl[0]
        b = _run_from(col == 4, t, g_igl)
        if b - t >= g_igl:
            it[i], ib[i] = t, b
    complete = (lead >= 0) & (it > eb) & (ib > it) & (ib < end - 2)
    return dict(lead=lead, eb=eb, it=it, ib=ib, end=end, complete=complete)


def offset(short, Z):
    """camera offset to subtract: 75 counts on the 10x 12-bit panoramas (measured outside the tissue of the three
    unmasked ones, s5_A / B / C: 70, 86, 72), x16 on the two stored stretched to 16 bit (s4_A, s4_C); unknown on
    the masked 20x panoramas (0)."""
    if genotype(short)[2] == "20x":
        return 0.0
    return 75.0 * (16 if max(int(Z[v].max()) for v in STAINS.values()) > 20000 else 1)


def section(short):
    Z = np.load(F.OUT2 / f"{short}.npz")
    R = np.load(FIX / f"{short}.npz")
    L, valid = R["L_new"], Z["valid"]
    use = R["ok_new"] & ~R["geo"]
    ks = um_scale(short)
    step = STEP * ks                                  # true um per strip sample (and between columns)
    C = column_layers(L, valid, step)
    cnt = use & C["complete"]
    N, Dn = L.shape
    off = offset(short, Z)
    S = {k: Z[v].astype(np.float32) - off for k, v in STAINS.items()}
    k = np.arange(Dn)[None, :]
    lead, eb, it, ib = (C[x][:, None] for x in ("lead", "eb", "it", "ib"))
    mid = (lead + eb) // 2
    bands = {"oEGL": (k >= lead) & (k < mid), "iEGL": (k >= mid) & (k < eb), "EGL": (k >= lead) & (k < eb),
             "ML": (k >= eb) & (k < it), "IGL": (k >= it) & (k < ib)}
    th = {"EGL": (C["eb"] - C["lead"]) * step, "ML": (C["it"] - C["eb"]) * step, "IGL": (C["ib"] - C["it"]) * step}
    animal, g, mag = genotype(short)
    row = dict(section=short, animal=animal, genotype=g, imaging=mag, genotype_project=genotype(short, "project")[1], offset=off,
               px_um=PX_TRUE[mag], um_scale=round(ks, 4),
               columns=int(N), line_mm=round(N * step / 1000, 2), covered_mm=round(use.sum() * step / 1000, 2),
               counted_mm=round(cnt.sum() * step / 1000, 2), counted_frac=round(float(cnt.sum() / max(use.sum(), 1)), 3))
    for ly, t in th.items():
        row[f"{ly}_thick_um"] = round(float(np.median(t[cnt])), 1)
        row[f"{ly}_thick_iqr_um"] = round(float(np.subtract(*np.percentile(t[cnt], [75, 25]))), 1)
    colmean = {}
    for st, A in S.items():
        for ly, m in bands.items():
            mm = m & valid & cnt[:, None]
            v = np.where(mm, A, 0).sum(1) / np.maximum(mm.sum(1), 1)
            colmean[(st, ly)] = np.where(mm.any(1), v, np.nan)
            row[f"{st}_{ly}"] = round(float(np.nanmedian(colmean[(st, ly)])), 1)
        for ly in LAYERS:
            row[f"{st}_{ly}_rel_ML"] = round(row[f"{st}_{ly}"] / max(row[f"{st}_ML"], 1e-6), 3)
            row[f"{st}_{ly}_ml_igl"] = round((row[f"{st}_{ly}"] - row[f"{st}_ML"]) / max(row[f"{st}_IGL"] - row[f"{st}_ML"], 1e-6), 3)
    # p27+ inner EGL: EGL samples whose p27 is past PIEGL of the way from the section's ML to its IGL level
    # (smoothed 2 um in depth); the outer, p27- EGL is the rest
    p27n = (ndi.gaussian_filter1d(S["p27"], 2.0 / step, axis=1) - row["p27_ML"]) / max(row["p27_IGL"] - row["p27_ML"], 1e-6)
    egl = bands["EGL"] & valid & cnt[:, None]
    piegl = {}
    for thr in PIEGL:
        pos = egl & (p27n >= thr)
        piegl[thr] = pos.sum(1) * step
        tag = "" if thr == PIEGL[0] else f"_{int(thr * 100)}"
        row[f"iEGLp27{tag}_thick_um"] = round(float(np.median(piegl[thr][cnt])), 1)
        row[f"oEGLp27{tag}_thick_um"] = round(float(np.median((th["EGL"] - piegl[thr])[cnt])), 1)
        row[f"iEGLp27{tag}_frac"] = round(float(np.median((piegl[thr] / np.maximum(th["EGL"], step))[cnt])), 3)
    # depth profiles from the surface (lead-aligned), on the section's samples, then on the common 2 um axis
    nd = int(DMAX_UM / step)
    depth_c = np.arange(int(DMAX_UM / STEP)) * STEP
    prof = {}
    rows_c = np.flatnonzero(cnt)
    for st, A in S.items():
        P = np.full((len(rows_c), nd), np.nan, np.float32)
        for j, i in enumerate(rows_c):
            a = C["lead"][i]
            seg = A[i, a:min(a + nd, C["end"][i])]
            P[j, :len(seg)] = seg
        enough = np.isfinite(P).sum(0) >= max(20, 0.2 * len(rows_c))
        pm = np.where(enough, np.nanmedian(P, 0), np.nan)
        dk = np.arange(nd) * step
        ok_ = np.isfinite(pm)
        prof[f"depth_{st}"] = np.where(depth_c <= (dk[ok_].max() if ok_.any() else -1), np.interp(depth_c, dk[ok_], pm[ok_]) if ok_.any() else np.nan, np.nan)
        # layer coordinates
        Q = np.full((len(rows_c), 3 * NBIN), np.nan, np.float32)
        for j, i in enumerate(rows_c):
            for q, (a, b) in enumerate(((C["lead"][i], C["eb"][i]), (C["eb"][i], C["it"][i]), (C["it"][i], C["ib"][i]))):
                if b - a < 2:
                    continue
                x = np.arange(a, b) + 0.5
                u = (np.arange(NBIN) + 0.5) / NBIN * (b - a) + a
                Q[j, q * NBIN:(q + 1) * NBIN] = np.interp(u, x, A[i, a:b])
        prof[f"layer_{st}"] = np.nanmedian(Q, 0)
    # areas on the final grid
    lab = tifffile.imread(FIX / "final" / f"{short}_labels.tif")
    for code, ly in ((2, "EGL"), (3, "ML"), (4, "IGL"), (5, "WM")):
        row[f"{ly}_area_mm2"] = round(float((lab == code).sum() * (G * ks) ** 2 / 1e6), 4)
    tot = sum(row[f"{ly}_area_mm2"] for ly in ("EGL", "ML", "IGL", "WM"))
    row["layers_area_mm2"] = round(tot, 4)
    for ly in ("EGL", "ML", "IGL", "WM"):
        row[f"{ly}_area_frac"] = round(row[f"{ly}_area_mm2"] / max(tot, 1e-9), 4)
        row[f"{ly}_area_per_mm"] = round(row[f"{ly}_area_mm2"] / max(row["covered_mm"], 1e-9) * 1000, 1)  # um2 per um = um
    # stain-only boundaries (a check on the labels): EGL bottom where DAPI drops past half way from the section's
    # EGL level to its ML level, IGL top where NeuN rises past half way from ML to IGL, each held >= 4 um, both
    # smoothed 2 um in depth; searched from the column start / the stain EGL bottom
    dn = (ndi.gaussian_filter1d(S["DAPI"], 2.0 / step, axis=1) - row["DAPI_ML"]) / max(row["DAPI_EGL"] - row["DAPI_ML"], 1e-6)
    nn = (ndi.gaussian_filter1d(S["NeuN"], 2.0 / step, axis=1) - row["NeuN_ML"]) / max(row["NeuN_IGL"] - row["NeuN_ML"], 1e-6)
    hold = max(1, round(4 / step))
    low = ndi.minimum_filter1d((dn < 0.5).astype(np.uint8), hold, axis=1, origin=-(hold // 2)) > 0
    high = ndi.minimum_filter1d((nn > 0.5).astype(np.uint8), hold, axis=1, origin=-(hold // 2)) > 0
    ebs = np.full(N, -1)
    its = np.full(N, -1)
    for i in np.flatnonzero(cnt):
        a, e = C["lead"][i], C["end"][i]
        w = np.flatnonzero(low[i, a:e])
        if len(w):
            ebs[i] = a + w[0]
            v = np.flatnonzero(high[i, ebs[i]:e])
            if len(v):
                its[i] = ebs[i] + v[0]
    okst = cnt & (ebs >= 0) & (its > ebs)
    row["stain_counted_frac"] = round(float(okst.sum() / max(cnt.sum(), 1)), 3)
    row["EGL_stain_thick_tm_um"] = round(float(trim_mean(((ebs - C["lead"]) * step)[okst], 0.1)), 2)
    row["ML_stain_thick_tm_um"] = round(float(trim_mean(((its - ebs) * step)[okst], 0.1)), 2)
    row["IGLtop_depth_tm_um"] = round(float(trim_mean(((C["it"] - C["lead"]) * step)[cnt], 0.1)), 2)
    row["IGLtop_stain_depth_tm_um"] = round(float(trim_mean(((its - C["lead"]) * step)[okst], 0.1)), 2)
    row["EGL_stain_label_r"] = round(float(np.corrcoef(((ebs - C["lead"]) * step)[okst], th["EGL"][okst])[0, 1]), 3)
    row["ML_stain_label_r"] = round(float(np.corrcoef(((its - ebs) * step)[okst], th["ML"][okst])[0, 1]), 3)
    # 10 % trimmed means beside the medians (column thickness comes in 2 um steps, so medians move in 2 um steps)
    tm = dict(EGL=th["EGL"], ML=th["ML"], IGL=th["IGL"], iEGLp27=piegl[PIEGL[0]], oEGLp27=th["EGL"] - piegl[PIEGL[0]])
    for ly, v in tm.items():
        row[f"{ly}_thick_tm_um"] = round(float(trim_mean(v[cnt], 0.1)), 2)
    cols = dict(thick_EGL=th["EGL"], thick_ML=th["ML"], thick_IGL=th["IGL"], counted=cnt, thick_iEGLp27=piegl[PIEGL[0]],
                **{f"{st}_{ly}": v for (st, ly), v in colmean.items()})
    return row, prof, cols


def main(shorts=None):
    OUT.mkdir(parents=True, exist_ok=True)
    rows, profs = [], {}
    for short in shorts or F.all_shorts():
        row, prof, cols = section(short)
        rows.append(row)
        profs.update({f"{short}__{k}": v for k, v in prof.items()})
        np.savez_compressed(OUT / f"columns_{short}.npz", **cols)
        print(json.dumps({k: row[k] for k in ("section", "animal", "genotype", "imaging", "counted_mm", "counted_frac",
                                               "EGL_thick_um", "ML_thick_um", "IGL_thick_um", "p27_iEGL", "p27_IGL")}), flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "sections.csv", index=False)
    np.savez_compressed(OUT / "profiles.npz", **profs, depth_um=np.arange(int(DMAX_UM / STEP)) * STEP, nbin=NBIN)
    (OUT / "genotype.json").write_text(json.dumps(dict(GENOTYPE, table={r["section"]: [r["animal"], r["genotype"], r["imaging"], r["genotype_project"]] for r in rows}), indent=1))
    return df


if __name__ == "__main__":
    main(sys.argv[1:] or None)
