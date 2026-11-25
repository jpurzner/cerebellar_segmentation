"""
Figures for the molecular-layer nuclei (2026-11-24), from scripts/ml_nuclei.py outputs.

    python scripts/ml_nuclei_figures.py      # -> figures/quant/Q4*.png|pdf|svg, outputs/quant/ml_stats.csv

Unit of replication = mouse (mean of its sections; each section = mean over its kept NeuN+ ML nuclei). 10x (2017
litter, 1.0239 um/px) and 20x (2018 litter) pooled; test: cKO - WT on mouse means with litter as a factor
(figstyle.tests). Q4s shows what the 10x pixel size alone does to the 20x numbers. The 2018 wild type exists only in
the old repo's folium crops, so the crops give a second pooled test (2017 sections + 2018 crops, 4 WT / 7 cKO mice),
written to ml_stats.csv as pooledcrops_*.
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=RuntimeWarning)
sys.path.insert(0, str(Path(__file__).parent))
import figstyle as S  # noqa: E402
import quant_layers as QL  # noqa: E402

plt = S.plt
OUT = QL.OUT
ELONG = 0.85                # "elongated": the nuclei whose axis means something (round ones come out uniformly oriented)
SECTIONS = [f.name.replace("ml_nuclei_", "").replace(".csv.gz", "") for f in sorted(OUT.glob("ml_nuclei_s*.csv.gz")) if "_down" not in f.name]


def nuclei(down=False):
    tag = "_down" if down else ""
    out = []
    for s in SECTIONS:
        p = OUT / f"ml_nuclei_{s}{tag}.csv.gz"
        if p.exists():
            out.append(pd.read_csv(p))
    return pd.concat(out, ignore_index=True)


def per_section(N, subset=None):
    rows = []
    for s, d in N.groupby("section"):
        pos = d[d.neun_positive]
        if subset is not None:
            pos = pos[subset(pos)]
        rows.append(dict(section=s, animal=d.animal.iloc[0], genotype=d.genotype.iloc[0], imaging=d.imaging.iloc[0],
                         n=len(d), n_pos=len(pos), neun_frac=d.neun_positive.mean(), ecc=pos.eccentricity.mean(),
                         angle=pos.angle_igl_deg.mean(), angle_egl=pos.angle_egl_deg.mean(), radial30=(pos.angle_igl_deg < 30).mean(),
                         area=pos.area_um2.median(), major=pos.major_um.median(), minor=pos.minor_um.median(),
                         aspect=(pos.major_um / pos.minor_um).median(),
                         elong_frac=(pos.eccentricity >= ELONG).mean(),
                         angle_elong=pos[pos.eccentricity >= ELONG].angle_igl_deg.mean(),
                         radial_elong=(pos[pos.eccentricity >= ELONG].angle_igl_deg < 30).mean()))
    return pd.DataFrame(rows)


def _mouse_curves(N, g, m, key, bins, val):
    """per mouse, the mean of val in bins of key (over NeuN+ nuclei of all its sections)."""
    out = []
    d = N[(N.genotype == g) & (N.imaging == m) & N.neun_positive]
    for a, x in d.groupby("animal"):
        idx = np.digitize(x[key], bins) - 1
        out.append([x[val][idx == k].mean() if (idx == k).sum() >= 10 else np.nan for k in range(len(bins) - 1)])
    return np.array(out)


def _mouse_hist(N, g, m, bins):
    out = []
    d = N[(N.genotype == g) & (N.imaging == m) & N.neun_positive]
    for a, x in d.groupby("animal"):
        h, _ = np.histogram(x.angle_igl_deg, bins=bins)
        out.append(h / h.sum() * (len(bins) - 1))            # 1 = uniform (no preferred direction)
    return np.array(out)


def _curve_panel(ax, xs, curves, ylab, title, ref=None, xlab=None):
    for (g, m), A in curves.items():
        if A is None or not len(A):
            continue
        mu = np.nanmean(A, 0)
        se = np.nanstd(A, 0, ddof=1) / np.sqrt(np.isfinite(A).sum(0))
        if m == "10x":
            ax.fill_between(xs, mu - se, mu + se, color=S.COL[g], alpha=0.18, lw=0)
            ax.plot(xs, mu, color=S.COL[g], lw=1.4, marker="o", ms=2.8, label=f"{g} 10x")
        else:
            ax.plot(xs, mu, color=S.COL[g], lw=1.0, ls=(0, (3, 1.6)), marker="^", ms=2.8, mfc="white", label=f"{g} 20x")
    if ref is not None:
        ax.axhline(ref, color=S.MUTED, lw=0.6, ls=(0, (2, 2)), zorder=0)
    ax.set_ylabel(ylab)
    if xlab:
        ax.set_xlabel(xlab)
    ax.set_title(title, loc="left", fontweight="bold")


def fig_main(N, P, rows):
    fig, axs = plt.subplots(2, 4, figsize=(7.1, 4.9))
    for ax, (col, yl, t, ref) in zip(axs[0], [("ecc", "eccentricity (NeuN+ ML nuclei)", "A  Eccentricity", None),
                                              ("elong_frac", f"share with eccentricity \u2265 {ELONG}", "B  Elongated nuclei", None),
                                              ("angle_elong", "angle to IGL perpendicular (\u00b0)", "C  Their orientation", 45),
                                              ("neun_frac", "NeuN+ share of ML nuclei", "D  NeuN+ nuclei", None)]):
        r = S.group_dots(ax, P, col, ylabel=yl, title=t)
        rows.append(dict(measure=col, **r))
        if ref is not None:
            y0, y1 = ax.get_ylim()
            if ref > y1 - 0.06 * (y1 - y0):
                ax.set_ylim(y0, ref + 0.1 * (ref - y0))
            ax.axhline(ref, color=S.MUTED, lw=0.6, ls=(0, (2, 2)), zorder=0)
            ax.text(1.53, ref, "random", fontsize=5.5, color=S.MUTED, ha="right", va="bottom")
        if col in ("ecc", "elong_frac", "neun_frac"):
            ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.2f}"))
    E = N[N.eccentricity >= ELONG]
    bins = np.arange(0, 91, 10)
    xs = bins[:-1] + 5
    _curve_panel(axs[1, 0], xs, {k: _mouse_hist(E, *k, bins) for k in S.GROUPS}, "density (1 = no preferred angle)",
                 "E  Angles, elongated", ref=1, xlab="angle to IGL perpendicular (\u00b0)")
    axs[1, 0].set_xticks([0, 30, 60, 90])
    db = np.linspace(0, 1, 6)
    _curve_panel(axs[1, 1], (db[:-1] + db[1:]) / 2, {k: _mouse_curves(E, *k, "ml_depth_frac", db, "angle_igl_deg") for k in S.GROUPS},
                 "angle to IGL perpendicular (\u00b0)", "F  Across the ML", ref=45, xlab="depth in ML (0 EGL side, 1 IGL)")
    _curve_panel(axs[1, 2], (db[:-1] + db[1:]) / 2, {k: _mouse_curves(N, *k, "ml_depth_frac", db, "eccentricity") for k in S.GROUPS},
                 "eccentricity (all NeuN+)", "G  Shape across the ML", xlab="depth in ML (0 EGL side, 1 IGL)")
    eb = np.array([0.0, 0.6, 0.7, 0.8, 0.9, 1.0])
    _curve_panel(axs[1, 3], (eb[:-1] + eb[1:]) / 2, {k: _mouse_curves(N, *k, "eccentricity", eb, "angle_igl_deg") for k in S.GROUPS},
                 "angle to IGL perpendicular (\u00b0)", "H  Angle by shape", ref=45, xlab="eccentricity")
    axs[1, 0].legend(loc="upper right", handlelength=1.6)
    S.legend_groups(fig)
    fig.text(0.5, 0.004, "NeuN+ DAPI nuclei inside the final ML label (\u2265 3 \u00b5m from EGL and IGL). Angle: major axis vs the line from the nucleus "
             f"perpendicular to the IGL\n(0\u00b0 radial, 90\u00b0 tangential); elongated = eccentricity \u2265 {ELONG}. A-D: dots = sections, ringed = mice, "
             "bar = mean \u00b1 SEM over mice; p: cKO vs WT on mouse\nmeans, litter as a factor. E-H: mean \u00b1 SEM over mice.",
             ha="center", va="bottom", fontsize=5.8, color=S.INK2)
    fig.subplots_adjust(left=0.075, right=0.985, top=0.885, bottom=0.17, wspace=0.62, hspace=0.68)
    return fig


def fig_supp(N, Nd, P, Pd, rows):
    """resolution test (20x native vs resampled to 10x pixels), geometry subset, IGL vs EGL perpendicular."""
    fig, axs = plt.subplots(2, 4, figsize=(7.1, 4.6))
    m = P.set_index("section").join(Pd.set_index("section"), rsuffix="_down", how="inner")
    for ax, (col, lab, t) in zip(axs[0], [("ecc", "eccentricity", "A  Eccentricity"), ("angle", "angle to IGL perp. (°)", "B  Angle"),
                                          ("area", "median area (µm²)", "C  Area"), ("neun_frac", "NeuN+ fraction", "D  NeuN+ fraction")]):
        x, y = m[col], m[f"{col}_down"]
        lo, hi = min(x.min(), y.min()), max(x.max(), y.max())
        pad = 0.08 * (hi - lo)
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color=S.AXIS, lw=0.6, zorder=0)
        ax.scatter(x, y, s=16, marker="^", facecolor="white", edgecolor=S.CKO, linewidth=0.9, zorder=3)
        ax.set_xlim(lo - pad, hi + pad)
        ax.set_ylim(lo - pad, hi + pad)
        ax.set_xlabel(f"20x native ({lab})")
        ax.set_ylabel("resampled to 1.0239 µm/px")
        ax.set_title(t, loc="left", fontweight="bold")
        ax.grid(axis="x")
        d = (y - x)
        ax.text(0.96, 0.04, f"mean shift {d.mean():+.3g}", transform=ax.transAxes, ha="right", va="bottom", fontsize=5.8, color=S.INK2)
        rows.append(dict(measure=f"resolution_shift_{col}", diff=float(d.mean()), lo=float(d.min()), hi=float(d.max()), p=np.nan, na=len(d), nb=len(d)))
    # geometry subset: nuclei where the IGL and EGL perpendiculars agree within 30 deg
    Pc = per_section(N, subset=lambda p: p.normal_dev_deg < 30)
    r = S.group_dots(axs[1, 0], Pc, "angle_elong", ylabel="angle, elongated (°)", title="E  Straight ML only")
    rows.append(dict(measure="angle_clean_geometry", **r))
    axs[1, 0].axhline(45, color=S.MUTED, lw=0.6, ls=(0, (2, 2)), zorder=0)
    r = S.group_dots(axs[1, 1], P, "angle_egl", ylabel="angle to EGL perp., all (°)", title="F  EGL reference")
    rows.append(dict(measure="angle_egl", **r))
    axs[1, 1].axhline(45, color=S.MUTED, lw=0.6, ls=(0, (2, 2)), zorder=0)
    r = S.group_dots(axs[1, 2], P, "angle", ylabel="angle, all NeuN+ (°)", title="G  All nuclei")
    rows.append(dict(measure="angle", **r))
    axs[1, 2].axhline(45, color=S.MUTED, lw=0.6, ls=(0, (2, 2)), zorder=0)
    P2 = P.assign(dens=P.n)
    r = S.group_dots(axs[1, 3], P, "minor", ylabel="median minor axis (µm)", title="H  Nucleus width")
    rows.append(dict(measure="minor_um", **r))
    S.legend_groups(fig)
    fig.text(0.5, 0.004, "A-D: the 20x cKO sections segmented again after resampling to the 10x pixel size (what pixel size alone does; the 10x optics are not simulated).\n"
             "E: elongated nuclei where the perpendiculars to the IGL and to the EGL agree within 30°. F, G: all NeuN+ nuclei. p: cKO vs WT, litter as a factor.",
             ha="center", va="bottom", fontsize=5.8, color=S.INK2)
    fig.subplots_adjust(left=0.075, right=0.985, top=0.885, bottom=0.15, wspace=0.62, hspace=0.72)
    return fig


def fig_crops(rows):
    """the old repo's 20x folium crops, layers from the stains: the 20x wild type (slide position 3) against the four
    20x cKO mice, and the crops against the whole labelled sections."""
    C = pd.read_csv(OUT / "ml_nuclei_crops.csv.gz")
    Pc = per_section(C.assign(imaging="20x"))
    P = pd.read_csv(OUT / "ml_nuclei_per_section.csv")
    fig, axs = plt.subplots(1, 5, figsize=(7.1, 2.75), gridspec_kw=dict(width_ratios=[1, 1, 1, 1, 1.25]))
    rng = np.random.default_rng(1)
    for ax, (col, yl, t) in zip(axs[:4], [("ecc", "eccentricity", "A  Eccentricity"), ("elong_frac", f"share ecc \u2265 {ELONG}", "B  Elongated"),
                                          ("angle_elong", "angle, elongated (\u00b0)", "C  Orientation"), ("neun_frac", "NeuN+ share", "D  NeuN+")]):
        am = Pc.groupby(["animal", "genotype"], as_index=False)[col].mean()
        for x, g in enumerate(("WT", "cKO")):
            sv = Pc[Pc.genotype == g][col].to_numpy(float)
            av = am[am.genotype == g][col].to_numpy(float)
            ax.scatter(x + rng.uniform(-0.15, 0.15, len(sv)), sv, s=9, marker="^", facecolor="none", edgecolor=S.COL[g], lw=0.6, alpha=0.6, zorder=3)
            ax.scatter(x + rng.uniform(-0.05, 0.05, len(av)), av, s=30, marker="^", facecolor="white", edgecolor=S.COL[g], lw=1.0, zorder=4)
            mu = av.mean()
            ax.plot([x - 0.28, x + 0.28], [mu, mu], color=S.INK, lw=1.0, zorder=5)
        ax.set_xticks([0, 1], ["WT\n1 mouse", "cKO\n4 mice"])
        ax.set_xlim(-0.6, 1.6)
        ax.tick_params(axis="x", length=0)
        ax.set_ylabel(yl)
        ax.set_title(t, loc="left", fontweight="bold")
        if col in ("ecc", "elong_frac", "neun_frac"):
            ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.2f}"))
        w = am[am.genotype == "WT"][col].iloc[0]
        k = am[am.genotype == "cKO"][col]
        rows.append(dict(measure=f"crops_{col}", diff=float(k.mean() - w), lo=float(k.min() - w), hi=float(k.max() - w), p=np.nan, na=1, nb=len(k)))
    ax = axs[4]
    m = Pc.set_index("section").join(P.set_index("section"), rsuffix="_full", how="inner")
    lo = min(m.angle_elong.min(), m.angle_elong_full.min()) - 2
    hi = max(m.angle_elong.max(), m.angle_elong_full.max()) + 2
    ax.plot([lo, hi], [lo, hi], color=S.AXIS, lw=0.6, zorder=0)
    ax.scatter(m.angle_elong_full, m.angle_elong, s=16, marker="^", facecolor="white", edgecolor=S.CKO, lw=0.9, zorder=3)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel("whole section, labels (\u00b0)")
    ax.set_ylabel("folium crops, stain layers (\u00b0)")
    ax.set_title("E  Crops vs sections", loc="left", fontweight="bold")
    ax.grid(axis="x")
    fig.text(0.5, 0.004, "The old repo's 43 hand-cut 20x folium crops (2018 litter), layers classified from DAPI and NeuN. The wild type is slide position 3 "
             "(s1_3, s2_3, s3_3), unlabelled in this project.\nSmall triangles = sections, ringed = mice, bar = mean. One wild-type mouse: no test. "
             "E: elongated-nucleus angle, the same 8 cKO sections both ways (crops are straight stretches).", ha="center", va="bottom", fontsize=5.8, color=S.INK2)
    fig.subplots_adjust(left=0.07, right=0.99, top=0.88, bottom=0.27, wspace=0.75)
    return fig


def main():
    N = nuclei()
    P = per_section(N)
    rows = []
    S.save(fig_main(N, P, rows), "Q4_ml_nuclei")
    Nd = nuclei(down=True)
    if len(Nd):
        Pd = per_section(Nd)
        S.save(fig_supp(N, Nd, P, Pd, rows), "Q4s_ml_nuclei_checks")
    if (OUT / "ml_nuclei_crops.csv.gz").exists():
        S.save(fig_crops(rows), "Q4c_crops_20x")
    if (OUT / "ml_nuclei_crops.csv.gz").exists():        # 2017 whole sections + 2018 crops (the 2018 litter's WT is only in the crops)
        Pc = per_section(pd.read_csv(OUT / "ml_nuclei_crops.csv.gz").assign(imaging="20x"))
        both = pd.concat([P[P.imaging == "10x"], Pc], ignore_index=True)
        for col in ("ecc", "elong_frac", "angle_elong", "neun_frac", "angle", "radial_elong"):
            rows.append(dict(measure=f"pooledcrops_{col}", **S.tests(S.animal_means(both, col), col)))
    st = pd.DataFrame(rows)
    st.to_csv(OUT / "ml_stats.csv", index=False)
    P.to_csv(OUT / "ml_nuclei_per_section.csv", index=False)
    pd.set_option("display.width", 200)
    print(st.round(4).to_string(index=False))
    print(P.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
