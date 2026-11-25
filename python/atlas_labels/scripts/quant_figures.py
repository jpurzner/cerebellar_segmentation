"""
Figures for layer size and stain intensity, WT vs Ezh2 cKO (2026-11-24), from scripts/quant_layers.py outputs.

    python scripts/quant_figures.py        # -> figures/quant/*.png|pdf|svg, outputs/quant/stats.csv

10x (2017 litter, 1.0239 um/px) and 20x (2018 litter, 0.5119 um/px) are pooled. Unit = mouse (mean of its
sections). Test: cKO - WT by least squares on mouse means with litter as a factor (figstyle.tests); the plain pooled
Welch test and the 2017-litter-only Welch test go to stats.csv beside it. The labelled 20x sections are slide
positions 2, 4, 5 (all cKO), so until the 20x wild type (position 3) is labelled the genotype difference comes from
the 2017 litter and the 2018 mice sharpen the error estimate.
"""
import sys
from pathlib import Path

import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=RuntimeWarning)
sys.path.insert(0, str(Path(__file__).parent))
import figstyle as S  # noqa: E402
import quant_layers as QL  # noqa: E402

plt = S.plt
OUT = QL.OUT


def load():
    df = pd.read_csv(OUT / "sections.csv")
    P = np.load(OUT / "profiles.npz")
    return df, P


def fig_size(df, stats_rows):
    fig, axs = plt.subplots(2, 4, figsize=(7.1, 4.6))
    panels = [("EGL_thick_um", "EGL thickness (µm)", "A  EGL"), ("ML_thick_um", "ML thickness (µm)", "B  Molecular layer"),
              ("IGL_thick_um", "IGL thickness (µm)", "C  IGL"), ("layers_area_mm2", "labelled area (mm²)", "D  Section size"),
              ("EGL_area_frac", "fraction of labelled area", "E  EGL area"), ("ML_area_frac", "fraction of labelled area", "F  ML area"),
              ("IGL_area_frac", "fraction of labelled area", "G  IGL area"), ("WM_area_frac", "fraction of labelled area", "H  White matter area")]
    for ax, (col, yl, t) in zip(axs.ravel(), panels):
        r = S.group_dots(ax, df, col, ylabel=yl, title=t)
        stats_rows.append(dict(measure=col, **r))
        if col.endswith("frac"):
            ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.2f}"))
    S.legend_groups(fig)
    fig.text(0.5, 0.005, "Thickness: median over the section's complete cortical columns (EGL at the surface, IGL ending inside the column). "
             "Areas: final labels. p: cKO vs WT on mouse means, litter as a factor.", ha="center", va="bottom", fontsize=5.8, color=S.INK2)
    fig.subplots_adjust(left=0.075, right=0.985, top=0.885, bottom=0.12, wspace=0.55, hspace=0.75)
    return fig


def _norm_prof(df, P, kind, st):
    """per section profile scaled 0 = its ML, 1 = its IGL (offset and gain cancel)."""
    out = {}
    for _, r in df.iterrows():
        p = P[f"{r.section}__{kind}_{st}"]
        out[r.section] = (p - r[f"{st}_ML"]) / max(r[f"{st}_IGL"] - r[f"{st}_ML"], 1e-6)
    return out


def _group_prof(df, prof, g, m):
    am = []
    for a, d in df[(df.genotype == g) & (df.imaging == m)].groupby("animal"):
        am.append(np.nanmean(np.stack([prof[s] for s in d.section]), 0))
    return np.stack(am) if am else None


def fig_profiles(df, P, stats_rows, kind="layer"):
    stains = ["p27", "NeuN", "DAPI"]
    fig, axs = plt.subplots(2, 3, figsize=(7.1, 4.4), gridspec_kw=dict(height_ratios=[2.2, 1]), sharex=True)
    nb = int(P["nbin"])
    if kind == "layer":
        x = (np.arange(3 * nb) + 0.5) / nb
        xl = "position in the layer (EGL | ML | IGL, each scaled to its own thickness)"
    else:
        x = P["depth_um"]
        xl = "depth below the surface (µm)"
    for j, st in enumerate(stains):
        prof = _norm_prof(df, P, kind, st)
        ax, axd = axs[0, j], axs[1, j]
        curves = {}
        for (g, m), ls, lw in ((("WT", "10x"), "-", 1.4), (("cKO", "10x"), "-", 1.4), (("cKO", "20x"), (0, (3, 1.6)), 1.0)):
            A = _group_prof(df, prof, g, m)
            mu = np.nanmean(A, 0)
            se = np.nanstd(A, 0, ddof=1) / np.sqrt(np.isfinite(A).sum(0))
            curves[(g, m)] = A
            if m == "10x":
                ax.fill_between(x, mu - se, mu + se, color=S.COL[g], alpha=0.18, lw=0)
            ax.plot(x, mu, color=S.COL[g], ls=ls, lw=lw, label=f"{g} {m}")
        ax.axhline(0, color=S.AXIS, lw=0.5, zorder=0)
        ax.axhline(1, color=S.AXIS, lw=0.5, zorder=0)
        ax.set_title(f"{'ABC'[j]}  {st}", loc="left", fontweight="bold")
        if j == 0:
            ax.set_ylabel("intensity (0 = section's ML, 1 = its IGL)")
        d, lo, hi = S.ols_curves(curves)
        axd.fill_between(x, lo, hi, color=S.CKO, alpha=0.18, lw=0)
        axd.plot(x, d, color=S.CKO, lw=1.2)
        axd.axhline(0, color=S.INK2, lw=0.6)
        if j == 0:
            axd.set_ylabel("cKO − WT, litter-adj.\n95 % CI")
        axd.set_title(f"{'DEF'[j]}  difference, {st}", loc="left", fontweight="bold")
        axd.set_xlabel("")
        if kind == "layer":
            for v in (1, 2):
                for a_ in (ax, axd):
                    a_.axvline(v, color=S.MUTED, lw=0.6, ls=(0, (2, 2)), zorder=0)
            y1 = ax.get_ylim()[1]
            for v, nm in ((0.5, "EGL"), (1.5, "ML"), (2.5, "IGL")):
                ax.text(v, 1.0, nm, transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=6, color=S.INK2)
            ax.set_xlim(0, 3)
            axd.set_xticks([0, 0.5, 1, 1.5, 2, 2.5, 3], ["", "", "", "", "", "", ""])
            ax.set_title(f"{'ABC'[j]}  {st}", loc="left", fontweight="bold", pad=10)
        else:
            ax.set_xlim(0, 240)
        for k_, (seg, nm) in enumerate(((slice(0, nb // 2), "oEGL"), (slice(nb // 2, nb), "iEGL"), (slice(nb, 2 * nb), "ML"), (slice(2 * nb, 3 * nb), "IGL"))) if kind == "layer" else []:
            am = pd.DataFrame([dict(genotype=g, imaging=m, v=v) for (g, m), A in curves.items() if A is not None for v in np.nanmean(A[:, seg], 1)])
            stats_rows.append(dict(measure=f"profile_{st}_{nm}_ml_igl_scale", **S.tests(am, "v")))
    fig.text(0.535, 0.055, xl, ha="center", va="bottom", fontsize=7, color=S.INK)
    ax = axs[0, 0]
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, 0.93), handlelength=1.6)
    fig.text(0.5, 0.005, "Median over each section's complete columns, mean over its sections per mouse; bands: ± SEM over mice (10x). "
             "D-F: cKO − WT on mouse means with litter as a factor, 95 % CI.", ha="center", va="bottom", fontsize=5.8, color=S.INK2)
    fig.subplots_adjust(left=0.085, right=0.985, top=0.92, bottom=0.15, wspace=0.28, hspace=0.42)
    return fig


def fig_intensity(df, stats_rows):
    stains = ["p27", "NeuN", "DAPI"]
    layers = [("oEGL", "outer EGL"), ("iEGL", "inner EGL"), ("IGL", "IGL")]
    fig, axs = plt.subplots(3, 3, figsize=(7.1, 6.4))
    for i, st in enumerate(stains):
        for j, (ly, nm) in enumerate(layers):
            ax = axs[i, j]
            col = f"{st}_{ly}_rel_ML"
            r = S.group_dots(ax, df, col, ylabel=f"{st}, {nm} / ML" if j == 0 else f"{nm} / ML", title=f"{'ABCDEFGHI'[3 * i + j]}  {st}, {nm}")
            ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.2f}"))
            stats_rows.append(dict(measure=col, **r))
    S.legend_groups(fig)
    fig.text(0.5, 0.004, "Mean stain in the layer over each complete column, median over the section, divided by the same for the section's ML\n"
             "(camera offset subtracted on the 10x set). Inner / outer EGL = deeper / outer half of each column's EGL. p: cKO vs WT on mouse means, litter as a factor.",
             ha="center", va="bottom", fontsize=5.8, color=S.INK2)
    fig.subplots_adjust(left=0.085, right=0.985, top=0.925, bottom=0.075, wspace=0.5, hspace=0.62)
    return fig


def fig_columns(df):
    """supplement: per-mouse distributions of column thickness, both litters (20x dashed)."""
    fig, axs = plt.subplots(1, 3, figsize=(7.1, 2.5))
    for ax, ly in zip(axs, ("EGL", "ML", "IGL")):
        allv, kk = [], 0
        for g in ("WT", "cKO"):
            for (a, m), d in df[df.genotype == g].groupby(["animal", "imaging"]):
                v = np.concatenate([(lambda C: C[f"thick_{ly}"][C["counted"]])(np.load(OUT / f"columns_{s}.npz")) for s in d.section])
                xs = np.sort(v)
                allv.append(xs)
                ax.plot(xs, np.linspace(0, 1, len(xs)), color=S.COL[g], lw=1.0 if m == "10x" else 0.9, alpha=0.9,
                        ls="-" if m == "10x" else (0, (3, 1.6)))
                f = 0.08 + 0.1 * kk
                kk += 1
                xi = xs[min(int(f * len(xs)), len(xs) - 1)]
                ax.text(xi + 1.5, f, a.split("-")[1], fontsize=5.5, color=S.INK2, ha="left", va="center")
        lo, hi = np.percentile(np.concatenate(allv), [0.5, 99.5])
        ax.set_xlim(max(0, lo - 4), hi + 6)
        ax.set_xlabel(f"{ly} thickness per column (µm)")
        ax.set_title(f"{'ABC'[['EGL', 'ML', 'IGL'].index(ly)]}  {ly}", loc="left", fontweight="bold")
        ax.grid(axis="x")
    axs[0].set_ylabel("fraction of columns")
    from matplotlib.lines import Line2D
    fig.legend(handles=[Line2D([], [], color=S.WT, lw=1.2, label="WT mouse, 2017 litter (C, F, G)"), Line2D([], [], color=S.CKO, lw=1.2, label="cKO mouse, 2017 litter (A, B, E)"),
                        Line2D([], [], color=S.CKO, lw=0.9, ls=(0, (3, 1.6)), label="cKO mouse, 2018 litter (slide positions 2, 4, 5)")],
               loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.0))
    fig.subplots_adjust(left=0.08, right=0.965, top=0.8, bottom=0.18, wspace=0.3)
    return fig


def fig_iegl(df, P, stats_rows):
    """the p27+ inner EGL and the p27- outer EGL (the paper's iEGL / oEGL), and thickness as trimmed means."""
    fig, axs = plt.subplots(2, 4, figsize=(7.1, 4.6))
    panels = [("iEGLp27_thick_tm_um", "p27+ inner EGL (\u00b5m)", "A  Inner EGL (p27+)"),
              ("oEGLp27_thick_tm_um", "p27\u2212 outer EGL (\u00b5m)", "B  Outer EGL (p27\u2212)"),
              ("iEGLp27_frac", "p27+ share of the EGL", "C  Inner share")]
    for ax, (col, yl, t) in zip(axs[0, :3], panels):
        r = S.group_dots(ax, df, col, ylabel=yl, title=t)
        stats_rows.append(dict(measure=col, **r))
        if col.endswith("frac"):
            ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.2f}"))
    # p27 across the EGL only (layer coordinate 0-1), ML = 0, IGL = 1
    ax = axs[0, 3]
    nb = int(P["nbin"])
    prof = _norm_prof(df, P, "layer", "p27")
    x = (np.arange(nb) + 0.5) / nb
    for (g, m), ls, lw in ((("WT", "10x"), "-", 1.4), (("cKO", "10x"), "-", 1.4), (("cKO", "20x"), (0, (3, 1.6)), 1.0)):
        A = _group_prof(df, prof, g, m)[:, :nb]
        mu = np.nanmean(A, 0)
        se = np.nanstd(A, 0, ddof=1) / np.sqrt(np.isfinite(A).sum(0))
        if m == "10x":
            ax.fill_between(x, mu - se, mu + se, color=S.COL[g], alpha=0.18, lw=0)
        ax.plot(x, mu, color=S.COL[g], ls=ls, lw=lw, label=f"{g} {m}")
    ax.axhline(QL.PIEGL[0], color=S.MUTED, lw=0.6, ls=(0, (2, 2)), zorder=0)
    ax.text(0.03, QL.PIEGL[0], "p27+ cut", fontsize=5.5, color=S.MUTED, ha="left", va="bottom")
    ax.set_xlim(0, 1)
    ax.set_xlabel("depth in EGL (0 top, 1 bottom)")
    ax.set_ylabel("p27 (0 = ML, 1 = IGL)")
    ax.set_title("D  p27 through the EGL", loc="left", fontweight="bold")
    ax.legend(loc="upper left", handlelength=1.6)
    for ax, (col, yl, t) in zip(axs[1], [("EGL_thick_tm_um", "EGL (\u00b5m)", "E  EGL, trimmed mean"),
                                         ("ML_thick_tm_um", "ML (\u00b5m)", "F  ML, trimmed mean"),
                                         ("IGL_thick_tm_um", "IGL (\u00b5m)", "G  IGL, trimmed mean"),
                                         ("iEGLp27_60_frac", "p27+ share, cut at 0.6", "H  Inner share, stricter")]):
        r = S.group_dots(ax, df, col, ylabel=yl, title=t)
        stats_rows.append(dict(measure=col, **r))
        if col.endswith("frac"):
            ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.2f}"))
    S.legend_groups(fig)
    fig.text(0.5, 0.004, "p27+ inner EGL: EGL samples whose p27 is past half way from the section's ML level to its IGL level (H: 0.6). "
             "A, B, E-G: 10 % trimmed means over columns\n(the medians of Figure 1 move in sample steps: 2 \u00b5m at 20x, 3.15 \u00b5m at 10x). p: cKO vs WT on mouse means, litter as a factor.",
             ha="center", va="bottom", fontsize=5.8, color=S.INK2)
    fig.subplots_adjust(left=0.075, right=0.985, top=0.885, bottom=0.15, wspace=0.6, hspace=0.75)
    return fig


def main():
    df, P = load()
    rows = []
    probs = {}
    probs["Q1_layer_size"] = S.save(fig_size(df, rows), "Q1_layer_size")
    probs["Q2_profiles_layer"] = S.save(fig_profiles(df, P, rows, "layer"), "Q2_profiles_layer")
    probs["Q2s_profiles_depth"] = S.save(fig_profiles(df, P, [], "depth"), "Q2s_profiles_depth")
    probs["Q3_layer_intensity"] = S.save(fig_intensity(df, rows), "Q3_layer_intensity")
    probs["Q1s_column_thickness"] = S.save(fig_columns(df), "Q1s_column_thickness")
    probs["Q5_inner_egl"] = S.save(fig_iegl(df, P, rows), "Q5_inner_egl")
    for col in ("EGL_stain_thick_tm_um", "ML_stain_thick_tm_um", "IGLtop_depth_tm_um", "IGLtop_stain_depth_tm_um"):   # label checks, no figure
        am = S.animal_means(df, col)
        rows.append(dict(measure=col, **S.tests(am, col)))
    st = pd.DataFrame(rows)
    st.to_csv(OUT / "stats.csv", index=False)
    pd.set_option("display.width", 200)
    print(st.round(4).to_string(index=False))
    return probs


if __name__ == "__main__":
    main()
