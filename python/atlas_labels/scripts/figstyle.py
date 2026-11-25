"""
Shared style for the quantification figures (2026-11-24): palette, fonts, group dot plots, statistics, and a
layout check run on every figure before it is saved (text overlaps, text outside the figure, PNG size).

Palette (dataviz reference palette, validated with validate_palette.py: CVD dE 24.7, normal 33.6, both pass):
WT = slot 1 blue, cKO = slot 2 orange. The litter / imaging set is a marker shape (2017 10x circle, 2018 20x triangle),
never a hue. 10x and 20x are pooled; the test is genotype adjusted for litter on mouse means (tests()).
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy import stats  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures/quant"

WT, CKO = "#2a78d6", "#eb6834"
INK, INK2, MUTED, GRID, AXIS = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
COL = {"WT": WT, "cKO": CKO}
MARK = {"10x": "o", "20x": "^"}
GROUPS = [("WT", "10x"), ("cKO", "10x"), ("cKO", "20x")]
MAX_PX = 1900

plt.rcParams.update({
    "font.family": "Arial", "font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7, "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5, "legend.fontsize": 6.5, "axes.edgecolor": AXIS, "axes.linewidth": 0.6,
    "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
    "axes.labelcolor": INK, "axes.titlecolor": INK, "text.color": INK, "axes.spines.top": False,
    "axes.spines.right": False, "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5,
    "ytick.major.size": 2.5, "axes.grid": True, "axes.grid.axis": "y", "grid.color": GRID, "grid.linewidth": 0.5,
    "axes.axisbelow": True, "legend.frameon": False, "savefig.facecolor": "white", "figure.facecolor": "white",
    "pdf.fonttype": 42, "svg.fonttype": "none", "lines.solid_capstyle": "round",
})


def welch(a, b):
    """difference of means (b - a), its 95 % CI (Welch) and the p value."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return dict(diff=np.nan, lo=np.nan, hi=np.nan, p=np.nan, na=len(a), nb=len(b))
    t = stats.ttest_ind(b, a, equal_var=False)
    va, vb = a.var(ddof=1) / len(a), b.var(ddof=1) / len(b)
    se = np.sqrt(va + vb)
    df = (va + vb) ** 2 / (va ** 2 / (len(a) - 1) + vb ** 2 / (len(b) - 1))
    q = stats.t.ppf(0.975, df)
    d = b.mean() - a.mean()
    return dict(diff=d, lo=d - q * se, hi=d + q * se, p=float(t.pvalue), na=len(a), nb=len(b))


def fmt_p(p):
    if not np.isfinite(p):
        return "n/a"
    return f"p = {p:.2g}" if p >= 0.001 else "p < 0.001"


def animal_means(df, col):
    """per animal: mean of its sections' values (sections are not independent)."""
    g = df.groupby(["animal", "genotype", "imaging"], as_index=False)[col].mean()
    return g


def ols_geno(y, cko, l2):
    """least squares y ~ 1 + cKO + litter(2018) on mouse values: the cKO coefficient, its 95 % CI and p. The litter
    term is dropped when only one litter is present."""
    y, cko, l2 = (np.asarray(v, float) for v in (y, cko, l2))
    ok = np.isfinite(y)
    y, cko, l2 = y[ok], cko[ok], l2[ok]
    cols = [np.ones_like(y), cko] + ([l2] if 0 < l2.sum() < len(l2) else [])
    X = np.column_stack(cols)
    n, k = X.shape
    if n - k < 1 or np.linalg.matrix_rank(X) < k:
        return dict(diff=np.nan, lo=np.nan, hi=np.nan, p=np.nan, df=0)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ beta
    s2 = res @ res / (n - k)
    se = np.sqrt(s2 * np.linalg.inv(X.T @ X)[1, 1])
    t = beta[1] / se if se > 0 else np.inf
    q = stats.t.ppf(0.975, n - k)
    return dict(diff=float(beta[1]), lo=float(beta[1] - q * se), hi=float(beta[1] + q * se), p=float(2 * stats.t.sf(abs(t), n - k)), df=n - k)


def tests(am, col):
    """genotype tests on mouse means. Main (diff, lo, hi, p): cKO - WT adjusted for litter, both litters pooled.
    Also: pooled without the litter term (_pooled) and the 2017 litter alone (_10x, Welch)."""
    wt, ko = am.genotype == "WT", am.genotype == "cKO"
    r = ols_geno(am[col], ko.astype(float), (am.imaging == "20x").astype(float))
    r.update(na=int((wt & am[col].notna()).sum()), nb=int((ko & am[col].notna()).sum()))
    w = welch(am[wt][col], am[ko][col])
    r.update(diff_pooled=w["diff"], lo_pooled=w["lo"], hi_pooled=w["hi"], p_pooled=w["p"])
    w = welch(am[wt & (am.imaging == "10x")][col], am[ko & (am.imaging == "10x")][col])
    r.update(diff_10x=w["diff"], lo_10x=w["lo"], hi_10x=w["hi"], p_10x=w["p"])
    return r


def ols_curves(curves):
    """per x position: cKO - WT adjusted for litter, from {(genotype, imaging): (mice, x) array}."""
    ys, ck, l2 = [], [], []
    for (g, m), A in curves.items():
        if A is None:
            continue
        ys.append(A)
        ck += [g == "cKO"] * len(A)
        l2 += [m == "20x"] * len(A)
    Y = np.concatenate(ys)
    out = [ols_geno(Y[:, k], ck, l2) for k in range(Y.shape[1])]
    return (np.array([o[f] for o in out]) for f in ("diff", "lo", "hi"))


OFF = {"10x": -0.15, "20x": 0.15}
LITTER = {"10x": "2017 litter, 10x", "20x": "2018 litter, 20x"}


def group_dots(ax, df, col, ylabel=None, title=None, test=True, ylim=None, fmt="{:.0f}", show_n=False):
    """WT | cKO, each with its 2017 (10x, filled circles, left) and 2018 (20x, hollow triangles, right) mice:
    sections small and light, mice large (mean of the mouse's sections); bar = mean +- SEM over all the genotype's
    mice. p: cKO vs WT adjusted for litter (tests())."""
    rng = np.random.default_rng(1)
    am = animal_means(df, col)
    labels = []
    for gi, g in enumerate(("WT", "cKO")):
        ag = am[am.genotype == g][col].to_numpy(float)
        ns = []
        for m in ("10x", "20x"):
            s = df[(df.genotype == g) & (df.imaging == m)][col].to_numpy(float)
            a = am[(am.genotype == g) & (am.imaging == m)][col].to_numpy(float)
            ns.append(len(a))
            if not len(s):
                continue
            x = gi + OFF[m]
            hollow = m == "20x"
            ax.scatter(x + rng.uniform(-0.07, 0.07, len(s)), s, s=8, marker=MARK[m], facecolor="none" if hollow else COL[g],
                       edgecolor=COL[g], linewidth=0.6, alpha=0.5, zorder=3)
            ax.scatter(x + rng.uniform(-0.025, 0.025, len(a)), a, s=28, marker=MARK[m], facecolor="white" if hollow else COL[g],
                       edgecolor=COL[g] if hollow else "white", linewidth=1.0 if hollow else 0.9, zorder=4)
        if len(ag):
            mu = np.nanmean(ag)
            sem = np.nanstd(ag, ddof=1) / np.sqrt(np.isfinite(ag).sum()) if len(ag) > 1 else 0
            ax.plot([gi - 0.33, gi + 0.33], [mu, mu], color=INK, lw=1.0, zorder=5)
            ax.plot([gi, gi], [mu - sem, mu + sem], color=INK, lw=0.8, zorder=5)
        labels.append(f"{g}\n{ns[0]} + {ns[1]} mice")
    ax.set_xticks([0, 1], labels)
    ax.set_xlim(-0.55, 1.55)
    ax.tick_params(axis="x", length=0)
    if ylim:
        ax.set_ylim(*ylim)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title, loc="left", fontweight="bold")
    res = None
    if test:
        res = tests(am, col)
        y0, y1 = ax.get_ylim()
        top = np.nanmax(np.r_[df[col].to_numpy(float), am[col].to_numpy(float)])
        yb = top + 0.06 * (y1 - y0)
        ax.plot([0, 0, 1, 1], [yb - 0.02 * (y1 - y0), yb, yb, yb - 0.02 * (y1 - y0)], color=MUTED, lw=0.6)
        ax.text(0.5, yb + 0.015 * (y1 - y0), fmt_p(res["p"]), ha="center", va="bottom", fontsize=6, color=INK2)
        ax.set_ylim(y0, max(y1, yb + 0.14 * (y1 - y0)))
    return res


def legend_groups(fig, y=0.995, x=0.5):
    """one shared legend: genotype by colour, litter / imaging set by shape, section vs mouse by size."""
    from matplotlib.lines import Line2D
    h = [Line2D([], [], marker="o", ls="", markerfacecolor=WT, markeredgecolor="white", markersize=5.5, label="WT"),
         Line2D([], [], marker="o", ls="", markerfacecolor=CKO, markeredgecolor="white", markersize=5.5, label="cKO"),
         Line2D([], [], marker="o", ls="", markerfacecolor=MUTED, markeredgecolor="white", markersize=5, label=LITTER["10x"]),
         Line2D([], [], marker="^", ls="", markerfacecolor="white", markeredgecolor=MUTED, markersize=5, label=LITTER["20x"]),
         Line2D([], [], marker="o", ls="", markerfacecolor=MUTED, alpha=0.5, markersize=3, label="section (large: mouse)"),
         Line2D([], [], color=INK, lw=1.0, label="mean ± SEM, all mice")]
    fig.legend(handles=h, loc="upper center", bbox_to_anchor=(x, y), ncol=6, handletextpad=0.3, columnspacing=1.1)


def qa(fig, name=""):
    """layout check: overlapping text, text outside the figure. Returns a list of problems."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    W, H = fig.canvas.get_width_height()
    hidden = set()                         # tick labels of ticks outside the view are never drawn
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            lo, hi = sorted(axis.get_view_interval())
            for tk in axis.get_major_ticks() + axis.get_minor_ticks():
                if not (lo - 1e-9 <= tk.get_loc() <= hi + 1e-9):
                    hidden |= {tk.label1, tk.label2}
    texts = []
    for t in fig.findobj(matplotlib.text.Text):
        if not t.get_visible() or not t.get_text().strip() or t in hidden:
            continue
        try:
            bb = t.get_window_extent(r)
        except Exception:
            continue
        if bb.width < 1 or bb.height < 1:
            continue
        texts.append((t, bb))
    probs = []
    for t, bb in texts:
        if bb.x0 < -1 or bb.y0 < -1 or bb.x1 > W + 1 or bb.y1 > H + 1:
            probs.append(f"outside figure: {t.get_text()[:40]!r}")
    for i in range(len(texts)):
        for j in range(i + 1, len(texts)):
            a, b = texts[i][1], texts[j][1]
            ix = min(a.x1, b.x1) - max(a.x0, b.x0)
            iy = min(a.y1, b.y1) - max(a.y0, b.y0)
            if ix > 1.5 and iy > 1.5:
                probs.append(f"overlap: {texts[i][0].get_text()[:30]!r} x {texts[j][0].get_text()[:30]!r} ({ix:.0f}x{iy:.0f} px)")
    # text (not tick labels) over data marks: scatter points and line vertices of any axes
    pts = []
    for ax in fig.axes:
        for c in ax.collections:
            if isinstance(c, matplotlib.collections.PathCollection) and len(c.get_offsets()):
                pts.append(c.get_offset_transform().transform(c.get_offsets()))
        for ln in ax.lines:
            xy = ln.get_xydata()
            if len(xy) and ln.get_marker() not in (None, "None", "", " "):
                pts.append(ln.get_transform().transform(xy))
    P = np.concatenate(pts) if pts else np.zeros((0, 2))
    ticks = set()
    for ax in fig.axes:
        ticks |= set(ax.get_xticklabels()) | set(ax.get_yticklabels())
    for t, bb in texts:
        if t in ticks or not len(P):
            continue
        inside = ((P[:, 0] > bb.x0 - 1) & (P[:, 0] < bb.x1 + 1) & (P[:, 1] > bb.y0 - 1) & (P[:, 1] < bb.y1 + 1)).sum()
        if inside:
            probs.append(f"text over {inside} data marks: {t.get_text()[:40]!r}")
    return probs


def save(fig, name, dpi=250):
    """PNG (long side <= MAX_PX), PDF and SVG in figures/quant; returns the layout problems."""
    FIG.mkdir(parents=True, exist_ok=True)
    w, h = fig.get_size_inches()
    dpi = min(dpi, int(MAX_PX / max(w, h)))
    fig.set_dpi(dpi)
    probs = qa(fig, name)
    fig.savefig(FIG / f"{name}.png", dpi=dpi)
    fig.savefig(FIG / f"{name}.pdf")
    fig.savefig(FIG / f"{name}.svg")
    from PIL import Image
    size = Image.open(FIG / f"{name}.png").size
    assert max(size) <= MAX_PX, size
    plt.close(fig)
    print(f"{name}: {size[0]}x{size[1]} px at {dpi} dpi; layout problems: {len(probs)}")
    for p in probs[:20]:
        print("   ", p)
    return probs
