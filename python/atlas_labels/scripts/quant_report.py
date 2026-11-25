"""
The review page for the layer and ML-nuclei quantification (2026-11-24): builds outputs/quant/report/index.html
from outputs/quant/*.csv and figures/quant/*.png (copied beside the page). Every number on the page is read from
those files.

    python scripts/quant_report.py
"""
import html
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import quant_layers as QL  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = QL.OUT
FIG = ROOT / "figures/quant"
REP = OUT / "report"

NAMES = {
    "EGL_stain_thick_tm_um": ("EGL, stain-only boundary", "µm"), "ML_stain_thick_tm_um": ("ML, stain-only boundaries", "µm"),
    "IGLtop_depth_tm_um": ("IGL top below surface, labels", "µm"), "IGLtop_stain_depth_tm_um": ("IGL top below surface, stains", "µm"),
    "EGL_thick_tm_um": ("EGL thickness, trimmed mean", "µm"), "ML_thick_tm_um": ("ML thickness, trimmed mean", "µm"),
    "IGL_thick_tm_um": ("IGL thickness, trimmed mean", "µm"), "iEGLp27_thick_tm_um": ("p27+ inner EGL", "µm"),
    "oEGLp27_thick_tm_um": ("p27− outer EGL", "µm"), "iEGLp27_frac": ("p27+ share of EGL (cut 0.5)", ""),
    "iEGLp27_60_frac": ("p27+ share of EGL (cut 0.6)", ""),
    "EGL_thick_um": ("EGL thickness", "µm"), "ML_thick_um": ("ML thickness", "µm"), "IGL_thick_um": ("IGL thickness", "µm"),
    "layers_area_mm2": ("labelled area", "mm²"), "EGL_area_frac": ("EGL share of area", ""), "ML_area_frac": ("ML share of area", ""),
    "IGL_area_frac": ("IGL share of area", ""), "WM_area_frac": ("white matter share", ""),
    "p27_oEGL_rel_ML": ("p27, outer EGL / ML", ""), "p27_iEGL_rel_ML": ("p27, inner EGL / ML", ""), "p27_IGL_rel_ML": ("p27, IGL / ML", ""),
    "NeuN_oEGL_rel_ML": ("NeuN, outer EGL / ML", ""), "NeuN_iEGL_rel_ML": ("NeuN, inner EGL / ML", ""), "NeuN_IGL_rel_ML": ("NeuN, IGL / ML", ""),
    "DAPI_oEGL_rel_ML": ("DAPI, outer EGL / ML", ""), "DAPI_iEGL_rel_ML": ("DAPI, inner EGL / ML", ""), "DAPI_IGL_rel_ML": ("DAPI, IGL / ML", ""),
    "ecc": ("eccentricity, NeuN+ ML nuclei", ""), "angle": ("angle to IGL perpendicular, all NeuN+", "°"), "radial30": ("share within 30° of radial", ""),
    "neun_frac": ("NeuN+ share of ML nuclei", ""), "angle_clean_geometry": ("angle, straight ML only", "°"),
    "angle_egl": ("angle to EGL perpendicular, all NeuN+", "°"), "pooledcrops_ecc": ("eccentricity, 2017 sections + 2018 crops", ""), "pooledcrops_elong_frac": ("share elongated, 2017 sections + 2018 crops", ""), "pooledcrops_angle_elong": ("angle of elongated nuclei, 2017 sections + 2018 crops", "°"), "pooledcrops_neun_frac": ("NeuN+ share, 2017 sections + 2018 crops", ""), "elong_frac": ("share elongated (ecc ≥ 0.85)", ""), "angle_elong": ("angle of elongated nuclei", "°"), "minor_um": ("nucleus width", "µm"),
}


def _f(v, unit=""):
    if v is None or not np.isfinite(v):
        return "–"
    a = abs(v)
    s = f"{v:.3f}" if a < 1 else f"{v:.2f}" if a < 10 else f"{v:.1f}"
    return s


def _p(v):
    return "–" if v is None or not np.isfinite(v) else "%.2g" % v


def group_means(df, col):
    """mouse means: WT (both litters), cKO 2017 (10x), cKO 2018 (20x)."""
    am = df.groupby(["animal", "genotype", "imaging"], as_index=False)[col].mean()
    g = lambda gg, m: am[(am.genotype == gg) & (am.imaging == m)][col]
    return am[am.genotype == "WT"][col].mean(), g("cKO", "10x").mean(), g("cKO", "20x").mean()


def stats_table(st, frames):
    rows = []
    for _, r in st.iterrows():
        nm = r.measure
        if nm not in NAMES:
            continue
        label, unit = NAMES[nm]
        src = frames.get(nm)
        wt = ck = c20 = np.nan
        if src is not None:
            df, col = src
            wt, ck, c20 = group_means(df, col)
        sig = r.p < 0.05 if np.isfinite(r.p) else False
        rows.append(f"<tr class=\"{'sig' if sig else ''}\"><th scope=\"row\">{html.escape(label)}{' (' + unit + ')' if unit else ''}</th>"
                    f"<td class=\"wt\">{_f(wt)}</td><td class=\"ko\">{_f(ck)}</td><td class=\"k20\">{_f(c20)}</td><td>{_f(r['diff'])}</td>"
                    f"<td>{_f(r.lo)} to {_f(r.hi)}</td><td>{_p(r.p)}</td><td class=\"muted\">{_p(r.get('p_pooled'))}</td><td class=\"muted\">{_p(r.get('p_10x'))}</td></tr>")
    return "\n".join(rows)


def main():
    REP.mkdir(parents=True, exist_ok=True)
    (REP / "img").mkdir(exist_ok=True)
    for p in FIG.glob("*.png"):
        shutil.copy(p, REP / "img" / p.name)
    for p in (FIG / "qc").glob("*.png"):
        shutil.copy(p, REP / "img" / p.name)
    sec = pd.read_csv(OUT / "sections.csv")
    st = pd.read_csv(OUT / "stats.csv")
    ml_st = pd.read_csv(OUT / "ml_stats.csv") if (OUT / "ml_stats.csv").exists() else pd.DataFrame()
    mlp = pd.read_csv(OUT / "ml_nuclei_per_section.csv") if (OUT / "ml_nuclei_per_section.csv").exists() else pd.DataFrame()
    frames = {c: (sec, c) for c in sec.columns}
    for c in ("ecc", "angle", "radial30", "neun_frac", "angle_egl", "minor", "elong_frac", "angle_elong", "radial_elong"):
        if c in mlp:
            frames[c] = (mlp, c)
    frames["minor_um"] = (mlp, "minor") if "minor" in mlp else None
    t_size = stats_table(st[st.measure.isin(["EGL_thick_um", "ML_thick_um", "IGL_thick_um", "layers_area_mm2", "EGL_area_frac", "ML_area_frac",
                                              "IGL_area_frac", "WM_area_frac"])], frames)
    t_iegl = stats_table(st[st.measure.isin(["iEGLp27_thick_tm_um", "oEGLp27_thick_tm_um", "iEGLp27_frac", "iEGLp27_60_frac", "EGL_thick_tm_um",
                                              "ML_thick_tm_um", "IGL_thick_tm_um"])], frames)
    t_int = stats_table(st[st.measure.str.endswith("rel_ML")], frames)
    t_ml = stats_table(ml_st[ml_st.measure.isin(["ecc", "elong_frac", "angle_elong", "angle", "angle_egl", "neun_frac", "minor_um"])], frames) if len(ml_st) else ""
    fr_c = {}
    if (OUT / "ml_nuclei_crops.csv.gz").exists() and len(mlp):     # 2017 whole sections + 2018 crops (the 2018 WT is only in the crops)
        import ml_nuclei_figures as MF
        both = pd.concat([mlp[mlp.imaging == "10x"], MF.per_section(pd.read_csv(OUT / "ml_nuclei_crops.csv.gz").assign(imaging="20x"))], ignore_index=True)
        fr_c = {f"pooledcrops_{c}": (both, c) for c in ("ecc", "elong_frac", "angle_elong", "neun_frac")}
        frames.update(fr_c)
    t_mlc = stats_table(ml_st[ml_st.measure.isin(list(fr_c))], fr_c) if len(ml_st) else ""
    geno = "\n".join(f"<tr><th scope=\"row\">{r.section}</th><td>{r.animal}</td><td class=\"{'wt' if r.genotype == 'WT' else 'ko'}\">{r.genotype}</td>"
                     f"<td>{r.imaging}</td><td>{r.px_um:.4f}</td><td>{r.counted_mm:.1f}</td></tr>" for r in sec.itertuples())
    t_check = stats_table(st[st.measure.isin(["EGL_thick_tm_um", "EGL_stain_thick_tm_um", "ML_thick_tm_um", "ML_stain_thick_tm_um",
                                               "IGLtop_depth_tm_um", "IGLtop_stain_depth_tm_um"])], frames)
    per = sec.merge(mlp[["section", "ecc", "elong_frac", "angle_elong", "neun_frac"]], on="section", how="left") if len(mlp) else sec
    t_sec = "\n".join(
        f"<tr><th scope=\"row\">{r.section}</th><td>{r.animal}</td><td class=\"{'wt' if r.genotype == 'WT' else 'ko'}\">{r.genotype}</td><td>{r.imaging}</td>"
        f"<td>{r.EGL_thick_tm_um:.1f}</td><td>{r.iEGLp27_thick_tm_um:.1f}</td><td>{r.ML_thick_tm_um:.1f}</td><td>{r.IGL_thick_tm_um:.1f}</td>"
        f"<td>{r.p27_iEGL_rel_ML:.2f}</td><td>{r.p27_IGL_rel_ML:.2f}</td><td>{r.ecc:.3f}</td><td>{r.angle_elong:.1f}</td><td>{r.neun_frac:.2f}</td></tr>"
        for r in per.itertuples())
    data = dict(t_size=t_size, t_iegl=t_iegl, t_int=t_int, t_ml=t_ml, t_mlc=t_mlc, geno=geno, t_check=t_check, t_sec=t_sec)
    tpl = (Path(__file__).parent / "quant_report_template.html").read_text()
    for k, v in data.items():
        tpl = tpl.replace("{{" + k + "}}", v)
    # numbers quoted in the text
    nums = {}
    pxs = pd.read_csv(ROOT / "outputs/px_scale/allen_fit_summary.csv") if (ROOT / "outputs/px_scale/allen_fit_summary.csv").exists() else None
    if pxs is not None:
        m = pxs.groupby(["animal", "imaging"]).scale_top5.mean().reset_index()
        nums["allen.s10"] = m[m.imaging == "10x"].scale_top5.mean()
        nums["allen.s20"] = m[m.imaging == "20x"].scale_top5.mean()
        nums["allen.rel10"] = QL.PX_TRUE["10x"] / QL.PX_BUILT["10x"] / nums["allen.s10"]
        nums["allen.rel20"] = 1 / nums["allen.s20"]
    for _, r in pd.concat([st, ml_st]).iterrows():
        for f in ("diff", "lo", "hi", "p", "diff_pooled", "p_pooled", "diff_10x", "p_10x"):
            if f in r and pd.notna(r[f]):
                nums[f"{r.measure}.{f}"] = r[f]
    for k, (df, col) in [(k, v) for k, v in frames.items() if v is not None]:
        try:
            wt, ck, c20 = group_means(df, col)
            nums[f"{k}.wt"], nums[f"{k}.ko"], nums[f"{k}.k20"] = wt, ck, c20
        except Exception:
            pass
    import re

    def sub(m):
        key, fmt = m.group(1), m.group(2) or ".1f"
        v = nums.get(key)
        if v is None or not np.isfinite(v):
            return "–"
        return format(v, fmt)
    if (OUT / "ml_nuclei_crops.csv.gz").exists():                 # the 20x crops: mouse means, WT (1) and cKO (4)
        import ml_nuclei_figures as MF
        Pc = MF.per_section(pd.read_csv(OUT / "ml_nuclei_crops.csv.gz").assign(imaging="20x"))
        for col in ("ecc", "elong_frac", "angle_elong", "radial_elong", "neun_frac"):
            am = Pc.groupby(["animal", "genotype"])[col].mean().reset_index()
            nums[f"crops.{col}.wt"] = am[am.genotype == "WT"][col].mean()
            ko = am[am.genotype == "cKO"][col]
            nums[f"crops.{col}.komin"], nums[f"crops.{col}.komax"] = ko.min(), ko.max()
    tpl = re.sub(r"\{\{n:([A-Za-z0-9_.]+)(?:\|([^}]+))?\}\}", sub, tpl)
    left = re.findall(r"\{\{[^}]+\}\}", tpl)
    assert not left, left
    (REP / "index.html").write_text(tpl)
    print("report", REP / "index.html", len(tpl), "chars; images", len(list((REP / "img").glob("*.png"))))


if __name__ == "__main__":
    main()
