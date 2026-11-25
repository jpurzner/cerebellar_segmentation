"""
Tune and score depth-model relabelling against objective metrics.

Scores every configuration on:
  - gold accuracy / per-layer IoU on s1_2 (the hand-painted slide; its IGL
    and ML labels are painted swapped, so 4<->5 are swapped back on load);
  - layer continuity (share of each layer's area in its largest connected
    piece, as in python/dino/test_segmentation_continuity.py) on the 8 good
    20x and 11 10x current-pipeline slides;
  - fraction of cortex pixels changed.

Band models (how per-layer depth distributions are estimated):
  pooled_s       per normalised-arclength bin, pooled over 10x slides
  pooled_global  per layer, pooled over 10x slides, ignoring arclength
  self           per arclength window from the target slide's own labels
                 (falls back to the slide's global per-layer stats)
Pooled models for a 10x target leave that slide out.

Everything runs at the unfold's downsampled resolution; nothing is written
except the results table.
"""

import sys
import argparse
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image
from scipy import ndimage as ndi

sys.path.insert(0, str(Path(__file__).parent))
from ribbon_unfold_laplace import (                                # noqa: E402
    LBL, SANDBOX_20X, SANDBOX_10X, OUT_ROOT, REPO_ROOT, load_labels,
)

Image.MAX_IMAGE_PIXELS = None

SLIDES_20X = [
    "2018_05_22_s1_2_p27-0021", "2018_05_22_s1_4_p27-0020",
    "2018_05_22_s1_5_p27-0018", "2018_05_22_s2_2_p27",
    "2018_05_22_s2_4_p27-0024", "2018_05_22_s2_5_p27-0026",
    "2018_05_22_s3_2_p27-0001", "2018_05_22_s3_5_p27-0006",
]
SLIDES_10X = [
    "s4_A_p27", "s4_C_p27", "s4_F_p27", "s5_A_p27", "s5_B_p27", "s5_C_p27",
    "s5_E_p27", "s5_F_p27", "s5_G_p27", "s6_C_p27", "s6_G_p27",
]
GOLD = {
    "2018_05_22_s1_2_p27-0021":
        REPO_ROOT / "python" / "labelled"
        / "2018_05_22_s1_2_p27-0021_labelled.tiff",
}
ROOT_20X = OUT_ROOT / "ribbon_unfold_20x"
ROOT_10X = OUT_ROOT / "ribbon_unfold_10x"

LAYERS = [LBL["oEGL"], LBL["iEGL"], LBL["ML"], LBL["PCL"], LBL["IGL"]]
EGL_PAIR = {LBL["iEGL"], LBL["oEGL"]}
CONTINUITY = {"EGL": [2, 3], "ML": [5], "IGL": [4], "DWL": [6]}
GOLD_SCORED = [2, 3, 4, 5, 6, 7]
STD_FLOOR = 0.01


# --------------------------------------------------------------------------- #
# Data
# --------------------------------------------------------------------------- #

def load_slide(sid, mag):
    root = ROOT_20X if mag == "20x" else ROOT_10X
    sandbox = SANDBOX_20X if mag == "20x" else SANDBOX_10X
    c = np.load(root / sid / f"{sid}_ribbon_coords.npz", allow_pickle=False)
    ds = int(c["downsample"])
    before = load_labels(sandbox / f"{sid}_native_labels.tif")[::ds, ::ds]
    before = before.astype(np.uint8)
    keep = (c["reached"] & np.isfinite(c["s"]) & np.isfinite(c["psi"])
            & np.isin(c["layer"], LAYERS))
    slide = dict(
        sid=sid, mag=mag, ds=ds, before=before,
        cy=c["cy"][keep], cx=c["cx"][keep], layer=c["layer"][keep],
        psi=c["psi"][keep].astype(np.float64),
        s=c["s"][keep].astype(np.float64),
        s_norm=(c["s"][keep] / float(c["s_max"])).astype(np.float64),
    )
    assert np.array_equal(before[slide["cy"], slide["cx"]], slide["layer"]), \
        f"{sid}: unfold labels do not match the current label file"
    if sid in GOLD:
        g = np.array(Image.open(GOLD[sid])).astype(np.int16)[::ds, ::ds]
        gs = g.copy()
        gs[g == 4] = 5
        gs[g == 5] = 4
        slide["gold"] = gs
    return slide


# --------------------------------------------------------------------------- #
# Band models: each returns (med, std) arrays of shape (n_pix, n_layers)
# --------------------------------------------------------------------------- #

def _robust(v):
    q25, q50, q75 = np.percentile(v, [25, 50, 75])
    return q50, max((q75 - q25) / 1.349, STD_FLOOR)


def bands_pooled_global(pool, slide):
    n = len(slide["psi"])
    med = np.empty((n, len(LAYERS)))
    std = np.empty((n, len(LAYERS)))
    for k, L in enumerate(LAYERS):
        m, s = _robust(pool["psi"][pool["layer"] == L])
        med[:, k], std[:, k] = m, s
    return med, std


def bands_pooled_s(pool, slide, n_bins=200, min_count=30):
    edges = np.linspace(0, 1, n_bins + 1)
    pb = np.clip(np.searchsorted(edges, pool["s_norm"], "right") - 1, 0, n_bins - 1)
    tb = np.clip(np.searchsorted(edges, slide["s_norm"], "right") - 1, 0, n_bins - 1)
    gmed, gstd = bands_pooled_global(pool, {"psi": np.zeros(1)})
    med = np.empty((len(slide["psi"]), len(LAYERS)))
    std = np.empty_like(med)
    for k, L in enumerate(LAYERS):
        sel = pool["layer"] == L
        bm = np.full(n_bins, gmed[0, k])
        bs = np.full(n_bins, gstd[0, k])
        psi_L, b_L = pool["psi"][sel], pb[sel]
        order = np.argsort(b_L, kind="stable")
        psi_L, b_L = psi_L[order], b_L[order]
        bounds = np.searchsorted(b_L, np.arange(n_bins + 1))
        for b in range(n_bins):
            v = psi_L[bounds[b]:bounds[b + 1]]
            if len(v) >= min_count:
                bm[b], bs[b] = _robust(v)
        med[:, k], std[:, k] = bm[tb], bs[tb]
    return med, std


def bands_self(slide, window_px=120.0, step_px=20.0, min_count=30):
    """Per-layer bands from the slide's own labels in a sliding s-window."""
    s, psi, lay = slide["s"], slide["psi"], slide["layer"]
    centres = np.arange(0.0, s.max() + step_px, step_px)
    tb = np.clip(np.rint(s / step_px).astype(int), 0, len(centres) - 1)
    med = np.empty((len(psi), len(LAYERS)))
    std = np.empty_like(med)
    for k, L in enumerate(LAYERS):
        sel = lay == L
        if sel.sum() < min_count:
            med[:, k], std[:, k] = np.nan, 1.0
            continue
        gm, gs = _robust(psi[sel])
        order = np.argsort(s[sel])
        s_L, psi_L = s[sel][order], psi[sel][order]
        lo = np.searchsorted(s_L, centres - window_px / 2)
        hi = np.searchsorted(s_L, centres + window_px / 2)
        bm = np.full(len(centres), gm)
        bs = np.full(len(centres), gs)
        for i in np.nonzero(hi - lo >= min_count)[0]:
            bm[i], bs[i] = _robust(psi_L[lo[i]:hi[i]])
        med[:, k], std[:, k] = bm[tb], bs[tb]
    return med, std


# --------------------------------------------------------------------------- #
# Decision + metrics
# --------------------------------------------------------------------------- #

def relabel(slide, med, std, conf, margin, prior_boost, allow_egl_swap=False):
    psi = slide["psi"][:, None]
    z = (psi - med) / std
    ll = -0.5 * z * z - np.log(std)
    ll = np.where(np.isfinite(ll), ll, -np.inf)
    lab2col = np.full(10, -1)
    lab2col[LAYERS] = np.arange(len(LAYERS))
    cur_col = lab2col[slide["layer"]]
    ll[np.arange(len(ll)), cur_col] += prior_boost
    ll -= ll.max(axis=1, keepdims=True)
    p = np.exp(ll)
    p /= p.sum(axis=1, keepdims=True)
    top_col = p.argmax(axis=1)
    top = p[np.arange(len(p)), top_col]
    second = np.partition(p, -2, axis=1)[:, -2]
    pred = np.array(LAYERS)[top_col]
    q = (pred != slide["layer"]) & (top >= conf) & (top - second >= margin)
    if not allow_egl_swap:
        q &= ~(np.isin(slide["layer"], list(EGL_PAIR))
               & np.isin(pred, list(EGL_PAIR)))
    after = slide["before"].copy()
    after[slide["cy"][q], slide["cx"][q]] = pred[q]
    return after, int(q.sum())


def continuity(labels, ds):
    out = {}
    big = 500 / ds ** 2
    st = np.ones((3, 3), dtype=bool)
    for name, ls in CONTINUITY.items():
        m = np.isin(labels, ls)
        tot = int(m.sum())
        if tot == 0:
            out[name] = (0.0, 0)
            continue
        cc, n = ndi.label(m, structure=st)
        areas = np.bincount(cc.ravel())[1:]
        out[name] = (100.0 * areas.max() / tot, int((areas >= big).sum()))
    return out


def gold_scores(labels, gold):
    m = np.isin(gold, GOLD_SCORED)
    acc = float((labels[m] == gold[m]).mean())
    ious = {}
    for L in GOLD_SCORED:
        a, b = (labels == L) & m, (gold == L) & m
        u = (a | b).sum()
        ious[L] = float((a & b).sum() / u) if u else float("nan")
    return acc, ious


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #

def pool_from(slides, exclude=None):
    ss = [s for s in slides if s["sid"] != exclude]
    return dict(
        psi=np.concatenate([s["psi"] for s in ss]),
        layer=np.concatenate([s["layer"] for s in ss]),
        s_norm=np.concatenate([s["s_norm"] for s in ss]),
    )


def evaluate_config(slides, slides_10x, cfg, bands_cache):
    rows = []
    for sl in slides:
        key = (sl["sid"], cfg["model"], cfg.get("window_px"))
        if key not in bands_cache:
            if cfg["model"] == "self":
                bands_cache[key] = bands_self(sl, window_px=cfg["window_px"])
            else:
                pool = pool_from(slides_10x,
                                 exclude=sl["sid"] if sl["mag"] == "10x" else None)
                fn = bands_pooled_s if cfg["model"] == "pooled_s" else bands_pooled_global
                bands_cache[key] = fn(pool, sl)
        med, std = bands_cache[key]
        after, n_changed = relabel(sl, med, std, cfg["conf"], cfg["margin"],
                                   cfg["prior_boost"])
        row = dict(sid=sl["sid"], mag=sl["mag"],
                   changed_pct=100.0 * n_changed / len(sl["psi"]))
        for name, (pct, nl) in continuity(after, sl["ds"]).items():
            row[f"{name}_largest_pct"] = pct
            row[f"{name}_n_large"] = nl
        if "gold" in sl:
            acc, ious = gold_scores(after, sl["gold"])
            row["gold_acc"] = acc
            for L, v in ious.items():
                row[f"iou_{L}"] = v
        rows.append(row)
    return pd.DataFrame(rows)


def summarise(df, cfg):
    g = df[df["sid"].isin(GOLD)].iloc[0]
    out = dict(cfg)
    out["gold_acc"] = g["gold_acc"]
    out["gold_mIoU"] = float(np.nanmean([g[f"iou_{L}"] for L in GOLD_SCORED]))
    for L, nm in [(2, "iEGL"), (3, "oEGL"), (4, "IGL"), (5, "ML"), (6, "DWL"), (7, "PCL")]:
        out[f"iou_{nm}"] = g[f"iou_{L}"]
    for mag in ("20x", "10x"):
        d = df[df["mag"] == mag]
        for name in ("EGL", "ML", "IGL"):
            out[f"{name}_cont_{mag}"] = d[f"{name}_largest_pct"].mean()
        out[f"changed_{mag}"] = d["changed_pct"].mean()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true",
                    help="small grid for a first look")
    args = ap.parse_args()

    print("loading slides ...", flush=True)
    s20 = [load_slide(s, "20x") for s in SLIDES_20X]
    s10 = [load_slide(s, "10x") for s in SLIDES_10X]
    slides = s20 + s10

    configs = [dict(model="none", conf=2.0, margin=1.0, prior_boost=0.0)]
    if args.quick:
        grid = itertools.product(
            [("pooled_s", None), ("pooled_global", None),
             ("self", 120.0), ("self", 400.0)],
            [0.7], [0.3], [1.5])
    else:
        grid = itertools.product(
            [("pooled_s", None), ("pooled_global", None), ("self", 60.0),
             ("self", 120.0), ("self", 250.0), ("self", 500.0)],
            [0.5, 0.6, 0.7, 0.8, 0.9], [0.0, 0.15, 0.3], [0.0, 0.75, 1.5])
    for (model, w), conf, margin, pb in grid:
        c = dict(model=model, conf=conf, margin=margin, prior_boost=pb)
        if w is not None:
            c["window_px"] = w
        configs.append(c)

    cache, summaries = {}, []
    for i, cfg in enumerate(configs):
        if cfg["model"] == "none":
            df = pd.DataFrame([dict(
                sid=sl["sid"], mag=sl["mag"], changed_pct=0.0,
                **{f"{k}_largest_pct": v[0] for k, v in continuity(sl["before"], sl["ds"]).items()},
                **({"gold_acc": gold_scores(sl["before"], sl["gold"])[0],
                    **{f"iou_{L}": v for L, v in gold_scores(sl["before"], sl["gold"])[1].items()}}
                   if "gold" in sl else {}),
            ) for sl in slides])
        else:
            df = evaluate_config(slides, s10, cfg, cache)
        summaries.append(summarise(df, cfg))
        if (i + 1) % 10 == 0 or args.quick:
            print(f"  {i + 1}/{len(configs)} configs", flush=True)

    res = pd.DataFrame(summaries)
    out = OUT_ROOT / ("tune_relabel_quick.csv" if args.quick else "tune_relabel.csv")
    res.to_csv(out, index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 40)
    cols = ["model", "window_px", "conf", "margin", "prior_boost", "gold_acc",
            "gold_mIoU", "iou_ML", "iou_PCL", "iou_IGL", "iou_iEGL",
            "ML_cont_20x", "IGL_cont_20x", "EGL_cont_20x", "ML_cont_10x",
            "IGL_cont_10x", "changed_20x", "changed_10x"]
    cols = [c for c in cols if c in res.columns]
    print("\nbaseline (no relabel):")
    print(res[cols].iloc[[0]].round(3).to_string(index=False))
    print("\ntop configs by gold accuracy:")
    print(res[cols].sort_values("gold_acc", ascending=False).head(15).round(3).to_string(index=False))
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
