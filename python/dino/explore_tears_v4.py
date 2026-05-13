"""V4 anomaly detector — final tuning before MATLAB port.

V4 = ground-truth-tuned tear detector:
  1. Morphological tissue mask (close + fill holes), eroded 50 px from edge
  2. Per-channel flat-field correction (handles s2_5 inhomogeneity)
  3. combined = mean of 3 ff-corrected channels
  4. Per-pixel: combined < 0.10 inside eroded tissue
  5. Connected components → shape stats
  6. Classify:
       TEAR:   area >= 1000 AND area <= 50000 AND (ecc > 0.95 OR axR > 4)
       BLEACH: area >= 2000 AND area <= 50000 AND ecc < 0.85 AND extent > 0.5

Validation: render per-slide panels for s1_2 (compare to gold), s1_5, s2_5, s3_5.
For s2_5 specifically also show before/after flat-field to demonstrate it kills
the inhomogeneity-induced fake dropouts.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage as ndi
from skimage import measure
from scipy.ndimage import binary_fill_holes
import tifffile

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from slide_manifest import manifest

ROOT = Path(__file__).resolve().parent
LABELLED_DIR = ROOT.parent / "labelled"
OUT = ROOT.parent / "anomaly_explore_v4"
OUT.mkdir(exist_ok=True)
PX_UM = 0.5119049

LABEL_COLORS = np.array([
    [0,   0,   0], [0,   0, 128], [0, 128, 255], [0, 255, 255],
    [255, 255, 0], [128, 255, 128],[255, 128, 0], [255,   0, 0],
], dtype=np.uint8)


def to_unit(im):
    im = im.astype(np.float32); lo, hi = im.min(), im.max()
    return (im - lo) / max(hi - lo, 1e-9)


def find_gold(slide_id):
    cands = sorted(
        list(LABELLED_DIR.glob(f"{slide_id}_corrected.tif*")) +
        list(LABELLED_DIR.glob(f"{slide_id}_labelled.tif*")) +
        list(LABELLED_DIR.glob(f"{slide_id}_labelld.tif*")),
        key=lambda p: p.stat().st_mtime, reverse=True)
    cands = [c for c in cands if c.suffix in (".tif", ".tiff")]
    return cands[0] if cands else None


def flat_field(channel, tissue, sigma=200):
    """Tissue-aware flat-field correction.
    Estimate slow-varying background only inside tissue, divide channel by it,
    rescale to mean. Removes illumination inhomogeneity."""
    ch_t = channel * tissue.astype(np.float32)
    bg_num = ndi.gaussian_filter(ch_t, sigma=sigma)
    bg_den = ndi.gaussian_filter(tissue.astype(np.float32), sigma=sigma)
    bg = bg_num / np.maximum(bg_den, 0.05)
    bg_mean = bg[tissue].mean() if tissue.sum() > 0 else 1.0
    out = channel / np.maximum(bg, 0.02) * bg_mean
    return np.clip(out, 0, 1), bg


def detect_anomalies(p27, neun, dapi, tissue_mask, edge_erode_px=50):
    """V4 detection. Returns dict of masks + stats."""
    # Flat-field each channel
    p27_ff, _ = flat_field(p27, tissue_mask)
    neun_ff, _ = flat_field(neun, tissue_mask)
    dapi_ff, dapi_bg = flat_field(dapi, tissue_mask)

    combined_ff = (p27_ff + neun_ff + dapi_ff) / 3

    # Erode tissue mask to exclude border artifacts
    tissue_inner = ndi.binary_erosion(tissue_mask,
                                       structure=np.ones((edge_erode_px, edge_erode_px)))

    # Per-pixel dropout
    drop = (combined_ff < 0.10) & tissue_inner
    drop = ndi.binary_closing(drop, structure=np.ones((3, 3)))

    # Connected components
    cc_lbl = measure.label(drop, connectivity=2)
    props = measure.regionprops(cc_lbl)

    tear_mask = np.zeros_like(drop)
    bleach_mask = np.zeros_like(drop)
    big_blobs = []
    n_tears = n_bleach = 0
    for p in props:
        if p.area < 500: continue
        major = p.axis_major_length
        minor = max(p.axis_minor_length, 0.5)
        ax = major / minor
        info = {"label": p.label, "area": p.area, "ecc": p.eccentricity,
                "ax": ax, "extent": p.extent, "bbox": p.bbox}
        # V4.1: drop upper area cap. Eccentricity alone discriminates real
        # tears (>= 0.95) from ML inter-cell networks (< 0.95).
        if p.area >= 1000 and (p.eccentricity >= 0.95 or ax >= 4):
            tear_mask |= (cc_lbl == p.label); n_tears += 1
        elif p.area >= 2000 and p.eccentricity < 0.85 and p.extent > 0.5:
            bleach_mask |= (cc_lbl == p.label); n_bleach += 1
        if p.area >= 1000:
            big_blobs.append(info)

    return {"tear": tear_mask, "bleach": bleach_mask,
            "drop_raw": drop, "combined_ff": combined_ff,
            "dapi_ff": dapi_ff, "dapi_bg": dapi_bg,
            "n_tears": n_tears, "n_bleach": n_bleach,
            "big_blobs": big_blobs, "tissue_inner": tissue_inner}


def explore(s):
    print(f"\n=== {s.slide_id} ===")
    stack = tifffile.imread(s.input_path)
    H, W = stack.shape[1], stack.shape[2]
    p27  = to_unit(stack[0].astype(np.float32))
    neun = to_unit(stack[1].astype(np.float32))
    dapi = to_unit(stack[2].astype(np.float32))
    rgb = np.stack([p27, neun, dapi], axis=-1)
    rgb /= max(rgb.max(), 1e-9)

    # Build morphological tissue mask
    combined_raw = (p27 + neun + dapi) / 3
    tissue = combined_raw > 0.15
    tissue = ndi.binary_closing(tissue, structure=np.ones((30, 30)))
    tissue = binary_fill_holes(tissue)

    res = detect_anomalies(p27, neun, dapi, tissue)
    print(f"  tears: {res['n_tears']}  bleach: {res['n_bleach']}")
    print(f"  tear pixels: {res['tear'].sum()} ({100*res['tear'].sum()/H/W:.3f}%)")
    print(f"  bleach pixels: {res['bleach'].sum()} ({100*res['bleach'].sum()/H/W:.3f}%)")
    print(f"  top 8 big blobs (area, ecc, axR, extent, bbox, classified):")
    for b in sorted(res["big_blobs"], key=lambda x: -x["area"])[:8]:
        cls = "TEAR" if b["area"]>=1000 and (b["ecc"]>=0.95 or b["ax"]>=4) else \
              ("BLEACH" if b["area"]>=2000 and b["ecc"]<0.85 and b["extent"]>0.5 else "?")
        print(f"    {b['area']:>7} ecc={b['ecc']:.2f} axR={b['ax']:.1f} ext={b['extent']:.2f} bbox={b['bbox']} {cls}")

    # Render
    gold_path = find_gold(s.slide_id)
    gold = tifffile.imread(gold_path).astype(np.uint8) if gold_path else None
    fig, axes = plt.subplots(2, 4, figsize=(28, 14))

    axes[0,0].imshow(rgb); axes[0,0].set_title("RGB", fontsize=11)
    axes[0,1].imshow(dapi, cmap="gray", vmin=0, vmax=0.5)
    axes[0,1].set_title("DAPI raw normalized", fontsize=11)
    axes[0,2].imshow(res["dapi_bg"], cmap="gray", vmin=0, vmax=0.3)
    axes[0,2].set_title("DAPI bg (tissue-aware sig=200) → inhomogeneity map", fontsize=10)
    axes[0,3].imshow(res["dapi_ff"], cmap="gray", vmin=0, vmax=0.5)
    axes[0,3].set_title("DAPI flat-field corrected", fontsize=11)

    axes[1,0].imshow(res["combined_ff"], cmap="gray", vmin=0, vmax=0.5)
    axes[1,0].set_title("combined ff", fontsize=11)
    axes[1,1].imshow(res["drop_raw"], cmap="hot")
    axes[1,1].set_title(f"drop mask combined_ff<0.10\n{100*res['drop_raw'].sum()/res['drop_raw'].size:.3f}%", fontsize=10)
    overlay = rgb.copy()
    overlay[res["tear"]] = 0.2 * overlay[res["tear"]] + 0.8 * np.array([1, 0, 1])
    overlay[res["bleach"]] = 0.2 * overlay[res["bleach"]] + 0.8 * np.array([1, 0.5, 0])
    axes[1,2].imshow(np.clip(overlay, 0, 1))
    axes[1,2].set_title(f"V4 detect: TEARS magenta(n={res['n_tears']}) BLEACH orange(n={res['n_bleach']})", fontsize=10)

    if gold is not None:
        Hg = min(gold.shape[0], H); Wg = min(gold.shape[1], W)
        axes[1,3].imshow(LABEL_COLORS[gold[:Hg,:Wg]])
        axes[1,3].set_title("GOLD", fontsize=11)
    else:
        axes[1,3].imshow(combined_raw, cmap="gray", vmin=0, vmax=0.5)
        axes[1,3].set_title("(no gold) raw combined", fontsize=11)

    for ax in axes.ravel(): ax.axis("off")
    fig.suptitle(s.slide_id, fontsize=14)
    fig.tight_layout()
    fig.savefig(OUT / f"{s.slide_id}_v4.png", dpi=80)
    plt.close(fig)

    # Zoom on detected tears
    if res['n_tears'] > 0:
        cc = measure.label(res['tear'])
        tprops = sorted(measure.regionprops(cc), key=lambda x: -x.area)[:6]
        n = len(tprops)
        fig, axes = plt.subplots(1, n, figsize=(5*n, 5))
        if n == 1: axes = [axes]
        for i, p in enumerate(tprops):
            y0,x0,y1,x1 = p.bbox
            pad = 80
            y0p = max(0, y0-pad); x0p = max(0, x0-pad)
            y1p = min(H, y1+pad); x1p = min(W, x1+pad)
            crop = rgb[y0p:y1p, x0p:x1p].copy()
            blob_crop = (cc[y0p:y1p, x0p:x1p] == p.label)
            crop[blob_crop] = 0.3 * crop[blob_crop] + 0.7 * np.array([1, 0, 1])
            axes[i].imshow(np.clip(crop, 0, 1))
            ax = p.axis_major_length / max(p.axis_minor_length, 0.5)
            axes[i].set_title(f"area={p.area} ecc={p.eccentricity:.2f} axR={ax:.1f}\n"
                              f"y={y0}-{y1} x={x0}-{x1}", fontsize=9)
            axes[i].axis("off")
        fig.suptitle(f"{s.slide_id} — V4 tear zooms", fontsize=12)
        fig.tight_layout()
        fig.savefig(OUT / f"{s.slide_id}_v4_tears.png", dpi=80)
        plt.close(fig)


def main():
    target_ids = ["2018_05_22_s1_2_p27-0021",
                  "2018_05_22_s1_5_p27-0018",
                  "2018_05_22_s2_5_p27-0026",
                  "2018_05_22_s3_5_p27-0006"]
    slides = manifest()
    for s in [s for s in slides if s.slide_id in target_ids]:
        explore(s)


if __name__ == "__main__":
    main()
