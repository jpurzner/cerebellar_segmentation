"""V3 anomaly detector. Two new tools:

1. FLAT-FIELD CORRECTION before any thresholding.
   dapi_ff = dapi / smooth(dapi, large_sigma) * mean
   Removes large-scale illumination inhomogeneity. Specifically targets s2_5
   where a huge swath of low DAPI from illumination problems was being
   detected as anomaly.

2. SMALL Gaussian smooth (sigma=5 px = 2.5 um) of the corrected DAPI before
   thresholding, so ML inter-cell gaps bridge through but real tears (>5 um
   wide gaps in tissue) survive.

Then: per-pixel threshold + connected components + shape stats.
TEAR: 500 < area < 20000 px AND (ecc > 0.93 OR axR > 3.5)
BLEACH: 1000 < area < 30000 px AND ecc < 0.85 AND extent > 0.5

Critical change vs V2: we don't try to use a giant local-max window. Instead
we smooth gently then threshold. This means tears down to ~5-10 um wide
should survive (Gaussian sigma=5 means signal from cells 10+ um away gets
diluted enough to not rescue the tear pixels).

Output: same diag panels + zoom-on-tears for each of s1_2, s1_5, s2_5, s3_5.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage as ndi
from skimage import measure
import tifffile

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from slide_manifest import manifest

ROOT = Path(__file__).resolve().parent
LABELLED_DIR = ROOT.parent / "labelled"
OUT = ROOT.parent / "anomaly_explore_v3"
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


def explore(s):
    print(f"\n=== {s.slide_id} ===")
    stack = tifffile.imread(s.input_path)
    H, W = stack.shape[1], stack.shape[2]

    p27  = to_unit(stack[0].astype(np.float32))
    neun = to_unit(stack[1].astype(np.float32))
    dapi = to_unit(stack[2].astype(np.float32))
    rgb = np.stack([p27, neun, dapi], axis=-1)
    rgb = rgb / max(rgb.max(), 1e-9)

    combined = (p27 + neun + dapi) / 3
    tissue = combined > 0.10

    # --- FLAT-FIELD CORRECTION ---
    # Estimate background (slow illumination variation) at sigma=200 px = ~100 um.
    # That's ~2x typical layer thickness so the bg follows illumination patterns
    # but not individual layers.
    # BUT: only compute bg from inside tissue (otherwise edges drag down the bg)
    bg = ndi.gaussian_filter(dapi * tissue.astype(np.float32), sigma=200)
    mask_blur = ndi.gaussian_filter(tissue.astype(np.float32), sigma=200)
    bg_local = bg / np.maximum(mask_blur, 0.05)   # tissue-aware mean
    # Mean of bg_local inside tissue (target intensity)
    bg_mean = bg_local[tissue].mean()
    # Apply correction
    dapi_ff = dapi / np.maximum(bg_local, 0.02) * bg_mean
    dapi_ff = np.clip(dapi_ff, 0, 1)

    # --- SMALL SMOOTHING for ML protection ---
    dapi_ff_s = ndi.gaussian_filter(dapi_ff, sigma=5)
    # Threshold: bottom-end signal in flat-field-corrected smoothed DAPI
    thr = 0.04
    drop = (dapi_ff_s < thr) & tissue

    # Tiny morphological closing to bridge 1-2 px noise gaps
    drop = ndi.binary_closing(drop, structure=np.ones((3, 3)))

    # Connected components + shape stats
    cc_lbl = measure.label(drop, connectivity=2)
    props = measure.regionprops(cc_lbl)

    tears, bleaches, big_blobs = [], [], []
    for p in props:
        if p.area < 200: continue
        major = p.axis_major_length
        minor = max(p.axis_minor_length, 0.5)
        ax = major / minor
        info = {"label": p.label, "area": p.area, "ecc": p.eccentricity,
                "ax": ax, "extent": p.extent, "euler": p.euler_number,
                "bbox": p.bbox}
        if 500 <= p.area <= 20000 and (p.eccentricity > 0.93 or ax > 3.5):
            tears.append(info)
        elif 1000 <= p.area <= 30000 and p.eccentricity < 0.85 and p.extent > 0.5:
            bleaches.append(info)
        if p.area >= 1000:
            big_blobs.append(info)

    print(f"  total CCs >=200: {len([p for p in props if p.area>=200])}")
    print(f"  TEARS: {len(tears)} (areas: {sorted([t['area'] for t in tears], reverse=True)[:5]})")
    print(f"  BLEACH: {len(bleaches)} (areas: {sorted([b['area'] for b in bleaches], reverse=True)[:5]})")
    if big_blobs:
        max_area = max(b['area'] for b in big_blobs)
        print(f"  largest blob: {max_area} px (ecc={max([b for b in big_blobs if b['area']==max_area], key=lambda x:0)['ecc']:.2f})")

    tear_mask = np.zeros_like(drop, dtype=bool)
    for t in tears: tear_mask |= (cc_lbl == t["label"])
    bleach_mask = np.zeros_like(drop, dtype=bool)
    for b in bleaches: bleach_mask |= (cc_lbl == b["label"])

    with open(OUT / f"{s.slide_id}_v3.txt", "w") as f:
        f.write(f"slide: {s.slide_id}\n")
        f.write(f"V3 detector: flat-field DAPI + Gaussian sig=5 + threshold 0.04\n\n")
        f.write(f"top blobs >=1000 px (sorted by area):\n")
        f.write(f"  {'area':>8} {'ecc':>6} {'axR':>6} {'extent':>7} {'euler':>7}  bbox  classified\n")
        for b in sorted(big_blobs, key=lambda x: -x["area"])[:15]:
            cls = "TEAR" if b in tears else ("BLEACH" if b in bleaches else "?")
            f.write(f"  {b['area']:>8} {b['ecc']:>6.3f} {b['ax']:>6.2f} "
                    f"{b['extent']:>7.3f} {b['euler']:>7d}  {b['bbox']}  {cls}\n")

    # Render: 2x4 panel
    gold_path = find_gold(s.slide_id)
    gold = tifffile.imread(gold_path).astype(np.uint8) if gold_path else None
    fig, axes = plt.subplots(2, 4, figsize=(28, 14))

    axes[0,0].imshow(rgb); axes[0,0].set_title("RGB", fontsize=11)
    axes[0,1].imshow(dapi, cmap="gray", vmin=0, vmax=0.5)
    axes[0,1].set_title("DAPI raw normalized", fontsize=11)
    axes[0,2].imshow(bg_local, cmap="gray", vmin=0, vmax=0.3)
    axes[0,2].set_title("DAPI background (tissue-aware sig=200 px)\n→ inhomogeneity map", fontsize=10)
    axes[0,3].imshow(dapi_ff, cmap="gray", vmin=0, vmax=0.5)
    axes[0,3].set_title("DAPI flat-field corrected", fontsize=11)

    axes[1,0].imshow(dapi_ff_s, cmap="gray", vmin=0, vmax=0.3)
    axes[1,0].set_title("DAPI ff + Gaussian sig=5 px", fontsize=11)
    axes[1,1].imshow(drop, cmap="hot")
    axes[1,1].set_title(f"V3 dropout mask (thr={thr}, {100*drop.sum()/drop.size:.3f}%)", fontsize=11)
    overlay = rgb.copy()
    overlay[tear_mask] = 0.2 * overlay[tear_mask] + 0.8 * np.array([1, 0, 1])
    overlay[bleach_mask] = 0.2 * overlay[bleach_mask] + 0.8 * np.array([1, 0.5, 0])
    axes[1,2].imshow(np.clip(overlay, 0, 1))
    axes[1,2].set_title(f"V3 classified: TEARS magenta (n={len(tears)}) BLEACH orange (n={len(bleaches)})", fontsize=10)

    if gold is not None:
        Hg = min(gold.shape[0], H); Wg = min(gold.shape[1], W)
        axes[1,3].imshow(LABEL_COLORS[gold[:Hg,:Wg]])
        axes[1,3].set_title(f"GOLD (does it show the tear too?)", fontsize=10)
    else:
        axes[1,3].imshow(combined, cmap="gray"); axes[1,3].set_title("(no gold)")

    for ax in axes.ravel(): ax.axis("off")
    fig.suptitle(s.slide_id, fontsize=14)
    fig.tight_layout()
    fig.savefig(OUT / f"{s.slide_id}_v3.png", dpi=80)
    plt.close(fig)

    # Zoom on the TEAR detections
    if tears:
        n = min(6, len(tears))
        fig, axes = plt.subplots(1, n, figsize=(5*n, 5))
        if n == 1: axes = [axes]
        for i, t in enumerate(sorted(tears, key=lambda x: -x["area"])[:n]):
            y0, x0, y1, x1 = t["bbox"]
            pad = 80
            y0p = max(0, y0-pad); x0p = max(0, x0-pad)
            y1p = min(H, y1+pad); x1p = min(W, x1+pad)
            crop = rgb[y0p:y1p, x0p:x1p].copy()
            blob_crop = (cc_lbl[y0p:y1p, x0p:x1p] == t["label"])
            crop[blob_crop] = 0.3 * crop[blob_crop] + 0.7 * np.array([1, 0, 1])
            axes[i].imshow(np.clip(crop, 0, 1))
            axes[i].set_title(f"TEAR area={t['area']} ecc={t['ecc']:.2f} axR={t['ax']:.1f}\n"
                              f"y={y0}-{y1} x={x0}-{x1}", fontsize=9)
            axes[i].axis("off")
        fig.suptitle(f"{s.slide_id} — V3 tear zooms", fontsize=12)
        fig.tight_layout()
        fig.savefig(OUT / f"{s.slide_id}_v3_tears.png", dpi=80)
        plt.close(fig)

    print(f"  saved {s.slide_id}_v3.png")


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
