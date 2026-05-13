"""V2 tear detector: proper PCA eccentricity, no opening (preserves thin shapes).

For each of s1_2, s1_5, s2_5, s3_5:
  1. Compute per-pixel DAPI<0.05 inside tissue
  2. Use skimage regionprops for proper eccentricity / major/minor / orientation
  3. Tear = eccentricity > 0.95 (linear) AND area > 500 px
  4. Bleach = solid (eccentricity < 0.85) AND area > 2000 px
  5. Render: RGB | per-pixel dropout | TEARS only | BLEACH only | gold (if exists)
  6. Plus: zoom panels for the 5 largest detected tears per slide

Also for s2_5 specifically: render flat-field corrected version + the
illumination-corrected dropout mask, to demonstrate that inhomogeneity
causes the biggest "dropouts" and flat-fielding kills them.
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
OUT = ROOT.parent / "anomaly_explore_v2"
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


def explore_slide(s):
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

    # Per-pixel dropout — NO opening (preserves thin tears)
    drop = (dapi < 0.05) & tissue
    # ONLY do a small closing to merge near-adjacent pixels (1-px gaps)
    drop = ndi.binary_closing(drop, structure=np.ones((3, 3)))

    # Connected components with proper shape stats
    cc_lbl = measure.label(drop, connectivity=2)
    props = measure.regionprops(cc_lbl)

    # Filter by minimum area (> 200 px²) and compute classification
    tears, bleaches, big_blobs = [], [], []
    for p in props:
        if p.area < 200: continue
        ecc = p.eccentricity
        # major/minor axis lengths from second moments
        major = p.major_axis_length
        minor = max(p.minor_axis_length, 0.5)
        axis_ratio = major / minor
        # bbox density = area / bbox_area
        bbox_area = (p.bbox[2] - p.bbox[0]) * (p.bbox[3] - p.bbox[1])
        density = p.area / max(bbox_area, 1)
        info = {"label": p.label, "area": p.area, "ecc": ecc,
                "major": major, "minor": minor, "axis_ratio": axis_ratio,
                "density": density, "bbox": p.bbox, "centroid": p.centroid}
        if p.area > 500 and (ecc > 0.95 or axis_ratio > 4):
            tears.append(info)
        elif p.area > 2000 and ecc < 0.85 and density > 0.5:
            bleaches.append(info)
        if p.area > 1000:
            big_blobs.append(info)

    print(f"  total CCs >=200 px: {len([p for p in props if p.area>=200])}")
    print(f"  tears: {len(tears)} (largest: {max([t['area'] for t in tears], default=0)} px²)")
    print(f"  bleaches: {len(bleaches)} (largest: {max([b['area'] for b in bleaches], default=0)} px²)")
    print(f"  big blobs >1000 px: {len(big_blobs)}")

    # Build masks
    tear_mask = np.zeros_like(drop, dtype=bool)
    for t in tears: tear_mask |= (cc_lbl == t["label"])
    bleach_mask = np.zeros_like(drop, dtype=bool)
    for b in bleaches: bleach_mask |= (cc_lbl == b["label"])

    # Print top blobs to characterize what's there
    sorted_blobs = sorted(big_blobs, key=lambda x: -x["area"])[:12]
    with open(OUT / f"{s.slide_id}_blobs.txt", "w") as f:
        f.write(f"slide: {s.slide_id}\n")
        f.write(f"top connected components (per-pixel DAPI<0.05 inside tissue, >=1000 px²):\n\n")
        f.write(f"  {'area':>7} {'ecc':>6} {'axR':>6} {'major':>6} {'minor':>6} {'dens':>5}  bbox(y0,x0,y1,x1)  classified\n")
        for b in sorted_blobs:
            cls = "tear" if b in tears else ("bleach" if b in bleaches else "?")
            f.write(f"  {b['area']:>7} {b['ecc']:>6.3f} {b['axis_ratio']:>6.2f} "
                    f"{b['major']:>6.0f} {b['minor']:>6.0f} {b['density']:>5.2f}  "
                    f"{b['bbox']}  {cls}\n")

    # Render diagnostic panel
    gold_path = find_gold(s.slide_id)
    gold = tifffile.imread(gold_path).astype(np.uint8) if gold_path else None

    fig, axes = plt.subplots(2, 3, figsize=(24, 16))

    axes[0,0].imshow(rgb); axes[0,0].set_title("RGB", fontsize=11)
    axes[0,1].imshow(dapi, cmap="gray", vmin=0, vmax=0.5)
    axes[0,1].set_title("DAPI normalized", fontsize=11)
    axes[0,2].imshow(drop, cmap="hot")
    axes[0,2].set_title(f"DAPI<0.05 raw mask ({100*drop.sum()/drop.size:.2f}%)", fontsize=11)

    # Overlay tears + bleaches on RGB
    overlay = rgb.copy()
    overlay[tear_mask] = 0.2 * overlay[tear_mask] + 0.8 * np.array([1, 0, 1])
    overlay[bleach_mask] = 0.2 * overlay[bleach_mask] + 0.8 * np.array([1, 0.5, 0])
    axes[1,0].imshow(np.clip(overlay, 0, 1))
    axes[1,0].set_title(
        f"TEARS (magenta, n={len(tears)}) + BLEACH (orange, n={len(bleaches)})",
        fontsize=11)

    # Big blobs visualization (numbered)
    blob_vis = rgb.copy()
    for b in sorted_blobs[:8]:
        y0, x0, y1, x1 = b["bbox"]
        # Draw bbox outline
        for yy in [y0, y1-1]:
            blob_vis[yy, x0:x1] = [1, 1, 0]
        for xx in [x0, x1-1]:
            blob_vis[y0:y1, xx] = [1, 1, 0]
    axes[1,1].imshow(np.clip(blob_vis, 0, 1))
    axes[1,1].set_title("Top 8 big blobs (yellow bboxes)\nsee {0}_blobs.txt".format(s.slide_id), fontsize=10)

    if gold is not None:
        Hg = min(gold.shape[0], H); Wg = min(gold.shape[1], W)
        axes[1,2].imshow(LABEL_COLORS[gold[:Hg,:Wg]])
        axes[1,2].set_title(f"GOLD: {gold_path.name}", fontsize=10)
    else:
        # Show flat-field corrected DAPI
        bg = ndi.gaussian_filter(dapi, sigma=200)
        ff = dapi / np.maximum(bg, 0.01) * bg.mean()
        ff = np.clip(ff, 0, 1)
        axes[1,2].imshow(ff, cmap="gray", vmin=0, vmax=0.5)
        axes[1,2].set_title("Flat-field DAPI (no gold)", fontsize=10)

    for ax in axes.ravel(): ax.axis("off")

    fig.suptitle(s.slide_id, fontsize=14)
    fig.tight_layout()
    fig.savefig(OUT / f"{s.slide_id}_tearfind.png", dpi=80)
    plt.close(fig)

    # Render zoom panels for the largest blobs
    if sorted_blobs:
        n = min(6, len(sorted_blobs))
        fig, axes = plt.subplots(1, n, figsize=(5*n, 5))
        if n == 1: axes = [axes]
        for i, b in enumerate(sorted_blobs[:n]):
            y0, x0, y1, x1 = b["bbox"]
            pad = 50
            y0p = max(0, y0 - pad); x0p = max(0, x0 - pad)
            y1p = min(H, y1 + pad); x1p = min(W, x1 + pad)
            crop = rgb[y0p:y1p, x0p:x1p].copy()
            # Highlight the blob in magenta
            blob_crop = (cc_lbl[y0p:y1p, x0p:x1p] == b["label"])
            crop[blob_crop] = 0.3 * crop[blob_crop] + 0.7 * np.array([1, 0, 1])
            axes[i].imshow(np.clip(crop, 0, 1))
            cls = "tear" if b in tears else ("bleach" if b in bleaches else "?")
            axes[i].set_title(
                f"area={b['area']} ecc={b['ecc']:.2f} axR={b['axis_ratio']:.1f}\n"
                f"y={y0}-{y1} x={x0}-{x1}  → {cls}", fontsize=9)
            axes[i].axis("off")
        fig.suptitle(f"{s.slide_id} — top blobs zoom", fontsize=12)
        fig.tight_layout()
        fig.savefig(OUT / f"{s.slide_id}_zoomblobs.png", dpi=80)
        plt.close(fig)

    print(f"  saved {s.slide_id}_tearfind.png, _zoomblobs.png, _blobs.txt")


def main():
    target_ids = ["2018_05_22_s1_2_p27-0021",
                  "2018_05_22_s1_5_p27-0018",
                  "2018_05_22_s2_5_p27-0026",
                  "2018_05_22_s3_5_p27-0006"]
    slides = manifest()
    targets = [s for s in slides if s.slide_id in target_ids]
    for s in targets:
        explore_slide(s)


if __name__ == "__main__":
    main()
