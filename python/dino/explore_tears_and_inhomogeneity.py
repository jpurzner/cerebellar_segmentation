"""Diagnostic explorer for the user-identified anomalies:
  - s1_2 (gold)   : small tear (labeled in gold standard)
  - s1_5          : small tear
  - s2_5          : signal inhomogeneity (main issue; no tear)
  - s3_5          : small tear

For each slide, render a 3-row × 3-col diagnostic panel:
  Row 1: RGB | DAPI raw normalized | DAPI <0.05 mask (per pixel, no dilation)
  Row 2: DAPI smoothed sigma=10px | DAPI smoothed sigma=30px | low-pass background
  Row 3: gold (if exists) | DAPI dropout candidates by area | flat-field corrected DAPI

Specifically for tear detection:
  - Per-pixel low-DAPI mask + size filter at 200, 1000, 5000 px (find thin elongated)
  - shape stats per connected component (aspect, area, eccentricity)

For inhomogeneity (s2_5):
  - Show large-scale DAPI background (Gaussian sigma=200 px)
  - Show flat-field-corrected DAPI (image / background)
  - Compare row/col mean intensity across image

Output:
  /python/anomaly_explore/<slide>_diag.png    per-slide panel
  /python/anomaly_explore/<slide>_tear_components.txt   stats on dropout components
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage as ndi
import tifffile

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from slide_manifest import manifest

ROOT = Path(__file__).resolve().parent
LABELLED_DIR = ROOT.parent / "labelled"
OUT = ROOT.parent / "anomaly_explore"
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
    print(f"  shape: {stack.shape}, pixel: {s.pixel_size_um:.4f} um")

    # Channel order at 20x: page 1=p27, page 2=NeuN, page 3=DAPI
    p27_raw  = stack[0].astype(np.float32)
    neun_raw = stack[1].astype(np.float32)
    dapi_raw = stack[2].astype(np.float32)

    p27  = to_unit(p27_raw)
    neun = to_unit(neun_raw)
    dapi = to_unit(dapi_raw)

    rgb = np.stack([p27, neun, dapi], axis=-1)
    rgb = rgb / max(rgb.max(), 1e-9)

    # Rough tissue mask
    combined = (p27 + neun + dapi) / 3
    tissue = combined > 0.10

    # Per-pixel DAPI dropout
    drop_px = (dapi < 0.05) & tissue

    # DAPI smoothed at multiple scales
    dapi_s10 = ndi.gaussian_filter(dapi, sigma=10)
    dapi_s30 = ndi.gaussian_filter(dapi, sigma=30)

    # Large-scale background (Gaussian sigma=200 px ≈ 100 um) for flat-field
    bg = ndi.gaussian_filter(dapi, sigma=200)
    flatfield = dapi / np.maximum(bg, 0.01) * bg.mean()
    flatfield = np.clip(flatfield, 0, 1)

    # Connected components of per-pixel dropout — find tear-like
    drop_filt = ndi.binary_opening(drop_px, structure=np.ones((3, 3)))   # remove single-px noise
    drop_filt = ndi.binary_closing(drop_filt, structure=np.ones((5, 5)))  # close tiny gaps
    cc_labels, n_cc = ndi.label(drop_filt)
    print(f"  {n_cc} connected components in per-pixel dropout mask")

    # Per-component stats
    stats = []
    for cid in range(1, n_cc + 1):
        mask = cc_labels == cid
        area = mask.sum()
        if area < 100: continue
        ys, xs = np.where(mask)
        h_extent = ys.max() - ys.min() + 1
        w_extent = xs.max() - xs.min() + 1
        major = max(h_extent, w_extent)
        minor = min(h_extent, w_extent)
        aspect = major / max(minor, 1)
        # bounding box density
        bbox_area = h_extent * w_extent
        density = area / bbox_area
        stats.append({"id": cid, "area": area, "aspect": aspect,
                      "major": major, "minor": minor, "density": density,
                      "y": (ys.min(), ys.max()), "x": (xs.min(), xs.max())})

    # Tears: large area + high aspect ratio + low density (curved/linear)
    tear_cands = [s for s in stats if s["area"] > 500 and s["aspect"] > 3]
    bleach_cands = [s for s in stats if s["area"] > 2000 and s["aspect"] <= 3 and s["density"] > 0.5]

    # Save stats
    with open(OUT / f"{s.slide_id}_tear_components.txt", "w") as f:
        f.write(f"Connected components in per-pixel DAPI<0.05 mask\n")
        f.write(f"slide: {s.slide_id}\n")
        f.write(f"total CCs >=100 px: {len(stats)}\n\n")
        f.write(f"  {'id':>5} {'area':>8} {'aspect':>7} {'major':>6} {'minor':>6} {'density':>7}  bbox\n")
        for st in sorted(stats, key=lambda x: -x["area"])[:30]:
            f.write(f"  {st['id']:>5} {st['area']:>8} {st['aspect']:>7.2f} "
                    f"{st['major']:>6} {st['minor']:>6} {st['density']:>7.3f}  "
                    f"y={st['y']} x={st['x']}\n")
    print(f"  top component areas: {[s['area'] for s in sorted(stats, key=lambda x: -x['area'])[:10]]}")
    print(f"  tear-like (area>500, aspect>3): {len(tear_cands)}")
    print(f"  bleach-like (area>2000, aspect<3, density>0.5): {len(bleach_cands)}")

    # Build masks for visualization
    tear_mask = np.zeros_like(drop_filt, dtype=bool)
    for st in tear_cands:
        tear_mask |= (cc_labels == st["id"])
    bleach_mask = np.zeros_like(drop_filt, dtype=bool)
    for st in bleach_cands:
        bleach_mask |= (cc_labels == st["id"])

    # Render diagnostic panel
    gold_path = find_gold(s.slide_id)
    gold = tifffile.imread(gold_path).astype(np.uint8) if gold_path else None
    if gold is not None:
        Hg = min(gold.shape[0], H); Wg = min(gold.shape[1], W)
    else:
        Hg, Wg = H, W

    fig, axes = plt.subplots(3, 3, figsize=(24, 24))

    # Row 1: RGB | DAPI raw | per-pixel dropout
    axes[0,0].imshow(rgb); axes[0,0].set_title("RGB (R=p27 G=NeuN B=DAPI)", fontsize=11)
    axes[0,1].imshow(dapi, cmap="gray", vmin=0, vmax=0.5)
    axes[0,1].set_title("DAPI normalized (0-0.5 range)", fontsize=11)
    axes[0,2].imshow(drop_px, cmap="hot")
    axes[0,2].set_title(f"DAPI<0.05 raw mask ({100*drop_px.sum()/drop_px.size:.2f}% of img)", fontsize=11)

    # Row 2: smoothed DAPI at two scales | low-pass background
    axes[1,0].imshow(dapi_s10, cmap="gray", vmin=0, vmax=0.3)
    axes[1,0].set_title("DAPI Gaussian smooth sig=10 px", fontsize=11)
    axes[1,1].imshow(dapi_s30, cmap="gray", vmin=0, vmax=0.3)
    axes[1,1].set_title("DAPI Gaussian smooth sig=30 px", fontsize=11)
    axes[1,2].imshow(bg, cmap="gray", vmin=0, vmax=0.3)
    axes[1,2].set_title("DAPI background sig=200 px (= ~100 um)\nshows inhomogeneity", fontsize=11)

    # Row 3: gold | tear vs bleach overlay | flat-field DAPI
    if gold is not None:
        axes[2,0].imshow(LABEL_COLORS[gold[:Hg,:Wg]])
        axes[2,0].set_title(f"GOLD: {gold_path.name}", fontsize=10)
    else:
        axes[2,0].imshow(combined, cmap="gray")
        axes[2,0].set_title("(no gold) combined channels", fontsize=11)

    # Tear/bleach overlay on RGB
    overlay = rgb.copy()
    overlay[tear_mask] = (0.3 * overlay[tear_mask] + 0.7 * np.array([1, 0, 1]))
    overlay[bleach_mask] = (0.3 * overlay[bleach_mask] + 0.7 * np.array([1, 0.5, 0]))
    axes[2,1].imshow(np.clip(overlay, 0, 1))
    axes[2,1].set_title(f"tears (magenta, n={len(tear_cands)}) + bleach (orange, n={len(bleach_cands)})", fontsize=11)

    axes[2,2].imshow(flatfield, cmap="gray", vmin=0, vmax=0.5)
    axes[2,2].set_title("DAPI / large-scale background (FLATTENED)", fontsize=11)

    for ax in axes.ravel():
        ax.axis("off")

    fig.suptitle(s.slide_id, fontsize=15)
    fig.tight_layout()
    out = OUT / f"{s.slide_id}_diag.png"
    fig.savefig(out, dpi=75)
    plt.close(fig)
    print(f"  saved {out.name}")


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
