"""Render EGL line diagrams for each 20x slide.

For each slide produces a 3-panel PNG:
  RGB | RGB + EGL boundary (yellow) + skeleton (red) | label segmentation

Useful for:
  - Reviewing the EGL ribbon's continuity at a glance
  - Spotting cut-surface, base-of-cerebellum, deep-invagination issues
  - Comparing EGL trace before/after pipeline changes

Outputs:
  /python/egl_lines/<slide>_egl_lines.png
  /python/egl_lines/_INDEX.png
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import tifffile
from skimage.morphology import skeletonize
from skimage.segmentation import find_boundaries

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from slide_manifest import manifest

ROOT = Path(__file__).resolve().parent
SAND = ROOT.parent / "matlab_anomaly" / "sandbox"
OUT = ROOT.parent / "egl_lines"
OUT.mkdir(exist_ok=True)

LABEL_COLORS = np.array([
    [0,   0,   0],   [0,   0, 128], [0, 128, 255], [0, 255, 255],
    [255, 255, 0], [128, 255, 128],[255, 128, 0], [255,   0, 0],
    [255, 0, 255],   # 8 DCN
    [80,  80,  80],  # 9 anomaly
], dtype=np.uint8)


def to_unit(im):
    im = im.astype(np.float32); lo, hi = im.min(), im.max()
    return (im - lo) / max(hi - lo, 1e-9)


def find_latest_labels(slide_id):
    """Find the most recently produced labels file for a slide."""
    candidates = sorted(SAND.glob(f"{slide_id}*_labels.tif"),
                        key=lambda p: p.stat().st_mtime, reverse=True)
    return candidates[0] if candidates else None


def make_line_overlay(rgb, labels):
    """Build an overlay showing EGL boundary (yellow) + EGL skeleton (red)."""
    egl_mask = (labels == 2) | (labels == 3)
    if egl_mask.sum() == 0:
        return rgb.copy(), 0, 0

    boundary = find_boundaries(egl_mask, mode='outer')
    skeleton = skeletonize(egl_mask)

    overlay = rgb.copy()
    # Boundary in yellow
    overlay[boundary] = [1.0, 1.0, 0.0]
    # Skeleton in red (overrides boundary)
    overlay[skeleton] = [1.0, 0.0, 0.0]

    return overlay, int(boundary.sum()), int(skeleton.sum())


def main():
    slides = manifest()
    twenty = [s for s in slides if s.mag == "20x"]
    print(f"rendering EGL line diagrams for {len(twenty)} 20x slides")

    rendered = []
    for s in twenty:
        labels_path = find_latest_labels(s.slide_id)
        if labels_path is None:
            print(f"  {s.slide_id}: no labels — skip"); continue

        labels = tifffile.imread(labels_path).astype(np.uint8)
        try:
            stack = tifffile.imread(s.input_path)
        except Exception as e:
            print(f"  {s.slide_id}: read failed — {e}"); continue

        # RGB display: R=p27, G=NeuN, B=DAPI (20x channel order: p27, NeuN, DAPI)
        rgb = np.stack([to_unit(stack[0]), to_unit(stack[1]), to_unit(stack[2])], axis=-1)
        rgb /= max(rgb.max(), 1e-9)
        H = min(labels.shape[0], rgb.shape[0])
        W = min(labels.shape[1], rgb.shape[1])
        rgb = rgb[:H, :W]; labels = labels[:H, :W]

        overlay, n_boundary, n_skel = make_line_overlay(rgb, labels)

        # Save 3-panel PNG
        fig, axes = plt.subplots(1, 3, figsize=(21, 7))
        axes[0].imshow(rgb); axes[0].axis("off")
        axes[0].set_title(f"{s.slide_id}\nRGB", fontsize=10)
        axes[1].imshow(overlay); axes[1].axis("off")
        axes[1].set_title(
            f"EGL: yellow=boundary, red=skeleton\n"
            f"boundary={n_boundary} px, skeleton={n_skel} px", fontsize=10)
        axes[2].imshow(LABEL_COLORS[labels]); axes[2].axis("off")
        axes[2].set_title("Segmentation", fontsize=10)
        fig.tight_layout()
        out_path = OUT / f"{s.slide_id}_egl_lines.png"
        fig.savefig(out_path, dpi=85)
        plt.close(fig)
        print(f"  {s.slide_id}: boundary={n_boundary}, skel={n_skel} -> {out_path.name}")
        rendered.append(out_path)

    # Contact sheet
    if rendered:
        n = len(rendered)
        fig, axes = plt.subplots(n, 1, figsize=(28, n*7))
        if n == 1: axes = [axes]
        for i, p in enumerate(rendered):
            img = plt.imread(p)
            axes[i].imshow(img); axes[i].axis("off")
            axes[i].set_title(p.stem, fontsize=11)
        fig.suptitle("EGL line diagrams — 20x slides "
                     "(yellow=outer boundary, red=skeleton/centerline)",
                     fontsize=14)
        fig.tight_layout()
        fig.savefig(OUT / "_INDEX.png", dpi=70)
        plt.close(fig)
        print(f"\ncontact sheet: {OUT / '_INDEX.png'}")


if __name__ == "__main__":
    main()
