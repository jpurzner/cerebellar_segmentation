"""Render side-by-side review images for all 19 slides using the gold-boosted
DINO predictions. Saves PNGs to /python/review_panels/ for quick browsing.
"""
from __future__ import annotations
import time
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import tifffile

from slide_manifest import manifest
from dataset import load_slide, matlab_segments_to_labels

ROOT = Path(__file__).resolve().parent
PRED_DIR = ROOT.parent / "predictions_for_correction"
LABELLED_DIR = ROOT.parent / "labelled"
OUT = ROOT.parent / "review_panels"   # /python/review_panels/
OUT.mkdir(exist_ok=True)

LABEL_COLORS = np.array([
    [0,   0,   0],   [0,   0, 128], [0, 128, 255], [0, 255, 255],
    [255, 255, 0], [128, 255, 128],[255, 128, 0], [255,   0, 0],
], dtype=np.uint8)


def to_unit(im):
    im = im.astype(np.float32); lo, hi = im.min(), im.max()
    return (im - lo) / max(hi - lo, 1e-9)


def main():
    slides = manifest()
    rows = []
    for slide in slides:
        pred_path = PRED_DIR / f"{slide.slide_id}_dino_pred.tif"
        if not pred_path.exists():
            print(f"  no pred for {slide.slide_id}, skip"); continue

        print(f"\n=== {slide.mag}  {slide.slide_id} ===")
        stack = load_slide(slide.input_path, channel_order=slide.channel_order)
        rgb = np.stack([to_unit(stack[0]), to_unit(stack[1]), to_unit(stack[2])], axis=-1)
        rgb /= max(rgb.max(), 1e-9)

        dino = tifffile.imread(pred_path).astype(np.uint8)
        gt_rgb = tifffile.imread(slide.gt_path)
        matlab = matlab_segments_to_labels(gt_rgb)

        # crop all to same shape
        H = min(dino.shape[0], matlab.shape[0], rgb.shape[0])
        W = min(dino.shape[1], matlab.shape[1], rgb.shape[1])
        rgb = rgb[:H, :W]
        dino = dino[:H, :W]
        matlab = matlab[:H, :W]

        # check if there's a hand-correction
        cands = sorted(
            list(LABELLED_DIR.glob(f"{slide.slide_id}_corrected.tif*")) +
            list(LABELLED_DIR.glob(f"{slide.slide_id}_labelled.tif*")) +
            list(LABELLED_DIR.glob(f"{slide.slide_id}_labelld.tif*")) +
            list(LABELLED_DIR.glob(f"{slide.slide_id}_labeld.tif*")),
            key=lambda p: p.stat().st_mtime, reverse=True)
        cands = [c for c in cands if c.suffix in (".tif", ".tiff")]
        if cands:
            human = tifffile.imread(cands[0]).astype(np.uint8)[:H, :W]
            n_panels = 4
            titles = ["RGB", "DINO (gold-boosted)", "MATLAB GT", f"Human ({cands[0].name})"]
            imgs = [rgb, LABEL_COLORS[dino], LABEL_COLORS[matlab], LABEL_COLORS[human]]
        else:
            n_panels = 3
            titles = ["RGB", "DINO (gold-boosted)", "MATLAB GT"]
            imgs = [rgb, LABEL_COLORS[dino], LABEL_COLORS[matlab]]

        # EGL IoU vs MATLAB
        de = (dino == 2) | (dino == 3); me = (matlab == 2) | (matlab == 3)
        iou_m = (de & me).sum() / max((de | me).sum(), 1)

        fig, ax = plt.subplots(1, n_panels, figsize=(n_panels * 7, 9))
        for k, (img, t) in enumerate(zip(imgs, titles)):
            ax[k].imshow(img); ax[k].axis("off")
            ax[k].set_title(t, fontsize=11)
        fig.suptitle(f"{slide.slide_id}  ({slide.mag})  EGL IoU vs MATLAB = {iou_m:.2f}",
                      fontsize=14)
        out_path = OUT / f"{slide.slide_id}_review.png"
        fig.tight_layout(); fig.savefig(out_path, dpi=110); plt.close(fig)
        print(f"  saved {out_path.name}  EGL IoU={iou_m:.3f}")
        rows.append({"slide": slide.slide_id, "mag": slide.mag, "egl_iou": iou_m,
                      "has_human_label": bool(cands)})

    # contact-sheet index
    print("\nbuilding index page...")
    sorted_rows = sorted(rows, key=lambda r: (r["mag"], r["slide"]))
    fig, axes = plt.subplots(5, 4, figsize=(20, 22))
    axes = axes.ravel()
    for i, r in enumerate(sorted_rows):
        if i >= 20: break
        img = plt.imread(OUT / f"{r['slide']}_review.png")
        axes[i].imshow(img); axes[i].axis("off")
        flag = "★" if r["has_human_label"] else ""
        axes[i].set_title(f"{r['mag']} {r['slide'][:25]}{flag}\nIoU={r['egl_iou']:.2f}",
                          fontsize=8)
    for i in range(len(sorted_rows), 20):
        axes[i].axis("off")
    fig.suptitle("DINO (gold-boosted) vs MATLAB GT — all 19 slides   ★=hand-corrected",
                  fontsize=14)
    fig.tight_layout(); fig.savefig(OUT / "_INDEX.png", dpi=80); plt.close(fig)
    print(f"\nindex: {OUT / '_INDEX.png'}")
    print(f"all panels: {OUT}")


if __name__ == "__main__":
    main()
