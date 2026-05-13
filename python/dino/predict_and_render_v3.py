"""Predict all 19 slides with DINOv3 head + render review panels.

Output:
  /python/predictions_v3/<slide>_dino_v3_pred.tif    (per-pixel labels)
  /python/review_panels_v3/<slide>_review.png         (RGB | DINOv3 | MATLAB | Human?)
"""
from __future__ import annotations
from pathlib import Path
import time
import numpy as np
import matplotlib.pyplot as plt
import torch
import tifffile
from skimage.transform import resize

from dataset import N_CLASSES, load_slide, matlab_segments_to_labels
from train_head import LinearHead
from predict import cerebellum_mask_from_stack
from slide_manifest import manifest

ROOT = Path(__file__).resolve().parent
CACHE = ROOT / "cache_v3_um0.5"
MODELS = ROOT / "models"
PREDS = ROOT.parent / "predictions_v3"
PANELS = ROOT.parent / "review_panels_v3"
LABELLED = ROOT.parent / "labelled"
PREDS.mkdir(exist_ok=True)
PANELS.mkdir(exist_ok=True)

LABEL_COLORS = np.array([
    [0,   0,   0],   [0,   0, 128], [0, 128, 255], [0, 255, 255],
    [255, 255, 0], [128, 255, 128],[255, 128, 0], [255,   0, 0],
], dtype=np.uint8)


def to_unit(im):
    im = im.astype(np.float32); lo, hi = im.min(), im.max()
    return (im - lo) / max(hi - lo, 1e-9)


def main():
    state = torch.load(MODELS / "head_v3_full_gold.pt", map_location="cpu",
                       weights_only=False)
    head = LinearHead(); head.load_state_dict(state["state_dict"]); head.eval()

    slides = manifest()
    rows = []
    for s in slides:
        cache = CACHE / f"{s.slide_id}_features.npz"
        if not cache.exists():
            print(f"  {s.slide_id}: no v3 cache, skip"); continue
        d = np.load(cache, allow_pickle=True)
        feats = d["features"].astype(np.float32)
        gy, gx = d["patch_grid"]
        H, W = int(d["H"]), int(d["W"])
        with torch.no_grad():
            pred = head(torch.from_numpy(feats)).argmax(dim=1).cpu().numpy().astype(np.uint8)
        label_grid = pred.reshape(gy, gx)
        label_full = resize(label_grid, (H, W), order=0, preserve_range=True,
                             anti_aliasing=False).astype(np.uint8)

        stack = load_slide(s.input_path, channel_order=s.channel_order)
        mask = cerebellum_mask_from_stack(stack)
        label_full = label_full * mask.astype(np.uint8)

        out_path = PREDS / f"{s.slide_id}_dino_v3_pred.tif"
        tifffile.imwrite(out_path, label_full, compression="zlib")

        # render review panel
        rgb = np.stack([to_unit(stack[0]), to_unit(stack[1]), to_unit(stack[2])], axis=-1)
        rgb /= max(rgb.max(), 1e-9)
        gt_rgb = tifffile.imread(s.gt_path)
        matlab = matlab_segments_to_labels(gt_rgb)
        Hc = min(label_full.shape[0], matlab.shape[0], rgb.shape[0])
        Wc = min(label_full.shape[1], matlab.shape[1], rgb.shape[1])
        rgb = rgb[:Hc, :Wc]; label_full = label_full[:Hc, :Wc]; matlab = matlab[:Hc, :Wc]

        cands = sorted(
            list(LABELLED.glob(f"{s.slide_id}_corrected.tif*")) +
            list(LABELLED.glob(f"{s.slide_id}_labelled.tif*")) +
            list(LABELLED.glob(f"{s.slide_id}_labelld.tif*")),
            key=lambda p: p.stat().st_mtime, reverse=True)
        cands = [c for c in cands if c.suffix in (".tif", ".tiff")]
        if cands:
            human = tifffile.imread(cands[0]).astype(np.uint8)[:Hc, :Wc]
            n_panels, titles = 4, ["RGB", "DINOv3 (gold-boosted)", "MATLAB GT",
                                     f"Human ({cands[0].name})"]
            imgs = [rgb, LABEL_COLORS[label_full], LABEL_COLORS[matlab], LABEL_COLORS[human]]
        else:
            n_panels, titles = 3, ["RGB", "DINOv3 (gold-boosted)", "MATLAB GT"]
            imgs = [rgb, LABEL_COLORS[label_full], LABEL_COLORS[matlab]]

        de = (label_full == 2) | (label_full == 3); me = (matlab == 2) | (matlab == 3)
        iou_m = (de & me).sum() / max((de | me).sum(), 1)

        fig, ax = plt.subplots(1, n_panels, figsize=(n_panels*7, 9))
        for k, (img, t) in enumerate(zip(imgs, titles)):
            ax[k].imshow(img); ax[k].axis("off"); ax[k].set_title(t, fontsize=11)
        fig.suptitle(f"{s.slide_id}  ({s.mag})  EGL IoU vs MATLAB = {iou_m:.2f}",
                      fontsize=14)
        out_panel = PANELS / f"{s.slide_id}_review.png"
        fig.tight_layout(); fig.savefig(out_panel, dpi=110); plt.close(fig)
        rows.append((s.slide_id, s.mag, iou_m, bool(cands)))
        print(f"  {s.mag} {s.slide_id:30s} IoU={iou_m:.3f} → {out_panel.name}")

    # contact sheet
    print("\nbuilding INDEX.png...")
    sorted_rows = sorted(rows, key=lambda r: (r[1], r[0]))
    fig, axes = plt.subplots(5, 4, figsize=(20, 22)); axes = axes.ravel()
    for i, (sid, mag, iou, has_h) in enumerate(sorted_rows):
        if i >= 20: break
        img = plt.imread(PANELS / f"{sid}_review.png")
        axes[i].imshow(img); axes[i].axis("off")
        flag = "★" if has_h else ""
        axes[i].set_title(f"{mag} {sid[:25]}{flag}\nIoU={iou:.2f}", fontsize=8)
    for i in range(len(sorted_rows), 20):
        axes[i].axis("off")
    fig.suptitle("DINOv3 ViT-S/16 (gold-boosted) vs MATLAB — all 19   ★=hand-corrected",
                  fontsize=14)
    fig.tight_layout(); fig.savefig(PANELS / "_INDEX.png", dpi=80); plt.close(fig)
    print(f"\nindex: {PANELS / '_INDEX.png'}")
    print(f"all panels: {PANELS}")


if __name__ == "__main__":
    main()
