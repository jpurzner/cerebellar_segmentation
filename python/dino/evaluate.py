"""Evaluate DINO predictions against MATLAB GT (and optionally vs port v11)."""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import tifffile

from dataset import N_CLASSES, matlab_segments_to_labels

PREDS = Path(__file__).resolve().parent / "preds"
SEG_BASE = Path("/Users/jpurzner/Dropbox/imaging_analysis/cerebellar_segmentation/_test_run/20x")
PX_UM = 0.5119049

LAYER_NAMES = {1: "all_cereb", 2: "iEGL", 3: "oEGL", 4: "IGL",
                5: "ML", 6: "DWL", 7: "PCL"}

LABEL_COLORS = np.array([
    [0,   0,   0],   [0,   0, 128], [0, 128, 255], [0, 255, 255],
    [255, 255, 0], [128, 255, 128],[255, 128, 0], [255,   0, 0],
], dtype=np.uint8)


def per_class_iou(pred, gt, n=N_CLASSES):
    out = {}
    for c in range(1, n):
        p = (pred == c); g = (gt == c)
        inter = (p & g).sum(); union = (p | g).sum()
        out[c] = inter / max(union, 1)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--slide", required=True)
    args = ap.parse_args()

    pred_files = list(PREDS.glob(f"*{args.slide}*_pred.npz"))
    if not pred_files:
        print(f"no prediction for {args.slide}"); return
    d = np.load(pred_files[0])
    pred = d["label_full"]
    print(f"prediction shape {pred.shape}")

    seg_path = SEG_BASE / args.slide / f"{args.slide}_fused_crop_segments.tif"
    if not seg_path.exists():
        print(f"missing GT {seg_path}"); return
    seg_rgb = tifffile.imread(seg_path)
    gt = matlab_segments_to_labels(seg_rgb)

    if pred.shape != gt.shape:
        print(f"shape mismatch pred {pred.shape} vs gt {gt.shape}")
        # crop to common
        H = min(pred.shape[0], gt.shape[0])
        W = min(pred.shape[1], gt.shape[1])
        pred = pred[:H, :W]; gt = gt[:H, :W]

    ious = per_class_iou(pred, gt)
    print("\nPer-class IoU:")
    for c, iou in ious.items():
        print(f"  {c}({LAYER_NAMES.get(c,'?'):10s})  IoU={iou:.3f}")

    # EGL combined (iEGL + oEGL)
    pred_egl = (pred == 2) | (pred == 3)
    gt_egl   = (gt == 2)   | (gt == 3)
    inter = (pred_egl & gt_egl).sum()
    union = (pred_egl | gt_egl).sum()
    egl_iou = inter / max(union, 1)
    egl_p   = inter / max(pred_egl.sum(), 1)
    egl_r   = inter / max(gt_egl.sum(), 1)
    print(f"\nEGL combined: P={egl_p:.3f}  R={egl_r:.3f}  IoU={egl_iou:.3f}")

    # save side-by-side comparison
    pred_rgb = LABEL_COLORS[pred]
    gt_rgb = LABEL_COLORS[gt]
    fig, ax = plt.subplots(1, 2, figsize=(16, 8))
    ax[0].imshow(pred_rgb); ax[0].set_title(f"DINO prediction ({args.slide})"); ax[0].axis("off")
    ax[1].imshow(gt_rgb);   ax[1].set_title("MATLAB GT"); ax[1].axis("off")
    out_path = PREDS / f"{args.slide}_dino_vs_matlab.png"
    fig.tight_layout(); fig.savefig(out_path, dpi=120); plt.close(fig)
    print(f"saved {out_path}")


if __name__ == "__main__":
    main()
