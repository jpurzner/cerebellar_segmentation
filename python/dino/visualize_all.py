"""Predict on all 14 slides + generate 4-panel comparison plots
    [RGB | DINO pred | MATLAB GT | EGL diff]
"""
from __future__ import annotations
import time
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import torch
import tifffile
from skimage.transform import resize

from dataset import N_CLASSES, load_slide, matlab_segments_to_labels
from train_head import LinearHead
from predict import cerebellum_mask_from_stack

CACHE = Path(__file__).resolve().parent / "cache"
MODELS = Path(__file__).resolve().parent / "models"
PREDS = Path(__file__).resolve().parent / "preds"
INPUT_BASE = Path("/Users/jpurzner/Dropbox/images/edu_repeat/p27")
SEG_BASE = Path("/Users/jpurzner/Dropbox/imaging_analysis/cerebellar_segmentation/_test_run/20x")
PORT_OUT = Path(__file__).resolve().parent.parent / "figs" / "v_all_14"
OUT = Path(__file__).resolve().parent / "viz_all"
OUT.mkdir(exist_ok=True)
PX_UM = 0.5119049

LABEL_COLORS = np.array([
    [0,   0,   0],   [0,   0, 128], [0, 128, 255], [0, 255, 255],
    [255, 255, 0], [128, 255, 128],[255, 128, 0], [255,   0, 0],
], dtype=np.uint8)


def main():
    head_path = MODELS / "head_v2.pt"
    state = torch.load(head_path, map_location="cpu", weights_only=False)
    head = LinearHead()
    head.load_state_dict(state["state_dict"])
    head.eval()

    cache_files = sorted(CACHE.glob("*_features.npz"))
    print(f"{'slide':40s} {'iEGL':>6} {'oEGL':>6} {'IGL':>6} {'ML':>6} {'DWL':>6} "
          f"{'EGL_P':>6} {'EGL_R':>6} {'EGL_IoU':>7}")
    rows = []
    for cache_path in cache_files:
        slide = cache_path.stem.replace("_features", "")
        d = np.load(cache_path, allow_pickle=True)
        feats = d["features"].astype(np.float32)
        gy, gx = d["patch_grid"]
        H, W = int(d["H"]), int(d["W"])

        with torch.no_grad():
            pred = head(torch.from_numpy(feats)).argmax(dim=1).cpu().numpy().astype(np.uint8)
        label_grid = pred.reshape(gy, gx)
        label_full = resize(label_grid, (H, W), order=0, preserve_range=True,
                             anti_aliasing=False).astype(np.uint8)

        # mask with cereb
        inp_path = INPUT_BASE / slide / f"{slide}_fused_crop.tif"
        stack = load_slide(inp_path, channel_order=("p27", "neun", "dapi"))
        mask = cerebellum_mask_from_stack(stack)
        label_full = label_full * mask.astype(np.uint8)

        # GT
        seg_path = SEG_BASE / slide / f"{slide}_fused_crop_segments.tif"
        gt_rgb = tifffile.imread(seg_path)
        gt = matlab_segments_to_labels(gt_rgb)

        if label_full.shape != gt.shape:
            Hc = min(label_full.shape[0], gt.shape[0])
            Wc = min(label_full.shape[1], gt.shape[1])
            label_full = label_full[:Hc, :Wc]; gt = gt[:Hc, :Wc]; mask = mask[:Hc, :Wc]
            stack = stack[:, :Hc, :Wc]

        # per-class IoU
        ious = {}
        for c in range(1, N_CLASSES):
            p = (label_full == c); g = (gt == c)
            ious[c] = (p & g).sum() / max((p | g).sum(), 1)
        # EGL combined
        pe = (label_full == 2) | (label_full == 3)
        ge = (gt == 2) | (gt == 3)
        inter = (pe & ge).sum(); union = (pe | ge).sum()
        eP = inter / max(pe.sum(), 1)
        eR = inter / max(ge.sum(), 1)
        eIoU = inter / max(union, 1)
        rows.append({"slide": slide, "iEGL": ious[2], "oEGL": ious[3],
                     "IGL": ious[4], "ML": ious[5], "DWL": ious[6],
                     "egl_P": eP, "egl_R": eR, "egl_IoU": eIoU})
        print(f"  {slide:40s} {ious[2]:6.3f} {ious[3]:6.3f} {ious[4]:6.3f} "
              f"{ious[5]:6.3f} {ious[6]:6.3f} {eP:6.3f} {eR:6.3f} {eIoU:7.3f}")

        # plot
        rgb = stack.transpose(1, 2, 0)
        rgb = rgb / max(rgb.max(), 1e-9)
        rgb_dim = rgb * 0.45
        tp = pe & ge; fn = ge & ~pe; fp = pe & ~ge
        diff = np.copy(rgb_dim)
        diff[..., 0] = np.clip(diff[..., 0] + fp.astype(np.float32)*0.9 + fn.astype(np.float32)*0.9, 0, 1)
        diff[..., 1] = np.clip(diff[..., 1] + tp.astype(np.float32)*0.9 + fp.astype(np.float32)*0.7, 0, 1)
        dino_seg = LABEL_COLORS[label_full]

        fig, ax = plt.subplots(1, 4, figsize=(28, 9))
        ax[0].imshow(rgb); ax[0].set_title(f"{slide}\nRGB"); ax[0].axis("off")
        ax[1].imshow(dino_seg); ax[1].set_title("DINO ViT-S/14 + linear head"); ax[1].axis("off")
        ax[2].imshow(gt_rgb); ax[2].set_title("MATLAB segments (GT)"); ax[2].axis("off")
        ax[3].imshow(diff); ax[3].set_title(
            f"EGL diff: green=TP red=FN yellow=FP\nP={eP:.2f} R={eR:.2f} IoU={eIoU:.2f}")
        ax[3].axis("off")
        fig.tight_layout(); fig.savefig(OUT / f"{slide}_dino_compare.png", dpi=130)
        plt.close(fig)

    # summary
    print("\n=== summary ===")
    iou_list = [r["egl_IoU"] for r in rows]
    print(f"  EGL IoU: mean={np.mean(iou_list):.3f}  min={min(iou_list):.3f}  max={max(iou_list):.3f}")
    for c, name in [(2,"iEGL"),(3,"oEGL"),(4,"IGL"),(5,"ML"),(6,"DWL")]:
        m = np.mean([list(r.values())[1+([(2,"iEGL"),(3,"oEGL"),(4,"IGL"),(5,"ML"),(6,"DWL")].index((c,name)))] for r in rows])
        # simpler:
        key = name if name in rows[0] else None
    # actually compute from dicts
    for k in ["iEGL","oEGL","IGL","ML","DWL"]:
        vals = [r[k] for r in rows]
        print(f"  {k:5s} IoU mean={np.mean(vals):.3f}")

    # save metrics CSV
    import csv
    with open(OUT / "metrics_dino_vs_matlab.csv", "w") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for r in rows: writer.writerow(r)
    print(f"\nsaved metrics to {OUT / 'metrics_dino_vs_matlab.csv'}")


if __name__ == "__main__":
    main()
