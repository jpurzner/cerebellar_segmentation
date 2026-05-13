"""Leave-one-out cross-validation: train 14 separate heads, each holding out one slide.

For each fold:
  - train head on 13 slides
  - predict on held-out slide
  - compute IoU
"""
from __future__ import annotations
import time
from pathlib import Path
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader

from dataset import N_CLASSES, load_slide, matlab_segments_to_labels
from train_head import LinearHead, PatchDataset, class_weights
from predict import cerebellum_mask_from_stack
from skimage.transform import resize
import tifffile

CACHE = Path(__file__).resolve().parent / "cache"
MODELS = Path(__file__).resolve().parent / "models"
PREDS = Path(__file__).resolve().parent / "preds"
INPUT_BASE = Path("/Users/jpurzner/Dropbox/images/edu_repeat/p27")
SEG_BASE = Path("/Users/jpurzner/Dropbox/imaging_analysis/cerebellar_segmentation/_test_run/20x")
OUT = Path(__file__).resolve().parent / "loo"
OUT.mkdir(exist_ok=True)


def train_fold(train_files, val_file, epochs=30, lr=1e-3, batch=256,
                samples_per_class=10000):
    train_ds = PatchDataset(train_files, samples_per_class=samples_per_class)
    val_ds = PatchDataset([val_file])
    weights = class_weights(train_ds.labels)
    device = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
    head = LinearHead().to(device)
    opt = torch.optim.Adam(head.parameters(), lr=lr)
    crit = nn.CrossEntropyLoss(weight=weights.to(device))
    train_loader = DataLoader(train_ds, batch_size=batch, shuffle=True, drop_last=True)
    val_loader   = DataLoader(val_ds, batch_size=batch)
    best_state = None; best_iou = -1
    for ep in range(epochs):
        head.train()
        for x, y in train_loader:
            x = x.to(device); y = y.to(device)
            loss = crit(head(x), y)
            opt.zero_grad(); loss.backward(); opt.step()
        head.eval()
        with torch.no_grad():
            cm = np.zeros((N_CLASSES, N_CLASSES), dtype=np.int64)
            for x, y in val_loader:
                x = x.to(device); y_np = y.numpy()
                pred = head(x).argmax(dim=1).cpu().numpy()
                for tt, pp in zip(y_np, pred):
                    cm[tt, pp] += 1
        ious = []
        for c in range(1, N_CLASSES):
            tp = cm[c, c]; fp = cm[:, c].sum() - tp; fn = cm[c, :].sum() - tp
            ious.append(tp / max(tp + fp + fn, 1))
        miou = float(np.mean(ious))
        if miou > best_iou:
            best_iou = miou
            best_state = {k: v.cpu().clone() for k, v in head.state_dict().items()}
    return best_state, best_iou


def evaluate_fold(state, val_slide):
    cache_path = CACHE / f"{val_slide}_features.npz"
    d = np.load(cache_path, allow_pickle=True)
    feats = d["features"].astype(np.float32)
    gy, gx = d["patch_grid"]
    H, W = int(d["H"]), int(d["W"])
    head = LinearHead()
    head.load_state_dict(state); head.eval()
    with torch.no_grad():
        pred = head(torch.from_numpy(feats)).argmax(dim=1).cpu().numpy().astype(np.uint8)
    label_grid = pred.reshape(gy, gx)
    label_full = resize(label_grid, (H, W), order=0, preserve_range=True,
                         anti_aliasing=False).astype(np.uint8)

    inp_path = INPUT_BASE / val_slide / f"{val_slide}_fused_crop.tif"
    stack = load_slide(inp_path, channel_order=("p27", "neun", "dapi"))
    mask = cerebellum_mask_from_stack(stack)
    label_full = label_full * mask.astype(np.uint8)

    seg_path = SEG_BASE / val_slide / f"{val_slide}_fused_crop_segments.tif"
    gt_rgb = tifffile.imread(seg_path)
    gt = matlab_segments_to_labels(gt_rgb)
    if label_full.shape != gt.shape:
        Hc = min(label_full.shape[0], gt.shape[0])
        Wc = min(label_full.shape[1], gt.shape[1])
        label_full = label_full[:Hc, :Wc]; gt = gt[:Hc, :Wc]

    ious = {}
    for c in range(1, N_CLASSES):
        p = (label_full == c); g = (gt == c)
        ious[c] = (p & g).sum() / max((p | g).sum(), 1)
    pe = (label_full == 2) | (label_full == 3)
    ge = (gt == 2) | (gt == 3)
    inter = (pe & ge).sum(); union = (pe | ge).sum()
    eP = inter / max(pe.sum(), 1); eR = inter / max(ge.sum(), 1)
    eIoU = inter / max(union, 1)
    return ious, eP, eR, eIoU


def main():
    cache_files = sorted(CACHE.glob("*_features.npz"))
    print(f"LOO CV across {len(cache_files)} slides")
    print(f"{'held_out':40s} {'iEGL':>6} {'oEGL':>6} {'IGL':>6} {'ML':>6} {'DWL':>6} "
          f"{'EGL_P':>6} {'EGL_R':>6} {'EGL_IoU':>7} {'time':>5}")
    rows = []
    for hold_out in cache_files:
        slide = hold_out.stem.replace("_features", "")
        train_files = [f for f in cache_files if f != hold_out]
        t0 = time.time()
        state, val_miou = train_fold(train_files, hold_out, epochs=20)
        elapsed = time.time() - t0
        ious, eP, eR, eIoU = evaluate_fold(state, slide)
        print(f"  {slide:40s} {ious[2]:6.3f} {ious[3]:6.3f} {ious[4]:6.3f} "
              f"{ious[5]:6.3f} {ious[6]:6.3f} {eP:6.3f} {eR:6.3f} {eIoU:7.3f} {elapsed:5.0f}s")
        rows.append({"held_out": slide, "iEGL": ious[2], "oEGL": ious[3],
                     "IGL": ious[4], "ML": ious[5], "DWL": ious[6],
                     "egl_P": eP, "egl_R": eR, "egl_IoU": eIoU})

    print("\n=== LOO summary ===")
    iou_list = [r["egl_IoU"] for r in rows]
    print(f"  EGL IoU: mean={np.mean(iou_list):.3f}  min={min(iou_list):.3f}  max={max(iou_list):.3f}")
    for k in ["iEGL","oEGL","IGL","ML","DWL"]:
        vals = [r[k] for r in rows]
        print(f"  {k:5s} IoU mean={np.mean(vals):.3f}")
    import csv
    with open(OUT / "loo_metrics.csv", "w") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for r in rows: writer.writerow(r)
    print(f"\nsaved {OUT / 'loo_metrics.csv'}")


if __name__ == "__main__":
    main()
