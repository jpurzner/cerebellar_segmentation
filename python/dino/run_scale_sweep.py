"""LOO CV at multiple patch scales. Compare to find optimum."""
from __future__ import annotations
import time, csv, argparse
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

from slide_manifest import manifest

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "scale_sweep"
OUT.mkdir(exist_ok=True)


def train_fold(train_files, val_file, epochs=15, lr=1e-3, batch=256,
                samples_per_class=10000):
    train_ds = PatchDataset(train_files, samples_per_class=samples_per_class)
    val_ds = PatchDataset([val_file])
    weights = class_weights(train_ds.labels)
    device = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
    head = LinearHead().to(device)
    opt = torch.optim.Adam(head.parameters(), lr=lr)
    crit = nn.CrossEntropyLoss(weight=weights.to(device))
    train_loader = DataLoader(train_ds, batch_size=batch, shuffle=True, drop_last=True)
    val_loader = DataLoader(val_ds, batch_size=batch)
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
                for tt, pp in zip(y_np, pred): cm[tt, pp] += 1
        ious = []
        for c in range(1, N_CLASSES):
            tp = cm[c, c]; fp = cm[:, c].sum() - tp; fn = cm[c, :].sum() - tp
            ious.append(tp / max(tp + fp + fn, 1))
        miou = float(np.mean(ious))
        if miou > best_iou:
            best_iou = miou
            best_state = {k: v.cpu().clone() for k, v in head.state_dict().items()}
    return best_state, best_iou


def evaluate_slide(state, slide, cache_dir):
    cache_path = cache_dir / f"{slide.slide_id}_features.npz"
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

    stack = load_slide(slide.input_path, channel_order=slide.channel_order)
    mask = cerebellum_mask_from_stack(stack)
    label_full = label_full * mask.astype(np.uint8)

    gt_rgb = tifffile.imread(slide.gt_path)
    gt = matlab_segments_to_labels(gt_rgb)
    if label_full.shape != gt.shape:
        Hc = min(label_full.shape[0], gt.shape[0])
        Wc = min(label_full.shape[1], gt.shape[1])
        label_full = label_full[:Hc, :Wc]; gt = gt[:Hc, :Wc]

    pe = (label_full == 2) | (label_full == 3)
    ge = (gt == 2) | (gt == 3)
    inter = (pe & ge).sum(); union = (pe | ge).sum()
    eP = inter / max(pe.sum(), 1); eR = inter / max(ge.sum(), 1)
    eIoU = inter / max(union, 1)
    return eP, eR, eIoU


def main():
    targets = [0.35, 0.5]   # μm/px → patch sizes 9.8, 19.6, 29.4, 42 μm
    slides = manifest()

    summary_rows = []
    for target in targets:
        cache_dir = ROOT / f"cache_um{target:g}"
        cache_files = {s.slide_id: cache_dir / f"{s.slide_id}_features.npz" for s in slides}
        missing = [sid for sid, p in cache_files.items() if not p.exists()]
        if missing:
            print(f"target {target}: MISSING cache, skipping"); continue
        print(f"\n========================================")
        print(f"target {target} um/px (patch ~{14*target:.1f} um)")
        print(f"========================================")
        rows = []
        for hold_out in slides:
            train_files = [cache_files[s.slide_id] for s in slides
                           if s.slide_id != hold_out.slide_id]
            t0 = time.time()
            state, _ = train_fold(train_files, cache_files[hold_out.slide_id])
            elapsed = time.time() - t0
            eP, eR, eIoU = evaluate_slide(state, hold_out, cache_dir)
            print(f"  {hold_out.mag} {hold_out.slide_id:30s} "
                  f"P={eP:.3f} R={eR:.3f} IoU={eIoU:.3f}  ({elapsed:.0f}s)")
            rows.append({"slide": hold_out.slide_id, "mag": hold_out.mag,
                          "egl_P": eP, "egl_R": eR, "egl_IoU": eIoU})
        ious_egl = [r["egl_IoU"] for r in rows]
        # exclude s5_E outlier for cleaner mean
        rows_clean = [r for r in rows if r["slide"] != "s5_E_p27"]
        ious_clean = [r["egl_IoU"] for r in rows_clean]
        mean_all = np.mean(ious_egl); mean_clean = np.mean(ious_clean)
        ious_10x = [r["egl_IoU"] for r in rows_clean if r["mag"] == "10x"]
        ious_20x = [r["egl_IoU"] for r in rows_clean if r["mag"] == "20x"]
        print(f"  EGL IoU all: {mean_all:.3f}  ex_s5E: {mean_clean:.3f}  "
              f"10x_clean: {np.mean(ious_10x):.3f}  20x: {np.mean(ious_20x):.3f}")
        summary_rows.append({
            "target_um_per_px": target, "patch_um": 14*target,
            "mean_all": mean_all, "mean_ex_s5E": mean_clean,
            "mean_10x_ex_s5E": np.mean(ious_10x), "mean_20x": np.mean(ious_20x),
        })
        with open(OUT / f"loo_um{target:g}.csv", "w") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            for r in rows: writer.writerow(r)

    print("\n=== sweep summary ===")
    print(f"{'target':>8} {'patch_um':>10} {'mean_all':>10} {'mean_ex_s5E':>12} "
          f"{'mean_10x':>10} {'mean_20x':>10}")
    for r in summary_rows:
        print(f"  {r['target_um_per_px']:>6.1f}  {r['patch_um']:>8.1f} "
              f"{r['mean_all']:>10.3f} {r['mean_ex_s5E']:>12.3f} "
              f"{r['mean_10x_ex_s5E']:>10.3f} {r['mean_20x']:>10.3f}")
    with open(OUT / "summary.csv", "w") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
        writer.writeheader()
        for r in summary_rows: writer.writerow(r)


if __name__ == "__main__":
    main()
