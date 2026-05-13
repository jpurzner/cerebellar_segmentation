"""LOO CV + train-with-val on the 19-slide good-only manifest (8 20x + 11 10x).

Two modes:
  --mode train_val   : train on N-2 slides, val on 2 (s1_4 20x + s5_A 10x)
  --mode loo         : leave-one-out across all 19 slides
"""
from __future__ import annotations
import argparse, time, csv
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

CACHE = Path(__file__).resolve().parent / "cache_unified"
MODELS = Path(__file__).resolve().parent / "models"
PREDS = Path(__file__).resolve().parent / "preds_unified"
PREDS.mkdir(exist_ok=True)
OUT = Path(__file__).resolve().parent / "unified_out"
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


def evaluate_slide(state, slide):
    cache_path = CACHE / f"{slide.slide_id}_features.npz"
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

    ious = {}
    for c in range(1, N_CLASSES):
        p = (label_full == c); g = (gt == c)
        ious[c] = (p & g).sum() / max((p | g).sum(), 1)
    pe = (label_full == 2) | (label_full == 3)
    ge = (gt == 2) | (gt == 3)
    inter = (pe & ge).sum(); union = (pe | ge).sum()
    eP = inter / max(pe.sum(), 1); eR = inter / max(ge.sum(), 1)
    eIoU = inter / max(union, 1)
    return ious, eP, eR, eIoU, label_full, mask


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["train_val", "loo"], default="loo")
    ap.add_argument("--epochs", type=int, default=20)
    args = ap.parse_args()

    slides = manifest()
    cache_files = {s.slide_id: CACHE / f"{s.slide_id}_features.npz" for s in slides}
    missing = [sid for sid, p in cache_files.items() if not p.exists()]
    if missing:
        print(f"MISSING cache for: {missing}"); return

    if args.mode == "train_val":
        val_ids = ["2018_05_22_s1_4_p27-0020", "s5_A_p27"]
        train_files = [cache_files[s.slide_id] for s in slides if s.slide_id not in val_ids]
        print(f"\ntrain mode: {len(train_files)} train slides")
        # train one model on all train slides, eval on each val slide
        # use first val for early-stopping
        val_slide = next(s for s in slides if s.slide_id == val_ids[0])
        state, miou = train_fold(train_files, cache_files[val_slide.slide_id],
                                  epochs=args.epochs)
        torch.save({"state_dict": state, "best_miou": miou}, MODELS / "head_unified.pt")
        print(f"\nbest val patch mIoU: {miou:.3f}")
        print(f"\n{'slide':30s} {'mag':5s} {'iEGL':>6} {'oEGL':>6} {'IGL':>6} {'ML':>6} "
              f"{'DWL':>6} {'EGL_P':>6} {'EGL_R':>6} {'EGL_IoU':>7}")
        for sid in val_ids:
            slide = next(s for s in slides if s.slide_id == sid)
            ious, eP, eR, eIoU, _, _ = evaluate_slide(state, slide)
            print(f"  {sid:30s} {slide.mag:5s} {ious[2]:6.3f} {ious[3]:6.3f} "
                  f"{ious[4]:6.3f} {ious[5]:6.3f} {ious[6]:6.3f} "
                  f"{eP:6.3f} {eR:6.3f} {eIoU:7.3f}")
    else:
        print(f"LOO across {len(slides)} slides")
        print(f"{'held_out':30s} {'mag':5s} {'iEGL':>6} {'oEGL':>6} {'IGL':>6} {'ML':>6} "
              f"{'DWL':>6} {'EGL_P':>6} {'EGL_R':>6} {'EGL_IoU':>7} {'time':>5}")
        rows = []
        for hold_out in slides:
            train_files = [cache_files[s.slide_id] for s in slides if s.slide_id != hold_out.slide_id]
            t0 = time.time()
            state, _ = train_fold(train_files, cache_files[hold_out.slide_id],
                                   epochs=args.epochs)
            elapsed = time.time() - t0
            ious, eP, eR, eIoU, _, _ = evaluate_slide(state, hold_out)
            print(f"  {hold_out.slide_id:30s} {hold_out.mag:5s} {ious[2]:6.3f} {ious[3]:6.3f} "
                  f"{ious[4]:6.3f} {ious[5]:6.3f} {ious[6]:6.3f} "
                  f"{eP:6.3f} {eR:6.3f} {eIoU:7.3f} {elapsed:5.0f}s")
            rows.append({"held_out": hold_out.slide_id, "mag": hold_out.mag,
                         "iEGL": ious[2], "oEGL": ious[3], "IGL": ious[4],
                         "ML": ious[5], "DWL": ious[6], "egl_P": eP,
                         "egl_R": eR, "egl_IoU": eIoU})

        print("\n=== summary ===")
        ious_egl = [r["egl_IoU"] for r in rows]
        print(f"  EGL IoU: mean={np.mean(ious_egl):.3f}  min={min(ious_egl):.3f}  max={max(ious_egl):.3f}")
        for k in ["iEGL","oEGL","IGL","ML","DWL"]:
            vals = [r[k] for r in rows]
            print(f"  {k:5s} IoU mean={np.mean(vals):.3f}")
        # split by mag
        for mag in ["10x", "20x"]:
            mag_rows = [r for r in rows if r["mag"] == mag]
            if mag_rows:
                ious_mag = [r["egl_IoU"] for r in mag_rows]
                print(f"  {mag} EGL IoU: n={len(ious_mag)}  mean={np.mean(ious_mag):.3f}")

        with open(OUT / "loo_unified_metrics.csv", "w") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            for r in rows: writer.writerow(r)
        print(f"\nsaved {OUT / 'loo_unified_metrics.csv'}")


if __name__ == "__main__":
    main()
