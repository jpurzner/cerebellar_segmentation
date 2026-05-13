"""LOO CV with gold-slide upweighting.

Compares: baseline (MATLAB labels for all) vs gold-augmented (corrected labels
for s1_2, with extra samples from it) on LOO across 19 slides.
"""
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
OUT = ROOT / "loo_gold"
OUT.mkdir(exist_ok=True)


def train_fold(train_files, val_file, gold_files=None,
                epochs=20, lr=1e-3, batch=256, samples_per_class=10000,
                gold_boost=4):
    """Train head; if gold_files provided, also include extra samples from those slides."""
    train_ds = PatchDataset(train_files, samples_per_class=samples_per_class)
    if gold_files:
        gold_ds = PatchDataset(gold_files, samples_per_class=samples_per_class * gold_boost)
        # combine
        from torch.utils.data import ConcatDataset
        train_ds_combined = ConcatDataset([train_ds, gold_ds])
        # for class_weights, combine labels arrays
        all_labels = np.concatenate([train_ds.labels, gold_ds.labels])
        weights = class_weights(all_labels)
    else:
        train_ds_combined = train_ds
        weights = class_weights(train_ds.labels)
    val_ds = PatchDataset([val_file])
    device = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
    head = LinearHead().to(device)
    opt = torch.optim.Adam(head.parameters(), lr=lr)
    crit = nn.CrossEntropyLoss(weight=weights.to(device))
    train_loader = DataLoader(train_ds_combined, batch_size=batch, shuffle=True, drop_last=True)
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


def evaluate_slide(state, slide, cache_dir, eval_against="matlab"):
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

    # Choose GT
    LABELLED = ROOT.parent / "labelled"
    gold_path = None
    for pat in [f"{slide.slide_id}_corrected.tif*", f"{slide.slide_id}_labelled.tif*",
                f"{slide.slide_id}_labelld.tif*", f"{slide.slide_id}_labeld.tif*"]:
        cands = sorted(LABELLED.glob(pat), key=lambda p: p.stat().st_mtime, reverse=True)
        cands = [c for c in cands if c.suffix in (".tif", ".tiff")]
        if cands: gold_path = cands[0]; break

    if eval_against == "gold" and gold_path is not None:
        gt = tifffile.imread(gold_path).astype(np.uint8)
        gt_source = gold_path.name
    else:
        gt_rgb = tifffile.imread(slide.gt_path)
        gt = matlab_segments_to_labels(gt_rgb)
        gt_source = "matlab"

    if label_full.shape != gt.shape:
        Hc = min(label_full.shape[0], gt.shape[0])
        Wc = min(label_full.shape[1], gt.shape[1])
        label_full = label_full[:Hc, :Wc]; gt = gt[:Hc, :Wc]

    pe = (label_full == 2) | (label_full == 3)
    ge = (gt == 2) | (gt == 3)
    inter = (pe & ge).sum(); union = (pe | ge).sum()
    eP = inter / max(pe.sum(), 1); eR = inter / max(ge.sum(), 1)
    eIoU = inter / max(union, 1)
    return eP, eR, eIoU, gt_source, label_full


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", type=Path,
                    default=Path(__file__).resolve().parent / "cache_um0.5")
    ap.add_argument("--gold_boost", type=int, default=4,
                    help="multiplier for gold-slide samples")
    ap.add_argument("--epochs", type=int, default=20)
    args = ap.parse_args()

    slides = manifest()
    cache_files = {s.slide_id: args.cache / f"{s.slide_id}_features.npz" for s in slides}

    # Find which slides have corrected labels (now in updated cache)
    LABELLED = ROOT.parent / "labelled"
    gold_slides = []
    for s in slides:
        for pat in [f"{s.slide_id}_corrected.tif*", f"{s.slide_id}_labelled.tif*",
                    f"{s.slide_id}_labelld.tif*"]:
            cands = [c for c in LABELLED.glob(pat) if c.suffix in (".tif", ".tiff")]
            if cands:
                gold_slides.append(s.slide_id); break
    print(f"gold slides ({len(gold_slides)}): {gold_slides}")
    print(f"gold sample boost factor: {args.gold_boost}")

    print(f"\n{'held_out':30s} {'mag':5s} {'P':>6} {'R':>6} {'IoU(M)':>7} {'IoU(G)':>7} "
          f"{'gt_used':>10} {'time':>5}")
    rows = []
    for hold_out in slides:
        train_files = []
        gold_files = []
        for s in slides:
            if s.slide_id == hold_out.slide_id: continue
            f = cache_files[s.slide_id]
            if not f.exists():
                print(f"  missing cache for {s.slide_id}"); continue
            if s.slide_id in gold_slides:
                gold_files.append(f)
            else:
                train_files.append(f)
        t0 = time.time()
        state, val_miou = train_fold(train_files, cache_files[hold_out.slide_id],
                                       gold_files=gold_files, gold_boost=args.gold_boost,
                                       epochs=args.epochs)
        elapsed = time.time() - t0
        eP_m, eR_m, eIoU_m, _, _ = evaluate_slide(state, hold_out, args.cache, "matlab")
        eP_g, eR_g, eIoU_g, gt_src, _ = evaluate_slide(state, hold_out, args.cache, "gold")
        gt_label = "gold" if gt_src != "matlab" else "matlab"
        print(f"  {hold_out.slide_id:30s} {hold_out.mag:5s} "
              f"{eP_m:6.3f} {eR_m:6.3f} {eIoU_m:7.3f} {eIoU_g:7.3f} {gt_label:>10} {elapsed:5.0f}s")
        rows.append({"slide": hold_out.slide_id, "mag": hold_out.mag,
                      "matlab_P": eP_m, "matlab_R": eR_m, "matlab_IoU": eIoU_m,
                      "gold_IoU": eIoU_g, "gt_used": gt_label})

    print("\n=== summary ===")
    iou_m = [r["matlab_IoU"] for r in rows]
    iou_g_only = [r["gold_IoU"] for r in rows if r["gt_used"] == "gold"]
    rows_clean = [r for r in rows if r["slide"] != "s5_E_p27"]
    iou_m_clean = [r["matlab_IoU"] for r in rows_clean]
    iou_10x = [r["matlab_IoU"] for r in rows_clean if r["mag"] == "10x"]
    iou_20x = [r["matlab_IoU"] for r in rows_clean if r["mag"] == "20x"]
    print(f"  EGL IoU vs MATLAB (all 19): mean={np.mean(iou_m):.3f}  "
          f"ex_s5E={np.mean(iou_m_clean):.3f}")
    print(f"  EGL IoU vs MATLAB (10x clean): {np.mean(iou_10x):.3f}")
    print(f"  EGL IoU vs MATLAB (20x):       {np.mean(iou_20x):.3f}")
    if iou_g_only:
        print(f"  EGL IoU vs GOLD ({len(iou_g_only)} slides): {np.mean(iou_g_only):.3f}")
    with open(OUT / f"loo_gold_boost{args.gold_boost}.csv", "w") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for r in rows: writer.writerow(r)
    print(f"\nsaved {OUT / f'loo_gold_boost{args.gold_boost}.csv'}")


if __name__ == "__main__":
    main()
