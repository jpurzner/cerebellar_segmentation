"""Train a small linear head on cached DINOv2 features.

Loads cache/{slide}_features.npz for each slide, splits train/val,
trains a 1-layer linear classifier (384 → 8) with class-weighted CE.

Outputs: head.pt
"""
from __future__ import annotations
import argparse
from pathlib import Path
import time
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader

from dataset import N_CLASSES

CACHE = Path(__file__).resolve().parent / "cache"
OUT = Path(__file__).resolve().parent / "models"
OUT.mkdir(exist_ok=True)


class PatchDataset(Dataset):
    def __init__(self, slide_files, samples_per_class=None):
        feats_list = []
        labels_list = []
        for f in slide_files:
            d = np.load(f, allow_pickle=True)
            if not bool(d["has_labels"]):
                continue
            feats = d["features"].astype(np.float32)   # (N, 384)
            lbls  = d["labels"].astype(np.int64).ravel()
            assert len(feats) == len(lbls), \
                f"{f.name}: {len(feats)} feats vs {len(lbls)} lbls"
            feats_list.append(feats)
            labels_list.append(lbls)
        self.feats  = np.concatenate(feats_list, axis=0)
        self.labels = np.concatenate(labels_list, axis=0)
        # drop background
        keep = self.labels != 0
        self.feats = self.feats[keep]
        self.labels = self.labels[keep]
        # per-class balanced subsample
        if samples_per_class is not None:
            keep_idx = []
            for c in range(1, N_CLASSES):
                idx = np.where(self.labels == c)[0]
                if len(idx) > samples_per_class:
                    idx = np.random.choice(idx, samples_per_class, replace=False)
                keep_idx.append(idx)
            keep_idx = np.concatenate(keep_idx)
            np.random.shuffle(keep_idx)
            self.feats = self.feats[keep_idx]
            self.labels = self.labels[keep_idx]
        print(f"  dataset: {len(self.labels)} patches  class counts: ",
              dict(zip(*np.unique(self.labels, return_counts=True))))

    def __len__(self): return len(self.labels)
    def __getitem__(self, i):
        return torch.from_numpy(self.feats[i]), int(self.labels[i])


def class_weights(labels, n=N_CLASSES):
    """Inverse-frequency class weights for classes that actually appear in
    the training set. Class 0 (background) gets weight 0 since it's filtered
    out of training data — without this, weight[0] would dominate the
    normalisation and crush all other classes' weights to ~0.
    """
    counts = np.bincount(labels, minlength=n).astype(np.float32)
    inv = np.zeros(n, dtype=np.float32)
    nonzero = counts > 0
    inv[nonzero] = 1.0 / counts[nonzero]
    n_active = nonzero.sum()
    if inv.sum() > 0:
        w = inv / inv.sum() * n_active   # mean weight = 1 over active classes
    else:
        w = np.ones(n, dtype=np.float32)
    return torch.from_numpy(w.astype(np.float32))


class LinearHead(nn.Module):
    def __init__(self, in_dim=384, hidden=64, n_classes=N_CLASSES):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(in_dim, hidden),
            nn.ReLU(inplace=True),
            nn.Dropout(0.2),
            nn.Linear(hidden, n_classes),
        )

    def forward(self, x): return self.net(x)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--val", nargs="+", required=True,
                    help="Slides to hold out for validation (substring match)")
    ap.add_argument("--epochs", type=int, default=50)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--batch", type=int, default=256)
    ap.add_argument("--samples_per_class", type=int, default=10000)
    ap.add_argument("--out", type=str, default="head.pt")
    args = ap.parse_args()

    np.random.seed(0); torch.manual_seed(0)

    cache_files = sorted(CACHE.glob("*_features.npz"))
    print(f"found {len(cache_files)} cached features")
    train_files = []
    val_files = []
    for f in cache_files:
        if any(v in f.name for v in args.val):
            val_files.append(f)
        else:
            train_files.append(f)
    print(f"train: {len(train_files)}  val: {len(val_files)}")
    if not train_files:
        print("no training files found!"); return

    print("loading train...")
    train_ds = PatchDataset(train_files, samples_per_class=args.samples_per_class)
    val_ds   = PatchDataset(val_files)

    weights = class_weights(train_ds.labels)
    print(f"class weights: {weights.numpy()}")

    device = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
    head = LinearHead().to(device)
    opt = torch.optim.Adam(head.parameters(), lr=args.lr)
    crit = nn.CrossEntropyLoss(weight=weights.to(device))

    train_loader = DataLoader(train_ds, batch_size=args.batch, shuffle=True,
                              num_workers=0, drop_last=True)
    val_loader   = DataLoader(val_ds,   batch_size=args.batch, shuffle=False,
                              num_workers=0)

    best_val = -1
    for ep in range(args.epochs):
        head.train()
        t0 = time.time()
        loss_sum = 0; n = 0
        for x, y in train_loader:
            x = x.to(device); y = y.to(device)
            logits = head(x)
            loss = crit(logits, y)
            opt.zero_grad(); loss.backward(); opt.step()
            loss_sum += loss.item() * len(y); n += len(y)
        train_loss = loss_sum / n

        head.eval()
        with torch.no_grad():
            correct = 0; total = 0
            cm = np.zeros((N_CLASSES, N_CLASSES), dtype=np.int64)
            for x, y in val_loader:
                x = x.to(device); y_np = y.numpy()
                pred = head(x).argmax(dim=1).cpu().numpy()
                for tt, pp in zip(y_np, pred):
                    cm[tt, pp] += 1
                correct += (pred == y_np).sum(); total += len(y_np)
        acc = correct / max(total, 1)
        # per-class IoU
        ious = []
        for c in range(1, N_CLASSES):
            tp = cm[c, c]
            fp = cm[:, c].sum() - tp
            fn = cm[c, :].sum() - tp
            iou = tp / max(tp + fp + fn, 1)
            ious.append(iou)
        mean_iou = float(np.mean(ious))
        print(f"  ep {ep:3d}  loss {train_loss:.3f}  val_acc {acc:.3f}  "
              f"mIoU {mean_iou:.3f}  per-class IoU {[f'{i:.2f}' for i in ious]}  "
              f"({time.time()-t0:.0f}s)")

        if mean_iou > best_val:
            best_val = mean_iou
            torch.save({"state_dict": head.state_dict(),
                        "weights": weights, "epoch": ep,
                        "val_iou": ious}, OUT / args.out)
    print(f"best val mIoU: {best_val:.3f}, saved to {OUT/args.out}")


if __name__ == "__main__":
    main()
