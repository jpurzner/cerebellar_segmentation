"""Train linear head on ViT-B/14 features (768-dim)."""
from __future__ import annotations
import argparse, time
from pathlib import Path
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader

from dataset import N_CLASSES
from train_head import PatchDataset, class_weights

CACHE = Path(__file__).resolve().parent / "cache_b"
MODELS = Path(__file__).resolve().parent / "models"
MODELS.mkdir(exist_ok=True)


class LinearHeadB(nn.Module):
    def __init__(self, in_dim=768, hidden=64, out=N_CLASSES):
        super().__init__()
        self.net = nn.Sequential(nn.Linear(in_dim, hidden), nn.ReLU(),
                                  nn.Linear(hidden, out))
    def forward(self, x): return self.net(x)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--epochs", type=int, default=50)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--batch", type=int, default=256)
    ap.add_argument("--samples_per_class", type=int, default=10000)
    ap.add_argument("--val", nargs="+", default=["s1_4", "s2_5"])
    ap.add_argument("--out", default="head_b.pt")
    args = ap.parse_args()

    cache_files = sorted(CACHE.glob("*_features.npz"))
    train_files = [f for f in cache_files if not any(v in f.stem for v in args.val)]
    val_files = [f for f in cache_files if any(v in f.stem for v in args.val)]
    print(f"train slides: {len(train_files)}, val slides: {len(val_files)}")

    train_ds = PatchDataset(train_files, samples_per_class=args.samples_per_class)
    val_ds   = PatchDataset(val_files)
    print(f"  train samples: {len(train_ds)}, val samples: {len(val_ds)}")

    weights = class_weights(train_ds.labels)
    print(f"class weights: {weights.numpy().round(3)}")

    device = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
    head = LinearHeadB().to(device)
    opt = torch.optim.Adam(head.parameters(), lr=args.lr)
    crit = nn.CrossEntropyLoss(weight=weights.to(device))

    train_loader = DataLoader(train_ds, batch_size=args.batch, shuffle=True, drop_last=True)
    val_loader = DataLoader(val_ds, batch_size=args.batch)
    best_iou = -1; best_state = None
    for ep in range(args.epochs):
        head.train()
        total_loss = 0; n_batches = 0
        t0 = time.time()
        for x, y in train_loader:
            x = x.to(device); y = y.to(device)
            loss = crit(head(x), y)
            opt.zero_grad(); loss.backward(); opt.step()
            total_loss += loss.item(); n_batches += 1
        head.eval()
        with torch.no_grad():
            cm = np.zeros((N_CLASSES, N_CLASSES), dtype=np.int64)
            correct = total = 0
            for x, y in val_loader:
                x = x.to(device); y_np = y.numpy()
                pred = head(x).argmax(dim=1).cpu().numpy()
                correct += (pred == y_np).sum(); total += len(y_np)
                for tt, pp in zip(y_np, pred):
                    cm[tt, pp] += 1
        ious = []
        for c in range(1, N_CLASSES):
            tp = cm[c, c]; fp = cm[:, c].sum() - tp; fn = cm[c, :].sum() - tp
            ious.append(tp / max(tp + fp + fn, 1))
        miou = float(np.mean(ious)); val_acc = correct / max(total, 1)
        elapsed = time.time() - t0
        print(f"  ep {ep:3d}  loss {total_loss/n_batches:.3f}  val_acc {val_acc:.3f}  "
              f"mIoU {miou:.3f}  per-class IoU {[f'{i:.2f}' for i in ious]}  ({elapsed:.0f}s)")
        if miou > best_iou:
            best_iou = miou
            best_state = {k: v.cpu().clone() for k, v in head.state_dict().items()}
    out_path = MODELS / args.out
    torch.save({"state_dict": best_state, "best_miou": best_iou}, out_path)
    print(f"best val mIoU: {best_iou:.3f}, saved to {out_path}")


if __name__ == "__main__":
    main()
