"""Train DINOv3 head on ALL 19 slides with gold-boost (no holdout)."""
from __future__ import annotations
import time
from pathlib import Path
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, ConcatDataset

from dataset import N_CLASSES
from train_head import LinearHead, PatchDataset, class_weights
from slide_manifest import manifest

ROOT = Path(__file__).resolve().parent
LABELLED = ROOT.parent / "labelled"
CACHE = ROOT / "cache_v3_um0.5"
MODELS = ROOT / "models"


def main():
    slides = manifest()
    cache_files = [CACHE / f"{s.slide_id}_features.npz" for s in slides]
    cache_files = [f for f in cache_files if f.exists()]
    print(f"all-slide training on {len(cache_files)} v3 caches")

    gold_files, weak_files = [], []
    for f in cache_files:
        sid = f.stem.replace("_features", "")
        cands = list(LABELLED.glob(f"{sid}_corrected.tif*")) + \
                list(LABELLED.glob(f"{sid}_labelled.tif*")) + \
                list(LABELLED.glob(f"{sid}_labelld.tif*"))
        cands = [c for c in cands if c.suffix in (".tif", ".tiff")]
        (gold_files if cands else weak_files).append(f)
    print(f"  gold: {[f.stem.replace('_features','') for f in gold_files]}")
    print(f"  weak: {len(weak_files)} slides")

    weak_ds = PatchDataset(weak_files, samples_per_class=10000)
    if gold_files:
        gold_ds = PatchDataset(gold_files, samples_per_class=40000)
        all_labels = np.concatenate([weak_ds.labels, gold_ds.labels])
        train_ds = ConcatDataset([weak_ds, gold_ds])
    else:
        all_labels = weak_ds.labels
        train_ds = weak_ds
    weights = class_weights(all_labels)

    device = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
    head = LinearHead().to(device)
    opt = torch.optim.Adam(head.parameters(), lr=1e-3)
    crit = nn.CrossEntropyLoss(weight=weights.to(device))
    train_loader = DataLoader(train_ds, batch_size=256, shuffle=True, drop_last=True)

    best_loss = float('inf'); last_state = None
    for ep in range(25):
        head.train()
        total = 0; nb = 0; t0 = time.time()
        for x, y in train_loader:
            x = x.to(device); y = y.to(device)
            loss = crit(head(x), y)
            opt.zero_grad(); loss.backward(); opt.step()
            total += loss.item(); nb += 1
        avg = total / max(nb, 1)
        if avg < best_loss:
            best_loss = avg
            last_state = {k: v.cpu().clone() for k, v in head.state_dict().items()}
        print(f"  ep {ep:3d} loss {avg:.3f} ({time.time()-t0:.0f}s)")

    out = MODELS / "head_v3_full_gold.pt"
    torch.save({"state_dict": last_state, "best_loss": best_loss}, out)
    print(f"\nsaved {out}")


if __name__ == "__main__":
    main()
