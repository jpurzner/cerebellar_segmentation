"""Predict all 19 slides with the gold-boosted full-train model.
Saves predictions to predictions_for_correction/ — these are starting points
for the next round of hand-corrections.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import torch
import tifffile
from skimage.transform import resize

from dataset import N_CLASSES, load_slide
from train_head import LinearHead
from predict import cerebellum_mask_from_stack
from slide_manifest import manifest

ROOT = Path(__file__).resolve().parent
CACHE = ROOT / "cache_um0.5"
MODELS = ROOT / "models"
OUT = ROOT.parent / "predictions_for_correction"   # /python/predictions_for_correction/
OUT.mkdir(exist_ok=True)


def main():
    head_path = MODELS / "head_full_gold.pt"
    state = torch.load(head_path, map_location="cpu", weights_only=False)
    head = LinearHead()
    head.load_state_dict(state["state_dict"])
    head.eval()

    slides = manifest()
    print(f"predicting {len(slides)} slides with {head_path.name}")
    for s in slides:
        cache = CACHE / f"{s.slide_id}_features.npz"
        if not cache.exists():
            print(f"  {s.slide_id}: no cache, skip"); continue
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

        out_path = OUT / f"{s.slide_id}_dino_pred.tif"
        tifffile.imwrite(out_path, label_full, compression="zlib")
        u, c = np.unique(label_full, return_counts=True)
        sz = label_full.size
        cls_pcts = ", ".join(f"{ui}:{100*ci/sz:.1f}%" for ui, ci in zip(u, c) if ui > 0)
        print(f"  {s.mag}  {s.slide_id:30s}  {label_full.shape}  classes: {cls_pcts}")


if __name__ == "__main__":
    main()
