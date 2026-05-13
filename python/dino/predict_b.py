"""Predict with ViT-B head + cereb mask."""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import torch
from skimage.transform import resize

from dataset import N_CLASSES, load_slide
from train_head_b import LinearHeadB
from predict import cerebellum_mask_from_stack

CACHE = Path(__file__).resolve().parent / "cache_b"
MODELS = Path(__file__).resolve().parent / "models"
PREDS = Path(__file__).resolve().parent / "preds_b"
PREDS.mkdir(exist_ok=True)
INPUT_BASE = Path("/Users/jpurzner/Dropbox/images/edu_repeat/p27")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--head", type=str, default="head_b.pt")
    ap.add_argument("--slide", type=str, required=True)
    args = ap.parse_args()

    cache_files = sorted(CACHE.glob(f"*{args.slide}*_features.npz"))
    if not cache_files:
        print(f"no cache for {args.slide}"); return
    cache_path = cache_files[0]
    print(f"loading {cache_path.name}")
    d = np.load(cache_path, allow_pickle=True)
    feats = d["features"].astype(np.float32)
    gy, gx = d["patch_grid"]
    H, W = int(d["H"]), int(d["W"])

    state = torch.load(MODELS / args.head, map_location="cpu", weights_only=False)
    head = LinearHeadB()
    head.load_state_dict(state["state_dict"])
    head.eval()

    print("predicting")
    with torch.no_grad():
        pred = head(torch.from_numpy(feats)).argmax(dim=1).cpu().numpy().astype(np.uint8)
    label_grid = pred.reshape(gy, gx)
    label_full = resize(label_grid, (H, W), order=0, preserve_range=True,
                         anti_aliasing=False).astype(np.uint8)

    slide_dir_name = cache_path.stem.replace("_features", "")
    inp_path = INPUT_BASE / slide_dir_name / f"{slide_dir_name}_fused_crop.tif"
    if inp_path.exists():
        stack = load_slide(inp_path, channel_order=("p27", "neun", "dapi"))
        mask = cerebellum_mask_from_stack(stack)
        label_full = label_full * mask.astype(np.uint8)
    out_path = PREDS / f"{slide_dir_name}_pred.npz"
    np.savez_compressed(out_path, label_grid=label_grid, label_full=label_full,
                        H=H, W=W, gy=gy, gx=gx)
    print(f"saved {out_path}")


if __name__ == "__main__":
    main()
