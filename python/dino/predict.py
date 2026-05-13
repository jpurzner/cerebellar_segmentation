"""Run trained head on a slide's features → per-pixel label map.

Masks predictions to inside-cerebellum only (background was dropped from
training so the model can't predict class 0).
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import torch
import tifffile
import scipy.ndimage as ndi
from skimage import exposure, morphology

from dataset import N_CLASSES, load_slide
from train_head import LinearHead

CACHE = Path(__file__).resolve().parent / "cache"
MODELS = Path(__file__).resolve().parent / "models"
PREDS = Path(__file__).resolve().parent / "preds"
PREDS.mkdir(exist_ok=True)
INPUT_BASE = Path("/Users/jpurzner/Dropbox/images/edu_repeat/p27")


def cerebellum_mask_from_stack(stack):
    """Quick cerebellum mask from p27+DAPI sum > threshold + cleanup."""
    p27 = stack[0]; dapi = stack[2]
    s = (p27 + dapi) / 2
    s = s / max(s.max(), 1e-9)
    raw = s > 0.10
    closed = ndi.binary_closing(raw, iterations=10)
    filled = ndi.binary_fill_holes(closed)
    lbl, n = ndi.label(filled)
    if n == 0: return filled
    sizes = np.bincount(lbl.ravel()); sizes[0] = 0
    return lbl == sizes.argmax()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--head", type=str, default="head.pt")
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

    print(f"loading head {args.head}")
    state = torch.load(MODELS / args.head, map_location="cpu", weights_only=False)
    head = LinearHead()
    head.load_state_dict(state["state_dict"])
    head.eval()

    print("predicting patch labels")
    with torch.no_grad():
        x = torch.from_numpy(feats)
        logits = head(x).cpu().numpy()
    pred = logits.argmax(axis=1).astype(np.uint8)
    label_grid = pred.reshape(gy, gx)

    from skimage.transform import resize
    label_full = resize(label_grid, (H, W), order=0, preserve_range=True,
                         anti_aliasing=False).astype(np.uint8)

    # mask out background using a quick cereb mask
    print("masking with cerebellum mask")
    slide_dir_name = cache_path.stem.replace("_features", "")
    inp_path = INPUT_BASE / slide_dir_name / f"{slide_dir_name}_fused_crop.tif"
    if inp_path.exists():
        stack = load_slide(inp_path, channel_order=("p27", "neun", "dapi"))
        mask = cerebellum_mask_from_stack(stack)
        label_full = label_full * mask.astype(np.uint8)
        print(f"  mask area = {mask.sum() * (0.5119049/1000)**2:.3f} mm^2")
    else:
        print(f"  WARN: input TIF not found at {inp_path}")

    out_path = PREDS / f"{slide_dir_name}_pred.npz"
    np.savez_compressed(out_path, label_grid=label_grid, label_full=label_full,
                        H=H, W=W, gy=gy, gx=gx)
    print(f"saved {out_path}")


if __name__ == "__main__":
    main()
