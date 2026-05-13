"""Extract DINO features using the slide_manifest. Skips already-cached slides
so we only process new ones (10x).

Cache: cache_unified/  (separate from cache/ to avoid confusion w/ ViT-S 20x-only)
"""
from __future__ import annotations
import time
import warnings
from pathlib import Path
import numpy as np
import torch
import tifffile
from skimage import exposure

from slide_manifest import manifest
from dataset import (load_slide, slide_to_dino_input, load_label_for_slide,
                     tile_iter, N_CLASSES)

CACHE = Path(__file__).resolve().parent / "cache_unified"
CACHE.mkdir(exist_ok=True)
PATCH = 14
TILE = 532
OVERLAP = 56
FEAT_DIM = 384  # ViT-S


def extract_one(slide, model, device):
    out_path = CACHE / f"{slide.slide_id}_features.npz"
    if out_path.exists():
        print(f"  cached: {out_path.name}"); return
    t0 = time.time()
    stack = load_slide(slide.input_path, channel_order=slide.channel_order)
    rgb = slide_to_dino_input(stack)
    H, W = rgb.shape[1], rgb.shape[2]
    print(f"    {H}x{W}, channel_order={slide.channel_order}, {time.time()-t0:.0f}s")

    labels = None
    if slide.gt_path.exists():
        labels = load_label_for_slide(slide.gt_path)
        if labels.shape != (H, W):
            print(f"    SHAPE MISMATCH: labels {labels.shape} vs image ({H},{W}) — cropping")
            Hc = min(labels.shape[0], H); Wc = min(labels.shape[1], W)
            labels = labels[:Hc, :Wc]
            rgb = rgb[:, :Hc, :Wc]
            H, W = Hc, Wc

    gy_total = H // PATCH
    gx_total = W // PATCH
    feat_grid = np.zeros((gy_total, gx_total, FEAT_DIM), dtype=np.float16)
    seen = np.zeros((gy_total, gx_total), dtype=np.uint8)

    t0 = time.time()
    n_tiles = 0
    for (y0, y1, x0, x1) in tile_iter(H, W, tile=TILE, overlap=OVERLAP, snap_to=PATCH):
        tile = rgb[:, y0:y1, x0:x1]
        h, w = tile.shape[1], tile.shape[2]
        h_p = (h // PATCH) * PATCH
        w_p = (w // PATCH) * PATCH
        if h_p == 0 or w_p == 0: continue
        tile = tile[:, :h_p, :w_p]
        x = torch.from_numpy(tile[None]).to(device, dtype=torch.float32)
        with torch.no_grad():
            out = model.forward_features(x)
            feats = out["x_norm_patchtokens"]
        gy = h_p // PATCH; gx = w_p // PATCH
        feats = feats.reshape(1, gy, gx, FEAT_DIM).cpu().numpy().astype(np.float16)[0]
        gy0, gx0 = y0 // PATCH, x0 // PATCH
        for i in range(gy):
            for j in range(gx):
                gi, gj = gy0 + i, gx0 + j
                if gi >= gy_total or gj >= gx_total: continue
                if seen[gi, gj] == 0:
                    feat_grid[gi, gj] = feats[i, j]
                else:
                    feat_grid[gi, gj] = ((feat_grid[gi, gj].astype(np.float32) +
                                          feats[i, j].astype(np.float32)) * 0.5).astype(np.float16)
                seen[gi, gj] += 1
        n_tiles += 1
    print(f"    {n_tiles} tiles in {time.time()-t0:.0f}s")

    label_grid = None
    if labels is not None:
        label_grid = np.zeros((gy_total, gx_total), dtype=np.uint8)
        for i in range(gy_total):
            for j in range(gx_total):
                cell = labels[i*PATCH:(i+1)*PATCH, j*PATCH:(j+1)*PATCH]
                counts = np.bincount(cell.ravel(), minlength=N_CLASSES)
                label_grid[i, j] = int(counts.argmax())

    np.savez_compressed(
        out_path,
        features=feat_grid.reshape(-1, FEAT_DIM),
        patch_grid=np.array([gy_total, gx_total]),
        H=H, W=W,
        labels=label_grid if label_grid is not None else np.zeros((1,), dtype=np.uint8),
        has_labels=label_grid is not None,
        slide_id=slide.slide_id,
        mag=slide.mag,
        pixel_size_um=slide.pixel_size_um,
    )
    print(f"  saved {out_path.name}  feat={feat_grid.nbytes/1e6:.1f}MB")


def main():
    device = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
    print(f"device: {device}")
    print("loading DINOv2 ViT-S/14...")
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        model = torch.hub.load('facebookresearch/dinov2', 'dinov2_vits14', verbose=False)
    model = model.to(device).eval()
    print("  done")

    slides = manifest()
    print(f"\n{len(slides)} slides in manifest, processing only un-cached:")
    todo = [s for s in slides if not (CACHE / f"{s.slide_id}_features.npz").exists()]
    cached = [s for s in slides if (CACHE / f"{s.slide_id}_features.npz").exists()]
    print(f"  already cached: {len(cached)}")
    print(f"  to extract:     {len(todo)}")
    for s in todo:
        print(f"\n=== {s.mag}  {s.slide_id} ===")
        extract_one(s, model, device)


if __name__ == "__main__":
    main()
