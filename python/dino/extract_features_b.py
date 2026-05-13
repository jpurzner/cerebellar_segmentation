"""Re-extract features with DINOv2 ViT-B/14 (87M params, 768-dim features).

Output cache: cache_b/ to avoid clobbering ViT-S cache.
"""
from __future__ import annotations
import argparse
import warnings
import time
from pathlib import Path
import numpy as np
import torch
import tifffile

from dataset import (load_slide, slide_to_dino_input, load_label_for_slide,
                     tile_iter, N_CLASSES)

INPUT_BASE = Path("/Users/jpurzner/Dropbox/images/edu_repeat/p27")
SEG_BASE   = Path("/Users/jpurzner/Dropbox/imaging_analysis/cerebellar_segmentation/_test_run/20x")
CACHE_BASE = Path(__file__).resolve().parent / "cache_b"
CACHE_BASE.mkdir(exist_ok=True)

PATCH = 14
TILE  = 532
OVERLAP = 56
FEAT_DIM = 768   # ViT-B output dim


def extract_one_slide(slide_dir: Path, model: torch.nn.Module,
                      device: torch.device, save_labels: bool = True):
    name = slide_dir.name
    inp_path = slide_dir / f"{name}_fused_crop.tif"
    seg_path = SEG_BASE / name / f"{name}_fused_crop_segments.tif"
    out_path = CACHE_BASE / f"{name}_features.npz"
    if out_path.exists():
        print(f"  already cached: {out_path.name}"); return
    if not inp_path.exists():
        print(f"  missing input: {inp_path}"); return

    print(f"  loading + CLAHE")
    t0 = time.time()
    stack = load_slide(inp_path, channel_order=("p27", "neun", "dapi"))
    rgb = slide_to_dino_input(stack)
    H, W = rgb.shape[1], rgb.shape[2]
    print(f"    {H}x{W}, {time.time()-t0:.0f}s")

    labels = None
    if save_labels and seg_path.exists():
        labels = load_label_for_slide(seg_path)
        if labels.shape != (H, W):
            print(f"    SHAPE MISMATCH labels {labels.shape} vs image ({H},{W})")
            labels = None

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

    np.savez_compressed(out_path,
                        features=feat_grid.reshape(-1, FEAT_DIM),
                        patch_grid=np.array([gy_total, gx_total]),
                        H=H, W=W,
                        labels=label_grid if label_grid is not None else np.zeros((1,), dtype=np.uint8),
                        has_labels=label_grid is not None)
    print(f"  saved {out_path.name}  feat={feat_grid.nbytes/1e6:.1f}MB")


def main():
    device = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
    print(f"device: {device}")
    print("loading DINOv2 ViT-B/14 (87M params)...")
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        model = torch.hub.load('facebookresearch/dinov2', 'dinov2_vitb14', verbose=False)
    model = model.to(device).eval()
    print(f"  done")

    slide_dirs = sorted([d for d in INPUT_BASE.iterdir() if d.is_dir()
                         and (d / f"{d.name}_fused_crop.tif").exists()])
    print(f"processing {len(slide_dirs)} slides")
    for d in slide_dirs:
        print(f"\n=== {d.name} ===")
        extract_one_slide(d, model, device)


if __name__ == "__main__":
    main()
