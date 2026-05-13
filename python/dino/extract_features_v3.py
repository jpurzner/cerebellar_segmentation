"""Extract DINOv3 ViT-S/16 features. Requires DINOv3 weights to be present at:
  ~/.cache/torch/hub/checkpoints/dinov3_vits16_pretrain_lvd1689m-*.pth

See DINOV3_INSTRUCTIONS.md for download steps.
"""
from __future__ import annotations
import argparse, time, warnings
from pathlib import Path
import numpy as np
import torch
import tifffile
from skimage import exposure
from scipy.ndimage import zoom

from slide_manifest import manifest
from dataset import (load_slide, load_label_for_slide, tile_iter, N_CLASSES)

PATCH = 16
TILE = 528  # multiple of PATCH (33 patches per tile side)
OVERLAP = 64
FEAT_DIM = 384


def slide_to_dino_input(stack: np.ndarray, target_h: int, target_w: int) -> np.ndarray:
    p27, ch2, ch3 = stack[0], stack[1], stack[2]
    def to_unit(im):
        im = im.astype(np.float32); lo, hi = im.min(), im.max()
        return (im - lo) / max(hi - lo, 1e-9)
    p27_n = exposure.equalize_adapthist(to_unit(p27), clip_limit=0.01)
    ch2_n = exposure.equalize_adapthist(to_unit(ch2), clip_limit=0.01)
    ch3_n = exposure.equalize_adapthist(to_unit(ch3), clip_limit=0.01)
    rgb = np.stack([p27_n, ch2_n, ch3_n], axis=0)
    H, W = rgb.shape[1], rgb.shape[2]
    if (H, W) != (target_h, target_w):
        sy = target_h / H; sx = target_w / W
        rgb = zoom(rgb, (1.0, sy, sx), order=1)
    mean = np.array([0.485, 0.456, 0.406], dtype=np.float32).reshape(3, 1, 1)
    std  = np.array([0.229, 0.224, 0.225], dtype=np.float32).reshape(3, 1, 1)
    rgb = (rgb - mean) / std
    return rgb.astype(np.float32)


def extract_one(slide, target_um_per_px, cache_dir, model, device):
    out_path = cache_dir / f"{slide.slide_id}_features.npz"
    if out_path.exists():
        print(f"  cached: {out_path.name}"); return
    t0 = time.time()
    stack = load_slide(slide.input_path, channel_order=slide.channel_order)
    H_orig, W_orig = stack.shape[1], stack.shape[2]
    scale = slide.pixel_size_um / target_um_per_px
    target_h = int(round(H_orig * scale))
    target_w = int(round(W_orig * scale))
    rgb = slide_to_dino_input(stack, target_h, target_w)
    H, W = rgb.shape[1], rgb.shape[2]
    print(f"    orig ({H_orig},{W_orig})  scale {scale:.3f}  → ({H},{W}), patch={PATCH*target_um_per_px:.1f}um")

    labels_orig = None
    if slide.gt_path.exists():
        labels_orig = load_label_for_slide(slide.gt_path)
        if labels_orig.shape != (H_orig, W_orig):
            Hc = min(labels_orig.shape[0], H_orig)
            Wc = min(labels_orig.shape[1], W_orig)
            labels_orig = labels_orig[:Hc, :Wc]
            stack = stack[:, :Hc, :Wc]

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
            # DINOv3 returns dict similar to v2
            if isinstance(out, dict):
                feats = out.get("x_norm_patchtokens", out.get("x_patchtokens", None))
                if feats is None:
                    # fallback: get last hidden state
                    feats = list(out.values())[0]
            else:
                feats = out
        gy = h_p // PATCH; gx = w_p // PATCH
        feats_np = feats.reshape(1, gy, gx, -1).cpu().numpy().astype(np.float16)[0]
        # adapt FEAT_DIM if model differs
        if feats_np.shape[-1] != FEAT_DIM:
            print(f"    WARN: feature dim {feats_np.shape[-1]} != expected {FEAT_DIM}")
        gy0, gx0 = y0 // PATCH, x0 // PATCH
        for i in range(gy):
            for j in range(gx):
                gi, gj = gy0 + i, gx0 + j
                if gi >= gy_total or gj >= gx_total: continue
                if seen[gi, gj] == 0: feat_grid[gi, gj] = feats_np[i, j]
                else: feat_grid[gi, gj] = ((feat_grid[gi, gj].astype(np.float32) +
                                             feats_np[i, j].astype(np.float32)) * 0.5).astype(np.float16)
                seen[gi, gj] += 1
        n_tiles += 1
    print(f"    {n_tiles} tiles in {time.time()-t0:.0f}s")

    label_grid = None
    if labels_orig is not None:
        label_grid = np.zeros((gy_total, gx_total), dtype=np.uint8)
        ylo = (np.arange(gy_total) * PATCH / scale).astype(int)
        yhi = ((np.arange(gy_total) + 1) * PATCH / scale).astype(int)
        xlo = (np.arange(gx_total) * PATCH / scale).astype(int)
        xhi = ((np.arange(gx_total) + 1) * PATCH / scale).astype(int)
        for i in range(gy_total):
            ya, yb = ylo[i], min(yhi[i], labels_orig.shape[0])
            if ya >= yb: continue
            for j in range(gx_total):
                xa, xb = xlo[j], min(xhi[j], labels_orig.shape[1])
                if xa >= xb: continue
                cell = labels_orig[ya:yb, xa:xb]
                if cell.size == 0: continue
                counts = np.bincount(cell.ravel(), minlength=N_CLASSES)
                label_grid[i, j] = int(counts.argmax())

    np.savez_compressed(out_path,
        features=feat_grid.reshape(-1, FEAT_DIM),
        patch_grid=np.array([gy_total, gx_total]),
        H=H_orig, W=W_orig, H_dino=H, W_dino=W,
        labels=label_grid if label_grid is not None else np.zeros((1,), dtype=np.uint8),
        has_labels=label_grid is not None,
        slide_id=slide.slide_id, mag=slide.mag,
        pixel_size_um=slide.pixel_size_um,
        target_um_per_px=target_um_per_px, scale=scale,
    )
    print(f"  saved {out_path.name}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", type=float, default=0.5)
    ap.add_argument("--model", default="dinov3_vits16",
                    choices=["dinov3_vits16", "dinov3_vits16plus", "dinov3_vitb16", "dinov3_vitl16"])
    args = ap.parse_args()
    cache_dir = Path(__file__).resolve().parent / f"cache_v3_um{args.target:g}"
    cache_dir.mkdir(exist_ok=True)

    device = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
    print(f"device: {device}, target {args.target} um/px (patch ~{PATCH*args.target:.1f}um)")
    print(f"loading {args.model}...")
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        try:
            model = torch.hub.load('facebookresearch/dinov3', args.model, verbose=False)
        except Exception as e:
            print(f"\nFAILED to load DINOv3: {e}")
            print("\nSee DINOV3_INSTRUCTIONS.md for manual weight download.")
            return
    model = model.to(device).eval()
    print("  done")

    slides = manifest()
    todo = [s for s in slides if not (cache_dir / f"{s.slide_id}_features.npz").exists()]
    print(f"\n{len(slides)} slides, {len(todo)} to extract")
    for s in todo:
        print(f"\n=== {s.mag}  {s.slide_id} ===")
        extract_one(s, args.target, cache_dir, model, device)


if __name__ == "__main__":
    main()
