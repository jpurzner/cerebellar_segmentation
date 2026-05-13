"""Test if downsampling/blurring helps the MATLAB-style port for s1_2.

Variants of s1_2:
  native        — original 20x (0.512 µm/px)
  ds2x          — downsampled 2× (matches 10x physical resolution, ~1.02 µm/px)
  blur_sig2     — Gaussian blur sigma=2 µm
  blur_sig4     — Gaussian blur sigma=4 µm
  blur_sig8     — Gaussian blur sigma=8 µm

For each, run Python-port (MATLAB-equivalent) and compare to:
  - MATLAB GT (the original MATLAB output)
  - GOLD standard (your hand-corrected labels)

Save panels to /python/variant_test/.
"""
from __future__ import annotations
import time
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage as ndi
import tifffile
from skimage import filters
from scipy.ndimage import zoom

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from cerebellum_pipeline import segment_cerebellum
sys.path.insert(0, str(Path(__file__).resolve().parent))
from dataset import matlab_segments_to_labels
from slide_manifest import manifest

ROOT = Path(__file__).resolve().parent
LABELLED_DIR = ROOT.parent / "labelled"
OUT = ROOT.parent / "variant_test"
OUT.mkdir(exist_ok=True)
PX_UM = 0.5119049

LABEL_COLORS = np.array([
    [0,   0,   0],   [0,   0, 128], [0, 128, 255], [0, 255, 255],
    [255, 255, 0], [128, 255, 128],[255, 128, 0], [255,   0, 0],
], dtype=np.uint8)


def to_unit(im):
    im = im.astype(np.float32); lo, hi = im.min(), im.max()
    return (im - lo) / max(hi - lo, 1e-9)


def make_variant(stack_3ch, kind):
    """Apply transformation to the channels-first (3, H, W) stack."""
    if kind == "native":
        return stack_3ch.copy(), PX_UM
    if kind == "ds2x":
        # downsample by 2 → physical res 1.024 µm/px (matches 10x)
        out = np.zeros((3, stack_3ch.shape[1] // 2, stack_3ch.shape[2] // 2),
                       dtype=stack_3ch.dtype)
        for c in range(3):
            out[c] = zoom(stack_3ch[c], 0.5, order=1)
        return out, PX_UM * 2
    if kind.startswith("blur_sig"):
        sig_um = float(kind.replace("blur_sig", ""))
        sig_px = sig_um / PX_UM
        out = np.zeros_like(stack_3ch)
        for c in range(3):
            out[c] = ndi.gaussian_filter(stack_3ch[c].astype(np.float32),
                                          sigma=sig_px).astype(stack_3ch.dtype)
        return out, PX_UM
    raise ValueError(kind)


def egl_iou(pred, gt):
    pe = (pred == 2) | (pred == 3); ge = (gt == 2) | (gt == 3)
    inter = (pe & ge).sum(); union = (pe | ge).sum()
    if pe.sum() == 0 or ge.sum() == 0:
        return 0.0, 0.0, 0.0
    return inter / max(pe.sum(), 1), inter / max(ge.sum(), 1), inter / max(union, 1)


def main():
    slides = manifest()
    s1_2 = next(s for s in slides if s.slide_id == "2018_05_22_s1_2_p27-0021")

    stack = tifffile.imread(s1_2.input_path)
    print(f"loaded s1_2: {stack.shape}")

    matlab_gt = matlab_segments_to_labels(tifffile.imread(s1_2.gt_path))

    # find gold standard
    cands = sorted(
        list(LABELLED_DIR.glob(f"{s1_2.slide_id}_corrected.tif*")) +
        list(LABELLED_DIR.glob(f"{s1_2.slide_id}_labelled.tif*")) +
        list(LABELLED_DIR.glob(f"{s1_2.slide_id}_labelld.tif*")),
        key=lambda p: p.stat().st_mtime, reverse=True)
    cands = [c for c in cands if c.suffix in (".tif", ".tiff")]
    gold = tifffile.imread(cands[0]).astype(np.uint8) if cands else None
    print(f"gold: {cands[0].name if cands else 'none'}")

    rgb = np.stack([to_unit(stack[0]), to_unit(stack[1]), to_unit(stack[2])], axis=-1)

    # Run each variant
    variants = ["native", "ds2x", "blur_sig2", "blur_sig4", "blur_sig8"]
    rows = []
    for kind in variants:
        print(f"\n--- variant: {kind} ---")
        t0 = time.time()
        stk, pixel_um = make_variant(stack, kind)
        try:
            res = segment_cerebellum(stk,
                                     channel_order=("p27", "neun", "dapi"),
                                     pixel_size_um=pixel_um)
        except Exception as e:
            print(f"  port failed: {e}"); continue
        elapsed = time.time() - t0
        port_seg = res.set_bin
        print(f"  port done in {elapsed:.0f}s, shape {port_seg.shape}")

        # upsample port output back to original resolution if needed (for ds2x)
        if port_seg.shape != matlab_gt.shape:
            from skimage.transform import resize
            port_seg_full = resize(port_seg, matlab_gt.shape, order=0,
                                    preserve_range=True, anti_aliasing=False).astype(np.uint8)
        else:
            port_seg_full = port_seg
        H = min(port_seg_full.shape[0], matlab_gt.shape[0])
        W = min(port_seg_full.shape[1], matlab_gt.shape[1])
        port_seg_full = port_seg_full[:H, :W]
        m_gt = matlab_gt[:H, :W]
        rgb_c = rgb[:H, :W]

        # IoU vs MATLAB GT
        Pm, Rm, IoUm = egl_iou(port_seg_full, m_gt)
        print(f"  EGL IoU vs MATLAB GT: {IoUm:.3f}")

        # IoU vs GOLD
        if gold is not None:
            g = gold[:H, :W]
            Pg, Rg, IoUg = egl_iou(port_seg_full, g)
            print(f"  EGL IoU vs GOLD:      {IoUg:.3f}")
        else:
            IoUg = float("nan")

        rows.append({"variant": kind, "elapsed_s": elapsed,
                     "iou_matlab": IoUm, "iou_gold": IoUg,
                     "pred": port_seg_full})

    # render combined panel: [variant rows × (RGB | port pred | matlab GT | gold)]
    print("\nrendering panels...")
    if gold is not None:
        n_cols = 4; col_titles = ["RGB", "Port output", "MATLAB GT", "GOLD"]
    else:
        n_cols = 3; col_titles = ["RGB", "Port output", "MATLAB GT"]
    fig, axes = plt.subplots(len(rows), n_cols, figsize=(n_cols*7, len(rows)*7))
    if len(rows) == 1: axes = axes[None, :]
    for r_i, r in enumerate(rows):
        H = r["pred"].shape[0]; W = r["pred"].shape[1]
        axes[r_i, 0].imshow(rgb[:H, :W]); axes[r_i, 0].axis("off")
        axes[r_i, 0].set_title(f"{r['variant']}\n{col_titles[0]}", fontsize=11)
        axes[r_i, 1].imshow(LABEL_COLORS[r["pred"]]); axes[r_i, 1].axis("off")
        axes[r_i, 1].set_title(f"{col_titles[1]}\nIoU(MATLAB)={r['iou_matlab']:.3f} "
                                f"IoU(GOLD)={r['iou_gold']:.3f}",
                                fontsize=10)
        axes[r_i, 2].imshow(LABEL_COLORS[matlab_gt[:H, :W]]); axes[r_i, 2].axis("off")
        axes[r_i, 2].set_title(col_titles[2], fontsize=11)
        if gold is not None:
            axes[r_i, 3].imshow(LABEL_COLORS[gold[:H, :W]]); axes[r_i, 3].axis("off")
            axes[r_i, 3].set_title(col_titles[3], fontsize=11)
    fig.suptitle("Port-on-variants of s1_2 vs MATLAB GT and your GOLD",
                  fontsize=14)
    fig.tight_layout()
    panel_path = OUT / "s1_2_variants_panel.png"
    fig.savefig(panel_path, dpi=85)
    plt.close(fig)
    print(f"\nsaved {panel_path}")

    # save individual variant pred TIFs
    for r in rows:
        out_t = OUT / f"s1_2_{r['variant']}_pred.tif"
        tifffile.imwrite(out_t, r["pred"], compression="zlib")

    print("\n=== summary ===")
    print(f"{'variant':12s} {'time':>6} {'IoU(MATLAB)':>12} {'IoU(GOLD)':>11}")
    for r in rows:
        print(f"  {r['variant']:10s} {r['elapsed_s']:>5.0f}s "
              f"{r['iou_matlab']:>12.3f} {r['iou_gold']:>11.3f}")


if __name__ == "__main__":
    main()
