"""Run python-port on 10x slides + compare to MATLAB output side-by-side.

Tests slides where MATLAB is known to do well — does our port reproduce that?
"""
from __future__ import annotations
import time
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import tifffile
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from cerebellum_pipeline import segment_cerebellum
sys.path.insert(0, str(Path(__file__).resolve().parent))
from dataset import matlab_segments_to_labels
from slide_manifest import manifest, GOOD_10X

OUT = Path(__file__).resolve().parent.parent / "port_vs_matlab_10x"
OUT.mkdir(exist_ok=True)

LABEL_COLORS = np.array([
    [0,   0,   0],   [0,   0, 128], [0, 128, 255], [0, 255, 255],
    [255, 255, 0], [128, 255, 128],[255, 128, 0], [255,   0, 0],
], dtype=np.uint8)


def to_unit(im):
    im = im.astype(np.float32); lo, hi = im.min(), im.max()
    return (im - lo) / max(hi - lo, 1e-9)


def egl_iou(pred, gt):
    pe = (pred == 2) | (pred == 3); ge = (gt == 2) | (gt == 3)
    inter = (pe & ge).sum(); union = (pe | ge).sum()
    return inter / max(pe.sum(), 1), inter / max(ge.sum(), 1), inter / max(union, 1)


def main():
    slides = manifest()
    # Pick 4 representative 10x slides where MATLAB is known to be good
    test_ids = ["s5_C_p27", "s5_F_p27", "s4_F_p27", "s5_G_p27"]
    test_slides = [s for s in slides if s.slide_id in test_ids]
    print(f"running port on {len(test_slides)} 10x slides")

    rows = []
    for slide in test_slides:
        print(f"\n=== {slide.slide_id} ({slide.mag}) ===")
        stack = tifffile.imread(slide.input_path)
        print(f"  shape {stack.shape}")
        print(f"  channel_order = {slide.channel_order}, pixel = {slide.pixel_size_um}um")

        t0 = time.time()
        try:
            res = segment_cerebellum(stack,
                                     channel_order=slide.channel_order,
                                     pixel_size_um=slide.pixel_size_um)
        except Exception as e:
            print(f"  port FAILED: {e}")
            continue
        elapsed = time.time() - t0
        print(f"  port done in {elapsed:.0f}s")
        port_seg = res.set_bin

        gt_rgb = tifffile.imread(slide.gt_path)
        matlab_gt = matlab_segments_to_labels(gt_rgb)

        H = min(port_seg.shape[0], matlab_gt.shape[0])
        W = min(port_seg.shape[1], matlab_gt.shape[1])
        port_seg = port_seg[:H, :W]; matlab_gt = matlab_gt[:H, :W]

        Pm, Rm, IoUm = egl_iou(port_seg, matlab_gt)
        print(f"  EGL IoU vs MATLAB: P={Pm:.3f} R={Rm:.3f} IoU={IoUm:.3f}")
        rows.append({"slide": slide.slide_id, "elapsed": elapsed,
                     "iou": IoUm, "P": Pm, "R": Rm,
                     "port": port_seg, "matlab": matlab_gt, "stack": stack[:, :H, :W]})

        # Per-class % (port vs matlab)
        names = {2:"iEGL",3:"oEGL",4:"IGL",5:"ML",6:"DWL",7:"PCL"}
        print(f"  class %  matlab     port")
        for c in [2,3,4,5,6,7]:
            mp = 100*(matlab_gt==c).sum()/matlab_gt.size
            pp = 100*(port_seg==c).sum()/port_seg.size
            print(f"    {names[c]:5s}  {mp:6.2f}    {pp:6.2f}")

    # Render side-by-side panels
    print("\nrendering panels...")
    fig, axes = plt.subplots(len(rows), 3, figsize=(21, len(rows)*7))
    if len(rows) == 1: axes = axes[None, :]
    for r_i, r in enumerate(rows):
        rgb = np.stack([to_unit(r["stack"][0]), to_unit(r["stack"][1]),
                         to_unit(r["stack"][2])], axis=-1)
        rgb /= max(rgb.max(), 1e-9)
        axes[r_i, 0].imshow(rgb); axes[r_i, 0].axis("off")
        axes[r_i, 0].set_title(f"{r['slide']}\nRGB", fontsize=10)
        axes[r_i, 1].imshow(LABEL_COLORS[r["port"]]); axes[r_i, 1].axis("off")
        axes[r_i, 1].set_title(f"Port ({r['elapsed']:.0f}s)\nIoU vs MATLAB={r['iou']:.3f}",
                                fontsize=10)
        axes[r_i, 2].imshow(LABEL_COLORS[r["matlab"]]); axes[r_i, 2].axis("off")
        axes[r_i, 2].set_title("MATLAB GT", fontsize=10)
    fig.suptitle("Python-port vs MATLAB on 10x slides where MATLAB is good",
                  fontsize=14)
    fig.tight_layout()
    out_path = OUT / "port_vs_matlab_10x.png"
    fig.savefig(out_path, dpi=85)
    plt.close(fig)
    print(f"\nsaved {out_path}")
    print(f"\n=== summary ===")
    print(f"{'slide':12s} {'time':>6} {'P':>6} {'R':>6} {'IoU':>6}")
    for r in rows:
        print(f"  {r['slide']:10s} {r['elapsed']:>5.0f}s {r['P']:>6.3f} "
              f"{r['R']:>6.3f} {r['iou']:>6.3f}")
    if rows:
        ious = [r['iou'] for r in rows]
        print(f"  mean IoU = {np.mean(ious):.3f}")


if __name__ == "__main__":
    main()
