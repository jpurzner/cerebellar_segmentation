"""Run port v11 (with phasesym EGL extension) on all 14 20x slides + save plots."""
from pathlib import Path
import time
import numpy as np
import matplotlib.pyplot as plt
import tifffile
from cerebellum_pipeline import segment_cerebellum

PX_UM = 0.5119049
INPUT_BASE = Path("/Users/jpurzner/Dropbox/images/edu_repeat/p27")
SEG_BASE = Path("/Users/jpurzner/Dropbox/imaging_analysis/cerebellar_segmentation/_test_run/20x")
OUT_BASE = Path(__file__).resolve().parent / "figs" / "v_all_14"
OUT_BASE.mkdir(parents=True, exist_ok=True)

SLIDES = [
    "2018_05_22_s1_2_p27-0021",
    "2018_05_22_s1_3_p27-0019",
    "2018_05_22_s1_4_p27-0020",
    "2018_05_22_s1_5_p27-0018",
    "2018_05_22_s2_1_p27-0023",
    "2018_05_22_s2_2_p27",
    "2018_05_22_s2_3_p27-0025",
    "2018_05_22_s2_4_p27-0024",
    "2018_05_22_s2_5_p27-0026",
    "2018_05_22_s3_1_p27-0002",
    "2018_05_22_s3_2_p27-0001",
    "2018_05_22_s3_3_p27-0005",
    "2018_05_22_s3_4_p27-0004",
    "2018_05_22_s3_5_p27-0006",
]

LABEL_COLORS = np.array([
    [0,   0,   0],   [0,   0, 128], [0, 128, 255], [0, 255, 255],
    [255, 255, 0], [128, 255, 128],[255, 128, 0], [255,   0, 0],
], dtype=np.uint8)


def to_unit(im):
    im = im.astype(np.float32); lo, hi = im.min(), im.max()
    return (im - lo) / max(hi - lo, 1e-9)


def extract_egl_gt(seg_rgb):
    s = seg_rgb.astype(int)
    def col(t): return np.all(np.abs(s - np.array(t)) < 16, axis=-1)
    return col([0, 128, 255]) | col([0, 255, 255])


def main():
    print(f"{'slide':40s} {'iEGL':>8} {'oEGL':>8} {'P':>6} {'R':>6} {'IoU':>6} {'time_s':>7}")
    print("-" * 100)
    rows = []
    for slide in SLIDES:
        inp_path = INPUT_BASE / slide / f"{slide}_fused_crop.tif"
        seg_path = SEG_BASE   / slide / f"{slide}_fused_crop_segments.tif"
        if not inp_path.exists() or not seg_path.exists():
            print(f"  {slide:40s}  MISSING file"); continue
        try:
            stack = tifffile.imread(inp_path)
            t0 = time.time()
            res = segment_cerebellum(stack, channel_order=("p27", "neun", "dapi"),
                                     pixel_size_um=PX_UM)
            elapsed = time.time() - t0
        except Exception as e:
            print(f"  {slide:40s}  FAILED: {e}"); continue

        H, W = res.mask.shape
        rgb = np.stack([to_unit(stack[0]), to_unit(stack[1]), to_unit(stack[2])], axis=-1)
        seg = tifffile.imread(seg_path)
        port_seg = LABEL_COLORS[res.set_bin]

        port_egl = res.a_final | res.b_final
        gt_egl = extract_egl_gt(seg)
        if port_egl.shape != gt_egl.shape:
            print(f"  {slide}: SHAPE MISMATCH"); continue
        inter = (port_egl & gt_egl).sum()
        P = inter / max(port_egl.sum(),1); R = inter / max(gt_egl.sum(),1)
        IoU = inter / max((port_egl|gt_egl).sum(),1)
        iegl_mm2 = res.a_final.sum() * (PX_UM/1000)**2
        oegl_mm2 = res.b_final.sum() * (PX_UM/1000)**2
        print(f"  {slide:40s} {iegl_mm2:8.3f} {oegl_mm2:8.3f} "
              f"{P:6.3f} {R:6.3f} {IoU:6.3f} {elapsed:7.0f}")
        rows.append((slide, iegl_mm2, oegl_mm2, P, R, IoU))

        # diff overlay
        rgb_dim = rgb * 0.45
        tp = port_egl & gt_egl; fn = gt_egl & ~port_egl; fp = port_egl & ~gt_egl
        diff = np.copy(rgb_dim)
        diff[..., 0] = np.clip(diff[..., 0] + fp.astype(np.float32)*0.9 + fn.astype(np.float32)*0.9, 0, 1)
        diff[..., 1] = np.clip(diff[..., 1] + tp.astype(np.float32)*0.9 + fp.astype(np.float32)*0.7, 0, 1)

        fig, ax = plt.subplots(1, 4, figsize=(28, 9))
        ax[0].imshow(rgb); ax[0].set_title(f"{slide}\nRGB input"); ax[0].axis("off")
        ax[1].imshow(port_seg); ax[1].set_title("Python port v11 (with phasesym)"); ax[1].axis("off")
        ax[2].imshow(seg); ax[2].set_title("MATLAB segments (reference)"); ax[2].axis("off")
        ax[3].imshow(diff); ax[3].set_title(
            f"EGL: green=TP red=FN yellow=FP\nP={P:.2f}  R={R:.2f}  IoU={IoU:.2f}")
        ax[3].axis("off")
        fig.tight_layout()
        fig.savefig(OUT_BASE / f"{slide}_compare.png", dpi=130)
        plt.close(fig)

    print()
    if rows:
        ious = [r[5] for r in rows]
        Ps = [r[3] for r in rows]; Rs = [r[4] for r in rows]
        print(f"  IoU: mean={np.mean(ious):.3f}  min={min(ious):.3f}  max={max(ious):.3f}")
        print(f"  P  : mean={np.mean(Ps):.3f}    R: mean={np.mean(Rs):.3f}")


if __name__ == "__main__":
    main()
