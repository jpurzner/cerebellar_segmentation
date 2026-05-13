"""Run port v8 on all 4 validation slides + save 3-panel comparison + diff overlay."""
from pathlib import Path
import time
import numpy as np
import matplotlib.pyplot as plt
import tifffile
from cerebellum_pipeline import segment_cerebellum

PX_UM = 0.5119049
INPUT_BASE = Path("/Users/jpurzner/Dropbox/images/edu_repeat/p27")
SEG_BASE = Path("/Users/jpurzner/Dropbox/imaging_analysis/cerebellar_segmentation/_test_run/20x")
OUT_BASE = Path(__file__).resolve().parent / "figs" / "v_visual_all"
OUT_BASE.mkdir(parents=True, exist_ok=True)

SLIDES = [
    "2018_05_22_s1_2_p27-0021",
    "2018_05_22_s1_3_p27-0019",
    "2018_05_22_s2_5_p27-0026",
    "2018_05_22_s3_2_p27-0001",
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
    print(f"{'slide':40s} {'P':>6} {'R':>6} {'IoU':>6}")
    for slide in SLIDES:
        print(f"\n--- {slide} ---")
        inp_path = INPUT_BASE / slide / f"{slide}_fused_crop.tif"
        seg_path = SEG_BASE   / slide / f"{slide}_fused_crop_segments.tif"
        if not inp_path.exists() or not seg_path.exists():
            print(f"  MISSING file"); continue

        stack = tifffile.imread(inp_path)
        t0 = time.time()
        res = segment_cerebellum(stack, channel_order=("p27", "neun", "dapi"),
                                 pixel_size_um=PX_UM)
        print(f"  segment in {time.time()-t0:.0f}s")

        # RGB
        H, W = res.mask.shape
        rgb = np.stack([to_unit(stack[0]), to_unit(stack[1]), to_unit(stack[2])], axis=-1)
        seg = tifffile.imread(seg_path)
        port_seg = LABEL_COLORS[res.set_bin]

        # EGL metrics
        port_egl = res.a_final | res.b_final
        gt_egl = extract_egl_gt(seg)
        if port_egl.shape != gt_egl.shape:
            print(f"  shape mismatch port {port_egl.shape} vs gt {gt_egl.shape}")
            continue
        inter = (port_egl & gt_egl).sum()
        P = inter / max(port_egl.sum(),1); R = inter / max(gt_egl.sum(),1)
        IoU = inter / max((port_egl|gt_egl).sum(),1)
        print(f"  EGL P={P:.3f} R={R:.3f} IoU={IoU:.3f}")

        # diff overlay
        rgb_dim = rgb * 0.45
        tp = port_egl & gt_egl; fn = gt_egl & ~port_egl; fp = port_egl & ~gt_egl
        diff = np.copy(rgb_dim)
        diff[..., 0] = np.clip(diff[..., 0] + fp.astype(np.float32)*0.9 + fn.astype(np.float32)*0.9, 0, 1)
        diff[..., 1] = np.clip(diff[..., 1] + tp.astype(np.float32)*0.9 + fp.astype(np.float32)*0.7, 0, 1)

        # 4-panel figure
        fig, ax = plt.subplots(1, 4, figsize=(28, 9))
        ax[0].imshow(rgb); ax[0].set_title(f"{slide}\nRGB input"); ax[0].axis("off")
        ax[1].imshow(port_seg); ax[1].set_title("Python port v8 (Otsu)"); ax[1].axis("off")
        ax[2].imshow(seg); ax[2].set_title("MATLAB segments (reference)"); ax[2].axis("off")
        ax[3].imshow(diff); ax[3].set_title(
            f"EGL diff: green=TP red=FN yellow=FP\nP={P:.2f}  R={R:.2f}  IoU={IoU:.2f}")
        ax[3].axis("off")
        fig.tight_layout()
        fig.savefig(OUT_BASE / f"{slide}_compare.png", dpi=130)
        plt.close(fig)
        print(f"  saved {OUT_BASE / (slide + '_compare.png')}")


if __name__ == "__main__":
    main()
