"""Run port v5 on multiple slides and report EGL IoU vs MATLAB GT."""
from pathlib import Path
import sys, time
import numpy as np
import tifffile
from cerebellum_pipeline import segment_cerebellum

PX_UM = 0.5119049

# slides to test: (name, input.tif, matlab_segments.tif)
INPUT_BASE = Path("/Users/jpurzner/Dropbox/images/edu_repeat/p27")
SEG_BASE = Path("/Users/jpurzner/Dropbox/imaging_analysis/cerebellar_segmentation/_test_run/20x")

# pick 3 representative slides (other than s1_2 which we already tested)
SLIDES = [
    "2018_05_22_s1_2_p27-0021",  # tuning slide ref
    "2018_05_22_s1_3_p27-0019",   # historically problematic in MATLAB
    "2018_05_22_s2_5_p27-0026",   # different mouse, clean folds
    "2018_05_22_s3_2_p27-0001",   # different mouse, different shape
]


def extract_egl_gt(seg_rgb):
    s = seg_rgb.astype(int)
    def col(t): return np.all(np.abs(s - np.array(t)) < 16, axis=-1)
    return col([0, 128, 255]) | col([0, 255, 255])   # iEGL azure | oEGL cyan


def main():
    print(f"{'slide':40s} {'iEGL_mm2':>10} {'oEGL_mm2':>10} {'P':>6} {'R':>6} {'IoU':>6} {'time_s':>7}")
    print("-" * 100)
    rows = []
    for slide in SLIDES:
        inp_path = INPUT_BASE / slide / f"{slide}_fused_crop.tif"
        seg_path = SEG_BASE   / slide / f"{slide}_fused_crop_segments.tif"
        if not inp_path.exists():
            print(f"  MISSING input: {inp_path}")
            continue
        if not seg_path.exists():
            print(f"  MISSING segments: {seg_path}")
            continue
        stack = tifffile.imread(inp_path)
        t0 = time.time()
        try:
            res = segment_cerebellum(stack, channel_order=("p27", "neun", "dapi"),
                                     pixel_size_um=PX_UM)
        except Exception as e:
            print(f"  {slide}: FAILED  {e}")
            continue
        elapsed = time.time() - t0
        port_egl = res.a_final | res.b_final
        seg = tifffile.imread(seg_path)
        gt_egl = extract_egl_gt(seg)
        # quick alignment check (shapes should match)
        if port_egl.shape != gt_egl.shape:
            print(f"  {slide}: SHAPE MISMATCH  port={port_egl.shape}  gt={gt_egl.shape}")
            continue
        inter = (port_egl & gt_egl).sum()
        union = (port_egl | gt_egl).sum()
        P = inter / max(port_egl.sum(), 1)
        R = inter / max(gt_egl.sum(), 1)
        IoU = inter / max(union, 1)
        iegl_mm2 = res.a_final.sum() * (PX_UM/1000)**2
        oegl_mm2 = res.b_final.sum() * (PX_UM/1000)**2
        print(f"  {slide:40s} {iegl_mm2:10.3f} {oegl_mm2:10.3f} "
              f"{P:6.3f} {R:6.3f} {IoU:6.3f} {elapsed:7.1f}")
        rows.append((slide, iegl_mm2, oegl_mm2, P, R, IoU))

    print()
    if rows:
        ious = [r[5] for r in rows]
        print(f"  IoU summary: mean={np.mean(ious):.3f}  min={min(ious):.3f}  max={max(ious):.3f}")
        print(f"  s1_2 reference: IoU=0.600 (this was our tuning slide)")


if __name__ == "__main__":
    main()
