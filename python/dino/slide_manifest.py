"""Slide manifest: good (non-catastrophic-MATLAB) slides for training/eval.

Excludes:
  - 20x catastrophic-MATLAB: s1_3, s2_1, s2_3, s3_1, s3_3, s3_4
  - 10x bad (user filter):   s4_B, s4_E, s4_G, s6_E

Channel orders (from MATLAB pipeline):
  - 10x files: ch1=p27, ch2=DAPI, ch3=NeuN
  - 20x files: ch1=p27, ch2=NeuN, ch3=DAPI
"""
from pathlib import Path
from dataclasses import dataclass


P10_INPUT = Path("/Users/jpurzner/Dropbox/images/2017_05_18_p27_p7_ezh2_cko")
P20_INPUT = Path("/Users/jpurzner/Dropbox/images/edu_repeat/p27")
P10_GT_BASE = P10_INPUT
P20_GT_BASE = Path("/Users/jpurzner/Dropbox/imaging_analysis/cerebellar_segmentation/_test_run/20x")


@dataclass
class Slide:
    slide_id: str
    input_path: Path
    gt_path: Path
    channel_order: tuple   # (R-channel role, G-channel role, B-channel role)
    pixel_size_um: float
    mag: str               # "10x" or "20x"


# 8 good 20x slides
GOOD_20X = [
    "2018_05_22_s1_2_p27-0021",
    "2018_05_22_s1_4_p27-0020",
    "2018_05_22_s1_5_p27-0018",
    "2018_05_22_s2_2_p27",
    "2018_05_22_s2_4_p27-0024",
    "2018_05_22_s2_5_p27-0026",
    "2018_05_22_s3_2_p27-0001",
    "2018_05_22_s3_5_p27-0006",
]


# 11 good 10x slides — input filenames vs slide IDs
GOOD_10X = [
    ("s4_A_p27", "a_p27_cor-prune.tif"),
    ("s4_C_p27", "2017_05_17_C_p27_dapi_10x_pano_red_cor_crop.tif"),
    ("s4_F_p27", "s4_F_p27_merge_crop.tif"),
    ("s5_A_p27", "s5_A_p27_cor_merge_crop.tif"),
    ("s5_B_p27", "s5_B_p27_cor_merge_crop.tif"),
    ("s5_C_p27", "s5_C_p27_cor_merge_crop.tif"),
    ("s5_E_p27", "s5_E_p27_cor_merge_crop.tif"),
    ("s5_F_p27", "s5_F_p27_cor_merge_crop.tif"),
    ("s5_G_p27", "s5_G_p27_cor_merge_crop.tif"),
    ("s6_C_p27", "s6_CC_p27_merge_crop.tif"),
    ("s6_G_p27", "s6_G_p27_merge_crop.tif"),
]


def manifest():
    slides = []
    # 20x
    for sid in GOOD_20X:
        inp = P20_INPUT / sid / f"{sid}_fused_crop.tif"
        gt = P20_GT_BASE / sid / f"{sid}_fused_crop_segments.tif"
        if inp.exists() and gt.exists():
            slides.append(Slide(
                slide_id=sid, input_path=inp, gt_path=gt,
                channel_order=("p27", "neun", "dapi"),
                pixel_size_um=0.5119049, mag="20x"))
        else:
            print(f"  20x MISSING: {sid}")
    # 10x
    for sid, fname in GOOD_10X:
        inp = P10_INPUT / fname
        gt_name = fname.replace(".tif", "_segments.tif")
        gt = P10_GT_BASE / gt_name
        if inp.exists() and gt.exists():
            slides.append(Slide(
                slide_id=sid, input_path=inp, gt_path=gt,
                channel_order=("p27", "dapi", "neun"),
                pixel_size_um=1.0239, mag="10x"))
        else:
            print(f"  10x MISSING: {sid} (inp={inp.exists()}, gt={gt.exists()})")
    return slides


if __name__ == "__main__":
    slides = manifest()
    print(f"\n{len(slides)} slides in manifest")
    for s in slides:
        print(f"  {s.mag}  {s.slide_id:30s}  pixel={s.pixel_size_um:.4f}  "
              f"chan={s.channel_order}")
