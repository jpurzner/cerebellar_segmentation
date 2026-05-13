"""Run the Python port on s1_2 and visually compare against MATLAB output."""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import tifffile
from cerebellum_pipeline import segment_cerebellum

INPUT = Path(
    "/Users/jpurzner/Dropbox/images/edu_repeat/p27/"
    "2018_05_22_s1_2_p27-0021/2018_05_22_s1_2_p27-0021_fused_crop.tif"
)
MATLAB_SEG = Path(
    "/Users/jpurzner/Dropbox/imaging_analysis/cerebellar_segmentation/_test_run/"
    "20x/2018_05_22_s1_2_p27-0021/"
    "2018_05_22_s1_2_p27-0021_fused_crop_segments.tif"
)
OUT = Path(__file__).resolve().parent / "figs" / "v_port"
OUT.mkdir(parents=True, exist_ok=True)
PX_UM = 0.5119049

# --- jet(8)-style colour map for our set_bin labels ---
LABEL_COLORS = np.array([
    [0,   0,   0],     # 0: bg
    [0,   0, 128],     # 1: all_cerebellum (dark blue)
    [0, 128, 255],     # 2: iEGL  (azure)
    [0, 255, 255],     # 3: oEGL  (cyan)
    [255, 255, 0],     # 4: IGL   (yellow)
    [128, 255, 128],   # 5: ML    (lt green)
    [255, 128, 0],     # 6: DWL   (orange)
    [255, 0, 0],       # 7: PCL   (red)
], dtype=np.uint8)


def label_to_rgb(set_bin):
    return LABEL_COLORS[set_bin]


def main():
    print(f"loading {INPUT.name}")
    stack = tifffile.imread(INPUT)
    print(f"  shape={stack.shape}  dtype={stack.dtype}")

    print("segmenting (Python port)…")
    import time; t0 = time.time()
    res = segment_cerebellum(stack, channel_order=("p27", "neun", "dapi"),
                             pixel_size_um=PX_UM)
    print(f"  done in {time.time()-t0:.1f}s")
    print(f"  cerebellum: {res.mask.sum()*(PX_UM/1000)**2:.3f} mm^2")
    print(f"  iEGL:       {res.a_final.sum()*(PX_UM/1000)**2:.3f} mm^2")
    print(f"  oEGL:       {res.b_final.sum()*(PX_UM/1000)**2:.3f} mm^2")
    print(f"  IGL:        {res.c_final.sum()*(PX_UM/1000)**2:.3f} mm^2")
    print(f"  ML:         {res.ml_final.sum()*(PX_UM/1000)**2:.3f} mm^2")
    print(f"  DWL:        {res.dwl_final.sum()*(PX_UM/1000)**2:.3f} mm^2")
    print(f"  PCL:        {res.pc_final.sum()*(PX_UM/1000)**2:.3f} mm^2")

    # render
    py_seg = label_to_rgb(res.set_bin)

    if MATLAB_SEG.exists():
        ml_seg = tifffile.imread(MATLAB_SEG)
    else:
        ml_seg = None

    fig, ax = plt.subplots(1, 3 if ml_seg is not None else 2, figsize=(20, 9))
    rgb = np.stack([res.p27_n, res.neun_n, res.dapi_n], axis=-1)
    rgb = rgb / max(rgb.max(), 1e-9)
    ax[0].imshow(rgb); ax[0].set_title("RGB input"); ax[0].axis("off")
    ax[1].imshow(py_seg); ax[1].set_title("Python port segmentation")
    ax[1].axis("off")
    if ml_seg is not None:
        ax[2].imshow(ml_seg); ax[2].set_title("MATLAB segmentation (reference)")
        ax[2].axis("off")
    fig.tight_layout(); fig.savefig(OUT / "01_port_vs_matlab.png", dpi=140)
    plt.close(fig)

    # overlay outlines
    fig, ax = plt.subplots(1, 1, figsize=(14, 12))
    ax.imshow(rgb * 0.55)
    contours_color = {
        "iEGL": (res.a_final, "magenta"),
        "oEGL": (res.b_final, "cyan"),
        "IGL":  (res.c_final, "yellow"),
        "DWL":  (res.dwl_final, "orange"),
        "PCL":  (res.pc_final, "red"),
    }
    for name, (m, c) in contours_color.items():
        ax.contour(m, levels=[0.5], colors=[c], linewidths=0.4)
    ax.set_title("Python port contours on RGB")
    ax.axis("off")
    fig.tight_layout(); fig.savefig(OUT / "02_contours_overlay.png", dpi=140)
    plt.close(fig)

    # save raw masks
    np.savez_compressed(OUT / "port_masks.npz",
                        set_bin=res.set_bin,
                        mask=res.mask,
                        a_final=res.a_final, b_final=res.b_final,
                        c_final=res.c_final, ml_final=res.ml_final,
                        dwl_final=res.dwl_final, pc_final=res.pc_final)
    print(f"figures + masks in {OUT}")


if __name__ == "__main__":
    main()
