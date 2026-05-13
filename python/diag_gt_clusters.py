"""Show each unique-color cluster from the s1_2 segments.tif as its own panel
so we can identify which is iEGL/oEGL/IGL/etc by spatial pattern."""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import tifffile

SEG_TIF = Path(
    "/Users/jpurzner/Dropbox/imaging_analysis/cerebellar_segmentation/_test_run/"
    "20x/2018_05_22_s1_2_p27-0021/"
    "2018_05_22_s1_2_p27-0021_fused_crop_segments.tif"
)
A_IEGL_TIF = SEG_TIF.parent / SEG_TIF.name.replace("_segments", "_a_iegl_only")

OUT = Path(__file__).resolve().parent / "figs" / "v_baseline_s1_2"
OUT.mkdir(parents=True, exist_ok=True)
PX_UM = 0.5119049


def main():
    seg = tifffile.imread(SEG_TIF)
    print("seg shape", seg.shape, "dtype", seg.dtype)
    # find unique colors (down-sample for speed)
    flat = seg.reshape(-1, 3)
    # quantise to nearest 8 per channel to merge near-duplicates
    q = (flat // 8) * 8
    uniq, counts = np.unique(q, axis=0, return_counts=True)
    order = np.argsort(-counts)
    uniq = uniq[order]; counts = counts[order]
    print(f"distinct quantised colors: {len(uniq)}")
    print("top 10 by area:")
    for i in range(min(10, len(uniq))):
        print(f"  {tuple(uniq[i].tolist())}  n={counts[i]}  "
              f"{counts[i]*(PX_UM/1000)**2:.3f} mm^2")

    # Pick the top 8 colors (other than black) as actual layer colors
    nonblack = [i for i in range(len(uniq)) if uniq[i].sum() > 10][:8]
    H, W = seg.shape[:2]
    fig, ax = plt.subplots(2, 4, figsize=(20, 11))
    for i, idx in enumerate(nonblack):
        c = uniq[idx]
        # find pixels close to this color in the FULL (non-quantised) image
        match = np.all(np.abs(seg.astype(int) - c.astype(int)) < 16, axis=-1)
        a = ax[i // 4, i % 4]
        a.imshow(match, cmap="gray")
        a.set_title(f"#{i}: rgb={tuple(int(v) for v in c)}  "
                    f"n={counts[idx]}  {counts[idx]*(PX_UM/1000)**2:.3f} mm^2")
        a.axis("off")
    fig.tight_layout(); fig.savefig(OUT / "00_cluster_per_panel.png", dpi=130)
    plt.close(fig)

    if A_IEGL_TIF.exists():
        a_iegl = tifffile.imread(A_IEGL_TIF)
        if a_iegl.ndim == 3: a_iegl = a_iegl.sum(-1)
        iegl_mask = a_iegl > 0
        print(f"\na_iegl_only mask area = {iegl_mask.sum()*(PX_UM/1000)**2:.3f} mm^2")
        # show on top of segments
        fig, ax = plt.subplots(1, 2, figsize=(16, 8))
        ax[0].imshow(seg)
        ax[0].set_title("segments.tif")
        ax[0].axis("off")
        ax[1].imshow(seg * 0.4)
        ax[1].contour(iegl_mask, levels=[0.5], colors=["lime"], linewidths=0.6)
        ax[1].set_title("a_iegl_only outline (lime) on dim segments")
        ax[1].axis("off")
        fig.tight_layout(); fig.savefig(OUT / "00_iegl_overlay.png", dpi=130)
        plt.close(fig)

    print(f"figures in {OUT}")


if __name__ == "__main__":
    main()
