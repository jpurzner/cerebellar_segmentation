"""Quick diagnostic: is the s1_3 NeuN signal much weaker in the upper half?
If so, that explains why MATLAB's global IGL threshold misses the upper IGL
entirely and the whole upper cerebellum gets called DWL by default.
"""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import tifffile
import scipy.ndimage as ndi
from skimage import exposure, morphology
from skimage.morphology import disk

REF = Path(
    "/Users/jpurzner/Dropbox/images/edu_repeat/p27/"
    "2018_05_22_s1_3_p27-0019/2018_05_22_s1_3_p27-0019_fused_crop.tif"
)
OUT = Path(__file__).resolve().parent / "figs" / "v_baseline_diag"
OUT.mkdir(parents=True, exist_ok=True)
PX_UM = 0.5119049
CH_P27, CH_NEUN, CH_DAPI = 0, 1, 2


def to_unit(img):
    img = img.astype(np.float32)
    lo, hi = img.min(), img.max()
    return (img - lo) / max(hi - lo, 1e-9)


def main():
    print(f"loading {REF.name}")
    stack = tifffile.imread(REF)
    if stack.shape[0] not in (3, 4) and stack.shape[-1] in (3, 4):
        stack = np.moveaxis(stack, -1, 0)
    p27_raw  = stack[CH_P27]
    neun_raw = stack[CH_NEUN]
    dapi_raw = stack[CH_DAPI]
    H, W = p27_raw.shape
    print(f"shape={p27_raw.shape} dtype={p27_raw.dtype}")

    # cerebellum mask using raw (no CLAHE) sums to avoid CLAHE artifacts
    s = (to_unit(p27_raw) + to_unit(dapi_raw)) / 2
    raw_mask = s > 0.10
    closed = morphology.closing(raw_mask, footprint=disk(int(round(5/PX_UM))))
    mask = ndi.binary_fill_holes(closed)
    lbl, _ = ndi.label(mask)
    sizes = np.bincount(lbl.ravel()); sizes[0] = 0
    mask = (lbl == sizes.argmax())

    # === plot 1: raw channel intensities + bg-subtracted local means ===
    print("plot 1: raw intensity heatmaps for each channel")
    fig, ax = plt.subplots(2, 3, figsize=(18, 12))
    for i, (raw, name) in enumerate([(p27_raw, "p27"), (neun_raw, "NeuN"),
                                       (dapi_raw, "DAPI")]):
        # raw image
        ax[0, i].imshow(np.where(mask, raw, np.nan), cmap="hot",
                          vmin=np.percentile(raw[mask], 5),
                          vmax=np.percentile(raw[mask], 99))
        ax[0, i].set_title(f"{name} raw  (5th..99th pct in mask)")
        ax[0, i].axis("off")
        # large-scale local mean (Gaussian blur with sigma 100 um) -- shows
        # regional bias / signal attenuation
        sigma_px = 100.0 / PX_UM
        smooth = ndi.gaussian_filter(raw.astype(np.float32) * mask, sigma_px)
        ax[1, i].imshow(np.where(mask, smooth, np.nan), cmap="viridis")
        ax[1, i].set_title(f"{name} smoothed sigma 100um  (regional bias)")
        ax[1, i].axis("off")
    fig.tight_layout(); fig.savefig(OUT / "01_channel_uniformity.png", dpi=130)
    plt.close(fig)

    # === plot 2: row-wise mean intensity (top-to-bottom signal trend) ===
    print("plot 2: row-wise mean intensity profiles")
    rows = np.arange(H) * PX_UM / 1000  # mm
    fig, ax = plt.subplots(1, 1, figsize=(12, 5))
    for raw, name, color in [(p27_raw, "p27", "red"),
                              (neun_raw, "NeuN", "green"),
                              (dapi_raw, "DAPI", "blue")]:
        # mean across columns where mask=true
        masked = raw.astype(np.float32) * mask
        n_per_row = mask.sum(axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            row_mean = masked.sum(axis=1) / np.maximum(n_per_row, 1)
        valid = n_per_row > 100
        ax.plot(rows[valid], row_mean[valid], color=color, label=name)
    ax.set_xlabel("row position (mm, top=0)")
    ax.set_ylabel("mean intensity inside mask (raw uint16 units)")
    ax.set_title("row-wise mean: is one half systematically dimmer?")
    ax.legend()
    fig.tight_layout(); fig.savefig(OUT / "02_row_means.png", dpi=120); plt.close(fig)

    # === stats ===
    H2 = H // 2
    upper = mask.copy(); upper[H2:] = False
    lower = mask.copy(); lower[:H2] = False
    print(f"\nupper-half mask area = {upper.sum() * (PX_UM/1000)**2:.3f} mm^2")
    print(f"lower-half mask area = {lower.sum() * (PX_UM/1000)**2:.3f} mm^2")
    for raw, name in [(p27_raw, "p27"), (neun_raw, "NeuN"), (dapi_raw, "DAPI")]:
        u_med = np.median(raw[upper])
        l_med = np.median(raw[lower])
        ratio = u_med / max(l_med, 1)
        print(f"  {name}: upper med={u_med:.0f}  lower med={l_med:.0f}  ratio={ratio:.2f}")

    print(f"done. figures in {OUT}")


if __name__ == "__main__":
    main()
