"""v2: spokes + 1D layer detection per spoke + back-projection to 2D.

Pipeline:
  1. Build spokes via radial_spokes.build_spokes (already done in v1)
  2. For each valid spoke, sample (p27, NeuN, DAPI) profile and detect layer
     boundaries
  3. Back-project: for each pixel inside the spoke "wedge", assign the layer
     label corresponding to its position along the spoke
  4. Compare to MATLAB s1_2 GT (extract iEGL+oEGL from segments tif)
"""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage as ndi
import tifffile
from skimage import exposure
from radial_spokes import build_spokes, PIAL_SAMPLE_UM
from spoke_profiles import sample_along_spoke, detect_layers_1d

INPUT = Path(
    "/Users/jpurzner/Dropbox/images/edu_repeat/p27/"
    "2018_05_22_s1_2_p27-0021/2018_05_22_s1_2_p27-0021_fused_crop.tif"
)
PORT_MASKS = Path(__file__).resolve().parent / "figs" / "v_port" / "port_masks.npz"
SEG_GT_TIF = Path(
    "/Users/jpurzner/Dropbox/imaging_analysis/cerebellar_segmentation/_test_run/"
    "20x/2018_05_22_s1_2_p27-0021/"
    "2018_05_22_s1_2_p27-0021_fused_crop_segments.tif"
)
OUT = Path(__file__).resolve().parent / "figs" / "v_radial_v2"
OUT.mkdir(parents=True, exist_ok=True)
PX_UM = 0.5119049


def to_unit(im):
    im = im.astype(np.float32); lo, hi = im.min(), im.max()
    return (im - lo) / max(hi - lo, 1e-9)


def clahe(im, clip=0.01):
    return exposure.equalize_adapthist(im, clip_limit=clip, nbins=256)


def extract_egl_gt(seg_rgb):
    s = seg_rgb.astype(int)
    def col(target): return np.all(np.abs(s - np.array(target)) < 16, axis=-1)
    iEGL = col([0, 128, 255]); oEGL = col([0, 255, 255])
    return iEGL | oEGL


def main():
    print(f"loading {INPUT.name}")
    stack = tifffile.imread(INPUT)
    p27_n  = clahe(to_unit(stack[0]))
    neun_n = clahe(to_unit(stack[1]))
    dapi_n = clahe(to_unit(stack[2]))
    chans = np.stack([p27_n, neun_n, dapi_n], axis=0)
    H, W = p27_n.shape
    rgb = np.stack([p27_n, neun_n, dapi_n], axis=-1); rgb /= max(rgb.max(), 1e-9)

    print(f"loading port masks")
    z = np.load(PORT_MASKS)
    cereb = z['mask'].astype(bool); igl = z['c_final'].astype(bool)

    print("building spokes")
    sl = build_spokes(cereb, igl, pixel_size_um=PX_UM, boundary_smooth_um=3.0)
    print(f"  {len(sl.pia_points)} spokes, {sl.valid.sum()} valid")

    print("sampling 1D profiles + detecting layer boundaries")
    spoke_layers = []
    sample_step_um = 1.0
    n = len(sl.pia_points)
    for i in range(n):
        if not sl.valid[i]:
            spoke_layers.append(None); continue
        prof = sample_along_spoke(chans, sl.pia_points[i], sl.igl_endpoints[i],
                                   sample_step_um, PX_UM)
        layers = detect_layers_1d(prof, sample_step_um,
                                   channels={"p27": 0, "neun": 1, "dapi": 2})
        spoke_layers.append(layers)
        if i < 5:
            print(f"  spoke {i}: pia={layers['pia_um']:.0f}um egl_end={layers['egl_end_um']:.0f}"
                  f" iegl=[{layers['iegl_start_um']:.0f},{layers['iegl_end_um']:.0f}] "
                  f"igl_start={layers['igl_start_um']:.0f} L={layers['L_um']:.0f}")

    # Back-project: build label map by walking each spoke and assigning each
    # discrete pixel to a layer based on its 1D position.
    print("back-projecting to 2D label map")
    set_bin = np.zeros((H, W), dtype=np.uint8)
    for i in range(n):
        if not sl.valid[i] or spoke_layers[i] is None:
            continue
        py, px = sl.pia_points[i]
        ye, xe = sl.igl_endpoints[i]
        L_px = np.hypot(ye - py, xe - px)
        L_um = L_px * PX_UM
        ns = max(2, int(L_um / sample_step_um) + 1)
        ys = np.linspace(py, ye, ns)
        xs = np.linspace(px, xe, ns)
        s_um = np.linspace(0, L_um, ns)
        layers = spoke_layers[i]
        # Stamp each step's pixel with its layer label
        for k in range(ns):
            yi = int(round(ys[k])); xi = int(round(xs[k]))
            if not (0 <= yi < H and 0 <= xi < W):
                continue
            s = s_um[k]
            if s < layers['pia_um']:
                lbl = 0  # background-ish (above pia)
            elif s < layers['iegl_start_um']:
                lbl = 3  # oEGL (between pia and iEGL)
            elif s <= layers['iegl_end_um']:
                lbl = 2  # iEGL
            elif s < layers['igl_start_um']:
                lbl = 5  # ML
            else:
                lbl = 4  # IGL
            # write only if not already set, OR if same priority (this gives last-spoke wins)
            set_bin[yi, xi] = lbl

    # Now we only have labels along spoke lines. We need to thicken them.
    # Use morphological dilation of each label to fill the spoke wedges.
    print("thickening labels to fill wedges (dilate by half spoke spacing)")
    # half spacing in pixels
    half_spacing_px = max(1, int(round(PIAL_SAMPLE_UM / PX_UM / 2)))
    # dilate each label separately and reassign
    set_bin_thick = np.zeros_like(set_bin)
    for lbl in [5, 4, 3, 2]:  # ML, IGL, oEGL, iEGL (assignment order)
        m = (set_bin == lbl)
        if not m.any(): continue
        m = ndi.binary_dilation(m, iterations=half_spacing_px)
        # don't overwrite existing higher-priority labels (lower number = higher priority)
        m = m & (set_bin_thick == 0)
        set_bin_thick[m] = lbl
    # restrict to cerebellum mask
    set_bin_thick = set_bin_thick * cereb

    # === GT comparison ===
    if SEG_GT_TIF.exists():
        seg_rgb = tifffile.imread(SEG_GT_TIF)
        egl_gt = extract_egl_gt(seg_rgb)
    else:
        egl_gt = None

    egl_pred = (set_bin_thick == 2) | (set_bin_thick == 3)
    if egl_gt is not None:
        inter = (egl_pred & egl_gt).sum()
        union = (egl_pred | egl_gt).sum()
        p = inter / max(egl_pred.sum(), 1); r = inter / max(egl_gt.sum(), 1)
        iou = inter / max(union, 1)
        print(f"\n=== EGL match vs MATLAB GT ===")
        print(f"  GT     = {egl_gt.sum()*(PX_UM/1000)**2:.3f} mm^2")
        print(f"  pred   = {egl_pred.sum()*(PX_UM/1000)**2:.3f} mm^2")
        print(f"  inter  = {inter*(PX_UM/1000)**2:.3f} mm^2")
        print(f"  P      = {p:.3f}  R = {r:.3f}  IoU = {iou:.3f}")

    # === plots ===
    LABEL_COLORS = np.array([
        [0,   0,   0],     [0,   0, 128],   [0, 128, 255],   [0, 255, 255],
        [255, 255, 0],   [128, 255, 128], [255, 128, 0],   [255, 0, 0],
    ], dtype=np.uint8)
    seg_color = LABEL_COLORS[set_bin_thick]

    fig, ax = plt.subplots(1, 3, figsize=(22, 9))
    ax[0].imshow(rgb); ax[0].set_title("RGB"); ax[0].axis("off")
    ax[1].imshow(seg_color); ax[1].set_title("v2 radial-spoke segmentation"); ax[1].axis("off")
    if egl_gt is not None:
        ax[2].imshow(rgb * 0.5)
        ax[2].contour(egl_gt,  levels=[0.5], colors=["lime"],   linewidths=0.4, alpha=0.8)
        ax[2].contour(egl_pred, levels=[0.5], colors=["yellow"], linewidths=0.4, alpha=0.8)
        ax[2].set_title(f"GT (lime) vs pred (yellow)  IoU={iou:.2f}"); ax[2].axis("off")
    fig.tight_layout(); fig.savefig(OUT / "01_seg_vs_gt.png", dpi=140); plt.close(fig)

    # 1D profile demo: show 6 representative spokes
    print("plotting 6 representative profiles")
    indices = np.linspace(0, n - 1, 8).astype(int)
    indices = [i for i in indices if sl.valid[i]][:6]
    fig, axes = plt.subplots(2, 3, figsize=(18, 10), sharex=True)
    axes = axes.ravel()
    for k, i in enumerate(indices):
        prof = sample_along_spoke(chans, sl.pia_points[i], sl.igl_endpoints[i],
                                   sample_step_um, PX_UM)
        s_um = np.arange(prof.shape[1]) * sample_step_um
        ax = axes[k]
        ax.plot(s_um, prof[0], color="red",   label="p27")
        ax.plot(s_um, prof[1], color="green", label="NeuN")
        ax.plot(s_um, prof[2], color="blue",  label="DAPI")
        layers = spoke_layers[i]
        for n_, c in [("pia_um", "k"), ("egl_end_um", "magenta"),
                      ("iegl_start_um", "orange"), ("iegl_end_um", "orange"),
                      ("igl_start_um", "cyan")]:
            ax.axvline(layers[n_], color=c, linestyle="--", lw=0.8, alpha=0.7)
        ax.set_xlabel("position along spoke (um)")
        ax.set_ylabel("intensity")
        ax.set_title(f"spoke #{i}  L={layers['L_um']:.0f}um")
        if k == 0: ax.legend(fontsize=8)
    fig.suptitle("1D profiles + detected boundaries (k=pia, magenta=EGL end, "
                 "orange=iEGL, cyan=IGL start)")
    fig.tight_layout(); fig.savefig(OUT / "02_profiles.png", dpi=120); plt.close(fig)

    print(f"figures saved to {OUT}")


if __name__ == "__main__":
    main()
