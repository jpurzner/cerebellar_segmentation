"""
Diagnostic for geodesic-depth approach to cerebellum layer segmentation.

Goal: on the reference image (s1_3, 20x WT), compute the geodesic distance
from the pia inward and check whether isodepth contours at 10/20/50/100/200 um
visually trace the layer boundaries.

If they do, the depth coordinate is a viable backbone for the new pipeline.
If they peel off the layers anywhere, we know what to fix.

Saves plots to figs/.
"""
from __future__ import annotations

import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as ndi
import skfmm
import tifffile
from skimage import exposure, filters, morphology

REF = Path(
    "/Users/jpurzner/Dropbox/images/edu_repeat/p27/"
    "2018_05_22_s1_3_p27-0019/2018_05_22_s1_3_p27-0019_fused_crop.tif"
)
OUT = Path(__file__).resolve().parent / "figs"
OUT.mkdir(exist_ok=True)

PX_UM = 0.5119049  # 20x pixel size in micrometers

# Channel order for 20x files: ch1=p27, ch2=NeuN, ch3=DAPI
# (from cerebellum_threshold_segment20x.m: a=ch1, c=ch2, b=ch3 with a=p27, b=DAPI, c=NeuN)
CH_P27, CH_NEUN, CH_DAPI = 0, 1, 2


def to_unit(img: np.ndarray) -> np.ndarray:
    """Per-channel rescale to [0, 1] like MATLAB's mat2gray."""
    img = img.astype(np.float32)
    lo, hi = img.min(), img.max()
    return (img - lo) / max(hi - lo, 1e-9)


def clahe(img: np.ndarray, clip: float = 0.01, tiles: tuple = (50, 50)) -> np.ndarray:
    """Approximate adapthisteq with rayleigh distribution; skimage's CLAHE is uniform."""
    return exposure.equalize_adapthist(img, kernel_size=None, clip_limit=clip, nbins=256)


def cerebellum_mask(p27: np.ndarray, dapi: np.ndarray) -> np.ndarray:
    """Port of all_cerebellum block from MATLAB.

    Sum p27 + DAPI background-subtract and threshold; close gaps; fill holes;
    keep largest connected component.
    """
    s = p27 + dapi
    s = s / s.max()
    raw = s > 0.10  # MATLAB: > 0.1
    # close gaps (50 um worth ~ 100 px at 20x)
    r_close = max(1, int(round(20 / PX_UM)))  # 20 um
    closed = morphology.binary_closing(raw, footprint=morphology.disk(r_close))
    filled = ndi.binary_fill_holes(closed)
    # keep largest piece
    lbl, n = ndi.label(filled)
    if n > 1:
        sizes = np.bincount(lbl.ravel())
        sizes[0] = 0
        keep = sizes.argmax()
        filled = lbl == keep
    return filled


def geodesic_depth(mask: np.ndarray) -> np.ndarray:
    """Geodesic distance from pia (mask boundary), within mask, in micrometers."""
    # Mark the boundary as zero level set; everything else = 1 inside the mask, masked out
    # outside.
    boundary = mask & ~ndi.binary_erosion(mask)
    phi = np.ones_like(mask, dtype=float)
    phi[boundary] = 0.0
    masked_phi = np.ma.masked_array(phi, mask=~mask)
    d_px = skfmm.distance(masked_phi).filled(0.0)
    return d_px * PX_UM


def main() -> None:
    print(f"loading {REF.name}")
    stack = tifffile.imread(REF)
    print(f"  shape={stack.shape}  dtype={stack.dtype}")
    # tifffile may return (frames, H, W) or (H, W, frames). Normalize to (frames, H, W)
    if stack.ndim != 3:
        raise SystemExit(f"unexpected stack shape {stack.shape}")
    if stack.shape[0] not in (3, 4) and stack.shape[-1] in (3, 4):
        stack = np.moveaxis(stack, -1, 0)
    print(f"  normalized to {stack.shape}  channels first")

    p27_raw, neun_raw, dapi_raw = stack[CH_P27], stack[CH_NEUN], stack[CH_DAPI]
    print(f"  p27  range=({p27_raw.min()}, {p27_raw.max()})")
    print(f"  neun range=({neun_raw.min()}, {neun_raw.max()})")
    print(f"  dapi range=({dapi_raw.min()}, {dapi_raw.max()})")

    # Histogram-equalize each channel so the mask threshold is on a normalized scale
    print("normalizing channels (CLAHE)")
    p27 = clahe(to_unit(p27_raw))
    neun = clahe(to_unit(neun_raw))
    dapi = clahe(to_unit(dapi_raw))

    print("computing cerebellum mask")
    mask = cerebellum_mask(p27, dapi)
    area_mm2 = mask.sum() * (PX_UM / 1000) ** 2
    print(f"  mask area = {area_mm2:.3f} mm^2  ({mask.sum()} px)")

    # --- save mask + RGB overlay ------------------------------------------------
    print("saving channel overview + mask")
    fig, ax = plt.subplots(2, 2, figsize=(12, 10))
    ax[0, 0].imshow(p27, cmap="Reds"); ax[0, 0].set_title("p27 (CLAHE)"); ax[0, 0].axis("off")
    ax[0, 1].imshow(neun, cmap="Greens"); ax[0, 1].set_title("NeuN (CLAHE)"); ax[0, 1].axis("off")
    ax[1, 0].imshow(dapi, cmap="Blues"); ax[1, 0].set_title("DAPI (CLAHE)"); ax[1, 0].axis("off")
    rgb = np.stack([p27, neun, dapi], axis=-1)
    rgb = rgb / rgb.max()
    ax[1, 1].imshow(rgb)
    ax[1, 1].contour(mask, levels=[0.5], colors=["yellow"], linewidths=1.0)
    ax[1, 1].set_title("RGB + cerebellum mask (yellow)")
    ax[1, 1].axis("off")
    fig.tight_layout()
    fig.savefig(OUT / "01_channels_and_mask.png", dpi=120)
    plt.close(fig)

    # --- geodesic depth -------------------------------------------------------
    print("computing geodesic distance from pia (this can take 30-60s)")
    depth_um = geodesic_depth(mask)
    print(f"  max depth = {depth_um.max():.1f} um")

    # --- depth heatmap with isodepth contours ---------------------------------
    print("saving depth heatmap")
    fig, ax = plt.subplots(1, 1, figsize=(12, 10))
    im = ax.imshow(np.where(mask, depth_um, np.nan), cmap="viridis", vmin=0, vmax=300)
    levels = [10, 20, 30, 50, 80, 120, 200]
    cs = ax.contour(depth_um, levels=levels, colors=["red", "orange", "yellow",
                                                      "lime", "cyan", "magenta", "white"],
                    linewidths=0.6)
    ax.clabel(cs, inline=True, fontsize=7)
    ax.set_title("Geodesic depth from pia (um) with isodepth contours")
    ax.axis("off")
    plt.colorbar(im, ax=ax, label="um", fraction=0.03)
    fig.tight_layout()
    fig.savefig(OUT / "02_depth_heatmap.png", dpi=140)
    plt.close(fig)

    # --- isodepth contours overlaid on RGB ------------------------------------
    print("saving isodepth-on-RGB")
    fig, ax = plt.subplots(1, 1, figsize=(14, 12))
    ax.imshow(rgb)
    cs = ax.contour(depth_um, levels=levels, colors=["red", "orange", "yellow",
                                                      "lime", "cyan", "magenta", "white"],
                    linewidths=0.6)
    ax.clabel(cs, inline=True, fontsize=7)
    ax.set_title("Isodepth contours over RGB (red=p27, green=NeuN, blue=DAPI)")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(OUT / "03_isodepth_on_rgb.png", dpi=160)
    plt.close(fig)

    # --- 1D profiles at 6 sample points ---------------------------------------
    print("sampling 1D profiles at 6 anchor points")
    # Sample along radii from the centroid to interesting boundary points.
    cy, cx = np.array(np.where(mask)).mean(axis=1).astype(int)
    # Pick 6 angles, find first cerebellum-mask-edge along each ray, walk inward.
    H, W = mask.shape
    angles = np.linspace(0, 2 * np.pi, 7)[:-1]  # 6 directions
    fig, axes = plt.subplots(2, 3, figsize=(16, 9), sharex=True, sharey=True)
    axes = axes.ravel()

    # Also remember each anchor for the overview plot.
    anchors = []
    for i, theta in enumerate(angles):
        # walk outward from centroid until we leave mask
        for r in np.arange(1, max(H, W), 2.0):
            y = int(cy + r * np.sin(theta))
            x = int(cx + r * np.cos(theta))
            if y < 0 or y >= H or x < 0 or x >= W or not mask[y, x]:
                # last in-mask pixel = pial boundary along this ray
                y0 = int(cy + (r - 2) * np.sin(theta))
                x0 = int(cx + (r - 2) * np.cos(theta))
                break
        else:
            continue
        anchors.append((y0, x0, theta))

        # Walk inward along the gradient of depth_um (toward higher depth)
        # Sample the 1D intensity profile at increasing depth using nearest-neighbor
        # in pixels along the inward ray.
        ys = np.linspace(y0, cy, 300).astype(int)
        xs = np.linspace(x0, cx, 300).astype(int)
        # Cull to mask
        valid = mask[ys, xs]
        ys = ys[valid]; xs = xs[valid]
        d_profile = depth_um[ys, xs]
        p_profile = p27[ys, xs]
        n_profile = neun[ys, xs]
        b_profile = dapi[ys, xs]
        ax = axes[i]
        ax.plot(d_profile, p_profile, color="red", label="p27")
        ax.plot(d_profile, n_profile, color="green", label="NeuN")
        ax.plot(d_profile, b_profile, color="blue", label="DAPI")
        ax.set_title(f"anchor {i}  (y={y0}, x={x0})")
        ax.set_xlim(0, 300)
        if i == 0:
            ax.legend(loc="upper right", fontsize=8)
        ax.set_xlabel("depth from pia (um)")
        ax.set_ylabel("intensity (CLAHE)")
    fig.suptitle("1D intensity profiles vs depth from pia at 6 anchor points")
    fig.tight_layout()
    fig.savefig(OUT / "04_profiles_at_anchors.png", dpi=120)
    plt.close(fig)

    # Show anchors on the RGB
    fig, ax = plt.subplots(1, 1, figsize=(14, 12))
    ax.imshow(rgb)
    cs = ax.contour(depth_um, levels=levels, colors="white", linewidths=0.4, alpha=0.6)
    for i, (y0, x0, _) in enumerate(anchors):
        ax.plot(x0, y0, "o", color="yellow", markersize=10, markerfacecolor="none",
                markeredgewidth=2)
        ax.text(x0 + 30, y0, str(i), color="yellow", fontsize=14, weight="bold")
    ax.set_title("Anchor points used for 1D profiles")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(OUT / "05_anchor_locations.png", dpi=140)
    plt.close(fig)

    print(f"done. figures in {OUT}")


if __name__ == "__main__":
    main()
