"""
v3 diagnostic: multivariate matched-filter EGL detector.

Cross-section template (perpendicular to the EGL ribbon, going from outside the
cerebellum inward):

    background  |  oEGL  |  iEGL  |  ML  | (deeper)
    DAPI:  low      HIGH    HIGH    low    low
    p27:   low      med     HIGH    low    low

We build:
  - dapi_template(d): zero-mean profile, positive on EGL, negative on flanks
  - p27_template(d):  zero-mean, peaked toward iEGL position (offset inward)

Then the 2D matched-filter response at orientation theta is

    R(x, y, theta) = sum_d [ DAPI(x + d sinT, y - d cosT) * dapi_template(d)
                           + p27 (x + d sinT, y - d cosT) * p27_template(d) ]

We compute this for N_THETA orientations (0 .. pi) and take the per-pixel max.
The negative lobes in the templates do the "flanking" enforcement: a uniform
bright region won't score high because the flanks cancel the center.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as ndi
import tifffile
from skimage import exposure, filters, morphology
from skimage.morphology import disk

REF = Path(
    "/Users/jpurzner/Dropbox/images/edu_repeat/p27/"
    "2018_05_22_s1_3_p27-0019/2018_05_22_s1_3_p27-0019_fused_crop.tif"
)
OUT = Path(__file__).resolve().parent / "figs" / "v3"
OUT.mkdir(parents=True, exist_ok=True)

PX_UM = 0.5119049
CH_P27, CH_NEUN, CH_DAPI = 0, 1, 2

# Cross-section template (in micrometers, perpendicular to ribbon)
# Coordinates: d=0 at the centre of the EGL; positive d = inward (toward ML)
# The template extends from d_min to d_max; lobes get carved out below.
TEMPLATE_HALFLEN_UM = 60.0    # template span: -60 .. +60 um (so ~120um window)
EGL_HALFWIDTH_UM    = 22.0    # EGL is roughly 30-50 um wide; positive lobe halfwidth
FLANK_WIDTH_UM      = 25.0    # negative-flank width on each side

P27_OFFSET_UM       = 8.0     # iEGL sits ~8 um inward of EGL centre
P27_PEAK_HALFWIDTH  = 12.0    # p27 peak is sharper

KERNEL_TICK_UM      = 1.0     # sample template every 1 um
KERNEL_THICK_UM     = 5.0     # along-ribbon kernel thickness
N_THETA             = 12      # orientations 0..pi

# weighting between DAPI and p27 in combined response
W_DAPI = 1.0
W_P27  = 1.5

EGL_THRESH_PCT      = 96      # take top 4% of inside-cerebellum response as EGL


# --------------------------------------------------------------------------
def to_unit(img):
    img = img.astype(np.float32)
    lo, hi = img.min(), img.max()
    return (img - lo) / max(hi - lo, 1e-9)


def clahe(img, clip=0.01):
    return exposure.equalize_adapthist(img, clip_limit=clip, nbins=256)


def cerebellum_mask(p27, dapi):
    s = (p27 + dapi) / 2
    s = s / s.max()
    raw = s > 0.10
    r_close = max(1, int(round(5 / PX_UM)))
    closed = morphology.closing(raw, footprint=disk(r_close))
    filled = ndi.binary_fill_holes(closed)
    lbl, n = ndi.label(filled)
    sizes = np.bincount(lbl.ravel()); sizes[0] = 0
    if n == 0:
        return filled
    return lbl == sizes.argmax()


def make_templates():
    """Return (d_um, dapi_template, p27_template, kernel_2d_dapi, kernel_2d_p27).

    Templates are zero-mean 1D profiles indexed by perpendicular position d in um.
    """
    half_n = int(TEMPLATE_HALFLEN_UM / KERNEL_TICK_UM)
    d_um = np.arange(-half_n, half_n + 1) * KERNEL_TICK_UM
    # DAPI: positive inside |d| < EGL_HALFWIDTH; negative on flanks beyond that
    dapi_t = np.zeros_like(d_um, dtype=np.float32)
    dapi_t[np.abs(d_um) <= EGL_HALFWIDTH_UM] = +1.0
    flank_zone = (np.abs(d_um) > EGL_HALFWIDTH_UM) & \
                 (np.abs(d_um) <= EGL_HALFWIDTH_UM + FLANK_WIDTH_UM)
    dapi_t[flank_zone] = -1.0
    # zero-mean
    dapi_t -= dapi_t.mean()
    # normalize so sum of squares = 1 (matched-filter convention)
    dapi_t /= max(np.sqrt((dapi_t ** 2).sum()), 1e-9)

    # p27: peaked at +P27_OFFSET (inward = iEGL side), narrower
    p27_t = np.zeros_like(d_um, dtype=np.float32)
    p27_peak = (d_um >= P27_OFFSET_UM - P27_PEAK_HALFWIDTH) & \
               (d_um <= P27_OFFSET_UM + P27_PEAK_HALFWIDTH)
    p27_t[p27_peak] = +1.0
    # negative flanks
    p27_neg = ((d_um < P27_OFFSET_UM - P27_PEAK_HALFWIDTH) & \
               (d_um > P27_OFFSET_UM - P27_PEAK_HALFWIDTH - FLANK_WIDTH_UM)) | \
              ((d_um > P27_OFFSET_UM + P27_PEAK_HALFWIDTH) & \
               (d_um < P27_OFFSET_UM + P27_PEAK_HALFWIDTH + FLANK_WIDTH_UM))
    p27_t[p27_neg] = -0.5
    p27_t -= p27_t.mean()
    p27_t /= max(np.sqrt((p27_t ** 2).sum()), 1e-9)

    # 2D kernel = template along the perpendicular axis, smoothed/replicated along
    # the parallel axis. We make a thin (KERNEL_THICK_UM wide) 2D version oriented
    # perpendicular = vertical (will be rotated later).
    thick_n = int(round(KERNEL_THICK_UM / KERNEL_TICK_UM))
    if thick_n < 1: thick_n = 1
    # Convert "1 um per template tick" to pixels at our PX_UM
    # We'll build the kernel in template-tick coordinates first, then resize to px.
    kernel_dapi_um = np.tile(dapi_t.reshape(-1, 1), (1, thick_n))
    kernel_p27_um  = np.tile(p27_t.reshape(-1, 1),  (1, thick_n))

    # Resample from um-grid to pixel-grid using zoom (1 um -> 1/PX_UM px)
    zoom = 1.0 / PX_UM
    k_dapi = ndi.zoom(kernel_dapi_um, (zoom, zoom), order=1)
    k_p27  = ndi.zoom(kernel_p27_um,  (zoom, zoom), order=1)
    # Force odd shape for nice rotation pivots
    def make_odd(K):
        h, w = K.shape
        if h % 2 == 0: K = K[:-1, :]
        if w % 2 == 0: K = K[:, :-1]
        return K
    k_dapi = make_odd(k_dapi)
    k_p27  = make_odd(k_p27)
    return d_um, dapi_t, p27_t, k_dapi, k_p27


def matched_filter_response(img, kernel, n_theta=N_THETA):
    """Apply oriented matched filter at n_theta angles, return per-pixel max.

    Uses scipy.signal.fftconvolve under the hood (~200x faster than direct
    correlate for our ~250x250 kernel).
    """
    from scipy.signal import fftconvolve
    img = img.astype(np.float32)
    H, W = img.shape
    best = np.full((H, W), -np.inf, dtype=np.float32)
    for k in range(n_theta):
        theta_deg = k * 180.0 / n_theta
        K = ndi.rotate(kernel, theta_deg, order=1, reshape=True, mode="constant",
                       cval=0.0).astype(np.float32)
        K = K - K.mean()
        # correlation = convolution with flipped kernel
        Kf = K[::-1, ::-1]
        r = fftconvolve(img, Kf, mode="same")
        np.maximum(best, r, out=best)
    return best


# --------------------------------------------------------------------------
def main():
    print(f"loading {REF.name}")
    stack = tifffile.imread(REF)
    if stack.shape[0] not in (3, 4) and stack.shape[-1] in (3, 4):
        stack = np.moveaxis(stack, -1, 0)
    p27 = clahe(to_unit(stack[CH_P27]))
    neun = clahe(to_unit(stack[CH_NEUN]))
    dapi = clahe(to_unit(stack[CH_DAPI]))
    rgb = np.stack([p27, neun, dapi], axis=-1); rgb = rgb / max(rgb.max(), 1e-9)

    print("cerebellum mask")
    mask = cerebellum_mask(p27, dapi)

    print("smoothing channels for kernel input (single-cell texture suppression)")
    sigma_um_pre = 4.0
    sigma_px_pre = sigma_um_pre / PX_UM
    dapi_smooth = filters.gaussian(dapi, sigma=sigma_px_pre, preserve_range=True)
    p27_smooth  = filters.gaussian(p27,  sigma=sigma_px_pre, preserve_range=True)
    dapi_smooth = dapi_smooth * mask
    p27_smooth  = p27_smooth  * mask

    print("building 1D + 2D templates")
    d_um, dapi_t, p27_t, k_dapi, k_p27 = make_templates()
    print(f"  template span = [{d_um.min()}, {d_um.max()}] um (n={len(d_um)})")
    print(f"  k_dapi shape  = {k_dapi.shape}  k_p27 shape = {k_p27.shape}")

    # save template plot
    fig, ax = plt.subplots(1, 2, figsize=(12, 4))
    ax[0].plot(d_um, dapi_t, color="blue", label="DAPI template")
    ax[0].plot(d_um, p27_t,  color="red",  label="p27 template")
    ax[0].axhline(0, color="gray", lw=0.5)
    ax[0].axvline(0, color="gray", lw=0.5)
    ax[0].set_xlabel("d (um, perpendicular to ribbon)  -- 0=EGL centre, +inward")
    ax[0].set_ylabel("template weight (zero-mean)")
    ax[0].legend(); ax[0].set_title("1D matched-filter templates")
    ax[1].imshow(np.stack([k_dapi, np.zeros_like(k_dapi), -k_dapi], axis=-1) * 0.5 + 0.5,
                 extent=[0, k_dapi.shape[1] * PX_UM, 0, k_dapi.shape[0] * PX_UM])
    ax[1].set_title(f"2D DAPI kernel (red=+, blue=-)  size {k_dapi.shape}")
    fig.tight_layout(); fig.savefig(OUT / "00_templates.png", dpi=120); plt.close(fig)

    print(f"computing matched-filter response at {N_THETA} orientations (slow)")
    print("  DAPI channel")
    r_dapi = matched_filter_response(dapi_smooth, k_dapi)
    print("  p27 channel")
    r_p27  = matched_filter_response(p27_smooth,  k_p27)
    response = W_DAPI * r_dapi + W_P27 * r_p27
    response = response * mask

    # Threshold inside cerebellum
    inside = response[mask]
    thr = np.percentile(inside, EGL_THRESH_PCT)
    egl = (response > thr) & mask
    # cleanup
    egl_clean = morphology.remove_small_objects(egl, min_size=int(50 / PX_UM ** 2))
    egl_clean = morphology.closing(egl_clean, disk(int(round(3 / PX_UM))))

    print(f"  response: thr={thr:.4f}, EGL-pixel fraction = "
          f"{egl.sum()/mask.sum():.3%}  (cleaned={egl_clean.sum()/mask.sum():.3%})")

    # ---------------- plots ----------------
    print("plot 01: kernel response heatmap")
    fig, ax = plt.subplots(1, 3, figsize=(20, 7))
    ax[0].imshow(rgb); ax[0].set_title("RGB"); ax[0].axis("off")
    ax[1].imshow(np.where(mask, response, np.nan), cmap="hot")
    ax[1].set_title(f"matched-filter response (max over {N_THETA} thetas, "
                    f"{W_DAPI}*DAPI + {W_P27}*p27)")
    ax[1].axis("off")
    ax[2].imshow(rgb * 0.4)
    ax[2].contour(egl_clean, levels=[0.5], colors=["yellow"], linewidths=0.6)
    ax[2].set_title(f"EGL mask (top {100-EGL_THRESH_PCT}% of response, cleaned)")
    ax[2].axis("off")
    fig.tight_layout(); fig.savefig(OUT / "01_response_and_mask.png", dpi=140)
    plt.close(fig)

    print("plot 02: response colored by best orientation (sanity check)")
    # rerun and remember argmax theta (FFT-convolve)
    from scipy.signal import fftconvolve
    H, W = dapi_smooth.shape
    best = np.full((H, W), -np.inf, dtype=np.float32)
    best_theta = np.zeros((H, W), dtype=np.float32)
    for k in range(N_THETA):
        theta_deg = k * 180.0 / N_THETA
        Kd = ndi.rotate(k_dapi, theta_deg, order=1, reshape=True, mode="constant",
                        cval=0.0).astype(np.float32)
        Kp = ndi.rotate(k_p27,  theta_deg, order=1, reshape=True, mode="constant",
                        cval=0.0).astype(np.float32)
        Kd -= Kd.mean(); Kp -= Kp.mean()
        rd = fftconvolve(dapi_smooth.astype(np.float32), Kd[::-1, ::-1], mode="same")
        rp = fftconvolve(p27_smooth.astype(np.float32),  Kp[::-1, ::-1], mode="same")
        r = W_DAPI * rd + W_P27 * rp
        upd = r > best
        best_theta = np.where(upd, theta_deg, best_theta)
        best = np.where(upd, r, best)
    fig, ax = plt.subplots(1, 1, figsize=(14, 12))
    ax.imshow(rgb * 0.4)
    show = np.where(egl_clean, best_theta, np.nan)
    im = ax.imshow(show, cmap="hsv", vmin=0, vmax=180, alpha=0.85)
    ax.set_title("Best orientation per EGL pixel (HSV cycle 0..180deg)")
    ax.axis("off")
    plt.colorbar(im, ax=ax, label="theta (deg)", fraction=0.03)
    fig.tight_layout(); fig.savefig(OUT / "02_best_orientation.png", dpi=140)
    plt.close(fig)

    print(f"done. figures in {OUT}")


if __name__ == "__main__":
    main()
