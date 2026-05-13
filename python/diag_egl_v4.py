"""
v4: deep-EGL-pair kernel.

Cross-section template of a paired EGL in a deep fold (perpendicular to the fold's
long axis):

    position d:    ML | EGL  | gap | EGL  | ML
    DAPI:          -    +     -     +     -
    p27:           -    +     -     +     -

Total template span ~120 um. Two positive lobes flanking a central negative gap,
with negative flanks beyond on both sides.

This kernel SHOULD fire only where the paired structure exists (deep folds), and
NOT on isolated EGL ribbons (gyral crowns) because there's no second positive
lobe to support the response.

We then union this with one of the single-ribbon detectors from v3.5 to get
full coverage.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as ndi
import tifffile
from scipy.signal import fftconvolve
from skimage import exposure, filters, morphology
from skimage.morphology import disk

REF = Path(
    "/Users/jpurzner/Dropbox/images/edu_repeat/p27/"
    "2018_05_22_s1_3_p27-0019/2018_05_22_s1_3_p27-0019_fused_crop.tif"
)
OUT = Path(__file__).resolve().parent / "figs" / "v4"
OUT.mkdir(parents=True, exist_ok=True)

PX_UM = 0.5119049
CH_P27, CH_NEUN, CH_DAPI = 0, 1, 2

# --- pair-template params ---
KERNEL_TICK_UM      = 1.0
PAIR_HALFLEN_UM     = 70.0   # template span: -70 .. +70 um
GAP_HALFWIDTH_UM    = 12.0   # central gap halfwidth (LCSF gap is ~5-15 um)
LOBE_HALFWIDTH_UM   = 18.0   # each EGL lobe half-width (so each lobe ~36 um wide)
ML_FLANK_WIDTH_UM   = 20.0   # outer ML negative flank width
KERNEL_THICK_UM     = 5.0    # along-ribbon kernel thickness
N_THETA             = 12

W_DAPI = 1.0
W_P27  = 1.5

PAIR_THRESH_PCT  = 95   # take top 5% inside cerebellum as pair-EGL
SINGLE_THRESH_PCT = 90  # for the single-ribbon detector (just blurred top-X%)


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
    if n == 0: return filled
    return lbl == sizes.argmax()


def make_pair_template():
    """Pair template: positive on lobes, negative on gap and outer flanks.

    Returns:
        d_um (1D)
        dapi_template (1D, zero-mean unit-norm)
        p27_template  (1D, zero-mean unit-norm)
        K_dapi (2D pixel-grid kernel)
        K_p27  (2D pixel-grid kernel)
    """
    half_n = int(PAIR_HALFLEN_UM / KERNEL_TICK_UM)
    d = np.arange(-half_n, half_n + 1) * KERNEL_TICK_UM

    def build_axis(lobe_pos, lobe_w, gap_w, flank_w):
        """Build template centred at d=0 with lobes at d=+/-lobe_pos."""
        t = np.zeros_like(d, dtype=np.float32)
        # central gap (negative)
        t[np.abs(d) <= gap_w] = -1.0
        # lobes at +/- lobe_pos (positive)
        for sign in (-1, +1):
            inner = sign * lobe_pos - lobe_w
            outer = sign * lobe_pos + lobe_w
            lo, hi = (inner, outer) if inner < outer else (outer, inner)
            t[(d >= lo) & (d <= hi)] = +1.0
        # outer flanks (negative)
        for sign in (-1, +1):
            base = sign * (lobe_pos + lobe_w)
            if sign > 0:
                t[(d > base) & (d <= base + flank_w)] = -0.7
            else:
                t[(d < base) & (d >= base - flank_w)] = -0.7
        # zero-mean + unit norm
        t -= t.mean()
        t /= max(np.sqrt((t ** 2).sum()), 1e-9)
        return t

    # Lobe centre at GAP_HALFWIDTH + LOBE_HALFWIDTH so each lobe sits just outside the gap
    lobe_pos = GAP_HALFWIDTH_UM + LOBE_HALFWIDTH_UM
    dapi_t = build_axis(lobe_pos, LOBE_HALFWIDTH_UM, GAP_HALFWIDTH_UM, ML_FLANK_WIDTH_UM)
    p27_t  = build_axis(lobe_pos, LOBE_HALFWIDTH_UM - 4, GAP_HALFWIDTH_UM,
                        ML_FLANK_WIDTH_UM)  # narrower lobes for p27 (iEGL only)

    thick = max(1, int(round(KERNEL_THICK_UM / KERNEL_TICK_UM)))

    def to_pix(t):
        K_um = np.tile(t.reshape(-1, 1), (1, thick))
        K = ndi.zoom(K_um, (1.0/PX_UM, 1.0/PX_UM), order=1).astype(np.float32)
        if K.shape[0] % 2 == 0: K = K[:-1, :]
        if K.shape[1] % 2 == 0: K = K[:, :-1]
        return K

    return d, dapi_t, p27_t, to_pix(dapi_t), to_pix(p27_t)


def matched_pair_response(dapi_smooth, p27_smooth, K_dapi, K_p27):
    H, W = dapi_smooth.shape
    best = np.full((H, W), -np.inf, dtype=np.float32)
    dapi32 = dapi_smooth.astype(np.float32)
    p2732  = p27_smooth.astype(np.float32)
    for i in range(N_THETA):
        theta_deg = i * 180.0 / N_THETA
        Kd = ndi.rotate(K_dapi, theta_deg, order=1, reshape=True,
                        mode="constant", cval=0.0).astype(np.float32)
        Kp = ndi.rotate(K_p27,  theta_deg, order=1, reshape=True,
                        mode="constant", cval=0.0).astype(np.float32)
        Kd -= Kd.mean(); Kp -= Kp.mean()
        rd = fftconvolve(dapi32, Kd[::-1, ::-1], mode="same")
        rp = fftconvolve(p2732,  Kp[::-1, ::-1], mode="same")
        r = W_DAPI * rd + W_P27 * rp
        np.maximum(best, r, out=best)
    return best


def single_ribbon_top_pct(dapi, p27, mask, pct):
    """Cheap baseline: top X% of (DAPI+p27) inside cerebellum."""
    blurred = filters.gaussian((dapi + p27) / 2, sigma=4.0/PX_UM, preserve_range=True)
    inside = blurred[mask]
    thr = np.percentile(inside, 100 - pct)
    egl = (blurred > thr) & mask
    egl = morphology.remove_small_objects(egl, min_size=int(50/PX_UM**2))
    return egl


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

    print("pre-smoothing channels")
    pre_um = 4.0
    dapi_sm = filters.gaussian(dapi, sigma=pre_um/PX_UM, preserve_range=True)
    p27_sm  = filters.gaussian(p27,  sigma=pre_um/PX_UM, preserve_range=True)

    print("building pair template")
    d_um, dapi_t, p27_t, K_dapi, K_p27 = make_pair_template()
    print(f"  template span [{d_um.min()}, {d_um.max()}]um, kernel size {K_dapi.shape}")

    print("plot 00: pair template")
    fig, ax = plt.subplots(1, 2, figsize=(14, 4))
    ax[0].plot(d_um, dapi_t, color="blue", label="DAPI pair template")
    ax[0].plot(d_um, p27_t,  color="red",  label="p27  pair template")
    ax[0].axhline(0, color="gray", lw=0.5)
    for d in [-(GAP_HALFWIDTH_UM+LOBE_HALFWIDTH_UM), -GAP_HALFWIDTH_UM,
               GAP_HALFWIDTH_UM,  GAP_HALFWIDTH_UM+LOBE_HALFWIDTH_UM]:
        ax[0].axvline(d, color="gray", lw=0.3, linestyle=":")
    ax[0].set_xlabel("d (um, perpendicular to fold)  -- 0=gap centre")
    ax[0].legend(); ax[0].set_title("Deep-EGL-pair 1D templates")
    K_disp = np.stack([np.maximum(K_dapi,0)*1.5,
                       np.zeros_like(K_dapi),
                       np.maximum(-K_dapi,0)*1.5], axis=-1)
    K_disp = np.clip(K_disp, 0, 1)
    ax[1].imshow(K_disp); ax[1].set_title(f"2D DAPI pair kernel (red=+, blue=-) {K_dapi.shape}")
    fig.tight_layout(); fig.savefig(OUT / "00_pair_template.png", dpi=120); plt.close(fig)

    print(f"computing pair response at {N_THETA} orientations (FFT)")
    pair_resp = matched_pair_response(dapi_sm * mask, p27_sm * mask, K_dapi, K_p27)
    pair_resp = pair_resp * mask

    print("thresholding pair response")
    inside = pair_resp[mask]
    thr_pair = np.percentile(inside, PAIR_THRESH_PCT)
    egl_pair = (pair_resp > thr_pair) & mask
    egl_pair = morphology.remove_small_objects(egl_pair, min_size=int(80/PX_UM**2))

    print("single-ribbon top-pct baseline (for union)")
    egl_single = single_ribbon_top_pct(dapi, p27, mask, SINGLE_THRESH_PCT)

    egl_union = egl_pair | egl_single

    # --- plots ---
    print("plot 01: full-image comparison")
    rgb_dim = rgb * 0.45
    fig, ax = plt.subplots(1, 4, figsize=(28, 8))
    ax[0].imshow(rgb); ax[0].set_title("RGB"); ax[0].axis("off")
    ax[1].imshow(np.where(mask, pair_resp, np.nan), cmap="hot")
    ax[1].set_title("pair-kernel response")
    ax[1].axis("off")
    ax[2].imshow(rgb_dim)
    ax[2].contour(egl_pair, levels=[0.5], colors=["yellow"], linewidths=0.5)
    ax[2].set_title(f"deep-pair EGL only (top {100-PAIR_THRESH_PCT}%)")
    ax[2].axis("off")
    ax[3].imshow(rgb_dim)
    ax[3].contour(egl_pair,   levels=[0.5], colors=["red"],    linewidths=0.5,
                   alpha=0.85)
    ax[3].contour(egl_single, levels=[0.5], colors=["yellow"], linewidths=0.4,
                   alpha=0.7)
    ax[3].set_title("union: red=pair, yellow=single (top 10% intensity)")
    ax[3].axis("off")
    fig.tight_layout(); fig.savefig(OUT / "01_pair_full.png", dpi=130); plt.close(fig)

    # --- zoom ---
    print("plot 02: zoom into deep folds")
    cy, cx = np.array(np.where(mask)).mean(axis=1).astype(int)
    half = 750
    y0, y1 = max(0, cy - half), min(rgb.shape[0], cy + half)
    x0, x1 = max(0, cx - half), min(rgb.shape[1], cx + half)
    fig, ax = plt.subplots(1, 4, figsize=(28, 8))
    ax[0].imshow(rgb[y0:y1, x0:x1]); ax[0].set_title("RGB (zoom)"); ax[0].axis("off")
    ax[1].imshow(np.where(mask[y0:y1, x0:x1], pair_resp[y0:y1, x0:x1], np.nan),
                  cmap="hot"); ax[1].set_title("pair response (zoom)")
    ax[1].axis("off")
    rgb_zd = rgb[y0:y1, x0:x1] * 0.45
    ax[2].imshow(rgb_zd)
    ax[2].contour(egl_pair[y0:y1, x0:x1], levels=[0.5], colors=["yellow"], linewidths=0.7)
    ax[2].set_title(f"deep-pair EGL (zoom)"); ax[2].axis("off")
    ax[3].imshow(rgb_zd)
    ax[3].contour(egl_pair[y0:y1, x0:x1],   levels=[0.5], colors=["red"],    linewidths=0.7)
    ax[3].contour(egl_single[y0:y1, x0:x1], levels=[0.5], colors=["yellow"], linewidths=0.5,
                   alpha=0.7)
    ax[3].set_title("union (zoom): red=pair, yellow=single"); ax[3].axis("off")
    fig.tight_layout(); fig.savefig(OUT / "02_pair_zoom.png", dpi=140); plt.close(fig)

    print(f"\npair-response stats: thr={thr_pair:.4f}")
    print(f"  egl_pair area:   {egl_pair.sum() * (PX_UM/1000)**2:.3f} mm^2")
    print(f"  egl_single area: {egl_single.sum() * (PX_UM/1000)**2:.3f} mm^2")
    print(f"  union area:      {egl_union.sum() * (PX_UM/1000)**2:.3f} mm^2")
    print(f"  pair only (no overlap with single): "
          f"{(egl_pair & ~egl_single).sum() * (PX_UM/1000)**2:.3f} mm^2")

    print(f"done. figures in {OUT}")


if __name__ == "__main__":
    main()
