"""
v3.5: try three EGL detectors side-by-side on the same image.

  (1) MATLAB baseline: top 15% of (DAPI+p27) inside cerebellum, masked away from IGL
  (2) Geometric: outer 50um band of the cerebellum mask (distance from background)
  (3) Asymmetric matched filter: like v3 but the negative flank is one-sided
      (one kernel for "negative flank inward (left)", one for "negative flank inward
      (right)"), max over both flank directions and N orientations

Then compare them visually overlaid on the RGB.
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
OUT = Path(__file__).resolve().parent / "figs" / "v35"
OUT.mkdir(parents=True, exist_ok=True)

PX_UM = 0.5119049
CH_P27, CH_NEUN, CH_DAPI = 0, 1, 2

# --- detector params ---
EGL_BAND_UM = 50.0      # geometric: outer band depth
TOP_PCT     = 15.0       # MATLAB baseline: top X% of (DAPI+p27) inside cerebellum

# matched-filter params (asymmetric)
KERNEL_TICK_UM      = 1.0
TEMPLATE_HALFLEN_UM = 60.0
EGL_HALFWIDTH_UM    = 22.0   # positive lobe halfwidth (um) -> ~44um wide
FLANK_WIDTH_UM      = 25.0   # one-sided negative flank width
KERNEL_THICK_UM     = 5.0
N_THETA             = 12
W_DAPI              = 1.0
W_P27               = 1.5
ASYM_THRESH_PCT     = 95     # take top 5% of asym response


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


def igl_mask(neun, mask):
    sigma_px = 8.0 / PX_UM
    blurred = filters.gaussian(neun, sigma=sigma_px, preserve_range=True) * mask
    inside = blurred[mask]
    thr = np.quantile(inside, 0.6)
    igl = (blurred > thr) & mask
    igl = morphology.closing(igl, disk(int(round(20 / PX_UM))))
    igl = ndi.binary_fill_holes(igl)
    lbl, _ = ndi.label(igl)
    if lbl.max() > 0:
        sizes = np.bincount(lbl.ravel()); sizes[0] = 0
        keep = np.where(sizes >= max(50000, sizes.max() * 0.05))[0]
        igl = np.isin(lbl, keep)
    return igl


def make_asym_templates():
    """Build two flank directions:
       template_L: positive lobe centred at 0, negative flank to the LEFT (d < 0)
       template_R: positive lobe centred at 0, negative flank to the RIGHT (d > 0)
       Each is converted to a 2D kernel of shape (perp_um, parallel_um) and zoomed
       to pixel grid.
    """
    half_n = int(TEMPLATE_HALFLEN_UM / KERNEL_TICK_UM)
    d = np.arange(-half_n, half_n + 1) * KERNEL_TICK_UM

    def build(side):
        # DAPI
        dt = np.zeros_like(d, dtype=np.float32)
        dt[np.abs(d) <= EGL_HALFWIDTH_UM] = +1.0
        if side == "L":
            zone = (d < -EGL_HALFWIDTH_UM) & (d >= -EGL_HALFWIDTH_UM - FLANK_WIDTH_UM)
        else:
            zone = (d > +EGL_HALFWIDTH_UM) & (d <= +EGL_HALFWIDTH_UM + FLANK_WIDTH_UM)
        dt[zone] = -1.5  # extra weight on the one flank since we have only one
        dt -= dt.mean(); dt /= max(np.sqrt((dt ** 2).sum()), 1e-9)
        # p27 (peaked at iEGL position) -- iEGL sits inward from oEGL; "inward" is
        # opposite of the flank direction. If side=L (flank on left), inward=right.
        offset = +8.0 if side == "L" else -8.0
        pt = np.zeros_like(d, dtype=np.float32)
        pt[(d >= offset - 12) & (d <= offset + 12)] = +1.0
        if side == "L":
            zone = (d < -EGL_HALFWIDTH_UM) & (d >= -EGL_HALFWIDTH_UM - FLANK_WIDTH_UM)
        else:
            zone = (d > +EGL_HALFWIDTH_UM) & (d <= +EGL_HALFWIDTH_UM + FLANK_WIDTH_UM)
        pt[zone] = -1.0
        pt -= pt.mean(); pt /= max(np.sqrt((pt ** 2).sum()), 1e-9)
        return dt, pt

    L_dapi, L_p27 = build("L")
    R_dapi, R_p27 = build("R")

    thick = max(1, int(round(KERNEL_THICK_UM / KERNEL_TICK_UM)))

    def to_pix(t):
        K_um = np.tile(t.reshape(-1, 1), (1, thick))
        K = ndi.zoom(K_um, (1.0 / PX_UM, 1.0 / PX_UM), order=1).astype(np.float32)
        # odd shape
        if K.shape[0] % 2 == 0: K = K[:-1, :]
        if K.shape[1] % 2 == 0: K = K[:, :-1]
        return K

    return d, (L_dapi, L_p27, R_dapi, R_p27), {
        "L_dapi": to_pix(L_dapi), "L_p27": to_pix(L_p27),
        "R_dapi": to_pix(R_dapi), "R_p27": to_pix(R_p27),
    }


def best_response(dapi_smooth, p27_smooth, ker):
    """Take max over orientations (0..pi) and flank-direction (L/R) of
       W_DAPI*conv(dapi, K_dapi) + W_P27*conv(p27, K_p27)."""
    H, W = dapi_smooth.shape
    best = np.full((H, W), -np.inf, dtype=np.float32)
    dapi32 = dapi_smooth.astype(np.float32)
    p2732  = p27_smooth.astype(np.float32)
    for side in ("L", "R"):
        Kd = ker[f"{side}_dapi"]
        Kp = ker[f"{side}_p27"]
        for i in range(N_THETA):
            theta_deg = i * 180.0 / N_THETA
            Kd_r = ndi.rotate(Kd, theta_deg, order=1, reshape=True,
                              mode="constant", cval=0.0).astype(np.float32)
            Kp_r = ndi.rotate(Kp, theta_deg, order=1, reshape=True,
                              mode="constant", cval=0.0).astype(np.float32)
            Kd_r -= Kd_r.mean(); Kp_r -= Kp_r.mean()
            rd = fftconvolve(dapi32, Kd_r[::-1, ::-1], mode="same")
            rp = fftconvolve(p2732,  Kp_r[::-1, ::-1], mode="same")
            r = W_DAPI * rd + W_P27 * rp
            np.maximum(best, r, out=best)
    return best


def main():
    print(f"loading {REF.name}")
    stack = tifffile.imread(REF)
    if stack.shape[0] not in (3, 4) and stack.shape[-1] in (3, 4):
        stack = np.moveaxis(stack, -1, 0)
    p27 = clahe(to_unit(stack[CH_P27]))
    neun = clahe(to_unit(stack[CH_NEUN]))
    dapi = clahe(to_unit(stack[CH_DAPI]))
    rgb = np.stack([p27, neun, dapi], axis=-1); rgb = rgb / max(rgb.max(), 1e-9)

    print("masks")
    mask = cerebellum_mask(p27, dapi)
    igl  = igl_mask(neun, mask)

    # --- detector 1: MATLAB-style top 15% of (DAPI+p27), masked away from IGL ---
    print("detector 1: MATLAB-style top {0}%".format(TOP_PCT))
    blurred = filters.gaussian((dapi + p27) / 2, sigma=4.0 / PX_UM,
                               preserve_range=True)
    egl_zone_1 = mask & ~ndi.binary_dilation(igl, structure=disk(int(20/PX_UM)))
    inside = blurred[egl_zone_1]
    thr = np.percentile(inside, 100 - TOP_PCT)
    egl1 = (blurred > thr) & egl_zone_1
    egl1 = morphology.remove_small_objects(egl1, min_size=int(50/PX_UM**2))

    # --- detector 2: geometric outer band ---
    print(f"detector 2: outer {EGL_BAND_UM}um band")
    bg = ~mask
    dist_to_bg_um = ndi.distance_transform_edt(mask) * PX_UM
    egl2 = mask & (dist_to_bg_um <= EGL_BAND_UM)

    # --- detector 3: asymmetric matched filter ---
    print("detector 3: asymmetric matched filter (24 kernels via FFT)")
    sigma_pre_um = 4.0
    dapi_sm = filters.gaussian(dapi, sigma=sigma_pre_um/PX_UM, preserve_range=True) * mask
    p27_sm  = filters.gaussian(p27,  sigma=sigma_pre_um/PX_UM, preserve_range=True) * mask
    d_um, ts1d, ker = make_asym_templates()
    response = best_response(dapi_sm, p27_sm, ker) * mask
    inside = response[mask]
    thr = np.percentile(inside, ASYM_THRESH_PCT)
    egl3 = (response > thr) & mask
    egl3 = morphology.remove_small_objects(egl3, min_size=int(50/PX_UM**2))
    egl3 = morphology.closing(egl3, disk(int(round(3/PX_UM))))

    # --- plots ---
    print("plot 00: 1D templates")
    Ld, Lp, Rd, Rp = ts1d
    fig, ax = plt.subplots(1, 2, figsize=(14, 4))
    ax[0].plot(d_um, Ld, color="blue", label="DAPI L (flank on left)")
    ax[0].plot(d_um, Lp, color="red",  label="p27  L (peak right of centre)")
    ax[0].axhline(0, color="gray", lw=0.5); ax[0].axvline(0, color="gray", lw=0.5)
    ax[0].legend(); ax[0].set_title("L-side template (negative flank d<0)")
    ax[1].plot(d_um, Rd, color="blue", label="DAPI R (flank on right)")
    ax[1].plot(d_um, Rp, color="red",  label="p27  R (peak left of centre)")
    ax[1].axhline(0, color="gray", lw=0.5); ax[1].axvline(0, color="gray", lw=0.5)
    ax[1].legend(); ax[1].set_title("R-side template (negative flank d>0)")
    fig.tight_layout(); fig.savefig(OUT / "00_templates_asym.png", dpi=120); plt.close(fig)

    print("plot 01: side-by-side comparison")
    fig, ax = plt.subplots(1, 4, figsize=(28, 8))
    rgb_dim = rgb * 0.45
    titles = [
        f"(1) MATLAB-style top {TOP_PCT}% on (DAPI+p27)",
        f"(2) outer {EGL_BAND_UM}um geometric band",
        f"(3) asymmetric matched filter (top {100-ASYM_THRESH_PCT}%)",
        "RGB reference",
    ]
    masks_ = [egl1, egl2, egl3, None]
    for i, (m, t) in enumerate(zip(masks_, titles)):
        a = ax[i]
        a.imshow(rgb_dim if m is not None else rgb)
        if m is not None:
            a.contour(m, levels=[0.5], colors=["yellow"], linewidths=0.5)
        a.contour(igl, levels=[0.5], colors=["cyan"], linewidths=0.4, alpha=0.7)
        a.set_title(t); a.axis("off")
    fig.tight_layout(); fig.savefig(OUT / "01_egl_comparison.png", dpi=130)
    plt.close(fig)

    # --- zoom into a deep fold ---
    print("plot 02: zoom into a deep fold for closer look")
    # Pick a window in the middle of the cerebellum (where folds are deepest).
    # Find centroid of mask, take a 1500 px window centred there.
    cy, cx = np.array(np.where(mask)).mean(axis=1).astype(int)
    half = 750
    y0, y1 = max(0, cy - half), min(rgb.shape[0], cy + half)
    x0, x1 = max(0, cx - half), min(rgb.shape[1], cx + half)
    fig, ax = plt.subplots(1, 4, figsize=(28, 8))
    rgb_z = rgb[y0:y1, x0:x1]
    rgb_zd = rgb_z * 0.45
    masks_z = [egl1[y0:y1, x0:x1], egl2[y0:y1, x0:x1], egl3[y0:y1, x0:x1], None]
    for i, (m, t) in enumerate(zip(masks_z, titles)):
        a = ax[i]
        a.imshow(rgb_zd if m is not None else rgb_z)
        if m is not None:
            a.contour(m, levels=[0.5], colors=["yellow"], linewidths=0.7)
        a.contour(igl[y0:y1, x0:x1], levels=[0.5], colors=["cyan"],
                   linewidths=0.5, alpha=0.7)
        a.set_title(t + " (zoom)"); a.axis("off")
    fig.tight_layout(); fig.savefig(OUT / "02_egl_zoom.png", dpi=140)
    plt.close(fig)

    print(f"done. figures in {OUT}")


if __name__ == "__main__":
    main()
