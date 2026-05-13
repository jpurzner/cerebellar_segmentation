"""
v2 diagnostic on the reference image.

Builds:
  - outer cerebellum mask (tighter than v1)
  - IGL mask (port of c_final logic from cerebellum_threshold_segment2.m)
  - pial outline, sampled at fixed arc-length spacing, with local tangent and inward normal
  - "best ray" at each sampled pial point: among +/- ANGLE_SEARCH degrees around the
    local inward normal, pick the direction whose perpendicular ray reaches the IGL with
    the shortest length (and stays inside the cerebellum mask)
  - EGL ridge response from Frangi filter on (DAPI + p27)

Saves figures to figs/v2/.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as ndi
import tifffile
from skimage import exposure, filters, morphology, measure
from skimage.morphology import disk

REF = Path(
    "/Users/jpurzner/Dropbox/images/edu_repeat/p27/"
    "2018_05_22_s1_3_p27-0019/2018_05_22_s1_3_p27-0019_fused_crop.tif"
)
OUT = Path(__file__).resolve().parent / "figs" / "v2"
OUT.mkdir(parents=True, exist_ok=True)

PX_UM = 0.5119049

# Channel order for 20x: ch1=p27, ch2=NeuN, ch3=DAPI (per cerebellum_threshold_segment20x.m)
CH_P27, CH_NEUN, CH_DAPI = 0, 1, 2

# Tuning
PIAL_SAMPLE_UM = 8.0       # spacing between sampled pial points
TANGENT_HALFWIN_UM = 25.0  # half-window for local tangent estimation (arc length)
ANGLE_SEARCH_DEG = 60      # +/- search range around inward normal (B); was 30, hit walls
ANGLE_SEARCH_N = 41        # number of angles to evaluate
RAY_MAX_UM = 600           # don't search rays longer than this
RAY_STEP_UM = 1.5

DRAW_EVERY = 12  # for visualization, draw 1 in N rays (was 25; ~100um spacing now)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def to_unit(img):
    img = img.astype(np.float32)
    lo, hi = img.min(), img.max()
    return (img - lo) / max(hi - lo, 1e-9)


def clahe(img, clip=0.01):
    return exposure.equalize_adapthist(img, clip_limit=clip, nbins=256)


def cerebellum_mask(p27, dapi):
    s = p27 + dapi
    s = s / s.max()
    raw = s > 0.10
    # only a small close to bridge sub-cell-size gaps -- much smaller than v1
    r_close = max(1, int(round(5 / PX_UM)))   # 5 um (was 20)
    closed = morphology.closing(raw, footprint=disk(r_close))
    filled = ndi.binary_fill_holes(closed)
    # keep only large pieces
    lbl, n = ndi.label(filled)
    sizes = np.bincount(lbl.ravel()); sizes[0] = 0
    if n == 0:
        return filled
    keep = sizes.argmax()
    return lbl == keep


def igl_mask(neun, mask):
    """Port of the IGL-detection backbone (NeuN-driven thresh_by_area + cleanup)."""
    # smooth NeuN within mask
    sigma_um = 8.0
    sigma_px = sigma_um / PX_UM
    blurred = filters.gaussian(neun, sigma=sigma_px, preserve_range=True)
    blurred = blurred * mask

    # thresh_by_area (~40% of cerebellum area): pick threshold so 40% of mask area
    # remains.
    inside = blurred[mask]
    thr = np.quantile(inside, 0.6)  # equivalent to "top 40%"
    igl = (blurred > thr) & mask

    # cleanup: close, fill, remove small
    r = max(1, int(round(20 / PX_UM)))
    igl = morphology.closing(igl, disk(r))
    igl = ndi.binary_fill_holes(igl)
    # keep regions larger than 50_000 px^2 (~13 mm^2) -- sized for big IGL only
    lbl, _ = ndi.label(igl)
    if lbl.max() > 0:
        sizes = np.bincount(lbl.ravel()); sizes[0] = 0
        keep_lbl = np.where(sizes >= max(50000, sizes.max() * 0.05))[0]
        igl = np.isin(lbl, keep_lbl)
    return igl


def smooth_outline(mask, smooth_sigma_um=10.0):
    """Smooth mask boundary, return ordered (y, x) outline as Nx2 float array."""
    sigma_px = smooth_sigma_um / PX_UM
    smooth_mask = filters.gaussian(mask.astype(float), sigma=sigma_px) > 0.5
    contours = measure.find_contours(smooth_mask.astype(float), 0.5)
    # take the longest contour
    contours.sort(key=len, reverse=True)
    outline = contours[0]  # (y, x)
    return outline, smooth_mask


def resample_outline(outline, spacing_um):
    """Resample an ordered outline to roughly uniform arc-length spacing in pixels."""
    spacing_px = spacing_um / PX_UM
    seg = np.diff(outline, axis=0)
    arc = np.concatenate([[0], np.cumsum(np.linalg.norm(seg, axis=1))])
    n = max(2, int(round(arc[-1] / spacing_px)))
    targets = np.linspace(0, arc[-1], n)
    ys = np.interp(targets, arc, outline[:, 0])
    xs = np.interp(targets, arc, outline[:, 1])
    return np.stack([ys, xs], axis=1)  # (N, 2) -> (y, x)


def local_tangents(outline, halfwin_um):
    """Local unit tangent per outline point via centered differences over a window."""
    halfwin_px = halfwin_um / PX_UM
    seg = np.diff(outline, axis=0)
    arc = np.concatenate([[0], np.cumsum(np.linalg.norm(seg, axis=1))])
    tangents = np.zeros_like(outline)
    for i in range(len(outline)):
        a = arc[i]
        # find neighbors within halfwin in arc length
        j_lo = np.searchsorted(arc, a - halfwin_px)
        j_hi = np.searchsorted(arc, a + halfwin_px)
        j_lo = max(0, j_lo); j_hi = min(len(outline) - 1, j_hi)
        v = outline[j_hi] - outline[j_lo]
        n = np.linalg.norm(v)
        if n > 0:
            tangents[i] = v / n
    return tangents  # (N, 2) -> (dy, dx)


def inward_normals(outline, tangents, mask):
    """For each outline point, compute the inward unit normal.

    Two perpendiculars exist: (-tx, ty)/(... ) and (tx, -ty)/(... ) wait -- 90deg rotation
    of (ty, tx) -> (-tx, ty) etc. Pick the one that steps INTO the mask.
    """
    H, W = mask.shape
    n_a = np.stack([-tangents[:, 1], tangents[:, 0]], axis=1)  # rotate +90
    n_b = -n_a
    # sample mask 3 px along each candidate
    out = np.zeros_like(outline)
    for i in range(len(outline)):
        py, px = outline[i]
        for n in (n_a[i], n_b[i]):
            ty = int(round(py + 3 * n[0]))
            tx = int(round(px + 3 * n[1]))
            if 0 <= ty < H and 0 <= tx < W and mask[ty, tx]:
                out[i] = n
                break
    return out


def best_ray(start, n0, mask, igl, max_um, step_um, search_deg, search_n):
    """At a pial point with inward normal n0, search +/-search_deg around n0 for the
    angle whose perpendicular ray reaches IGL with the shortest distance, while staying
    in mask. Returns (best_endpoint_y, best_endpoint_x, best_length_um, best_angle_deg)
    or (None, None, np.inf, 0.0) if no angle reaches IGL within max_um.
    """
    H, W = mask.shape
    max_steps = int(max_um / step_um)
    step_px = step_um / PX_UM
    angles = np.deg2rad(np.linspace(-search_deg, search_deg, search_n))
    best = (None, None, np.inf, 0.0)
    for a in angles:
        ca, sa = np.cos(a), np.sin(a)
        # rotate n0 by a
        ny = ca * n0[0] - sa * n0[1]
        nx = sa * n0[0] + ca * n0[1]
        for k in range(1, max_steps + 1):
            y = int(round(start[0] + k * step_px * ny))
            x = int(round(start[1] + k * step_px * nx))
            if not (0 <= y < H and 0 <= x < W) or not mask[y, x]:
                break
            if igl[y, x]:
                length_um = k * step_um
                if length_um < best[2]:
                    best = (y, x, length_um, np.rad2deg(a))
                break
    return best


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    print(f"loading {REF.name}")
    stack = tifffile.imread(REF)
    if stack.ndim != 3:
        raise SystemExit("unexpected stack shape")
    if stack.shape[0] not in (3, 4) and stack.shape[-1] in (3, 4):
        stack = np.moveaxis(stack, -1, 0)
    p27 = clahe(to_unit(stack[CH_P27]))
    neun = clahe(to_unit(stack[CH_NEUN]))
    dapi = clahe(to_unit(stack[CH_DAPI]))
    H, W = p27.shape
    rgb = np.stack([p27, neun, dapi], axis=-1)
    rgb = rgb / max(rgb.max(), 1e-9)
    print(f"  shape={p27.shape}")

    print("cerebellum mask")
    mask = cerebellum_mask(p27, dapi)
    print(f"  mask area = {mask.sum() * (PX_UM/1000)**2:.3f} mm^2")

    print("IGL mask")
    igl = igl_mask(neun, mask)
    print(f"  IGL area  = {igl.sum() * (PX_UM/1000)**2:.3f} mm^2")

    print("pial outline (smoothed) + sample pial points")
    outline, smooth_mask = smooth_outline(mask, smooth_sigma_um=10.0)
    pts = resample_outline(outline, PIAL_SAMPLE_UM)
    print(f"  {len(pts)} sampled pial points (spacing ~{PIAL_SAMPLE_UM} um)")

    print("local tangents + inward normals")
    tangs = local_tangents(pts, TANGENT_HALFWIN_UM)
    norms = inward_normals(pts, tangs, mask)

    print(f"best rays (angle search +/-{ANGLE_SEARCH_DEG}deg, {ANGLE_SEARCH_N} angles)")
    rays = []
    for i, p in enumerate(pts):
        n0 = norms[i]
        if np.linalg.norm(n0) < 0.5:
            rays.append(None)
            continue
        ye, xe, L, a = best_ray(p, n0, mask, igl, RAY_MAX_UM, RAY_STEP_UM,
                                 ANGLE_SEARCH_DEG, ANGLE_SEARCH_N)
        rays.append((ye, xe, L, a))
    hits = [r for r in rays if r is not None and r[0] is not None]
    print(f"  rays that hit IGL: {len(hits)} / {len(rays)}")
    if hits:
        Ls = np.array([r[2] for r in hits])
        print(f"  ray length um: median={np.median(Ls):.0f}  p10={np.percentile(Ls,10):.0f}  p90={np.percentile(Ls,90):.0f}")

    # ---------------- plots ----------------
    print("plot 01: cerebellum mask (tight)")
    fig, ax = plt.subplots(1, 2, figsize=(16, 8))
    ax[0].imshow(rgb)
    ax[0].contour(mask, levels=[0.5], colors=["yellow"], linewidths=0.7)
    ax[0].set_title(f"v2 mask ({mask.sum()*(PX_UM/1000)**2:.2f} mm^2)")
    ax[0].axis("off")
    ax[1].imshow(rgb)
    ax[1].contour(igl, levels=[0.5], colors=["cyan"], linewidths=0.8)
    ax[1].set_title(f"IGL ({igl.sum()*(PX_UM/1000)**2:.2f} mm^2)")
    ax[1].axis("off")
    fig.tight_layout(); fig.savefig(OUT / "01_mask_and_igl.png", dpi=140); plt.close(fig)

    print("plot 02: pial points + tangents + normals")
    fig, ax = plt.subplots(1, 1, figsize=(14, 12))
    ax.imshow(rgb)
    ax.contour(mask, levels=[0.5], colors=["yellow"], linewidths=0.4, alpha=0.6)
    ax.contour(igl,  levels=[0.5], colors=["cyan"],   linewidths=0.6, alpha=0.8)
    every = DRAW_EVERY
    L_PX = 20 / PX_UM
    for i in range(0, len(pts), every):
        py, px = pts[i]
        ny, nx = norms[i]
        ax.plot([px, px + L_PX * nx], [py, py + L_PX * ny], "-", color="white", lw=0.6)
        ax.plot(px, py, ".", color="orange", markersize=2)
    ax.set_title(f"Sampled pial points (every {PIAL_SAMPLE_UM*every:.0f} um) + inward normals")
    ax.axis("off")
    fig.tight_layout(); fig.savefig(OUT / "02_pial_normals.png", dpi=160); plt.close(fig)

    print("plot 03: best perpendicular rays to IGL (dim background + bright rays)")
    rgb_dim = rgb * 0.35
    fig, ax = plt.subplots(1, 1, figsize=(16, 14))
    ax.imshow(rgb_dim)
    ax.contour(igl, levels=[0.5], colors=["cyan"], linewidths=0.7, alpha=0.9)
    ax.contour(mask, levels=[0.5], colors=["yellow"], linewidths=0.4, alpha=0.6)
    for i in range(0, len(pts), every):
        r = rays[i]
        if r is None or r[0] is None:
            continue
        py, px = pts[i]
        ye, xe, L, a = r
        ax.plot([px, xe], [py, ye], "-", color="magenta", lw=1.0, alpha=0.95)
        ax.plot(px, py, ".", color="orange", markersize=3)
        ax.plot(xe, ye, ".", color="lime",   markersize=3)
    ax.set_title(f"Best perpendicular ray to IGL (every {PIAL_SAMPLE_UM*every:.0f} um), "
                 f"angle search +/-{ANGLE_SEARCH_DEG}deg | "
                 f"orange=pia, lime=IGL hit")
    ax.axis("off")
    fig.tight_layout(); fig.savefig(OUT / "03_rays_to_igl.png", dpi=160); plt.close(fig)

    print("plot 03b: where do rays fail / saturate the angle search?")
    fail_pts = np.array([pts[i] for i in range(len(pts))
                         if rays[i] is None or rays[i][0] is None])
    sat_pts  = np.array([pts[i] for i in range(len(pts))
                         if rays[i] is not None and rays[i][0] is not None
                         and abs(rays[i][3]) > ANGLE_SEARCH_DEG - 3])
    fig, ax = plt.subplots(1, 1, figsize=(16, 14))
    ax.imshow(rgb_dim)
    ax.contour(igl, levels=[0.5], colors=["cyan"], linewidths=0.6, alpha=0.7)
    if len(fail_pts):
        ax.plot(fail_pts[:, 1], fail_pts[:, 0], "x", color="red",   markersize=6,
                label=f"no IGL within {RAY_MAX_UM}um (n={len(fail_pts)})")
    if len(sat_pts):
        ax.plot(sat_pts[:, 1],  sat_pts[:, 0],  "o", color="yellow",
                markerfacecolor="none", markersize=8,
                label=f"angle saturated (|a|>{ANGLE_SEARCH_DEG-3}, n={len(sat_pts)})")
    ax.legend(loc="lower right")
    ax.set_title("Where do rays fail or saturate?")
    ax.axis("off")
    fig.tight_layout(); fig.savefig(OUT / "03b_ray_failures.png", dpi=160); plt.close(fig)

    print("plot 04: ray-length distribution + angle-deflection histogram")
    if hits:
        fig, ax = plt.subplots(1, 2, figsize=(14, 5))
        Ls = np.array([r[2] for r in hits])
        As = np.array([r[3] for r in hits])
        ax[0].hist(Ls, bins=60); ax[0].set_xlabel("ray length to IGL (um)")
        ax[0].set_title(f"n={len(Ls)}; median={np.median(Ls):.0f}um")
        ax[1].hist(As, bins=np.linspace(-ANGLE_SEARCH_DEG, ANGLE_SEARCH_DEG, 31))
        ax[1].set_xlabel("angle deflection from local normal (deg)")
        ax[1].set_title("how often was the angle search needed?")
        fig.tight_layout(); fig.savefig(OUT / "04_ray_stats.png", dpi=120); plt.close(fig)

    print("plot 05: Frangi ridge response on heavily-blurred (DAPI+p27) (EGL detector)")
    # The EGL is ~30-50 um wide; we want a ridge filter that responds at THAT scale,
    # not the cell scale. Smooth with sigma >> cell size first to wipe out cellular
    # texture, then run ridge filter at the ribbon half-width range.
    pre_sigma_um = 8.0
    egl_input = filters.gaussian((dapi + p27) / 2,
                                  sigma=pre_sigma_um / PX_UM,
                                  preserve_range=True)
    egl_input = egl_input * mask
    sigmas_um = [10, 15, 20, 25]  # ridge half-widths in um (target EGL ribbon)
    sigmas_px = [s / PX_UM for s in sigmas_um]
    ridge = filters.frangi(egl_input, sigmas=sigmas_px,
                            black_ridges=False, alpha=0.5, beta=0.5)
    ridge = ridge * mask
    ridge_norm = ridge / max(ridge.max(), 1e-9)
    fig, ax = plt.subplots(1, 3, figsize=(20, 7))
    ax[0].imshow(rgb); ax[0].set_title("RGB"); ax[0].axis("off")
    ax[1].imshow(egl_input, cmap="gray")
    ax[1].set_title(f"smoothed (DAPI+p27)/2 sigma {pre_sigma_um}um")
    ax[1].axis("off")
    ax[2].imshow(rgb_dim)
    ridge_thr = ridge_norm > 0.05
    ax[2].contour(ridge_thr, levels=[0.5], colors=["yellow"], linewidths=0.5, alpha=0.9)
    ax[2].set_title(f"Frangi ridges (sigmas {sigmas_um}um) > 0.05, on RGB")
    ax[2].axis("off")
    fig.tight_layout(); fig.savefig(OUT / "05_frangi_ridges.png", dpi=140); plt.close(fig)

    print(f"done. figures in {OUT}")


if __name__ == "__main__":
    main()
