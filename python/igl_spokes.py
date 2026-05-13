"""IGL-anchored spokes: walk OUTWARD from the IGL boundary toward the pia.

Symmetric to radial_spokes (which goes pia → IGL). The IGL is a robust anchor
in deep folds where the outer pial boundary is hard to detect (folds touch each
other across thin LCSF gaps), so IGL-anchored spokes catch what pia-anchored
spokes miss.

For each IGL-boundary point:
  1. compute the local OUTWARD normal (pointing into the cortex)
  2. search +/- ANGLE_SEARCH_DEG for direction of shortest ray to the pia
     (= cerebellum_mask boundary, OR background)
  3. smooth the angle field along arc length

Layer order along these spokes (s increasing from IGL outward):
  IGL -> ML -> iEGL -> oEGL -> pia/background
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple

import numpy as np
import scipy.ndimage as ndi
from skimage import filters, measure


# tuning -----------------------------------------------------------------
PIAL_SAMPLE_UM = 8.0
TANGENT_HALFWIN_UM = 25.0
ANGLE_SEARCH_DEG = 60
ANGLE_SEARCH_N = 41
RAY_MAX_UM = 250.0
RAY_STEP_UM = 1.5
SMOOTH_SIGMA_UM = 30.0


def _resample_outline(outline, spacing_um, px_per_um):
    seg = np.diff(outline, axis=0)
    arc = np.concatenate([[0], np.cumsum(np.linalg.norm(seg, axis=1))])
    spacing_px = spacing_um * px_per_um
    n = max(2, int(round(arc[-1] / spacing_px)))
    targets = np.linspace(0, arc[-1], n, endpoint=False)
    ys = np.interp(targets, arc, outline[:, 0])
    xs = np.interp(targets, arc, outline[:, 1])
    return np.stack([ys, xs], axis=1)


def _local_tangents(pts, halfwin_um, px_per_um):
    halfwin_px = halfwin_um * px_per_um
    seg = np.diff(pts, axis=0)
    arc = np.concatenate([[0], np.cumsum(np.linalg.norm(seg, axis=1))])
    n = len(pts)
    tangents = np.zeros_like(pts)
    for i in range(n):
        a = arc[i]
        j_lo = max(0, np.searchsorted(arc, a - halfwin_px))
        j_hi = min(n - 1, np.searchsorted(arc, a + halfwin_px))
        v = pts[j_hi] - pts[j_lo]
        m = np.linalg.norm(v)
        if m > 0: tangents[i] = v / m
    return tangents


def _outward_normals(pts, tangents, igl_mask, cereb_mask):
    """Pick the perpendicular that steps OUTWARD (out of IGL but into cereb)."""
    H, W = igl_mask.shape
    n_a = np.stack([-tangents[:, 1], tangents[:, 0]], axis=1)
    n_b = -n_a
    out = np.zeros_like(pts)
    for i in range(len(pts)):
        py, px = pts[i]
        for n in (n_a[i], n_b[i]):
            ty = int(round(py + 5 * n[0]))
            tx = int(round(px + 5 * n[1]))
            if (0 <= ty < H and 0 <= tx < W and not igl_mask[ty, tx]
                and cereb_mask[ty, tx]):
                out[i] = n
                break
    return out


def _ray_to_target(start, n, mask_in, target_mask, max_um, step_um, px_per_um):
    """Walk along n from start while staying in mask_in.
    Stop on hitting target_mask. Returns (ye, xe, length_um) or None."""
    H, W = mask_in.shape
    step_px = step_um * px_per_um
    max_steps = int(max_um / step_um)
    for k in range(1, max_steps + 1):
        y = int(round(start[0] + k * step_px * n[0]))
        x = int(round(start[1] + k * step_px * n[1]))
        if not (0 <= y < H and 0 <= x < W):
            return None, None, np.inf
        if target_mask[y, x]:
            return y, x, k * step_um
        if not mask_in[y, x]:
            # exited the cereb mask: this IS the pia boundary
            return y, x, k * step_um
    return None, None, np.inf


def _best_ray_to_pia(start, n0, cereb_mask, target_mask, max_um, step_um,
                     px_per_um, search_deg, search_n):
    """Search +/- search_deg around n0 for the shortest ray to target_mask boundary."""
    angles = np.deg2rad(np.linspace(-search_deg, search_deg, search_n))
    best = (None, None, np.inf, 0.0, n0)
    for a in angles:
        ca, sa = np.cos(a), np.sin(a)
        ny = ca * n0[0] - sa * n0[1]
        nx = sa * n0[0] + ca * n0[1]
        n = np.array([ny, nx])
        ye, xe, L = _ray_to_target(start, n, cereb_mask, target_mask,
                                    max_um, step_um, px_per_um)
        if L < best[2]:
            best = (ye, xe, L, np.rad2deg(a), n)
    return best


def _smooth_directions_along_arc(dirs, sigma_um, px_per_um, spacing_um):
    sigma_pts = sigma_um / spacing_um
    sm_y = ndi.gaussian_filter1d(dirs[:, 0], sigma=sigma_pts, mode='wrap')
    sm_x = ndi.gaussian_filter1d(dirs[:, 1], sigma=sigma_pts, mode='wrap')
    sm = np.stack([sm_y, sm_x], axis=1)
    norms = np.linalg.norm(sm, axis=1, keepdims=True)
    return sm / np.maximum(norms, 1e-9)


@dataclass
class IGLSpokeLayout:
    igl_anchor: np.ndarray      # (N, 2) (y, x) on IGL boundary
    spoke_dirs: np.ndarray      # (N, 2) outward direction
    pia_endpoint: np.ndarray    # (N, 2) (y, x) where ray exits cereb / hits pia
    lengths_um: np.ndarray      # (N,) NaN if invalid
    valid: np.ndarray
    angle_deflection_deg: np.ndarray


def build_igl_spokes(cereb_mask, igl_mask, pixel_size_um=0.5119049,
                     boundary_smooth_um=3.0,
                     pial_sample_um=PIAL_SAMPLE_UM,
                     tangent_halfwin_um=TANGENT_HALFWIN_UM,
                     angle_search_deg=ANGLE_SEARCH_DEG,
                     angle_search_n=ANGLE_SEARCH_N,
                     ray_max_um=RAY_MAX_UM,
                     ray_step_um=RAY_STEP_UM,
                     smooth_sigma_um=SMOOTH_SIGMA_UM,
                     max_dist_to_pia_um=180.0):
    """Build spokes anchored on the cortical-side IGL boundary, going outward.

    `max_dist_to_pia_um` filters out IGL boundary points farther than this
    from the cerebellum mask boundary (i.e., the deep IGL-DWL edge whose
    'outward' direction would cross the WM tract).
    """
    px_per_um = 1.0 / pixel_size_um

    # outer boundary of IGL = IGL pixels adjacent to non-IGL within cereb
    sm_igl = filters.gaussian(igl_mask.astype(float),
                               sigma=boundary_smooth_um*px_per_um) > 0.5
    contours = measure.find_contours(sm_igl.astype(float), 0.5)
    if not contours:
        raise RuntimeError("no IGL boundary found")
    contours.sort(key=len, reverse=True)

    # combine ALL contours (each IGL piece contributes its own boundary)
    # rather than only taking the longest, so we get IGL boundary even in
    # disjoint folds. We'll filter by distance-to-pia next.
    outline = np.concatenate(contours, axis=0)

    pts = _resample_outline(outline, pial_sample_um, px_per_um)
    # Filter: only keep points within max_dist_to_pia_um of the cerebellum
    # boundary (i.e., the cortical-side IGL edge).
    bg = ~cereb_mask
    dist_to_bg_px = ndi.distance_transform_edt(cereb_mask)
    dist_to_bg_um = dist_to_bg_px * pixel_size_um
    yi = np.clip(pts[:, 0].astype(int), 0, cereb_mask.shape[0] - 1)
    xi = np.clip(pts[:, 1].astype(int), 0, cereb_mask.shape[1] - 1)
    keep_mask = dist_to_bg_um[yi, xi] <= max_dist_to_pia_um
    pts = pts[keep_mask]
    n = len(pts)

    tangs = _local_tangents(pts, tangent_halfwin_um, px_per_um)
    norms = _outward_normals(pts, tangs, igl_mask, cereb_mask)

    # target = pia (i.e., outside the cerebellum mask)
    pia_target = ~cereb_mask

    end_y = np.full(n, np.nan)
    end_x = np.full(n, np.nan)
    L_um  = np.full(n, np.nan)
    deflect = np.zeros(n)
    chosen_dirs = norms.copy()
    for i, (p, n0) in enumerate(zip(pts, norms)):
        if np.linalg.norm(n0) < 0.5:
            continue
        ye, xe, L, ang, nbest = _best_ray_to_pia(p, n0, cereb_mask, pia_target,
                                                  ray_max_um, ray_step_um, px_per_um,
                                                  angle_search_deg, angle_search_n)
        if ye is not None:
            end_y[i] = ye; end_x[i] = xe; L_um[i] = L
            deflect[i] = ang; chosen_dirs[i] = nbest

    smooth_dirs = _smooth_directions_along_arc(chosen_dirs, smooth_sigma_um,
                                                px_per_um, pial_sample_um)

    # re-cast each ray with smoothed direction
    end_y_s = np.full(n, np.nan); end_x_s = np.full(n, np.nan); L_s = np.full(n, np.nan)
    for i in range(n):
        if np.linalg.norm(smooth_dirs[i]) < 0.5: continue
        ye, xe, L = _ray_to_target(pts[i], smooth_dirs[i], cereb_mask, pia_target,
                                    ray_max_um, ray_step_um, px_per_um)
        if ye is not None:
            end_y_s[i] = ye; end_x_s[i] = xe; L_s[i] = L

    valid = ~np.isnan(L_s)
    return IGLSpokeLayout(
        igl_anchor=pts,
        spoke_dirs=smooth_dirs,
        pia_endpoint=np.stack([end_y_s, end_x_s], axis=1),
        lengths_um=L_s,
        valid=valid,
        angle_deflection_deg=deflect,
    )
