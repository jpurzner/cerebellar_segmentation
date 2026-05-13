"""Smooth-angle radial detector — first iteration.

Algorithm:
  1. Take the cerebellum mask + IGL mask (from the Python port).
  2. Sample pia points uniformly along the outer mask boundary.
  3. At each pia point, compute the local pial tangent (smoothed over arc length)
     and the inward normal.
  4. Search angles +/- ANGLE_SEARCH_DEG around the inward normal for the
     direction whose perpendicular ray reaches the IGL with the shortest length
     (constraint: must stay inside the cerebellum mask).
  5. Smooth the resulting spoke directions along arc length so adjacent spokes
     don't suddenly jump.
  6. Walk each spoke from pia to its IGL endpoint; record its (y,x) trace and
     length.

Outputs a dict containing:
    pia_points (N,2)     -- pia (y,x) anchor for each spoke
    spoke_dirs (N,2)     -- unit (dy,dx) inward direction
    igl_endpoints (N,2)  -- (y,x) where the spoke first hits IGL
    lengths_um (N,)      -- physical length of each spoke
    valid (N,)           -- whether the spoke reaches IGL within RAY_MAX_UM
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np
import scipy.ndimage as ndi
from skimage import filters, measure
from skimage.morphology import disk


# tuning ------------------------------------------------------------------
PIAL_SAMPLE_UM = 8.0           # spoke spacing along pia
TANGENT_HALFWIN_UM = 25.0      # window for local tangent estimation
ANGLE_SEARCH_DEG = 60          # +/- angle search around local normal
ANGLE_SEARCH_N = 41            # number of angles tested
RAY_MAX_UM = 600               # don't search rays longer than this
RAY_STEP_UM = 1.5              # sampling step along ray
SMOOTH_SIGMA_UM = 30.0         # arc-length sigma for spoke direction smoothing


# ------------------------------------------------------------------------

def _resample_outline(outline: np.ndarray, spacing_um: float, px_per_um: float):
    """Resample a closed outline (Nx2 of (y,x)) to roughly uniform arc-length."""
    seg = np.diff(outline, axis=0)
    arc = np.concatenate([[0], np.cumsum(np.linalg.norm(seg, axis=1))])
    spacing_px = spacing_um * px_per_um
    n = max(2, int(round(arc[-1] / spacing_px)))
    targets = np.linspace(0, arc[-1], n, endpoint=False)
    ys = np.interp(targets, arc, outline[:, 0])
    xs = np.interp(targets, arc, outline[:, 1])
    return np.stack([ys, xs], axis=1)


def _local_tangents(pts: np.ndarray, halfwin_um: float, px_per_um: float) -> np.ndarray:
    """Local unit tangent per point via centered diff over a window."""
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
        if m > 0:
            tangents[i] = v / m
    return tangents


def _inward_normals(pts: np.ndarray, tangents: np.ndarray,
                    mask: np.ndarray) -> np.ndarray:
    """Pick inward normal (perpendicular to tangent that steps INTO mask)."""
    H, W = mask.shape
    n_a = np.stack([-tangents[:, 1], tangents[:, 0]], axis=1)
    n_b = -n_a
    out = np.zeros_like(pts)
    for i in range(len(pts)):
        py, px = pts[i]
        for n in (n_a[i], n_b[i]):
            ty = int(round(py + 3 * n[0]))
            tx = int(round(px + 3 * n[1]))
            if 0 <= ty < H and 0 <= tx < W and mask[ty, tx]:
                out[i] = n
                break
    return out


def _ray_to_igl(start: Tuple[float, float], n: np.ndarray,
                mask: np.ndarray, igl: np.ndarray,
                max_um: float, step_um: float, px_per_um: float):
    """Walk along direction n from start until hitting IGL or leaving mask.
    Returns (end_y, end_x, length_um) or (None, None, np.inf) if never hits.
    """
    H, W = mask.shape
    max_steps = int(max_um * px_per_um / (step_um * px_per_um))
    step_px = step_um * px_per_um
    for k in range(1, max_steps + 1):
        y = int(round(start[0] + k * step_px * n[0]))
        x = int(round(start[1] + k * step_px * n[1]))
        if not (0 <= y < H and 0 <= x < W) or not mask[y, x]:
            return None, None, np.inf
        if igl[y, x]:
            return y, x, k * step_um
    return None, None, np.inf


def _best_ray(start: Tuple[float, float], n0: np.ndarray,
              mask: np.ndarray, igl: np.ndarray,
              max_um: float, step_um: float, px_per_um: float,
              search_deg: float, search_n: int):
    """Search +/-search_deg around n0 for the direction with shortest ray to IGL."""
    angles = np.deg2rad(np.linspace(-search_deg, search_deg, search_n))
    best = (None, None, np.inf, 0.0, n0)
    for a in angles:
        ca, sa = np.cos(a), np.sin(a)
        ny = ca * n0[0] - sa * n0[1]
        nx = sa * n0[0] + ca * n0[1]
        n = np.array([ny, nx])
        ye, xe, L = _ray_to_igl(start, n, mask, igl, max_um, step_um, px_per_um)
        if L < best[2]:
            best = (ye, xe, L, np.rad2deg(a), n)
    return best


def _smooth_directions_along_arc(dirs: np.ndarray, sigma_um: float,
                                  px_per_um: float, spacing_um: float) -> np.ndarray:
    """Gaussian-smooth (dy,dx) along arc length, then renormalize.

    Doing this in (dy,dx) space rather than angle space avoids wrap-around
    issues at +/-180.  The closed-loop nature of the cerebellum boundary is
    handled by `mode='wrap'` in scipy.ndimage.gaussian_filter1d.
    """
    sigma_pts = sigma_um / spacing_um
    sm_y = ndi.gaussian_filter1d(dirs[:, 0], sigma=sigma_pts, mode='wrap')
    sm_x = ndi.gaussian_filter1d(dirs[:, 1], sigma=sigma_pts, mode='wrap')
    sm = np.stack([sm_y, sm_x], axis=1)
    norms = np.linalg.norm(sm, axis=1, keepdims=True)
    return sm / np.maximum(norms, 1e-9)


@dataclass
class SpokeLayout:
    pia_points: np.ndarray       # (N, 2) (y, x)
    spoke_dirs: np.ndarray       # (N, 2) unit (dy, dx) AFTER smoothing
    igl_endpoints: np.ndarray    # (N, 2) (y, x) NaN if invalid
    lengths_um: np.ndarray       # (N,) NaN if invalid
    valid: np.ndarray            # (N,) bool
    angle_deflection_deg: np.ndarray  # (N,) angle picked by best_ray search
    tangents: np.ndarray         # (N, 2) for diagnostics
    initial_normals: np.ndarray  # (N, 2) for diagnostics


def build_spokes(cereb_mask: np.ndarray, igl_mask: np.ndarray,
                 pixel_size_um: float = 0.5119049,
                 *,
                 boundary_smooth_um: float = 10.0,
                 pial_sample_um: float = PIAL_SAMPLE_UM,
                 tangent_halfwin_um: float = TANGENT_HALFWIN_UM,
                 angle_search_deg: float = ANGLE_SEARCH_DEG,
                 angle_search_n: int = ANGLE_SEARCH_N,
                 ray_max_um: float = RAY_MAX_UM,
                 ray_step_um: float = RAY_STEP_UM,
                 smooth_sigma_um: float = SMOOTH_SIGMA_UM) -> SpokeLayout:
    px_per_um = 1.0 / pixel_size_um

    # 1. smoothed pia outline
    sm_mask = filters.gaussian(cereb_mask.astype(float),
                                sigma=boundary_smooth_um * px_per_um) > 0.5
    contours = measure.find_contours(sm_mask.astype(float), 0.5)
    if not contours:
        raise RuntimeError("no cerebellum boundary found")
    contours.sort(key=len, reverse=True)
    outline = contours[0]    # (M, 2) (y, x)

    pts = _resample_outline(outline, pial_sample_um, px_per_um)
    n_spokes = len(pts)

    # 2. local tangents + inward normals
    tangs = _local_tangents(pts, tangent_halfwin_um, px_per_um)
    norms_inward = _inward_normals(pts, tangs, cereb_mask)

    # 3. for each pia point, search +/- angle_search_deg for shortest ray to IGL
    end_y = np.full(n_spokes, np.nan)
    end_x = np.full(n_spokes, np.nan)
    L_um  = np.full(n_spokes, np.nan)
    deflect = np.zeros(n_spokes)
    chosen_dirs = norms_inward.copy()
    for i, (p, n0) in enumerate(zip(pts, norms_inward)):
        if np.linalg.norm(n0) < 0.5:
            continue
        ye, xe, L, ang, nbest = _best_ray(p, n0, cereb_mask, igl_mask,
                                           ray_max_um, ray_step_um, px_per_um,
                                           angle_search_deg, angle_search_n)
        if ye is not None:
            end_y[i] = ye
            end_x[i] = xe
            L_um[i] = L
            deflect[i] = ang
            chosen_dirs[i] = nbest

    # 4. smooth chosen direction field along arc length
    smooth_dirs = _smooth_directions_along_arc(
        chosen_dirs, smooth_sigma_um, px_per_um, pial_sample_um)

    # 5. RE-cast each spoke along its smoothed direction (so endpoints reflect smoothing)
    end_y_sm = np.full(n_spokes, np.nan)
    end_x_sm = np.full(n_spokes, np.nan)
    L_sm = np.full(n_spokes, np.nan)
    for i in range(n_spokes):
        if np.linalg.norm(smooth_dirs[i]) < 0.5:
            continue
        ye, xe, L = _ray_to_igl(pts[i], smooth_dirs[i], cereb_mask, igl_mask,
                                 ray_max_um, ray_step_um, px_per_um)
        if ye is not None:
            end_y_sm[i] = ye; end_x_sm[i] = xe; L_sm[i] = L

    valid = ~np.isnan(L_sm)
    return SpokeLayout(
        pia_points=pts,
        spoke_dirs=smooth_dirs,
        igl_endpoints=np.stack([end_y_sm, end_x_sm], axis=1),
        lengths_um=L_sm,
        valid=valid,
        angle_deflection_deg=deflect,
        tangents=tangs,
        initial_normals=norms_inward,
    )
