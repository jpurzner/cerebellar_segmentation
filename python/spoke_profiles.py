"""1D layer detection along each spoke (layer_quant style).

For each spoke:
  - sample (p27, DAPI, NeuN) at every step_um along the spoke pia->IGL
  - find pia anchor via pia_prunner trick: peaks in d(p27)/d(s) * d(NeuN)/d(s)
  - identify EGL: characteristic high-DAPI ribbon
  - identify iEGL within EGL: high p27 sub-region
  - identify IGL: high NeuN at the deep end
  - ML lies between EGL and IGL

Key signal-processing recipes from layer_quant/pia_prunner:
  pial_finder = diff(p27) * diff(NeuN)        # peaks at pial boundary
  EGL: high (DAPI + p27) within outer ~50 um
  iEGL: peak p27 inside EGL
  IGL: high NeuN, plateau

Returns per-spoke layer boundaries in the 1D coordinate (um from pia).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np
import scipy.ndimage as ndi
from scipy.signal import find_peaks


# ---------------------------------------------------------------------------
SAMPLE_STEP_UM = 1.0          # physical step along spoke
SMOOTH_PROFILE_UM = 4.0       # gaussian smoothing on 1D profile
EGL_MAX_DEPTH_UM = 60.0       # search EGL only in the outer 60 um


def sample_along_spoke(image_stack: np.ndarray,
                        pia_pt: Tuple[float, float],
                        igl_pt: Tuple[float, float],
                        step_um: float,
                        pixel_size_um: float) -> np.ndarray:
    """Sample (n_channels, n_steps) intensity profile along a spoke.

    image_stack shape: (n_channels, H, W). pia_pt and igl_pt are (y, x) floats.
    """
    n_ch = image_stack.shape[0]
    px_per_um = 1.0 / pixel_size_um
    L_px = np.hypot(igl_pt[0] - pia_pt[0], igl_pt[1] - pia_pt[1])
    L_um = L_px / px_per_um
    n = max(2, int(round(L_um / step_um)) + 1)
    ys = np.linspace(pia_pt[0], igl_pt[0], n)
    xs = np.linspace(pia_pt[1], igl_pt[1], n)
    out = np.zeros((n_ch, n), dtype=np.float32)
    for c in range(n_ch):
        out[c] = ndi.map_coordinates(image_stack[c], [ys, xs],
                                      order=1, mode='nearest')
    return out


def detect_layers_1d(profile: np.ndarray, step_um: float,
                     channels: dict = {"p27": 0, "neun": 1, "dapi": 2},
                     egl_max_depth_um: float = EGL_MAX_DEPTH_UM,
                     smooth_um: float = SMOOTH_PROFILE_UM) -> dict:
    """Detect layer boundaries on a 1D profile from a spoke.

    Spoke starts at the pial surface (s=0) by construction. Layer order
    along the spoke (s increasing inward):

        EGL (oEGL near pia, iEGL deeper) -> ML -> IGL

    Strategy:
      - EGL_end: where DAPI drops below half its peak (going inward), within
                 the first egl_max_depth_um.
      - IGL_start: where NeuN rises above half its peak (going inward), after
                 the EGL_end.
      - iEGL_start, iEGL_end: where p27 within EGL is above half its EGL-peak.
    """
    p27_idx  = channels["p27"]
    dapi_idx = channels["dapi"]
    neun_idx = channels["neun"]
    n = profile.shape[1]
    s_um = np.arange(n) * step_um

    sigma = max(0.5, smooth_um / step_um)
    p27  = ndi.gaussian_filter1d(profile[p27_idx],  sigma=sigma)
    dapi = ndi.gaussian_filter1d(profile[dapi_idx], sigma=sigma)
    neun = ndi.gaussian_filter1d(profile[neun_idx], sigma=sigma)

    pia_step = 0
    pia_um = 0.0

    # --- EGL_end: location of maximum NEGATIVE gradient of DAPI in outer
    #    egl_max_depth_um. The strongest "edge from bright to dark" along the
    #    spoke is the EGL->ML transition, regardless of absolute intensity. ---
    cap = min(n - 1, int(egl_max_depth_um / step_um))
    if cap > 3:
        # Smoothed gradient via Sobel-like central difference
        grad_dapi = np.gradient(dapi)
        # restrict search to the outer 1.5x egl_max_depth_um (allow EGL to be a bit
        # deeper than the cap but never beyond)
        search_end = min(n, int(egl_max_depth_um * 1.5 / step_um))
        seg = grad_dapi[:search_end]
        # find min of gradient (most negative = sharp drop)
        # but skip the very first 5um (image-edge artifacts) and require a
        # non-trivially negative value
        min_skip = max(2, int(5.0 / step_um))
        seg_search = seg[min_skip:]
        if len(seg_search) > 0 and seg_search.min() < -0.005:
            egl_end_step = min_skip + int(np.argmin(seg_search))
        else:
            egl_end_step = cap // 2
    else:
        egl_end_step = max(1, cap)

    # --- IGL_start: location of maximum POSITIVE gradient of NeuN, AFTER
    #    egl_end_step. ---
    if egl_end_step + 3 < n:
        grad_neun = np.gradient(neun)
        seg_after = grad_neun[egl_end_step:]
        if seg_after.max() > 0.005:
            igl_start_step = egl_end_step + int(np.argmax(seg_after))
        else:
            # NeuN never sharply rises -> IGL touches EGL with no ML
            igl_start_step = egl_end_step
    else:
        igl_start_step = max(egl_end_step, n - 1)

    # --- iEGL: p27 above half its EGL-peak, within EGL ---
    iegl_start_step = egl_end_step
    iegl_end_step = egl_end_step
    if egl_end_step > pia_step + 2:
        p27_seg = p27[pia_step:egl_end_step + 1]
        peak = p27_seg.max()
        thr = max(0.20, peak * 0.55)
        above = p27_seg >= thr
        idx = np.where(above)[0]
        if len(idx) > 0:
            iegl_start_step = pia_step + idx[0]
            iegl_end_step   = pia_step + idx[-1]

    return {
        "pia_um":         pia_um,
        "egl_end_um":     s_um[egl_end_step],
        "iegl_start_um":  s_um[iegl_start_step],
        "iegl_end_um":    s_um[iegl_end_step],
        "igl_start_um":   s_um[igl_start_step],
        "igl_end_um":     s_um[-1],
        "n_steps":        n,
        "L_um":           s_um[-1] if len(s_um) else 0,
    }
