"""
Re-create synthetic sections from their seeds and draw sections (synthetic or real) as SVG
paths for the report.
"""
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi
from skimage import measure

sys.path.insert(0, str(Path(__file__).parent))
import igl_synth_gen as GEN  # noqa: E402

LAYERS = ("tissue", "wm", "igl", "egl")


def regenerate(fam, seed, params=None):
    rng = np.random.default_rng(seed)
    prior, gen = GEN.FAMILIES[fam]
    p = prior(rng) if params is None else params
    return gen(p, np.random.default_rng(seed + 10 ** 7)), p


def paths(mask, tol=1.0):
    m = np.pad(mask.astype(float), 1)
    out = []
    for c in measure.find_contours(m, 0.5):
        c = measure.approximate_polygon(c - 1, tol)
        if len(c) < 3:
            continue
        out.append("M" + "L".join(f"{x:.0f},{y:.0f}" for y, x in c) + "Z")
    return "".join(out)


def svg_data(g, max_px=700):
    """masks -> {n: [h, w], <layer>: path}, downsampled so the longer side is <= max_px."""
    h, w = g["tissue"].shape
    f = max(1, int(np.ceil(max(h, w) / max_px)))
    out = dict(n=[int(np.ceil(h / f)), int(np.ceil(w / f))], um_per_unit=4.0 * f)
    for k in LAYERS:
        m = g[k][::f, ::f] if f > 1 else g[k]
        out[k] = paths(m)
    return out
