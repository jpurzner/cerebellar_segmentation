"""
Find photobleached rectangles (40x imaging fields) in a 10x IF section: boxes of the 40x
field size (default 226 x 169 um, both orientations, axis-aligned in the native image)
whose interior is dimmer than a surrounding ring in ALL channels (p27, DAPI, NeuN).

    from bleach_boxes import find_boxes
    boxes = find_boxes(raw (C,H,W), tissue (H,W) bool, px_um)

Contrast per channel = mean log-intensity inside minus mean in the ring (tissue pixels
only). Score = the least negative channel contrast (all channels must drop). Boxes are
reported best first, non-overlapping, with score <= max_score.
"""
import numpy as np
from scipy import ndimage as ndi


def _boxsum(ii, h, w):
    """sum over every h x w box (top-left anchored) from an integral image ii (H+1, W+1)."""
    return ii[h:, w:] - ii[:-h, w:] - ii[h:, :-w] + ii[:-h, :-w]


def find_boxes(raw, tissue, px, field_um=(226.0, 169.0), ring_um=40.0, grid_um=2.6, max_score=-0.12,
               min_tissue=0.6, top=8):
    ds = max(1, int(round(grid_um / px)))
    t = tissue[::ds, ::ds].astype(np.float64)
    L = [np.log(np.maximum(ndi.uniform_filter(c[::ds, ::ds].astype(np.float64), 3), 1.0)) * t for c in raw]
    g = px * ds
    r = int(round(ring_um / g))
    ii_t = np.pad(t.cumsum(0).cumsum(1), ((1, 0), (1, 0)))
    ii_L = [np.pad(l.cumsum(0).cumsum(1), ((1, 0), (1, 0))) for l in L]
    found = []
    for (a, b) in (field_um, field_um[::-1]):
        h, w = int(round(b / g)), int(round(a / g))
        H2, W2 = h + 2 * r, w + 2 * r
        n_in = _boxsum(ii_t, h, w)[r:, r:]
        n_out = _boxsum(ii_t, H2, W2)
        Hc, Wc = n_out.shape
        n_in = n_in[:Hc, :Wc]
        ring_n = n_out - n_in
        ok = (n_in >= min_tissue * h * w) & (ring_n >= min_tissue * (H2 * W2 - h * w))
        sc = np.full(n_out.shape, np.inf)
        cons = []
        for ii in ii_L:
            s_in = _boxsum(ii, h, w)[r:, r:][:Hc, :Wc]
            s_out = _boxsum(ii, H2, W2)
            cons.append(np.where(ok, s_in / np.maximum(n_in, 1) - (s_out - s_in) / np.maximum(ring_n, 1), np.inf))
        sc = np.max(np.stack(cons), 0)          # least negative channel
        for _ in range(top):
            k = np.unravel_index(np.argmin(sc), sc.shape)
            if not np.isfinite(sc[k]) or sc[k] > max_score:
                break
            y, x = (k[0] + r) * ds, (k[1] + r) * ds
            found.append(dict(score=float(sc[k]), contrasts=[float(c[k]) for c in cons],
                              row0=int(y), col0=int(x), rows=int(h * ds), cols=int(w * ds),
                              size_um=[round(a), round(b)]))
            y0, x0 = max(k[0] - h, 0), max(k[1] - w, 0)
            sc[y0:k[0] + h, x0:k[1] + w] = np.inf    # suppress overlapping boxes
    found.sort(key=lambda f: f["score"])
    keep = []
    for f in found:
        if all(not (f["row0"] < q["row0"] + q["rows"] and q["row0"] < f["row0"] + f["rows"] and
                    f["col0"] < q["col0"] + q["cols"] and q["col0"] < f["col0"] + f["cols"]) for q in keep):
            keep.append(f)
    return keep
