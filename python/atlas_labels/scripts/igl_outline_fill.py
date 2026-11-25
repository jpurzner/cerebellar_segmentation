"""
Outline-driven IGL gap filler.

The synthetic-IGL study (igl_synth_*.py) found that the IGL is best described as the room
left over: everything below the molecular layer that is not white matter, where the white
matter is a core along each folium's centre line, present where the folium's half-width h
exceeds h_min, with radius rho * (h - m). This module turns that rule into a predictor
for missing IGL.

    import igl_outline_fill as OF
    geo = OF.geometry(g)                                   # g: 4 um label grid (igl_ribbon.grid_labels)
    prm = OF.fit(geo, igl_known, g["wm"], mode)            # m (+ rho, h_min) fitted on the section's own IGL
    P   = OF.predict(geo, prm, g["wm"], mode)              # predicted IGL over the whole section
    F   = OF.fill_between(P, igl_known, cand, a, b, wa, wb)    # the stretch of P joining two ribbon ends

Geometry (4 um grid):
  surface  EGL, plus background within 60 um of the EGL (pial surface where the EGL label
           has a short hole; the base cut, far from any EGL, does not count)
  dS       depth below the surface
  inner    tissue minus EGL; its medial axis = the centre line of every folium and lobule,
           h = distance from the axis to the EGL / background (half-width inside the EGL)
Rule, for a pixel x of `inner`:
  IGL  <=>  dS(x) >= m  and  x is not white matter
White matter by mode:
  label    the white-matter label
  outline  within rho * (h_a - h0) of a, the nearest axis point with h_a >= h0
           (no labels used: survives white matter that is mislabelled or missing). The
           layer measurements give rho ~0.84, h0 ~95 um for the IF controls; in pixels
           this core matches the labelled white matter only moderately (Dice 0.5-0.7).
  gated    the white-matter label, but only within `margin` of the outline core
           (labelled white matter far from where a core can be is ignored)
  union    label or outline
  depth    (predict_depth) no white matter at all: the IGL is the band m <= depth <= m + w,
           m and w measured on the ribbon either side of the gap
Fit: grid search maximising the Dice overlap of the rule's IGL with the known IGL of the
same section (pixels subsampled 2x2). The gap itself is never seen by the fit: the caller
passes the IGL with the gap removed.
Gap fill: candidates = predicted IGL that is not already IGL (1 px clearance), inside the
allowed tissue. Geodesic distances inside the candidates from each end's face; the fill is
the geodesic "ellipse" d_a + d_b <= shortest + slack (slack = the flank width + 20 um), so
only the corridor joining the two ends is taken, never the whole predicted ribbon. If the
ends are not connected inside the candidates, nothing is filled (status "disconnected").
`band` instead takes the predicted IGL inside a tube around a given centre line (e.g. the
ribbon model's), the ribbon model choosing the route and the rule setting the edges.
"""
import numpy as np
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
from skimage.graph import MCP_Geometric
from skimage.morphology import medial_axis

G = 4.0
SURF_BG_UM = 60.0
H_MIN_GRID = np.arange(40.0, 301.0, 20.0)        # h0: inner half-width (um) where the core starts
RHO_GRID = np.arange(0.2, 1.41, 0.1)
MARGIN_GRID = (16.0, 32.0, 48.0)
M_GRID = np.arange(20.0, 101.0, 4.0)
M_GRID_FINE = np.arange(16.0, 121.0, 2.0)
SLACK_UM = 20.0


def geometry(g):
    tissue = g["tissue"] | g.get("anom", np.zeros_like(g["tissue"]))
    egl = g["egl"] & tissue
    dEgl = ndi.distance_transform_edt(~egl) * G
    surface = egl | (~tissue & (dEgl <= SURF_BG_UM))
    dS = ndi.distance_transform_edt(~surface) * G
    inner = tissue & ~egl
    ax, dist = medial_axis(inner, return_distance=True)
    h = dist * G
    cores = {}
    for hm in H_MIN_GRID:
        ok = ax & (h >= hm)
        if not ok.any():
            cores[float(hm)] = None
            continue
        dax, (iy, ix) = ndi.distance_transform_edt(~ok, return_indices=True)
        cores[float(hm)] = (dax * G, h[iy, ix])
    return dict(tissue=tissue, egl=egl, inner=inner, dS=dS, ax=ax, h=h, cores=cores)


def _core_for(geo, h_min):
    hm = float(H_MIN_GRID[np.argmin(np.abs(H_MIN_GRID - h_min))])
    return geo["cores"][hm]


def wm_outline(geo, prm, m=None):
    m = prm["m"] if m is None else m
    c = _core_for(geo, prm["h_min"])
    if c is None:
        return np.zeros_like(geo["inner"])
    dax, hn = c
    return geo["inner"] & (dax <= prm["rho"] * (hn - prm["h_min"]))


def white_matter(geo, prm, wm_label, mode, m=None):
    if mode == "label":
        return wm_label
    core = wm_outline(geo, prm, m)
    if mode == "outline":
        return core
    if mode == "union":
        return wm_label | core
    if mode == "gated":
        c = _core_for(geo, prm["h_min"])
        if c is None:
            return np.zeros_like(wm_label)
        dax, hn = c
        return wm_label & (dax <= prm["rho"] * (hn - prm["h_min"]) + prm["margin"])
    raise ValueError(mode)


def predict(geo, prm, wm_label, mode, m=None):
    m = prm["m"] if m is None else m
    return geo["inner"] & (geo["dS"] >= m) & ~white_matter(geo, prm, wm_label, mode, m)


def predict_depth(geo, m, w, slack=8.0):
    """depth band: IGL = m <= depth <= m + w + slack, with m and w taken from the ribbon on
    either side of the gap. Uses no white-matter labels at all."""
    return geo["inner"] & (geo["dS"] >= m) & (geo["dS"] <= m + w + slack)


def _dice_counts(P, T):
    return 2.0 * np.count_nonzero(P & T) / max(np.count_nonzero(P) + np.count_nonzero(T), 1)


def fit(geo, igl_known, wm_label, mode, sub=2):
    """grid search on a 2x2-subsampled pixel set of the section's inner tissue."""
    sl = (slice(None, None, sub), slice(None, None, sub))
    inner = geo["inner"][sl]
    idx = np.nonzero(inner.ravel())[0]
    T = igl_known[sl].ravel()[idx]
    dS = geo["dS"][sl].ravel()[idx]
    WL = wm_label[sl].ravel()[idx]
    best = (-1.0, None)
    if mode == "label":
        for m in M_GRID_FINE:
            d = _dice_counts((dS >= m) & ~WL, T)
            if d > best[0]:
                best = (d, dict(m=float(m), rho=np.nan, h_min=np.nan))
        return dict(best[1], dice=round(best[0], 4))
    for hm in H_MIN_GRID:
        c = geo["cores"][float(hm)]
        if c is None:
            continue
        dax = c[0][sl].ravel()[idx]
        hn = c[1][sl].ravel()[idx]
        room = hn - hm
        for m in M_GRID:
            below = dS >= m
            for rho in RHO_GRID:
                core = dax <= rho * room
                for mg in (MARGIN_GRID if mode == "gated" else (0.0,)):
                    if mode == "outline":
                        W = core
                    elif mode == "union":
                        W = core | WL
                    else:  # gated
                        W = WL & (dax <= rho * room + mg)
                    d = _dice_counts(below & ~W, T)
                    if d > best[0]:
                        best = (d, dict(m=float(m), rho=float(rho), h_min=float(hm), margin=float(mg)))
    return dict(best[1], dice=round(best[0], 4))


def _geodesic(cand, seeds):
    cost = np.where(cand, 1.0, np.inf)
    mcp = MCP_Geometric(cost)
    D, _ = mcp.find_costs([tuple(s) for s in seeds])
    return D


def _seeds(cand, p, r_px):
    y0, x0 = np.round(p).astype(int)
    r = int(np.ceil(r_px))
    ys, xs = np.mgrid[max(y0 - r, 0):min(y0 + r + 1, cand.shape[0]), max(x0 - r, 0):min(x0 + r + 1, cand.shape[1])]
    sel = ((ys - p[0]) ** 2 + (xs - p[1]) ** 2 <= r_px ** 2) & cand[ys, xs]
    return np.stack([ys[sel], xs[sel]], 1)


def fill_between(P, igl_known, allowed, a, b, wa, wb):
    """the stretch of predicted IGL P joining ribbon ends a and b (px coords); wa, wb flank
    widths (um). Returns (fill mask, status)."""
    cand = P & allowed & ~ndi.binary_dilation(igl_known, iterations=1)
    w = max(wa, wb)
    ra, rb = (wa / 2 + 12.0) / G, (wb / 2 + 12.0) / G
    sa, sb = _seeds(cand, a, ra), _seeds(cand, b, rb)
    if len(sa) == 0 or len(sb) == 0:
        return np.zeros_like(P), "no_seed"
    # work in a crop around the two ends for speed
    pad = int(max(np.hypot(*(np.asarray(a) - b)) * 1.5, 60))
    y0 = max(int(min(a[0], b[0])) - pad, 0)
    x0 = max(int(min(a[1], b[1])) - pad, 0)
    y1 = min(int(max(a[0], b[0])) + pad + 1, P.shape[0])
    x1 = min(int(max(a[1], b[1])) + pad + 1, P.shape[1])
    c = cand[y0:y1, x0:x1]
    off = np.array([y0, x0])
    Da = _geodesic(c, sa - off)
    Db = _geodesic(c, sb - off)
    S = Da + Db
    smin = np.min(S)
    if not np.isfinite(smin):
        return np.zeros_like(P), "disconnected"
    keep = np.isfinite(S) & (S <= smin + (w + SLACK_UM) / G)
    out = np.zeros_like(P)
    out[y0:y1, x0:x1] = keep
    # close the 1 px seam left against the known IGL
    out |= ndi.binary_dilation(out, iterations=1) & P & allowed & ~igl_known
    return out, "ok"


def band(P, igl_known, allowed, C, wa, wb, extra_um=16.0):
    """predicted IGL P inside a tube around centre line C (px), radius = half the larger
    flank width + extra."""
    C = np.asarray(C, float)
    L = np.concatenate([[0], np.cumsum(np.hypot(*np.diff(C, axis=0).T))])
    n = max(int(L[-1] / 0.5), 2)
    t = np.linspace(0, L[-1], n)
    C = np.stack([np.interp(t, L, C[:, k]) for k in range(2)], 1)
    r = (max(wa, wb) / 2 + extra_um) / G
    y0, x0 = np.maximum(np.floor(C.min(0) - r - 1).astype(int), 0)
    y1, x1 = np.minimum(np.ceil(C.max(0) + r + 2).astype(int), P.shape)
    yy, xx = np.mgrid[y0:y1, x0:x1]
    d, _ = cKDTree(C).query(np.stack([yy.ravel(), xx.ravel()], 1))
    tube = np.zeros(P.shape, bool)
    sel = d <= r
    tube[yy.ravel()[sel], xx.ravel()[sel]] = True
    return tube & P & allowed & ~igl_known
