"""
A mathematical model of the IGL ribbon, learned from intact ribbons, for predicting
missing stretches.

The ribbon is a centre-line curve gamma(s) (s = arclength) carrying a width w(s); along
it we also track the depth below the EGL d(s). On intact ribbons of the good-IGL controls
we measure curvature kappa(s) = dtheta/ds, w(s), d(s), their variability and how far
along s each stays correlated (correlation length = lag where the autocorrelation of the
detrended signal falls to 1/e).

Gap completion: given the two flank ends (position, direction, width, depth), predict the
missing centre line by one of
  line     - straight segment
  hermite  - cubic Hermite curve matching both end directions
  elastic  - minimum-energy curve: sum of squared curvature (elastica) + depth term
             (stay at the flank depth below the EGL, sigma from the statistics) + barrier
             (stay out of EGL, white matter, background) + even spacing; ends and end
             directions fixed; initialised from the hermite curve
and the width by cubic Hermite interpolation of the flank widths.
"""
import numpy as np
from scipy import ndimage as ndi
from scipy import optimize

import igl_ribbon as RB

GRID = RB.GRID


def smooth_path(path, sigma_px=2.0, step_px=1.0):
    """smooth an ordered pixel path and resample it at an even spacing."""
    p = path.astype(float)
    if len(p) < 5:
        return p
    p = np.stack([ndi.gaussian_filter1d(p[:, k], sigma_px, mode="nearest") for k in range(2)], 1)
    seg = np.hypot(*np.diff(p, axis=0).T)
    s = np.concatenate([[0], np.cumsum(seg)])
    n = max(int(s[-1] / step_px), 2)
    t = np.linspace(0, s[-1], n)
    return np.stack([np.interp(t, s, p[:, k]) for k in range(2)], 1)


def features(curve, width_map, depth_map):
    d = np.diff(curve, axis=0)
    theta = np.unwrap(np.arctan2(d[:, 0], d[:, 1]))
    kappa = np.gradient(theta) / GRID                          # 1/um, per step of 1 px = 4 um
    idx = np.clip(np.round(curve).astype(int), 0, np.array(width_map.shape) - 1)
    w = width_map[idx[:, 0], idx[:, 1]]
    dep = depth_map[idx[:, 0], idx[:, 1]]
    return theta, kappa, w, dep


def corr_length(x, max_lag):
    x = x - ndi.uniform_filter1d(x, 101, mode="nearest") if len(x) > 120 else x - x.mean()
    x = x - x.mean()
    v = (x * x).mean()
    if v <= 0:
        return np.nan
    for lag in range(1, min(max_lag, len(x) - 1)):
        if (x[:-lag] * x[lag:]).mean() / v < 1 / np.e:
            return lag * GRID
    return np.nan


# --------------------------------------------------------------------------- #
# completion models
# --------------------------------------------------------------------------- #
def hermite(a, ta, b, tb_in, n):
    """cubic Hermite from a (direction ta) to b (arriving with direction tb_in)."""
    L = np.hypot(*(b - a))
    s = np.linspace(0, 1, n)[:, None]
    h00, h10 = 2 * s ** 3 - 3 * s ** 2 + 1, s ** 3 - 2 * s ** 2 + s
    h01, h11 = -2 * s ** 3 + 3 * s ** 2, s ** 3 - s ** 2
    return h00 * a + h10 * (L * ta) + h01 * b + h11 * (L * tb_in)


def line(a, b, n):
    return a + np.linspace(0, 1, n)[:, None] * (b - a)


class Field:
    """smooth fields of a section on the 4 um grid used by the elastic model: depth below
    the EGL (um) and distance (px) to forbidden tissue (EGL, white matter, background),
    with their gradients."""

    def __init__(self, depth_map, forbid_mask, smooth_px=1.0):
        self.forbid = forbid_mask
        self.dep_raw = depth_map
        self.dep = ndi.gaussian_filter(depth_map.astype(float), smooth_px)
        self.dgy, self.dgx = np.gradient(self.dep)
        self.fd = ndi.gaussian_filter(ndi.distance_transform_edt(~forbid_mask), smooth_px)
        self.fgy, self.fgx = np.gradient(self.fd)

    @staticmethod
    def at(a, Q):
        return ndi.map_coordinates(a, Q.T, order=1, mode="nearest")


def elastic(a, ta, b, tb_in, field, d_a, d_b, prm=None, init=None, return_energy=False):
    """minimum-energy completion between the flank ends a and b (px, end directions ta
    and tb_in). Energy, in continuous form (s = arclength in um):
        int kappa^2 ds / (2 q)                           curvature as white noise of intensity q
      + int (d(s) - dhat(s))^2 ds / (2 sigma_d^2 l_d)    depth below the EGL near the flank-interpolated depth
      + w_barrier * sum max(0, margin - dist_forbidden)^2
      + w_spacing * sum (segment - h)^2 / h^2            even spacing
    Ends fixed; end directions fixed by ghost points. Analytic gradient, L-BFGS-B."""
    prm = {**DEF, **(prm or {})}
    if init is None:
        L0 = np.hypot(*(b - a))
        init = hermite(a, ta, b, tb_in, 50)
        lh = np.hypot(*np.diff(init, axis=0).T).sum()
        n = int(np.clip(round(lh / prm["spacing_px"]), 6, 120))
        init = hermite(a, ta, b, tb_in, n)
    n = len(init)
    h = np.hypot(*np.diff(init, axis=0).T).sum() / (n - 1)
    pre, post = a - ta * h, b + tb_in * h
    dt = np.linspace(d_a, d_b, n)
    q = prm["q_curv"] * prm["bend_scale"]
    c_b = 1.0 / (2 * q * h ** 3 * GRID)
    sd = prm["sigma_depth_um"]
    c_d = 0.0 if not np.isfinite(sd) else h * GRID / (2 * sd ** 2 * prm["corr_depth_um"])
    c_s = prm["w_spacing"] / h ** 2
    c_r, m = prm["w_barrier"], prm["barrier_px"]

    def full(z):
        return np.vstack([pre, a, z.reshape(-1, 2), b, post])

    def fg(z):
        X = full(z)
        G = np.zeros_like(X)
        D = X[2:] - 2 * X[1:-1] + X[:-2]
        e = c_b * (D ** 2).sum()
        G[:-2] += 2 * c_b * D
        G[1:-1] -= 4 * c_b * D
        G[2:] += 2 * c_b * D
        C = X[1:-1]
        dv = np.diff(C, axis=0)
        seg = np.maximum(np.hypot(*dv.T), 1e-9)
        r = seg - h
        e += c_s * (r ** 2).sum()
        u = dv / seg[:, None] * (2 * c_s * r)[:, None]
        G[2:-1] += u
        G[1:-2] -= u
        if c_d:
            dep = Field.at(field.dep, C)
            res = dep - dt
            e += c_d * (res ** 2).sum()
            G[1:-1] += (2 * c_d * res)[:, None] * np.stack([Field.at(field.dgy, C), Field.at(field.dgx, C)], 1)
        fd = Field.at(field.fd, C)
        pen = np.clip(m - fd, 0, None)
        e += c_r * (pen ** 2).sum()
        G[1:-1] -= (2 * c_r * pen)[:, None] * np.stack([Field.at(field.fgy, C), Field.at(field.fgx, C)], 1)
        return e, G[2:-2].ravel()

    r = optimize.minimize(fg, init[1:-1].ravel(), jac=True, method="L-BFGS-B", options=dict(maxiter=prm["maxiter"]))
    if return_energy:
        return full(r.x)[1:-1], float(r.fun)
    return full(r.x)[1:-1]


def energy_terms(C, a, ta, b, tb_in, field, d_a, d_b, prm=None):
    """the energy of a finished curve, term by term (same formulas as elastic)."""
    prm = {**DEF, **(prm or {})}
    C = np.asarray(C, float)
    n = len(C)
    h = np.hypot(*np.diff(C, axis=0).T).mean()
    X = np.vstack([a - ta * h, C, b + tb_in * h])
    D = X[2:] - 2 * X[1:-1] + X[:-2]
    q = prm["q_curv"] * prm["bend_scale"]
    bend = float((D ** 2).sum() / (2 * q * h ** 3 * GRID))
    dep = Field.at(field.dep, C)
    res = dep - np.linspace(d_a, d_b, n)
    sd = prm["sigma_depth_um"]
    depth = 0.0 if not np.isfinite(sd) else float(h * GRID * (res ** 2).sum() / (2 * sd ** 2 * prm["corr_depth_um"]))
    fd = Field.at(field.fd, C)
    barrier = float(prm["w_barrier"] * (np.clip(prm["barrier_px"] - fd, 0, None) ** 2).sum())
    return dict(bend=bend, depth=depth, barrier=barrier, depth_rms_um=float(np.sqrt((res ** 2).mean())),
                length_um=float(h * (n - 1) * GRID))


def resample_n(C, n):
    C = np.asarray(C, float)
    s = np.concatenate([[0], np.cumsum(np.hypot(*np.diff(C, axis=0).T))])
    t = np.linspace(0, s[-1], n)
    return np.stack([np.interp(t, s, C[:, k]) for k in range(2)], 1)


def geodesic(a, b, field, d_hat, prm=None, avoid=None):
    """least-cost route from a to b through tissue that is not EGL, white matter or
    background, where each step costs 1 + lam * ((depth below EGL - d_hat) / sigma)^2:
    the route of the ribbon's own depth. avoid: pixels the route may not use (e.g. the
    rest of the IGL, so the route cannot run along the ribbon it is meant to complete);
    returns None when no route exists inside the window."""
    from skimage.graph import MCP_Geometric
    prm = {**DEF, **(prm or {})}
    H, W = field.dep.shape
    r = int(prm["geo_window_px"])
    c = np.round((a + b) / 2).astype(int)
    y0, x0 = max(c[0] - r, 0), max(c[1] - r, 0)
    y1, x1 = min(c[0] + r + 1, H), min(c[1] + r + 1, W)
    dep = field.dep_raw[y0:y1, x0:x1]
    cost = 1.0 + prm["geo_lambda"] * ((dep - d_hat) / prm["sigma_depth_um"]) ** 2
    bad = field.forbid[y0:y1, x0:x1].copy()
    if avoid is not None:
        bad |= avoid[y0:y1, x0:x1]
    ia = tuple(np.round(a).astype(int) - [y0, x0])
    ib = tuple(np.round(b).astype(int) - [y0, x0])
    if not all(0 <= v < s for p in (ia, ib) for v, s in zip(p, bad.shape)):
        return None
    yy, xx = np.ogrid[:bad.shape[0], :bad.shape[1]]
    for p in (ia, ib):
        bad &= (yy - p[0]) ** 2 + (xx - p[1]) ** 2 > 4
    cost[bad] = np.inf
    m = MCP_Geometric(cost)
    cum, _ = m.find_costs([ia], [ib])
    if not np.isfinite(cum[ib]):
        return None
    path = np.array(m.traceback(ib), float) + [y0, x0]
    return path


def complete(a, ta, b, tb_in, field, d_a, d_b, prm=None, avoid=None):
    """the model's prediction of a missing stretch: elastic refinement started from the
    smooth (Hermite) curve and from the depth route (geodesic); the lower-energy result
    wins. Returns (curve, info)."""
    prm = {**DEF, **(prm or {})}
    out = []
    h0 = hermite(a, ta, b, tb_in, 50)
    starts = [("hermite", h0)]
    g = geodesic(a, b, field, (d_a + d_b) / 2, prm, avoid)
    if g is not None and len(g) >= 4:
        starts.append(("geodesic", g))
    for name, C0 in starts:
        lc = np.hypot(*np.diff(C0, axis=0).T).sum()
        n = int(np.clip(round(lc / prm["spacing_px"]), 6, 160))
        init = resample_n(C0, n)
        if name == "geodesic":
            init = np.vstack([a, ndi.gaussian_filter1d(init, 2.0, axis=0, mode="nearest")[1:-1], b])
        C, e = elastic(a, ta, b, tb_in, field, d_a, d_b, prm, init=init, return_energy=True)
        out.append((e, name, C))
    e, name, C = min(out, key=lambda t: t[0])
    return C, dict(start=name, energy=e, energies={n: round(v, 2) for v, n, _ in out})


def hermite_width(wa, wb, n):
    s = np.linspace(0, 1, n)
    return wa + (wb - wa) * (3 * s ** 2 - 2 * s ** 3)


# q_curv: curvature intensity from the intact ribbons (sd of kappa ~0.012 /um, correlation
# length ~12 um -> 0.012^2 * 12); corr_depth_um: correlation length of the depth below the EGL.
DEF = dict(spacing_px=2.0, q_curv=1.7e-3, bend_scale=1.0, sigma_depth_um=15.0, corr_depth_um=55.0, w_spacing=10.0,
           w_barrier=50.0, barrier_px=2.0, maxiter=300, geo_lambda=1.0, geo_window_px=130)
