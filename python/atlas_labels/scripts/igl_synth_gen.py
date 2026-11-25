"""
Synthetic P7 cerebellum sections with an IGL, from three families of generative models.

    python scripts/igl_synth_gen.py run <family> <n> [workers]   # sample the prior, write outputs/igl_synth/samples_<family>.jsonl
    python scripts/igl_synth_gen.py refine <family> <IF|Allen> <n>  # local search around the sensible parameter sets (after igl_synth_abc.py)
    python scripts/igl_synth_gen.py demo                          # one default section per family, with its statistics

Families
  fold     folded-surface laminar model. A dome-shaped outline (ellipse, base cut) has
           fissures cut into it, running from the surface toward a core near the base
           (main fissures deep, minor ones shallow, each with a bend), open as a wedge
           of width 'open' at the surface and closed at the bottom.
           Depth below the pial surface (fissures included) sets the layers:
             EGL  0 .. e,  ML  e .. e+m,  IGL  e+m .. e+m+g,  white matter deeper.
           e, m, g vary as smooth random fields (CV, correlation length), with the local
           folding (k_m, k_g: thicker at fissure bases if > 0, at folium tips if < 0), and
           with cell-scale edge roughness ('rough': relative amplitude, 12 um correlation,
           applied to each layer's thickness separately, which can pinch the IGL), and in
           the later families (fold2, finger) with coherent depth jitter ('jit', um): all
           layer boundaries shift together, so the IGL wiggles without pinching.
  arbor    white-matter tree model. A tree of white matter (lobule branches from a core,
           some splitting into folia, tapering) is surrounded by tissue; fissures run
           between neighbouring branches down to a set depth. EGL and ML are laminae
           under the pial surface; the IGL fills everything between the ML and the
           white-matter tree (its thickness is not set directly).
  fold2    fold, plus two structural options: fissures aim at an inner core ellipse of
           relative size core_r (0 = one core point, so deep fissures converge), and in
           narrow folia EGL + ML take at most (1 - q_core) of the local half-width (surface
           to medial axis), so every folium keeps an IGL core.
  foldma   folded surface (as fold, with the core ellipse) but the white matter is a tree
           along the folia's medial axis: present where the half-width h >= h_min, radius
           rho * (h - e - m); the IGL fills everything between the ML and that tree, so its
           thickness follows the room left in each folium instead of being a parameter.
  finger   lobules as parallel-sided fingers (capsules) on a central core, some splitting
           into sub-folia at the tip; fissures are the gaps between neighbouring fingers;
           layers by depth below the pial surface as in fold.
  fold3    fold2 with two corrections: the half-width is the local thickness (radius of the
           largest disk inside the tissue that contains the pixel), not the medial-axis
           distance, whose spurs squeezed the EGL away near the surface; and the depth jitter
           moves only the ML/IGL/WM boundaries, so the EGL stays on the surface.
  fourier  smooth-surface laminar model. The outline is a circle deformed by random
           Fourier modes (power-law spectrum), with no fissures; layers by depth as in fold.
Every family renders masks on the 4 um grid (tissue, egl, ml, igl, wm), measured with
igl_synth_stats.summarize.
"""
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi
from skimage.draw import line as draw_line

sys.path.insert(0, str(Path(__file__).parent))
import igl_synth_stats as S  # noqa: E402

G = S.G
OUT = S.OUT
CANVAS = 1000                    # 4 mm


# --------------------------------------------------------------------------- #
def rand_field(shape, ell_um, rng):
    f = ndi.gaussian_filter(rng.standard_normal(shape), max(ell_um / G, 1.0))
    return f / max(f.std(), 1e-9)


def local_thickness(mask, step=1.2):
    """radius (px) of the largest disk inside mask that contains each pixel"""
    D = ndi.distance_transform_edt(mask)
    T = np.zeros(mask.shape)
    r = 1.5
    while r <= D.max():
        reach = ndi.distance_transform_edt(~(D >= r)) <= r
        T[reach & mask] = r
        r *= step
    return T


def layers(tissue, pial_bg, p, rng, wm_override=None):
    """EGL / ML / IGL / WM by depth below the pial surface, with smooth variation."""
    dS = ndi.distance_transform_edt(~pial_bg) * G
    dSj = dS
    if p.get("jit", 0.0) > 0:
        # coherent depth jitter: the layer boundaries shift together by a smooth random
        # amount (24 um correlation), so the IGL wiggles without pinching. fold2/finger shift
        # the EGL boundary too (which can strip the EGL off the surface); jit_mode "deep"
        # (fold3) keeps the EGL on the surface and shifts only the deeper boundaries.
        dSj = dS + p["jit"] * rand_field(tissue.shape, 24.0, rng)
        if p.get("jit_mode") != "deep":
            dS = dSj
    loc = ndi.uniform_filter(tissue.astype(float), int(2 * 150 / G) + 1)
    z = np.clip((loc - 0.6) / 0.2, -2, 2)                     # > 0 fissure base (concave), < 0 folium tip
    sh = tissue.shape
    rough = p.get("rough", 0.0)
    rf = lambda: rough * rand_field(sh, 12.0, rng) if rough > 0 else 0.0     # cell-scale edge roughness
    e = p["e"] * np.clip(1 + p["cv_e"] * rand_field(sh, p["ell"], rng) + rf(), 0.3, 3)
    m = p["m"] * np.clip(1 + p["cv_m"] * rand_field(sh, p["ell"], rng) + p["k_m"] * z + rf(), 0.3, 3)
    q = p.get("q_core", 0.0)
    if q > 0:
        # narrow folia: EGL + ML take at most (1 - q) of the local half-width h (surface to
        # medial axis), so the IGL keeps at least q of it
        if p.get("hw_mode") == "thick":
            hw = local_thickness(tissue) * G          # fold3: largest inscribed disk containing the pixel
        else:
            # fold2: half-width at the nearest medial-axis point; the axis has short spurs
            # toward every bump of the outline, so near the surface this is often far too
            # small and the EGL and ML get squeezed away
            from skimage.morphology import medial_axis
            ax, dist = medial_axis(tissue, return_distance=True)
            _, (iy, ix) = ndi.distance_transform_edt(~ax, return_indices=True)
            hw = dist[iy, ix] * G
        lam = np.clip((1 - q) * hw / np.maximum(e + m, 1e-6), 0, 1)
        e, m = e * lam, m * lam
    egl = tissue & (dS < e)
    ml = tissue & ~egl & (dSj < e + m)
    if p.get("wm_mode") == "medial":
        # white matter along the folia's medial axis: only where the half-width h >= h_min,
        # with radius rho * (h - e - m); the IGL fills everything between the ML and it
        from skimage.morphology import medial_axis
        ax, dist = medial_axis(tissue, return_distance=True)
        hw = dist * G
        ax_ok = ax & (hw >= p["h_min"])
        if ax_ok.any():
            rw = np.where(ax_ok, np.clip(p["rho"] * (hw - e - m), 0, None), 0.0)
            dax, (iy, ix) = ndi.distance_transform_edt(~ax_ok, return_indices=True)
            wm = tissue & ~egl & ~ml & (dax * G <= rw[iy, ix])
        else:
            wm = np.zeros_like(tissue)
        igl = tissue & ~egl & ~ml & ~wm
        z0 = np.zeros(sh, bool)
        return dict(tissue=tissue, egl=egl, ml=ml, igl=igl, wm=wm, anom=z0, unass=z0.copy())
    if wm_override is None:
        gg = p["g"] * np.clip(1 + p["cv_g"] * rand_field(sh, p["ell"], rng) + p["k_g"] * z + rf(), 0.3, 3)
        igl = tissue & ~egl & ~ml & (dSj < e + m + gg)
        wm = tissue & ~egl & ~ml & ~igl
    else:
        wm = tissue & wm_override & ~egl & ~ml
        igl = tissue & ~egl & ~ml & ~wm
    z0 = np.zeros(sh, bool)
    return dict(tissue=tissue, egl=egl, ml=ml, igl=igl, wm=wm, anom=z0, unass=z0.copy())


def crop(g, margin=10):
    ys, xs = np.nonzero(g["tissue"])
    if not len(ys):
        return g
    sl = (slice(max(ys.min() - margin, 0), ys.max() + margin + 1), slice(max(xs.min() - margin, 0), xs.max() + margin + 1))
    return {k: v[sl] for k, v in g.items()}


def raster_curve(mask, P):
    for p0, p1 in zip(P[:-1], P[1:]):
        rr, cc = draw_line(int(round(p0[0])), int(round(p0[1])), int(round(p1[0])), int(round(p1[1])))
        ok = (rr >= 0) & (rr < mask.shape[0]) & (cc >= 0) & (cc < mask.shape[1])
        mask[rr[ok], cc[ok]] = True


def raster_wedge(mask, P, open_px):
    """fissure that is open_px wide at the surface (P[0]) and closes to 1 px at its bottom."""
    from skimage.draw import disk
    raster_curve(mask, P)
    n = len(P)
    for i, q in enumerate(P):
        r = open_px / 2 * (1 - i / max(n - 1, 1)) ** 1.5
        if r >= 0.75:
            rr, cc = disk((q[0], q[1]), r, shape=mask.shape)
            mask[rr, cc] = True


def bezier(a, b, bend, n=60):
    d = b - a
    nrm = np.array([-d[1], d[0]])
    c = (a + b) / 2 + bend * nrm / 2
    t = np.linspace(0, 1, n)[:, None]
    return (1 - t) ** 2 * a + 2 * (1 - t) * t * c + t ** 2 * b


# --------------------------------------------------------------------------- #
# fold
def prior_fold(rng):
    return dict(A=rng.uniform(800, 1700), B=rng.uniform(700, 1500), cut=rng.uniform(-0.1, 0.6), core_h=rng.uniform(0.05, 0.45),
                th_a=rng.uniform(160, 240), th_p=rng.uniform(-60, 20), n_main=int(rng.integers(2, 15)), n_minor=int(rng.integers(0, 31)),
                f_main=rng.uniform(0.3, 0.95), f_minor=rng.uniform(0.03, 0.4), f_sd=rng.uniform(0, 0.25), jitter=rng.uniform(0, 0.45),
                bend=rng.uniform(0, 0.5), open=rng.uniform(4, 160),
                e=rng.uniform(8, 50), m=rng.uniform(20, 100), g=rng.uniform(20, 200),
                cv_e=rng.uniform(0, 0.4), cv_m=rng.uniform(0, 0.5), cv_g=rng.uniform(0, 0.5), ell=rng.uniform(30, 400),
                k_m=rng.uniform(-0.6, 0.6), k_g=rng.uniform(-0.6, 0.6), rough=rng.uniform(0, 0.4))


def gen_fold(p, rng):
    H = W = CANVAS
    cy, cx = H / 2 + 100, W / 2
    A, B = p["A"] / G, p["B"] / G
    yy, xx = np.mgrid[0:H, 0:W]
    ell = ((xx - cx) / A) ** 2 + ((yy - cy) / B) ** 2 <= 1
    y_cut = cy + p["cut"] * B
    top = cy - B
    core = np.array([y_cut - p["core_h"] * (y_cut - top), cx])
    n = p["n_main"] + p["n_minor"]
    base = np.linspace(np.radians(p["th_a"]), np.radians(p["th_p"]), n + 2)[1:-1]
    sp = abs(base[1] - base[0]) if n > 1 else 0.3
    ang = np.sort(base + rng.normal(0, p["jitter"] * sp, n))
    is_main = np.zeros(n, bool)
    is_main[rng.choice(n, size=min(p["n_main"], n), replace=False)] = True
    slits = np.zeros((H, W), bool)
    for th, mn in zip(ang, is_main):
        b = np.array([cy - B * np.sin(th), cx + A * np.cos(th)])
        if b[0] >= y_cut - 5:
            continue
        cr = p.get("core_r", 0.0)                       # fissures aim at an inner core ellipse, not one point
        tgt = core + cr * np.array([-B * np.sin(th), A * np.cos(th)])
        L = np.hypot(*(tgt - b))
        f = np.clip(rng.normal(p["f_main"] if mn else p["f_minor"], p["f_sd"]), 0.02, 0.97)
        u = (tgt - b) / L
        end = b + u * f * L
        P = bezier(b - u * 8, end, rng.uniform(-p["bend"], p["bend"]))
        raster_wedge(slits, P, p["open"] / G * (1.0 if mn else 0.6))
    tissue = ell & ~slits & (yy <= y_cut)
    pial_bg = ~ell | slits
    return crop(layers(tissue, pial_bg, p, rng))


# --------------------------------------------------------------------------- #
# fourier
def prior_fourier(rng):
    return dict(R=rng.uniform(800, 1400), amp=rng.uniform(0.0, 0.25), beta=rng.uniform(0.5, 3.0), kmax=int(rng.integers(6, 40)),
                cut=rng.uniform(-0.1, 0.6), e=rng.uniform(8, 50), m=rng.uniform(20, 100), g=rng.uniform(20, 200),
                cv_e=rng.uniform(0, 0.4), cv_m=rng.uniform(0, 0.5), cv_g=rng.uniform(0, 0.5), ell=rng.uniform(30, 400),
                k_m=rng.uniform(-0.6, 0.6), k_g=rng.uniform(-0.6, 0.6), rough=rng.uniform(0, 0.4))


def gen_fourier(p, rng):
    H = W = CANVAS
    cy, cx = H / 2 + 100, W / 2
    yy, xx = np.mgrid[0:H, 0:W]
    th = np.arctan2(cy - yy, xx - cx)
    r = np.hypot(yy - cy, xx - cx)
    ks = np.arange(2, p["kmax"] + 1)
    a = rng.normal(0, 1, len(ks)) * ks ** (-p["beta"] / 2)
    a = a / max(np.sqrt((a ** 2).sum()), 1e-9) * p["amp"]
    ph = rng.uniform(0, 2 * np.pi, len(ks))
    Rth = p["R"] / G * (1 + (a[:, None, None] * np.cos(ks[:, None, None] * th[None] + ph[:, None, None])).sum(0))
    inside = r <= Rth
    y_cut = cy + p["cut"] * p["R"] / G
    tissue = inside & (yy <= y_cut)
    return crop(layers(tissue, ~inside, p, rng))


# --------------------------------------------------------------------------- #
# arbor
def prior_arbor(rng):
    return dict(n_lob=int(rng.integers(4, 12)), th_a=rng.uniform(160, 240), th_p=rng.uniform(-60, 20), jitter=rng.uniform(0, 0.4),
                L=rng.uniform(500, 1500), L_sd=rng.uniform(0, 0.4), p_split=rng.uniform(0, 0.9), split_at=rng.uniform(0.3, 0.8),
                split_ang=rng.uniform(10, 45), r0=rng.uniform(30, 150), r_tip=rng.uniform(4, 40), cover=rng.uniform(80, 300),
                fiss=rng.uniform(0.2, 1.0), open=rng.uniform(4, 160), cut=rng.uniform(0.0, 0.5),
                e=rng.uniform(8, 50), m=rng.uniform(20, 100), g=0.0,
                cv_e=rng.uniform(0, 0.4), cv_m=rng.uniform(0, 0.5), cv_g=0.0, ell=rng.uniform(30, 400), k_m=rng.uniform(-0.6, 0.6), k_g=0.0,
                rough=rng.uniform(0, 0.4))


def gen_arbor(p, rng):
    H = W = CANVAS
    root = np.array([H / 2 + 250, W / 2])
    segs = []                                  # (a, b, r_a, r_b, leaf_id)
    n = p["n_lob"]
    base = np.linspace(np.radians(p["th_a"]), np.radians(p["th_p"]), n + 2)[1:-1]
    sp = abs(base[1] - base[0]) if n > 1 else 0.3
    ang = base + rng.normal(0, p["jitter"] * sp, n)
    leaf = 0
    for th in ang:
        L = p["L"] / G * max(0.3, 1 + rng.normal(0, p["L_sd"]))
        u = np.array([-np.sin(th), np.cos(th)])
        r_of = lambda s: (p["r0"] + (p["r_tip"] - p["r0"]) * s) / G
        if rng.random() < p["p_split"]:
            s0 = p["split_at"]
            mid = root + u * L * s0
            segs.append((root, mid, r_of(0), r_of(s0), None))
            for sgn in (1, -1):
                t = th + sgn * np.radians(p["split_ang"])
                v = np.array([-np.sin(t), np.cos(t)])
                segs.append((mid, mid + v * L * (1 - s0), r_of(s0), r_of(1), leaf))
                leaf += 1
        else:
            segs.append((root, root + u * L, r_of(0), r_of(1), leaf))
            leaf += 1
    # skeleton raster with radius and leaf identity
    sk = np.zeros((H, W), bool)
    rad = np.zeros((H, W))
    lid = np.full((H, W), -1)
    for a, b, ra, rb, lf in segs:
        k = int(np.hypot(*(b - a))) + 2
        for t in np.linspace(0, 1, k):
            q = np.round(a + (b - a) * t).astype(int)
            if 0 <= q[0] < H and 0 <= q[1] < W:
                sk[q[0], q[1]] = True
                rad[q[0], q[1]] = ra + (rb - ra) * t
                lid[q[0], q[1]] = lf if lf is not None else -2
    d, (iy, ix) = ndi.distance_transform_edt(~sk, return_indices=True)
    r_near = rad[iy, ix]
    wm = d <= r_near
    body = d <= r_near + p["cover"] / G
    # fissures: boundaries between different leaves (and leaves vs trunk), from the surface
    # down to a fraction 'fiss' of the cover thickness
    L_near = lid[iy, ix]
    bnd = np.zeros((H, W), bool)
    for ax in (0, 1):
        diff = np.diff(L_near, axis=ax) != 0
        sl = (slice(0, -1), slice(None)) if ax == 0 else (slice(None), slice(0, -1))
        bnd[sl] |= diff
    depth_from_wm = d - r_near
    slits = bnd & body & (depth_from_wm > (1 - p["fiss"]) * p["cover"] / G) & (L_near != -2)
    r_open = p["open"] / 2 / G
    if r_open >= 1:
        near_out = depth_from_wm > p["cover"] / G - 2 * r_open
        slits |= ndi.binary_dilation(slits, iterations=int(round(r_open))) & near_out & body
    y_cut = root[0] - p["cut"] * 250
    yy = np.arange(H)[:, None]
    tissue = body & ~slits & (yy <= y_cut)
    pial_bg = ~body | slits
    return crop(layers(tissue, pial_bg, p, rng, wm_override=wm))


def prior_fold2(rng):
    p = prior_fold(rng)
    p["q_core"] = rng.uniform(0.0, 0.7)
    p["core_r"] = rng.uniform(0.0, 0.5)
    p["rough"] = rng.uniform(0.0, 0.1)
    p["jit"] = rng.uniform(0.0, 25.0)
    return p


def prior_fold3(rng):
    p = prior_fold2(rng)
    p["hw_mode"] = "thick"
    p["jit_mode"] = "deep"
    return p


# finger: lobules as parallel-sided capsules on a central core
def prior_finger(rng):
    return dict(Rc=rng.uniform(250, 800), n_lob=int(rng.integers(5, 13)), th_a=rng.uniform(170, 240), th_p=rng.uniform(-60, 10),
                jitter=rng.uniform(0, 0.35), L=rng.uniform(250, 1300), L_sd=rng.uniform(0, 0.5), W=rng.uniform(90, 350),
                W_sd=rng.uniform(0, 0.4), p_sub=rng.uniform(0, 0.9), n_sub=int(rng.integers(2, 4)), sub_len=rng.uniform(0.15, 0.6),
                sub_ang=rng.uniform(5, 40), gap=rng.uniform(4, 40), cut=rng.uniform(0.1, 0.9),
                e=rng.uniform(8, 50), m=rng.uniform(20, 100), g=rng.uniform(20, 200),
                cv_e=rng.uniform(0, 0.4), cv_m=rng.uniform(0, 0.5), cv_g=rng.uniform(0, 0.5), ell=rng.uniform(30, 400),
                k_m=rng.uniform(-0.6, 0.6), k_g=rng.uniform(-0.6, 0.6), rough=rng.uniform(0, 0.1), jit=rng.uniform(0, 25))


def _capsule(shape, a, b, r):
    yy, xx = np.ogrid[:shape[0], :shape[1]]
    d = b - a
    L2 = max(float(d @ d), 1e-9)
    t = np.clip(((yy - a[0]) * d[0] + (xx - a[1]) * d[1]) / L2, 0, 1)
    py, px = a[0] + t * d[0], a[1] + t * d[1]
    return (yy - py) ** 2 + (xx - px) ** 2 <= r * r


def gen_finger(p, rng):
    """lobule k is a capsule of half-width W_k from the core outward (length L_k); with
    probability p_sub its tip splits into n_sub narrower sub-folia fanning by sub_ang.
    Each finger is its own piece of tissue; where neighbours do not overlap, the space
    between them is a fissure (pial surface); where they overlap they are kept apart by a
    fissure line of width 'gap' along their Voronoi boundary, down to the core circle."""
    H = W = CANVAS
    c = np.array([H / 2 + 150, W / 2])
    Rc = p["Rc"] / G
    n = p["n_lob"]
    base = np.linspace(np.radians(p["th_a"]), np.radians(p["th_p"]), n + 2)[1:-1]
    sp = abs(base[1] - base[0]) if n > 1 else 0.3
    ang = base + rng.normal(0, p["jitter"] * sp, n)
    owner = np.full((H, W), -1)
    best_d = np.full((H, W), np.inf)
    yy, xx = np.mgrid[0:H, 0:W]
    pieces = []
    for k, th in enumerate(ang):
        u = np.array([-np.sin(th), np.cos(th)])
        L = p["L"] / G * max(0.25, 1 + rng.normal(0, p["L_sd"]))
        Wk = p["W"] / G * max(0.4, 1 + rng.normal(0, p["W_sd"]))
        a = c + u * Rc * 0.6
        b = c + u * (Rc + L)
        segs = [(a, b, Wk)]
        if rng.random() < p["p_sub"]:
            ns = int(round(p["n_sub"]))
            tip_len = p["sub_len"] * L
            b0 = c + u * (Rc + L - tip_len)
            segs = [(a, b0, Wk)]
            for j in range(ns):
                off = (j - (ns - 1) / 2)
                t2 = th + np.radians(p["sub_ang"]) * off
                v = np.array([-np.sin(t2), np.cos(t2)])
                segs.append((b0 + np.array([-v[1], v[0]]) * off * Wk / ns, b0 + v * tip_len, Wk / ns * 1.15))
        m = np.zeros((H, W), bool)
        for (sa, sb, r) in segs:
            m |= _capsule((H, W), sa, sb, r)
        # distance to this lobule's axis (for the Voronoi split between overlapping lobules)
        dk = np.full((H, W), np.inf)
        for (sa, sb, r) in segs:
            d = sb - sa
            L2 = max(float(d @ d), 1e-9)
            t = np.clip(((yy - sa[0]) * d[0] + (xx - sa[1]) * d[1]) / L2, 0, 1)
            dk = np.minimum(dk, np.hypot(yy - (sa[0] + t * d[0]), xx - (sa[1] + t * d[1])) - r)
        upd = m & (dk < best_d)
        owner[upd] = k
        best_d[upd] = dk[upd]
        pieces.append(m)
    core = (yy - c[0]) ** 2 + (xx - c[1]) ** 2 <= Rc ** 2
    body = core | np.any(pieces, axis=0)
    # fissures: owner changes between lobules outside the core
    bnd = np.zeros((H, W), bool)
    o = owner
    bnd[:-1] |= (o[:-1] != o[1:]) & (o[:-1] >= 0) & (o[1:] >= 0)
    bnd[:, :-1] |= (o[:, :-1] != o[:, 1:]) & (o[:, :-1] >= 0) & (o[:, 1:] >= 0)
    rr = np.hypot(yy - c[0], xx - c[1])
    bnd &= rr > Rc * 0.98
    gp = p["gap"] / G
    slits = ndi.binary_dilation(bnd, iterations=max(int(round(gp / 2)), 0)) if gp >= 2 else bnd
    y_cut = c[0] + Rc * (1 - 2 * p["cut"]) * 0.9
    tissue = body & ~slits & (yy <= y_cut)
    pial_bg = ~body | slits
    return crop(layers(tissue, pial_bg, p, rng))


def prior_foldma(rng):
    p = prior_fold(rng)
    for k in ("g", "cv_g", "k_g"):
        p[k] = 0.0
    p.update(wm_mode="medial", rho=rng.uniform(0.0, 0.9), h_min=rng.uniform(40, 300), core_r=rng.uniform(0.0, 0.5))
    return p


FAMILIES = {"foldma": (prior_foldma, gen_fold), "fold": (prior_fold, gen_fold), "fourier": (prior_fourier, gen_fourier), "arbor": (prior_arbor, gen_arbor),
            "fold2": (prior_fold2, gen_fold), "finger": (prior_finger, gen_finger), "fold3": (prior_fold3, gen_fold)}


def sample(args):
    fam, seed = args
    rng = np.random.default_rng(seed)
    prior, gen = FAMILIES[fam]
    p = prior(rng)
    try:
        g = gen(p, np.random.default_rng(seed + 10 ** 7))
        s = S.summarize(g)
    except Exception as ex:                       # degenerate geometry
        s = None
    return dict(family=fam, seed=int(seed), params=p, stats=s)


def sample_given(args):
    fam, seed, p, parent = args
    try:
        g = FAMILIES[fam][1](p, np.random.default_rng(seed + 10 ** 7))
        st = S.summarize(g)
    except Exception:
        st = None
    return dict(family=fam, seed=int(seed), params=p, stats=st, parent=parent)


INT_PARAMS = {"n_main", "n_minor", "kmax", "n_lob", "n_sub"}


def refine(fam, target, n, workers=9):
    """ABC-SMC style local search: perturb the parameter sets that were sensible for a target
    (half their spread per parameter, kept inside the prior range) and sample again."""
    from multiprocessing import Pool
    abc = json.loads((OUT / "abc.json").read_text())
    ft = abc["families"][fam]["targets"][target]
    thr = abc["targets"][target]["threshold"]
    d = {int(k): v for k, v in ft["distances"].items()}
    rows = {r["seed"]: r for r in (json.loads(l) for l in (OUT / f"samples_{fam}.jsonl").read_text().splitlines() if l.strip())}
    acc = [s_ for s_, v in d.items() if v <= thr]
    if len(acc) < 30:
        acc = sorted(d, key=d.get)[:30]
    keys = [k for k in rows[acc[0]]["params"] if isinstance(rows[acc[0]]["params"][k], (int, float))]
    Pm = np.array([[rows[s_]["params"][k] for k in keys] for s_ in acc], float)
    sd = 0.5 * Pm.std(0) + 1e-9
    bounds = {k: ft["posterior"][k]["prior"] for k in keys}
    rng = np.random.default_rng(777 + len(target))
    jobs = []
    base_seed = 2_000_000 + (0 if target == "IF" else 1_000_000)
    for i in range(n):
        par = acc[rng.integers(len(acc))]
        p = dict(rows[par]["params"])
        for j, k in enumerate(keys):
            v = p[k] + rng.normal(0, sd[j])
            lo, hi = bounds[k]
            v = float(np.clip(v, lo, hi))
            p[k] = int(round(v)) if k in INT_PARAMS else v
        jobs.append((fam, base_seed + i, p, int(par)))
    path = OUT / f"samples_{fam}_refine_{target}.jsonl"
    t0 = time.time()
    (OUT / f"refine_{target}_family.txt").write_text(fam)
    with Pool(int(workers)) as pool, open(path, "w") as f:
        for i, r in enumerate(pool.imap_unordered(sample_given, jobs, chunksize=4), 1):
            f.write(json.dumps(r, default=float) + "\n")
            if i % 200 == 0:
                f.flush()
                print(f"refine {fam}/{target}: {i}/{n} in {time.time() - t0:.0f}s", flush=True)
    print(f"refine {fam}/{target}: done {n} from {len(acc)} parents in {time.time() - t0:.0f}s", flush=True)


def refine_iter(fam, target, n, src_tags, out_tag, workers=9, keep=30, shrink=0.5):
    """another local-search round: restart from the closest sets found so far (in the prior
    samples and the earlier refinement files), with steps of 'shrink' x their spread."""
    from multiprocessing import Pool
    import igl_synth_abc as ABC
    real = json.loads((OUT / "real_stats.json").read_text())
    med, sc = ABC.target_model([r for r in real if r["group"] in ABC.TARGETS[target]])
    rows = []
    for tag in src_tags:
        p = OUT / (f"samples_{fam}.jsonl" if tag == "prior" else f"samples_{fam}_{tag}_{target}.jsonl")
        if p.exists():
            rows += [json.loads(l) for l in p.read_text().splitlines() if l.strip()]
    rows = [r for r in rows if r["stats"]]
    X = np.array([[r["stats"][k] for k in ABC.USE] for r in rows], float)
    d, _ = ABC.dist(X, med, sc)
    order = np.argsort(d)[:keep]
    parents = [rows[i] for i in order]
    keys = [k for k in parents[0]["params"] if isinstance(parents[0]["params"][k], (int, float)) and not isinstance(parents[0]["params"][k], bool)]
    prior_rows = [json.loads(l) for l in (OUT / f"samples_{fam}.jsonl").read_text().splitlines() if l.strip()]
    lo = {k: min(r["params"][k] for r in prior_rows) for k in keys}
    hi = {k: max(r["params"][k] for r in prior_rows) for k in keys}
    Pm = np.array([[r["params"][k] for k in keys] for r in parents], float)
    sd = shrink * Pm.std(0) + 1e-9
    rnd = int(out_tag[-1]) if out_tag[-1].isdigit() else 0      # refine2 -> 4.2M, refine3 -> 4.3M ...
    rng = np.random.default_rng(4242 + 10 * rnd + len(target))
    base_seed = 4_000_000 + 100_000 * rnd + {"IF": 0, "Allen": 50_000}.get(target, 75_000)
    jobs = []
    for i in range(n):
        par = parents[rng.integers(len(parents))]
        p = dict(par["params"])
        for j, k in enumerate(keys):
            v = float(np.clip(p[k] + rng.normal(0, sd[j]), lo[k], hi[k]))
            p[k] = int(round(v)) if k in INT_PARAMS else v
        jobs.append((fam, base_seed + i, p, int(par["seed"])))
    print(f"refine_iter {fam}/{target}/{out_tag}: {len(rows)} candidates, parents' distances {d[order].min():.2f}-{d[order].max():.2f}", flush=True)
    path = OUT / f"samples_{fam}_{out_tag}_{target}.jsonl"
    t0 = time.time()
    with Pool(int(workers)) as pool, open(path, "w") as f:
        for i, r in enumerate(pool.imap_unordered(sample_given, jobs, chunksize=4), 1):
            f.write(json.dumps(r, default=float) + "\n")
            if i % 500 == 0:
                f.flush()
    rows2 = [json.loads(l) for l in path.read_text().splitlines() if l.strip()]
    rows2 = [r for r in rows2 if r["stats"]]
    X2 = np.array([[r["stats"][k] for k in ABC.USE] for r in rows2], float)
    d2, _ = ABC.dist(X2, med, sc)
    print(f"refine_iter {fam}/{target}/{out_tag}: done in {time.time() - t0:.0f}s; best {d2.min():.2f}, median {np.median(d2):.2f}", flush=True)


def run(fam, n, workers=9):
    from multiprocessing import Pool
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / f"samples_{fam}.jsonl"
    done = set()
    if path.exists():
        done = {json.loads(l)["seed"] for l in path.read_text().splitlines() if l.strip()}
    seeds = [s for s in range(1, n + 1) if s not in done]
    t0 = time.time()
    with Pool(int(workers)) as pool, open(path, "a") as f:
        for i, r in enumerate(pool.imap_unordered(sample, [(fam, s) for s in seeds], chunksize=4), 1):
            f.write(json.dumps(r, default=float) + "\n")
            if i % 200 == 0:
                f.flush()
                print(f"{fam}: {i}/{len(seeds)} in {time.time() - t0:.0f}s", flush=True)
    print(f"{fam}: done {len(seeds)} in {time.time() - t0:.0f}s", flush=True)


def demo():
    for fam, (prior, gen) in FAMILIES.items():
        t = time.time()
        r = sample((fam, 12345))
        s = r["stats"]
        print(fam, f"{time.time() - t:.1f}s", None if s is None else {k: round(v, 3) for k, v in s.items()})


if __name__ == "__main__":
    if sys.argv[1] == "run":
        run(sys.argv[2], int(sys.argv[3]), int(sys.argv[4]) if len(sys.argv) > 4 else 9)
    elif sys.argv[1] == "refine_iter":
        refine_iter(sys.argv[2], sys.argv[3], int(sys.argv[4]), sys.argv[5].split(","), sys.argv[6])
    elif sys.argv[1] == "refine":
        refine(sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]) if len(sys.argv) > 5 else 9)
    else:
        demo()
