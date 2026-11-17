"""
Laplace-based ribbon unfolding of a P7 cerebellum section.

Domain Ω = cortical tissue (labels 2,3,4,5,7) plus a 1-pixel ring of
white-matter / core tissue (labels 1, 6) used as the inner Dirichlet BC.

Boundary conditions:
  ψ = 0 on the pial ring   (cortex pixels with a background 4-neighbour)
  ψ = 1 on the inner ring  (label 1 / 6 pixels with a cortex 4-neighbour)
  Neumann (no flux) on any other domain boundary edge.

Solve ∇²ψ = 0  →  thickness-normalised depth field, immune to the
opposing-pia problem because gradient lines stay inside the tissue.

Arclength s is assigned by tracing each cortical pixel back along −∇ψ
until it hits the pial ring, then looking up the arclength of the
nearest point on the smoothed pial contour.
"""

from pathlib import Path
import numpy as np
import cv2
from scipy import ndimage as ndi
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve
from scipy.spatial import cKDTree
from PIL import Image
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

LBL = dict(BG=0, CEREB=1, iEGL=2, oEGL=3, IGL=4, ML=5, DWL=6, PCL=7,
           DCN=8, ANOMALY=9)
# Canonical anatomical convention (sandbox files): 4 = IGL, 5 = ML.
# Expected depth order from pia → core: oEGL < iEGL < ML(5) < PCL(7) < IGL(4)
#                                       < DWL(6) ≈ DCN(8).
CORTEX_LABELS = (LBL["iEGL"], LBL["oEGL"], LBL["IGL"], LBL["ML"], LBL["PCL"])
INNER_LABELS = (LBL["CEREB"], LBL["DWL"], LBL["DCN"])
# Tissue = anything except background and anomaly (tears/bleaches are holes,
# not tissue, and would distort the Laplace solve if included).
TISSUE_LABELS = (1, 2, 3, 4, 5, 6, 7, 8)

COLORS = {
    0: (0.00, 0.00, 0.00),   # bg              — black
    1: (0.10, 0.10, 0.40),   # cereb           — dark navy
    2: (0.40, 0.70, 0.95),   # iEGL            — light blue
    3: (0.00, 0.80, 0.90),   # oEGL            — cyan
    4: (1.00, 0.95, 0.20),   # IGL             — yellow
    5: (0.55, 0.90, 0.55),   # ML              — pale green
    6: (1.00, 0.55, 0.10),   # DWL             — orange
    7: (0.95, 0.10, 0.10),   # PCL             — red
    8: (0.60, 0.30, 0.80),   # DCN             — purple
    9: (0.55, 0.55, 0.55),   # anomaly         — grey (excluded from domain)
}


# --------------------------------------------------------------------------- #
# I/O & geometry
# --------------------------------------------------------------------------- #

def load_labels(path: Path) -> np.ndarray:
    return np.array(Image.open(path)).astype(np.int32)


def downsample(arr: np.ndarray, factor: int) -> np.ndarray:
    return arr[::factor, ::factor]


def tissue_mask(lbl: np.ndarray) -> np.ndarray:
    m = np.isin(lbl, TISSUE_LABELS).astype(np.uint8)
    n, comp, stats, _ = cv2.connectedComponentsWithStats(m, 8)
    biggest = 1 + np.argmax(stats[1:, cv2.CC_STAT_AREA])
    return comp == biggest


def pial_contour_xys(tissue: np.ndarray, smooth_win: int = 21):
    """Outer hull of the tissue (used only for visualisation overlay)."""
    contours, _ = cv2.findContours(
        tissue.astype(np.uint8), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE
    )
    contour = max(contours, key=lambda c: c.shape[0])[:, 0, :].astype(np.float32)
    k = np.ones(smooth_win) / smooth_win
    pad = np.concatenate([contour[-smooth_win:], contour, contour[:smooth_win]])
    sx = np.convolve(pad[:, 0], k, mode="same")[smooth_win:-smooth_win]
    sy = np.convolve(pad[:, 1], k, mode="same")[smooth_win:-smooth_win]
    xy = np.stack([sx, sy], axis=1)
    seg = np.linalg.norm(np.diff(xy, axis=0, append=xy[:1]), axis=1)
    s = np.concatenate([[0.0], np.cumsum(seg)[:-1]])
    return xy, s


def oegl_arclength(oegl_mask: np.ndarray, min_cc_size: int = 200,
                   smooth_win: int = 9, gap_between_components: float = 12.0):
    """Decompose the oEGL into connected components and parameterise each
    by arclength along its outer contour. Every oEGL pixel inside a CC is
    then assigned the arclength of its nearest outer-contour point — so
    the oEGL functions as a thick "ψ=0 marker" with a well-defined s.

    Returns:
      s_field   (H, W) float — s value at every oEGL pixel, -1 elsewhere
      contour_xy   (M, 2)   xy of all contour points (across components)
      contour_s    (M,)     cumulative arclength
      contour_cid  (M,)     component id
      total_s   float
      comp_ranges  list of (cid, s_lo, s_hi, cc_size_px)
    """
    H, W = oegl_mask.shape
    s_field = np.full((H, W), -1.0, dtype=np.float32)

    n, ccl, st, _ = cv2.connectedComponentsWithStats(
        oegl_mask.astype(np.uint8), 8
    )
    sizes = st[1:, cv2.CC_STAT_AREA]
    big_ids = (np.where(sizes >= min_cc_size)[0] + 1)
    order = big_ids[np.argsort(-sizes[big_ids - 1])]

    xys, ss, cids, ranges = [], [], [], []
    s_offset = 0.0
    new_cid = 0

    for ci in order:
        cc_mask = (ccl == ci).astype(np.uint8)
        contours, _ = cv2.findContours(
            cc_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE
        )
        if not contours:
            continue
        contour = max(contours, key=lambda c: c.shape[0])[:, 0, :].astype(np.float32)
        if len(contour) < 4:
            continue

        if len(contour) >= 2 * smooth_win:
            k = np.ones(smooth_win) / smooth_win
            pad = np.concatenate(
                [contour[-smooth_win:], contour, contour[:smooth_win]]
            )
            sx = np.convolve(pad[:, 0], k, mode="same")[smooth_win:-smooth_win]
            sy = np.convolve(pad[:, 1], k, mode="same")[smooth_win:-smooth_win]
            contour = np.stack([sx, sy], axis=1)

        seg = np.linalg.norm(np.diff(contour, axis=0, append=contour[:1]), axis=1)
        s_local = np.concatenate([[0.0], np.cumsum(seg)[:-1]])
        s_global = s_local + s_offset

        new_cid += 1
        xys.append(contour)
        ss.append(s_global)
        cids.append(np.full(len(contour), new_cid, dtype=np.int32))
        s_lo = s_offset
        s_hi = s_offset + s_local[-1] + seg[-1]
        ranges.append((new_cid, s_lo, s_hi, int(cc_mask.sum())))

        # Transfer s from contour to every oEGL pixel in this CC via nearest
        # contour point (Euclidean inside this CC only).
        cc_ys, cc_xs = np.where(cc_mask.astype(bool))
        kdt_local = cKDTree(contour[:, ::-1])    # (y, x) order
        _, nn = kdt_local.query(np.stack([cc_ys, cc_xs], axis=1), k=1)
        s_field[cc_ys, cc_xs] = s_global[nn].astype(np.float32)

        s_offset = s_hi + gap_between_components

    xy_all = np.concatenate(xys, axis=0)
    s_all = np.concatenate(ss, axis=0)
    cid_all = np.concatenate(cids, axis=0)
    return s_field, xy_all, s_all, cid_all, s_offset, ranges


def label_rgb(label_img: np.ndarray) -> np.ndarray:
    rgb = np.zeros((*label_img.shape, 3), dtype=np.float32)
    for v, c in COLORS.items():
        rgb[label_img == v] = c
    return rgb


# --------------------------------------------------------------------------- #
# Laplace depth field
# --------------------------------------------------------------------------- #

def _shift_neighbours(mask, dy, dx):
    """Return mask shifted by (dy, dx) with zero padding."""
    out = np.zeros_like(mask)
    H, W = mask.shape
    y0_src = max(0, -dy); y1_src = H - max(0, dy)
    x0_src = max(0, -dx); x1_src = W - max(0, dx)
    y0_dst = max(0, dy);  y1_dst = H - max(0, -dy)
    x0_dst = max(0, dx);  x1_dst = W - max(0, -dx)
    out[y0_dst:y1_dst, x0_dst:x1_dst] = mask[y0_src:y1_src, x0_src:x1_src]
    return out


def adjacent_mask(target: np.ndarray, source: np.ndarray) -> np.ndarray:
    """target pixels with any 4-neighbour in source."""
    nb = (_shift_neighbours(source, 1, 0) | _shift_neighbours(source, -1, 0)
          | _shift_neighbours(source, 0, 1) | _shift_neighbours(source, 0, -1))
    return target & nb


def solve_laplace_depth(lbl: np.ndarray):
    """Returns psi (NaN outside domain), domain mask, outer_bc, inner_bc.

    The true pia (tissue-adjacent-to-bg) is unreliable — segmentation often
    drops the meningeal edge. We use the oEGL (label 3), but only its
    **outer pial surface** as the ψ=0 boundary — the oEGL pixels that
    have a bg or label-1 (unclassified-cereb) 4-neighbour. ψ then varies
    smoothly THROUGH the oEGL band, so the streamline path length from
    any oEGL pixel back to ψ=0 gives that pixel's local oEGL depth, and
    the 95th-percentile depth across all oEGL pixels at a given arclength
    gives the local oEGL thickness.

      domain    = ALL tissue
      outer_bc  = oEGL ∩ adjacent(bg ∪ label-1)     → ψ = 0
      inner_bc  = ALL DWL pixels (label 6)          → ψ = 1
    """
    H, W = lbl.shape
    tissue = np.isin(lbl, TISSUE_LABELS)
    oegl = lbl == LBL["oEGL"]
    dwl = lbl == LBL["DWL"]

    # Drop tissue connected components (4-conn matches the stencil) that
    # contain no oEGL — otherwise their sub-block of the Laplace matrix
    # has no Dirichlet anchor and the whole system is singular.
    # An anomaly (label 9) is a tear/bleach hole — for BC purposes it is an
    # exposed surface just like bg, so we include it as "outside cortex".
    outside_cortex = ((lbl == LBL["BG"])
                      | (lbl == LBL["CEREB"])
                      | (lbl == LBL["ANOMALY"]))
    provisional_outer = adjacent_mask(oegl, outside_cortex) & tissue
    provisional_inner = dwl & tissue

    # Keep only sufficiently-large tissue connected components that contain
    # at least one BC pixel of EITHER kind. Tiny stray fragments (≤ ~200 px)
    # produce degenerate matrix rows even when nominally well-posed; their
    # contribution to the cortex unfolding is negligible.
    MIN_CC_SIZE = 200
    n_cc, ccl = cv2.connectedComponents(
        tissue.astype(np.uint8), connectivity=4
    )
    has_bc = provisional_outer | provisional_inner
    keep_ids = []
    for c in np.unique(ccl):
        if c == 0:
            continue
        m = ccl == c
        if m.sum() >= MIN_CC_SIZE and has_bc[m].any():
            keep_ids.append(int(c))
    valid_tissue = np.isin(ccl, keep_ids) if keep_ids else np.zeros_like(tissue)

    outer_bc = provisional_outer & valid_tissue
    inner_bc = provisional_inner & valid_tissue

    domain = valid_tissue.copy()

    n = int(domain.sum())
    idx = -np.ones((H, W), dtype=np.int64)
    idx[domain] = np.arange(n)

    is_bc  = outer_bc | inner_bc
    is_int = domain & ~is_bc

    rows_list, cols_list, vals_list = [], [], []
    b = np.zeros(n, dtype=np.float64)

    # Dirichlet rows
    bc_idx = idx[is_bc]
    rows_list.append(bc_idx)
    cols_list.append(bc_idx)
    vals_list.append(np.ones(len(bc_idx)))
    b[idx[outer_bc]] = 0.0
    b[idx[inner_bc]] = 1.0

    # Interior rows. First pass: discover which interior pixels actually
    # have ≥1 in-domain 4-neighbour. Anything with K_valid=0 is an isolated
    # outlier (image-edge corner case) and would make the matrix singular —
    # pin those to ψ=0 as Dirichlet instead.
    iy, ix = np.where(is_int)
    int_idx = idx[iy, ix]

    K_valid = np.zeros(len(int_idx), dtype=np.float64)
    neighbour_packets = []
    for dy, dx in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
        ny, nx = iy + dy, ix + dx
        ok = (ny >= 0) & (ny < H) & (nx >= 0) & (nx < W)
        in_dom = np.zeros(len(ok), dtype=bool)
        in_dom[ok] = domain[ny[ok], nx[ok]]
        K_valid += in_dom.astype(np.float64)
        neighbour_packets.append((in_dom, ny, nx))

    isolated = K_valid == 0
    if isolated.any():
        bad_idx = int_idx[isolated]
        rows_list.append(bad_idx)
        cols_list.append(bad_idx)
        vals_list.append(np.ones(len(bad_idx)))
        # b[bad_idx] left at 0 — fine.

    keep = ~isolated
    int_idx_k = int_idx[keep]

    for in_dom, ny, nx in neighbour_packets:
        # Only emit neighbour edges for kept rows.
        sel = in_dom & keep
        rows_list.append(int_idx[sel])
        cols_list.append(idx[ny[sel], nx[sel]])
        vals_list.append(np.ones(int(sel.sum())))

    rows_list.append(int_idx_k)
    cols_list.append(int_idx_k)
    vals_list.append(-K_valid[keep])

    rows = np.concatenate(rows_list)
    cols = np.concatenate(cols_list)
    vals = np.concatenate(vals_list)
    A = csr_matrix((vals, (rows, cols)), shape=(n, n))
    # Tiny diagonal regularisation — guards against rank deficiency from
    # degenerate sub-topologies without measurably perturbing the solution.
    from scipy.sparse import eye as _sp_eye
    A = A + _sp_eye(n) * 1e-9

    psi_flat = spsolve(A, b)
    psi = np.full((H, W), np.nan, dtype=np.float64)
    psi[domain] = psi_flat
    if not np.isfinite(psi[domain]).any():
        raise RuntimeError(
            "Laplace solve returned all-NaN — system was singular even after "
            "BC-anchored CC filtering. Likely a malformed segmentation."
        )
    return psi, domain, outer_bc, inner_bc


# --------------------------------------------------------------------------- #
# Streamline trace for arclength
# --------------------------------------------------------------------------- #

def _bilinear(field, y, x):
    H, W = field.shape
    y0 = np.clip(np.floor(y).astype(int), 0, H - 2)
    x0 = np.clip(np.floor(x).astype(int), 0, W - 2)
    fy = y - y0
    fx = x - x0
    return (field[y0, x0] * (1 - fy) * (1 - fx)
            + field[y0 + 1, x0] * fy * (1 - fx)
            + field[y0, x0 + 1] * (1 - fy) * fx
            + field[y0 + 1, x0 + 1] * fy * fx)


def trace_to_pia(psi, cortex_mask, step=0.4, max_steps=1500,
                 pia_thresh=0.005):
    """For each cortex pixel, integrate −∇ψ back to ψ ≤ pia_thresh.

    Returns:
      landing_y, landing_x  (n_cortex,) final positions
      reached               (n_cortex,) bool — did it actually land?
      path_len              (n_cortex,) cumulative streamline length (downsampled
                            pixels) from the pixel's start to its landing.
                            For an oEGL pixel near the pial surface this is its
                            local oEGL thickness.
    """
    psi_f = np.where(np.isnan(psi), 1.0, psi).astype(np.float64)
    gy, gx = np.gradient(psi_f)

    cy, cx = np.where(cortex_mask)
    n = len(cy)
    y = cy.astype(np.float64).copy()
    x = cx.astype(np.float64).copy()
    reached = np.zeros(n, dtype=bool)
    path_len = np.zeros(n, dtype=np.float64)

    for _ in range(max_steps):
        active_idx = np.where(~reached)[0]
        if len(active_idx) == 0:
            break

        ya = y[active_idx]
        xa = x[active_idx]
        psi_v = _bilinear(psi_f, ya, xa)
        landed_mask = psi_v <= pia_thresh
        if landed_mask.any():
            reached[active_idx[landed_mask]] = True

        # Step only those still active and NOT just landed.
        move = active_idx[~landed_mask]
        if len(move) == 0:
            continue
        ym = y[move]
        xm = x[move]

        gy_v = _bilinear(gy, ym, xm)
        gx_v = _bilinear(gx, ym, xm)
        mag = np.sqrt(gy_v ** 2 + gx_v ** 2) + 1e-9
        dy_step = -step * gy_v / mag
        dx_step = -step * gx_v / mag

        new_y = np.clip(ym + dy_step, 0, psi.shape[0] - 1)
        new_x = np.clip(xm + dx_step, 0, psi.shape[1] - 1)

        actual = np.sqrt((new_y - ym) ** 2 + (new_x - xm) ** 2)
        path_len[move] += actual

        y[move] = new_y
        x[move] = new_x

    return y, x, reached, path_len


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #

REPO_ROOT = Path(__file__).resolve().parents[2]
SANDBOX_20X = REPO_ROOT / "python" / "matlab_anomaly" / "sandbox"
SANDBOX_10X = REPO_ROOT / "python" / "matlab_anomaly_10x" / "sandbox"
OUT_ROOT = REPO_ROOT / "python" / "depth_model" / "out"

DEFAULT_SRC = SANDBOX_20X / "2018_05_22_s1_2_p27-0021_iter1deep_native_labels.tif"
OUTPUT_ROOT = OUT_ROOT / "ribbon_unfold"


def slide_id_from_path(src: Path) -> str:
    stem = src.stem
    for suffix in ("_native_labels", "_labelled"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    return stem


def main(src: Path = DEFAULT_SRC, ds: int = 4):
    slide_id = slide_id_from_path(src)
    out_dir = OUTPUT_ROOT / slide_id
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"slide: {slide_id}")
    print(f"output dir: {out_dir}")

    print(f"loading {src.name}")
    lbl_full = load_labels(src)
    DS = ds
    lbl = downsample(lbl_full, DS)
    print(f"  downsampled {DS}x to {lbl.shape}")

    # Show label histogram so we know what we're working with.
    vals, cnts = np.unique(lbl, return_counts=True)
    label_names = {0: "bg", 1: "cereb", 2: "iEGL", 3: "oEGL", 4: "IGL",
                   5: "ML", 6: "DWL", 7: "PCL", 8: "DCN", 9: "anomaly"}
    print("label histogram (downsampled):")
    for v, c in zip(vals.tolist(), cnts.tolist()):
        print(f"  {v} ({label_names.get(v, '?')}): {c} px")

    tissue = tissue_mask(lbl)
    print(f"tissue pixels (excludes bg + anomaly): {tissue.sum()}")
    print(f"anomaly pixels (excluded from domain): "
          f"{int((lbl == LBL['ANOMALY']).sum())}")

    pia_xy_outer, _ = pial_contour_xys(tissue, smooth_win=21)
    print(f"outer hull contour (for overlay only): {pia_xy_outer.shape[0]} pts")

    print("solving Laplace depth field ...")
    psi, domain, outer_bc, inner_bc = solve_laplace_depth(lbl)
    print(f"  domain n={domain.sum()}, outer_bc (oEGL) n={outer_bc.sum()}, "
          f"inner_bc (DWL) n={inner_bc.sum()}")
    print(f"  psi range: {np.nanmin(psi):.3f} .. {np.nanmax(psi):.3f}")

    s_field, pia_xy, pia_s, pia_cid, total_s, comp_ranges = oegl_arclength(
        outer_bc, min_cc_size=200, smooth_win=9, gap_between_components=12.0
    )
    print(f"oEGL components used: {len(comp_ranges)}  "
          f"(total arclength {total_s:.1f})")
    for cid, lo, hi, sz in comp_ranges[:10]:
        print(f"  comp {cid}: s=[{lo:.1f}, {hi:.1f}]  len={hi-lo:.1f}  oEGL px={sz}")

    # Include the oEGL itself in the unrolled cortex so both EGL bands are
    # visible (oEGL sits at ψ≈0, iEGL just below). The streamline tracer
    # already terminates immediately for ψ≤thresh pixels, so oEGL pixels
    # simply land on themselves and get s from the local s_field.
    cortex_mask = np.isin(lbl, CORTEX_LABELS)
    cy, cx = np.where(cortex_mask)
    psi_at_cortex = psi[cy, cx]

    print("psi per layer:")
    for v in CORTEX_LABELS:
        sel = lbl[cy, cx] == v
        if sel.sum() == 0:
            continue
        p = psi_at_cortex[sel]
        n_nan = int(np.isnan(p).sum())
        print(f"  label {v}: n={sel.sum()} (nan={n_nan}), "
              f"mean={np.nanmean(p):.3f}, p10={np.nanpercentile(p,10):.3f}, "
              f"p50={np.nanpercentile(p,50):.3f}, p90={np.nanpercentile(p,90):.3f}")

    print("tracing streamlines to pial surface (ψ ≤ 0.005) ...")
    land_y, land_x, reached, path_len = trace_to_pia(
        psi, cortex_mask, step=0.4, max_steps=1500, pia_thresh=0.005
    )
    print(f"  reached pia: {reached.sum()}/{len(reached)} "
          f"({100*reached.sum()/len(reached):.1f}%)")
    # Per-layer streamline length (≈ physical depth from pial surface, in
    # downsampled pixels). Multiply by DS for original-resolution px.
    print("streamline path length from pial surface (downsampled px):")
    layer_per_pixel = lbl[cy, cx]
    for v in CORTEX_LABELS:
        m = (layer_per_pixel == v) & reached
        if m.sum() == 0:
            continue
        pl = path_len[m]
        print(f"  label {v}: n={m.sum()}, "
              f"mean={pl.mean():.2f}, p10={np.percentile(pl,10):.2f}, "
              f"p50={np.percentile(pl,50):.2f}, p90={np.percentile(pl,90):.2f}")

    # Build KDTree on every oEGL pixel that has an s assigned, then look up
    # the s of the nearest oEGL pixel to each streamline's landing point.
    oegl_ys, oegl_xs = np.where(s_field >= 0)
    oegl_s_vals = s_field[oegl_ys, oegl_xs]
    # Component id at each oEGL pixel — use nearest contour point cid.
    kdt_contour = cKDTree(pia_xy[:, ::-1])
    _, nn_c = kdt_contour.query(np.stack([oegl_ys, oegl_xs], axis=1), k=1)
    oegl_cid_vals = pia_cid[nn_c]

    kdt_oegl = cKDTree(np.stack([oegl_ys, oegl_xs], axis=1).astype(np.float32))
    land_pts = np.stack([land_y, land_x], axis=1)
    _, nn = kdt_oegl.query(land_pts, k=1)
    s_per_pixel = oegl_s_vals[nn].astype(np.float64)
    comp_per_pixel = oegl_cid_vals[nn]
    s_per_pixel[~reached] = np.nan
    comp_per_pixel[~reached] = -1

    d_per_pixel = psi_at_cortex

    valid = reached & np.isfinite(s_per_pixel) & np.isfinite(d_per_pixel)
    print(f"valid (s,d) pixels: {valid.sum()}/{len(valid)}")

    # --- oEGL thickness profile vs arclength s -------------------------
    # For every oEGL pixel, path_len is the streamline distance from the
    # pixel back to the outer pial surface — i.e. its depth INSIDE the
    # oEGL. The deepest oEGL pixel at a given s is the inner edge of the
    # oEGL band there, so its path_len = local oEGL thickness.
    is_oegl = (layer_per_pixel == LBL["oEGL"]) & reached
    n_thick_bins = 600
    edges = np.linspace(0.0, total_s, n_thick_bins + 1)
    bin_idx = np.clip(
        np.searchsorted(edges, s_per_pixel[is_oegl], side="right") - 1,
        0, n_thick_bins - 1,
    )
    oegl_pl = path_len[is_oegl]
    thickness_ds_px = np.full(n_thick_bins, np.nan, dtype=np.float64)
    thickness_med_px = np.full(n_thick_bins, np.nan, dtype=np.float64)
    counts_per_bin = np.zeros(n_thick_bins, dtype=np.int64)
    for b in range(n_thick_bins):
        sel = bin_idx == b
        if sel.sum() >= 4:
            counts_per_bin[b] = sel.sum()
            thickness_ds_px[b] = np.percentile(oegl_pl[sel], 95)
            thickness_med_px[b] = np.median(oegl_pl[sel])
    bin_centers = 0.5 * (edges[:-1] + edges[1:])

    # Mask out bin centres that fall inside the gaps between oEGL components.
    inside_component = np.zeros(n_thick_bins, dtype=bool)
    for cid, lo, hi, sz in comp_ranges:
        inside_component |= (bin_centers >= lo) & (bin_centers <= hi)
    thickness_ds_px[~inside_component] = np.nan
    thickness_med_px[~inside_component] = np.nan

    valid_thk = np.isfinite(thickness_ds_px)
    if valid_thk.any():
        # Multiply by downsample factor for original-resolution pixels.
        thk95_orig = thickness_ds_px * DS
        thkmed_orig = thickness_med_px * DS
        print(f"oEGL thickness (orig-res px): "
              f"median across s = {np.nanmedian(thk95_orig):.1f}, "
              f"p10–p90 across s = "
              f"{np.nanpercentile(thk95_orig, 10):.1f}–"
              f"{np.nanpercentile(thk95_orig, 90):.1f}")

    # --- rasterise unrolled ribbon ---
    s_use = s_per_pixel[valid]
    d_use = d_per_pixel[valid]
    lay_use = layer_per_pixel[valid]

    s_bins, d_bins = 2400, 220
    s_max = total_s
    s_i = np.clip((s_use / s_max * (s_bins - 1)).astype(int), 0, s_bins - 1)
    d_i = np.clip((d_use * (d_bins - 1)).astype(int), 0, d_bins - 1)
    flat = d_i * s_bins + s_i
    order = np.argsort(flat, kind="stable")
    fi = flat[order]; lo = lay_use[order]
    unrolled = np.zeros((d_bins, s_bins), dtype=np.int32)
    if len(fi) > 0:
        boundaries = np.concatenate(
            [[0], np.where(np.diff(fi) != 0)[0] + 1, [len(fi)]]
        )
        for a, b_ in zip(boundaries[:-1], boundaries[1:]):
            vals, cnts = np.unique(lo[a:b_], return_counts=True)
            if len(vals):
                unrolled.flat[fi[a]] = vals[cnts.argmax()]
    else:
        print("  WARNING: no valid (s, ψ) pixels — unrolled raster left empty")

    # Fill rasterisation-undersampling gaps with the nearest non-zero label,
    # but only out to a small distance — that way the inter-component gaps
    # and the "no cortex past this ψ" gaps stay black.
    FILL_DIST = 2.0
    zero_mask = unrolled == 0
    if zero_mask.any() and (~zero_mask).any():
        dist, (yi, xi) = ndi.distance_transform_edt(zero_mask, return_indices=True)
        fill_here = zero_mask & (dist <= FILL_DIST)
        unrolled_filled = unrolled.copy()
        unrolled_filled[fill_here] = unrolled[yi[fill_here], xi[fill_here]]
        n_filled = int(fill_here.sum())
        n_real_gap = int(zero_mask.sum() - n_filled)
        print(f"  rasterisation gaps filled: {n_filled}  "
              f"(remaining real gaps: {n_real_gap})")
        unrolled = unrolled_filled

    # --- figure ---
    fig = plt.figure(figsize=(14, 17), constrained_layout=True)
    gs = fig.add_gridspec(5, 1, height_ratios=[5, 4, 4, 2, 0.45])
    axes = [fig.add_subplot(gs[i]) for i in range(3)]
    ax_thick = fig.add_subplot(gs[3], sharex=axes[2])
    ax_legend = fig.add_subplot(gs[4])

    ax = axes[0]
    ax.imshow(label_rgb(lbl))
    ax.plot(pia_xy_outer[:, 0], pia_xy_outer[:, 1], "w-", lw=0.3, alpha=0.4)
    palette = plt.get_cmap("tab20")
    for cid, lo, hi, sz in comp_ranges:
        m = pia_cid == cid
        col = palette((cid - 1) % 20)
        ax.scatter(pia_xy[m, 0], pia_xy[m, 1], c=[col], s=0.6)
    ax.set_title(f"native space — oEGL outer contour, "
                 f"{len(comp_ranges)} component(s)")
    ax.set_axis_off()

    ax = axes[1]
    im = ax.imshow(psi, cmap="viridis", vmin=0, vmax=1)
    ax.contour(np.where(np.isnan(psi), 1.0, psi),
               levels=[0.2, 0.4, 0.6, 0.8], colors="white",
               linewidths=0.4, alpha=0.7)
    ax.set_title("Laplace depth ψ (0 = oEGL, 1 = DWL)")
    ax.set_axis_off()
    fig.colorbar(im, ax=ax, fraction=0.025)

    ax = axes[2]
    ax.imshow(label_rgb(unrolled),
              aspect="auto", extent=[0, s_max, 1.0, 0.0],
              interpolation="nearest")
    # Draw a vertical line at the boundary of every pia component so the
    # outer hull and the individual sulci are visually separated.
    for cid, lo, hi, sz in comp_ranges:
        ax.axvline(hi, color="white", lw=0.6, alpha=0.6)
        ax.text(0.5 * (lo + hi), 1.04, f"c{cid}\n{sz}px",
                ha="center", va="top", fontsize=6, color="black",
                transform=ax.get_xaxis_transform())
    ax.set_xlabel("")  # the thickness pane below shares this x-axis
    ax.set_ylabel("Laplace depth ψ (0 = oEGL, 1 = DWL)")
    ax.set_title("unrolled ribbon — depth from oEGL × arclength along oEGL")

    # --- oEGL thickness profile (shares x-axis with the unrolled ribbon) ---
    ax_thick.plot(bin_centers, thickness_ds_px * DS,
                  color="#1565c0", lw=1.2, label="p95 (oEGL thickness)")
    ax_thick.plot(bin_centers, thickness_med_px * DS,
                  color="#43a047", lw=0.8, alpha=0.7, label="median")
    for cid, lo, hi, sz in comp_ranges:
        ax_thick.axvline(hi, color="gray", lw=0.5, alpha=0.6)
    ax_thick.set_xlim(0, s_max)
    ax_thick.set_ylim(bottom=0)
    ax_thick.set_xlabel("concatenated arclength s along oEGL components (px)")
    ax_thick.set_ylabel("oEGL thickness\n(orig-res px)")
    ax_thick.set_title("oEGL local thickness — streamline length from inner "
                       "oEGL pixels back to pial surface")
    ax_thick.legend(loc="upper right", fontsize=8, framealpha=0.9)
    ax_thick.grid(alpha=0.3)

    # Colour-key strip (canonical sandbox convention: 4=IGL, 5=ML).
    legend_items = [
        (3, "oEGL (ψ=0 BC)"),
        (2, "iEGL"),
        (5, "ML"),
        (7, "PCL"),
        (4, "IGL"),
        (6, "DWL (ψ=1 BC)"),
        (8, "DCN"),
        (1, "cereb / unclassified"),
        (9, "anomaly (excluded)"),
        (0, "background"),
    ]
    ax_legend.set_xlim(0, len(legend_items))
    ax_legend.set_ylim(0, 1)
    for i, (lbl_val, name) in enumerate(legend_items):
        ax_legend.add_patch(
            plt.Rectangle((i + 0.05, 0.2), 0.18, 0.6,
                          facecolor=COLORS[lbl_val], edgecolor="black", lw=0.4)
        )
        ax_legend.text(i + 0.28, 0.5, f"{lbl_val} — {name}",
                       va="center", ha="left", fontsize=7.5)
    ax_legend.set_axis_off()
    ax_legend.set_title("colour key (panes 1 + 3)", fontsize=9, loc="left")

    fig.suptitle(f"ribbon unfolding — {slide_id}", fontsize=11, fontweight="bold")
    out_png = out_dir / f"{slide_id}_ribbon_unfold.png"
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    print(f"wrote {out_png}")

    np.savez_compressed(
        out_dir / f"{slide_id}_ribbon_coords.npz",
        cy=cy.astype(np.int32),
        cx=cx.astype(np.int32),
        layer=layer_per_pixel.astype(np.uint8),
        psi=d_per_pixel.astype(np.float32),
        s=s_per_pixel.astype(np.float32),
        component=comp_per_pixel.astype(np.int32),
        reached=reached,
        path_len_ds_px=path_len.astype(np.float32),
        pia_xy=pia_xy.astype(np.float32),
        pia_s=pia_s.astype(np.float32),
        pia_cid=pia_cid.astype(np.int32),
        comp_ranges=np.array(comp_ranges, dtype=np.float32),
        s_max=float(s_max),
        downsample=DS,
        oegl_thickness_p95_ds_px=thickness_ds_px.astype(np.float32),
        oegl_thickness_med_ds_px=thickness_med_px.astype(np.float32),
        oegl_thickness_bin_centers=bin_centers.astype(np.float32),
    )
    print(f"wrote {out_dir / (slide_id + '_ribbon_coords.npz')}")


if __name__ == "__main__":
    import sys
    src_arg = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_SRC
    ds_arg = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    main(src_arg, ds_arg)
