"""
Anatomical relabelling rules for the current pipeline's label maps.

R1 and R2 use the same cue: whether IGL lies on the straight line from a
pixel to its nearest EGL. In a folium, white matter and DCN sit under the
IGL, while ML sits between EGL and IGL.

  R1  DWL with no IGL between it and the EGL, within 125 um of the EGL
      -> ML. On the gold slide (s1_2) these "white matter" chunks are
      mostly ML.
  R2  DCN with no IGL between it and the EGL, within 100 um of the EGL
      -> nearest cortex label (EGL, ML, PCL or IGL). DCN at the surface is
      mostly gold EGL.
  R3  DWL pieces cut off from the main white-matter tree -> the cortex
      label that surrounds them. Pieces that reconnect to the tree through
      anomaly (tear) pixels are left alone.
  R4  Tissue enclosed by IGL, after closing gaps of up to 25 um in the
      IGL + white-matter ring -> DWL. Recovers folium white-matter cores the
      pipeline missed because their IGL ring was broken (its own imfill
      step then never fills them). Only unassigned, EGL, ML and PCL pixels
      are changed.

Decisions are made at the downsampled resolution and applied at native
resolution only to pixels that carry the rule's source label, so native
boundaries elsewhere are untouched.

Outputs (under python/depth_model/out/anatomy_fixes/):
  <slide>_fixed_labels.tif    native-resolution labels with the fixes
  summary.csv                 per-slide changes, continuity, gold scores
"""

import sys
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image
from scipy import ndimage as ndi

sys.path.insert(0, str(Path(__file__).parent))
from ribbon_unfold_laplace import (                                # noqa: E402
    LBL, SANDBOX_20X, SANDBOX_10X, OUT_ROOT, load_labels,
)
from tune_relabel import (                                         # noqa: E402
    SLIDES_20X, SLIDES_10X, GOLD, continuity, gold_scores,
)

Image.MAX_IMAGE_PIXELS = None

UM_PER_PX = {"20x": 0.5119049, "10x": 1.0239}
CORTEX = [LBL["iEGL"], LBL["oEGL"], LBL["IGL"], LBL["ML"], LBL["PCL"]]
EGL = [LBL["iEGL"], LBL["oEGL"]]
OUT_DIR = OUT_ROOT / "anatomy_fixes"


def igl_fraction_on_path(lab, py, px, nearest_y, nearest_x, n_samples=40):
    t = np.linspace(0.0, 1.0, n_samples)[None, :]
    qy, qx = nearest_y[py, px], nearest_x[py, px]
    sy = np.rint(py[:, None] + t * (qy - py)[:, None]).astype(np.intp)
    sx = np.rint(px[:, None] + t * (qx - px)[:, None]).astype(np.intp)
    return (lab[sy, sx] == LBL["IGL"]).mean(axis=1)


def drop_bridges(dwl, cand, min_piece_frac=0.02, min_piece_px=200):
    """Remove candidate chunks whose conversion would split the DWL tree.

    Chunks are tested largest first against the DWL as already trimmed, and
    kept only if no remaining piece of significant size gets cut off from
    the DWL component the chunk belonged to.
    """
    st = np.ones((3, 3), dtype=bool)
    chunks, n = ndi.label(cand, structure=st)
    if n == 0:
        return cand
    order = np.argsort(-np.bincount(chunks.ravel())[1:]) + 1
    cur = dwl.copy()
    keep = np.zeros(n + 1, dtype=bool)
    for k in order:
        ck = chunks == k
        parent_lab, _ = ndi.label(cur, structure=st)
        parent_ids = np.unique(parent_lab[ck])
        parent_ids = parent_ids[parent_ids > 0]
        parent = np.isin(parent_lab, parent_ids)
        rest = parent & ~ck
        pieces, m = ndi.label(rest, structure=st)
        if m > 1:
            sizes = np.sort(np.bincount(pieces.ravel())[1:])[::-1]
            if sizes[1] >= max(min_piece_px, min_piece_frac * parent.sum()):
                continue
        keep[k] = True
        cur &= ~ck
    return keep[chunks]


def decide(lab, um_per_px, dwl_cap_um=125.0, dcn_cap_um=100.0, close_um=25.0):
    """Return the relabelled image at the resolution of `lab`."""
    egl = np.isin(lab, EGL)
    dist, (ny, nx) = ndi.distance_transform_edt(~egl, return_indices=True)
    depth_um = dist * um_per_px
    out = lab.copy()

    py, px = np.nonzero(lab == LBL["DWL"])
    f = igl_fraction_on_path(lab, py, px, ny, nx)
    r1 = (f == 0) & (depth_um[py, px] < dwl_cap_um)
    cand = np.zeros(lab.shape, dtype=bool)
    cand[py[r1], px[r1]] = True
    cand = drop_bridges(lab == LBL["DWL"], cand)
    r1 = cand[py, px]
    out[py[r1], px[r1]] = LBL["ML"]

    py, px = np.nonzero(lab == LBL["DCN"])
    f = igl_fraction_on_path(lab, py, px, ny, nx)
    r2 = (f == 0) & (depth_um[py, px] < dcn_cap_um)
    if r2.any():
        _, (cy, cx) = ndi.distance_transform_edt(~np.isin(lab, CORTEX),
                                                 return_indices=True)
        ty, tx = py[r2], px[r2]
        out[ty, tx] = lab[cy[ty, tx], cx[ty, tx]]

    n3 = relabel_dwl_islands(out)
    n4 = fill_igl_enclosed(out, um_per_px, close_um)
    return out, int(r1.sum()), int(r2.sum()), n3, n4


def fill_igl_enclosed(lab, um_per_px, close_um=25.0):
    """R4: tissue enclosed by IGL (gaps up to close_um closed) -> DWL.

    Modifies `lab` in place; returns the number of pixels changed.
    """
    r = max(1, int(round(close_um / um_per_px)))
    yy, xx = np.ogrid[-r:r + 1, -r:r + 1]
    disk = xx * xx + yy * yy <= r * r
    ring = ndi.binary_closing(np.isin(lab, [LBL["IGL"], LBL["DWL"]]), structure=disk)
    inside = ndi.binary_fill_holes(ring) & ~ring
    target = inside & np.isin(lab, [LBL["CEREB"], LBL["iEGL"], LBL["oEGL"],
                                    LBL["ML"], LBL["PCL"]])
    lab[target] = LBL["DWL"]
    return int(target.sum())


def relabel_dwl_islands(lab, ring_px=2):
    """R3: DWL pieces cut off from the main tree take their surrounding label.

    Pieces that reconnect to the tree through anomaly pixels are left alone:
    those are white matter separated from the tree by a tear.
    Modifies `lab` in place; returns the number of pixels changed.
    """
    st = np.ones((3, 3), dtype=bool)
    dwl = lab == LBL["DWL"]
    cc, n = ndi.label(dwl, structure=st)
    if n < 2:
        return 0
    main_id = np.bincount(cc.ravel())[1:].argmax() + 1
    via, _ = ndi.label(dwl | (lab == LBL["ANOMALY"]), structure=st)
    tree_ids = np.unique(via[cc == main_id])
    changed = 0
    for k, sl in enumerate(ndi.find_objects(cc), start=1):
        if k == main_id:
            continue
        pad = ring_px + 1
        ys = slice(max(0, sl[0].start - pad), sl[0].stop + pad)
        xs = slice(max(0, sl[1].start - pad), sl[1].stop + pad)
        comp = cc[ys, xs] == k
        if np.isin(via[ys, xs][comp], tree_ids).any():
            continue
        ring = ndi.binary_dilation(comp, st, iterations=ring_px) & ~comp
        neigh = lab[ys, xs][ring]
        neigh = neigh[np.isin(neigh, CORTEX)]
        if neigh.size == 0:
            continue
        target = np.bincount(neigh).argmax()
        lab[ys, xs][comp] = target
        changed += int(comp.sum())
    return changed


def apply_native(native, ds, decided_ds):
    """Stamp ds-level changes onto native pixels that had the source label."""
    H, W = native.shape
    ys = np.minimum(np.arange(H) // ds, decided_ds.shape[0] - 1)
    xs = np.minimum(np.arange(W) // ds, decided_ds.shape[1] - 1)
    orig_up = native[::ds, ::ds][np.ix_(ys, xs)]
    new_up = decided_ds[np.ix_(ys, xs)]
    change = (new_up != orig_up) & (native == orig_up)
    out = native.copy()
    out[change] = new_up[change]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ds", type=int, default=4)
    ap.add_argument("--dwl_cap_um", type=float, default=125.0)
    ap.add_argument("--dcn_cap_um", type=float, default=100.0)
    ap.add_argument("--igl_close_um", type=float, default=25.0)
    ap.add_argument("--no_write", action="store_true",
                    help="score only; do not write fixed label TIFFs")
    args = ap.parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    rows = []
    jobs = [(s, "20x", SANDBOX_20X) for s in SLIDES_20X] + \
           [(s, "10x", SANDBOX_10X) for s in SLIDES_10X]
    for sid, mag, sandbox in jobs:
        native = load_labels(sandbox / f"{sid}_native_labels.tif").astype(np.uint8)
        ds = args.ds
        decided, n1, n2, n3, n4 = decide(native[::ds, ::ds], UM_PER_PX[mag] * ds,
                                         args.dwl_cap_um, args.dcn_cap_um,
                                         args.igl_close_um)
        fixed = apply_native(native, ds, decided)
        changed = fixed != native
        cortex_tissue = np.isin(native, CORTEX + [LBL["DWL"], LBL["DCN"]])
        row = dict(sid=sid, mag=mag,
                   changed_pct=100.0 * changed.sum() / cortex_tissue.sum(),
                   dwl_to_ml_px=int(((native == LBL["DWL"]) & changed).sum()),
                   dcn_fixed_px=int(((native == LBL["DCN"]) & changed).sum()))
        before_c = continuity(native[::ds, ::ds], ds)
        after_c = continuity(fixed[::ds, ::ds], ds)
        for name in ("EGL", "ML", "IGL", "DWL"):
            row[f"{name}_cont_before"] = before_c[name][0]
            row[f"{name}_cont_after"] = after_c[name][0]
        if sid in GOLD:
            g = np.array(Image.open(GOLD[sid])).astype(np.int16)
            gs = g.copy()
            gs[g == 4] = 5
            gs[g == 5] = 4
            acc0, iou0 = gold_scores(native, gs)
            acc1, iou1 = gold_scores(fixed, gs)
            row.update(gold_acc_before=acc0, gold_acc_after=acc1,
                       ML_iou_before=iou0[5], ML_iou_after=iou1[5],
                       DWL_iou_before=iou0[6], DWL_iou_after=iou1[6],
                       EGL_iou_before=(iou0[2] + iou0[3]) / 2,
                       EGL_iou_after=(iou1[2] + iou1[3]) / 2)
        if not args.no_write:
            Image.fromarray(fixed).save(OUT_DIR / f"{sid}_fixed_labels.tif")
        rows.append(row)
        print(f"{sid:28s} {mag}  changed {row['changed_pct']:5.2f}%  "
              f"DWL changed {row['dwl_to_ml_px']:8d}  DCN {row['dcn_fixed_px']:7d}  "
              f"ML cont {row['ML_cont_before']:5.1f} -> {row['ML_cont_after']:5.1f}",
              flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT_DIR / "summary.csv", index=False)
    g = df.dropna(subset=["gold_acc_before"])
    for _, r in g.iterrows():
        print(f"\ngold ({r['sid']}, native resolution):")
        print(f"  accuracy {r['gold_acc_before']:.4f} -> {r['gold_acc_after']:.4f} "
              f"({100 * (r['gold_acc_after'] - r['gold_acc_before']):+.2f} points)")
        for nm in ("ML", "DWL", "EGL"):
            print(f"  {nm} IoU {r[f'{nm}_iou_before']:.3f} -> {r[f'{nm}_iou_after']:.3f}")
    for mag in ("20x", "10x"):
        d = df[df["mag"] == mag]
        print(f"\n{mag} mean: changed {d['changed_pct'].mean():.2f}% of cortex+WM; "
              + "  ".join(f"{n} cont {d[f'{n}_cont_before'].mean():.1f}->{d[f'{n}_cont_after'].mean():.1f}"
                          for n in ("EGL", "ML", "IGL", "DWL")))
    print(f"\nwrote {OUT_DIR / 'summary.csv'}")


if __name__ == "__main__":
    main()
