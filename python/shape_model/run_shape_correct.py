"""Run shape-based correction on every good slide (20x + 10x) and write review panels.

    python run_shape_correct.py            # all slides in dino/slide_manifest.py
    python run_shape_correct.py s2_5 s4_A  # substring filter

Outputs (python/shape_model/out/, gitignored):
    <slide>_shape_labels.tif   corrected labels, native resolution, same index as MATLAB
    <slide>_overview.png       RGB | MATLAB | corrected | changed pixels
    <slide>_zoom.png           4 largest changed regions, RGB | MATLAB | corrected
    _SUMMARY.txt               per-slide change counts and top transitions
"""
from __future__ import annotations
import os
import sys
from pathlib import Path
import numpy as np
import tifffile
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import ndimage as ndi

HERE = Path(__file__).resolve().parent
PY = HERE.parent
# Data (label maps, gold) is gitignored and lives next to the code in the main checkout;
# set CEREB_PY_DATA to point elsewhere (e.g. when running from a git worktree).
DATA = Path(os.environ.get("CEREB_PY_DATA", PY))
if not (DATA / "matlab_anomaly").exists():
    DATA = Path("/Users/jpurzner/Dropbox/imaging_analysis/cerebellar_segmentation/python")
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(PY / "dino"))
import shape_correct as sc
from slide_manifest import manifest

OUT = HERE / "out"
LABELS = {"20x": DATA / "matlab_anomaly" / "sandbox", "10x": DATA / "matlab_anomaly_10x" / "sandbox"}
WORK_UM = 1.0          # work at ~1 um/px: 20x -> ds 2, 10x -> ds 1
COLORS = np.array([[0, 0, 0], [30, 30, 90], [100, 180, 250], [0, 210, 230], [255, 235, 40],
                   [140, 230, 140], [255, 140, 20], [230, 20, 20], [140, 70, 200],
                   [130, 130, 130]], np.uint8)
NAMES = {0: "bg", 1: "cereb", 2: "iEGL", 3: "oEGL", 4: "IGL", 5: "ML", 6: "DWL", 7: "PCL",
         8: "DCN", 9: "anomaly"}


def load(slide):
    lab = tifffile.imread(LABELS[slide.mag] / f"{slide.slide_id}_native_labels.tif").astype(np.uint8)
    img = tifffile.imread(slide.input_path)
    if img.ndim == 3 and img.shape[-1] == 3 and img.shape[0] != 3:
        img = np.moveaxis(img, -1, 0)
    img = img[:3]
    H, W = min(lab.shape[0], img.shape[1]), min(lab.shape[1], img.shape[2])
    return lab[:H, :W], img[:, :H, :W]


def panels(sid, rgb, before, after, out):
    chg = before != after
    dim = COLORS[after].copy()
    dim[~chg] = (dim[~chg] * 0.25).astype(np.uint8)
    fig, ax = plt.subplots(1, 4, figsize=(28, 7.5))
    for a, im, t in zip(ax, [rgb, COLORS[before], COLORS[after], dim],
                        ["RGB (p27, NeuN, DAPI)", "MATLAB iter-6", "shape-corrected",
                         "changed pixels (rest dimmed)"]):
        a.imshow(im)
        a.set_title(f"{sid}  {t}")
        a.axis("off")
    plt.tight_layout()
    plt.savefig(out / f"{sid}_overview.png", dpi=60)
    plt.close()

    c = ndi.binary_closing(chg, iterations=4)
    lb, n = ndi.label(c)
    if n == 0:
        return
    sz = ndi.sum(c, lb, range(1, n + 1))
    order = np.argsort(sz)[::-1][:4]
    objs = ndi.find_objects(lb)
    H, W = after.shape
    fig, ax = plt.subplots(len(order), 3, figsize=(18, 6 * len(order)), squeeze=False)
    for r, k in enumerate(order):
        sl = objs[k]
        pad = 80
        r0, r1 = max(sl[0].start - pad, 0), min(sl[0].stop + pad, H)
        c0, c1 = max(sl[1].start - pad, 0), min(sl[1].stop + pad, W)
        for j, (im, t) in enumerate([(rgb[r0:r1, c0:c1], "RGB"),
                                     (COLORS[before[r0:r1, c0:c1]], "MATLAB iter-6"),
                                     (COLORS[after[r0:r1, c0:c1]], "shape-corrected")]):
            ax[r, j].imshow(im)
            ax[r, j].set_title(f"{sid} {t}  rows {r0}-{r1} cols {c0}-{c1} (working res)")
            ax[r, j].axis("off")
        ax[r, 2].contour(chg[r0:r1, c0:c1], levels=[0.5], colors="magenta", linewidths=0.6)
    plt.tight_layout()
    plt.savefig(out / f"{sid}_zoom.png", dpi=50)
    plt.close()


def main(filters):
    OUT.mkdir(exist_ok=True)
    rows = []
    for s in manifest():
        if filters and not any(f in s.slide_id for f in filters):
            continue
        lab_full, img = load(s)
        ds = max(int(round(WORK_UM / s.pixel_size_um)), 1)
        um = s.pixel_size_um * ds
        small = lab_full[::ds, ::ds]
        ch = sc.smooth_channels(img, ds, order=s.channel_order)
        new, log = sc.correct(small, um, ch)
        full_new = sc.apply_at_full_res(lab_full, small, new, ds)
        tifffile.imwrite(OUT / f"{s.slide_id}_shape_labels.tif", full_new, compression="zlib")

        rgb = np.stack([ch["p27"], ch["neun"], ch["dapi"]], -1)
        rgb = np.clip(rgb / np.percentile(rgb, 99.5, axis=(0, 1)), 0, 1)
        panels(s.slide_id, rgb, small, new, OUT)

        chg = small != new
        tr = {}
        for a, b in zip(small[chg], new[chg]):
            tr[(a, b)] = tr.get((a, b), 0) + 1
        top = sorted(tr.items(), key=lambda x: -x[1])[:8]
        pct = 100 * chg.sum() / max((small > 0).sum(), 1)
        line = (f"{s.slide_id} ({s.mag}): changed {chg.sum()} px at {um:.2f} um/px "
                f"({pct:.2f}% of tissue)  log={log}\n    " +
                ", ".join(f"{NAMES[a]}->{NAMES[b]} {n}" for (a, b), n in top))
        print(line, flush=True)
        rows.append(line)
    (OUT / "_SUMMARY.txt").write_text("\n".join(rows) + "\n")


if __name__ == "__main__":
    main(sys.argv[1:])
