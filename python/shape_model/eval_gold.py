"""Score the shape correction against the s1_2 gold standard (read-only; gold is never written).

Gold uses 4=ML, 5=IGL (swapped vs the canonical index) and has no DCN, so gold is
swapped to canonical and DCN is merged into DWL before IoU. Evaluated at ds=2 (~1 um/px).
"""
from __future__ import annotations
import os
import sys
from pathlib import Path
import numpy as np
import tifffile

HERE = Path(__file__).resolve().parent
PY = HERE.parent
# Data (label maps, gold) is gitignored and lives next to the code in the main checkout;
# set CEREB_PY_DATA to point elsewhere (e.g. when running from a git worktree).
DATA = Path(os.environ.get("CEREB_PY_DATA", PY))
if not (DATA / "matlab_anomaly").exists():
    DATA = Path("/Users/jpurzner/Dropbox/imaging_analysis/cerebellar_segmentation/python")
sys.path.insert(0, str(HERE))
import shape_correct as sc

SID = "2018_05_22_s1_2_p27-0021"
GOLD = DATA / "labelled" / f"{SID}_labelled.tiff"
PRED = DATA / "matlab_anomaly" / "sandbox" / f"{SID}_native_labels.tif"
IMG = Path("/Users/jpurzner/Dropbox/images/edu_repeat/p27") / SID / f"{SID}_fused_crop.tif"
DS, UM = 2, 0.5119049 * 2


def iou_table(lab, g):
    o = {}
    for L, n in [(5, "ML"), (7, "PCL"), (4, "IGL"), (2, "iEGL"), (3, "oEGL")]:
        a, b = lab == L, g == L
        o[n] = (a & b).sum() / max((a | b).sum(), 1)
    for n, Ls in [("DWL", (6, 8)), ("EGL", (2, 3))]:
        a, b = np.isin(lab, Ls), np.isin(g, Ls)
        o[n] = (a & b).sum() / max((a | b).sum(), 1)
    return o


def main():
    gold = tifffile.imread(GOLD).astype(np.uint8)
    g = gold.copy()
    g[gold == 4], g[gold == 5] = 5, 4
    p = tifffile.imread(PRED).astype(np.uint8)
    img = tifffile.imread(IMG)
    H = min(g.shape[0], p.shape[0], img.shape[1])
    W = min(g.shape[1], p.shape[1], img.shape[2])
    g, p = g[:H:DS, :W:DS], p[:H:DS, :W:DS]
    ch = sc.smooth_channels(img[:, :H, :W], DS)
    new, log = sc.correct(p, UM, ch)
    K = ["ML", "PCL", "DWL", "EGL", "iEGL", "oEGL", "IGL"]
    b, a = iou_table(p, g), iou_table(new, g)
    chg = new != p
    print("IoU vs gold (s1_2, ds2)")
    print("          " + "  ".join(f"{k:>6s}" for k in K))
    print("MATLAB    " + "  ".join(f"{b[k]:6.3f}" for k in K))
    print("shape     " + "  ".join(f"{a[k]:6.3f}" for k in K))
    print(f"changed {chg.sum()} px; now correct {(new[chg] == g[chg]).mean():.2f}, "
          f"were correct {(p[chg] == g[chg]).mean():.2f}")
    print(f"gold ML labelled DWL/DCN: {((g == 5) & np.isin(p, (6, 8))).sum()} -> "
          f"{((g == 5) & np.isin(new, (6, 8))).sum()}")
    print(f"gold DWL labelled ML/PCL: {((g == 6) & np.isin(p, (5, 7))).sum()} -> "
          f"{((g == 6) & np.isin(new, (5, 7))).sum()}")
    print("log:", log)


if __name__ == "__main__":
    main()
