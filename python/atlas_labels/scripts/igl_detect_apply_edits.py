"""
Write the user's paint / erase edits from the New IGL Step page (round 1, 2026-11-20; saved in
outputs/igl_detect/edits_round1/igl_detect/<section>.json) into corrected label maps for the
EZH2 cKO sections they edited.

On the 4 um grid: IGL = the proposal the page showed (MATLAB's IGL plus the new step's additions)
with the painted cells added and the erased cells removed (igl_detect_edits). At native
resolution: every tissue pixel in a cell that became IGL is set to IGL (4); an IGL pixel in a
cell that lost it takes the layer the model rates most likely among the others (EGL -> iEGL 2,
the inner EGL being what is pulled into the IGL; ML/PCL -> ML 5; white matter -> DWL 6).
Controls are left to outputs/final_labels/controls (their edits started from MATLAB's labels,
not the repaired ones, so they are kept as test data only).

    python scripts/igl_detect_apply_edits.py          # -> outputs/igl_detect/corrected_round1/
"""
import glob
import json
import sys
from pathlib import Path

import numpy as np
import tifffile

sys.path.insert(0, str(Path(__file__).parent))
import igl_detect as D  # noqa: E402
import igl_detect_edits as E  # noqa: E402
import igl_detect_page as P  # noqa: E402
import igl_gapfill_apply as GA  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/igl_detect/corrected_round1"
EDITS = ROOT / "outputs/igl_detect/edits_round1/igl_detect"
CODE_OF = {0: 2, 1: 5, 3: 6}          # model class index (EGL, ML/PCL, IGL, WM) -> label code


def upsample(grid, shape):
    """grid cells to native pixels, the same corner-aligned mapping as igl_ribbon.grid_labels."""
    H, W = shape
    gh, gw = grid.shape
    ys = np.clip(np.round(np.arange(H) * (gh - 1) / max(H - 1, 1)).astype(int), 0, gh - 1)
    xs = np.clip(np.round(np.arange(W) * (gw - 1) / max(W - 1, 1)).astype(int), 0, gw - 1)
    return grid[np.ix_(ys, xs)]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for f in sorted(glob.glob(str(EDITS / "*.json"))):
        d = json.loads(Path(f).read_text())
        doc = d.get("data", d)
        short = Path(f).stem
        if short in D.controls():
            continue
        sid = next(s for s in GA.section_list() if GA.short(s) == short)
        px = GA.info(sid)[0]
        _, _, auto, tissue = D.load(short)
        Z = np.load(ROOT / f"outputs/igl_detect/apply2/{short}.npz")
        Pall = Z["P"].astype(np.float32)
        ch = P.rule(Z["p_igl"].astype(np.float32), auto, tissue)
        prop = (auto == 4) | ch["add"]
        new = E.corrected_igl(prop, doc.get("strokes", []))
        other = Pall.copy()
        other[2] = -1                                   # best class other than IGL
        alt = np.vectorize(CODE_OF.get)(np.argmax(other, 0)).astype(np.uint8)
        native = tifffile.imread(GA.RE.FX / f"{sid}_fixed_labels.tif")
        gain, loss = upsample(new & ~(auto == 4), native.shape), upsample(~new & (auto == 4), native.shape)
        out = native.copy()
        out[gain & (native > 0)] = 4
        lost = loss & (native == 4)
        out[lost] = upsample(alt, native.shape)[lost]
        tifffile.imwrite(OUT / f"{short}_labels.tif", out, compression="zlib", resolution=(1 / px, 1 / px))
        a2 = px * px / 1e6
        r = dict(section=short, strokes=len(doc.get("strokes", [])), igl_mm2_matlab=round(float((native == 4).sum() * a2), 4),
                 igl_mm2_corrected=round(float((out == 4).sum() * a2), 4),
                 igl_removed_to={n: round(float(((native == 4) & (out == c)).sum() * a2), 4) for n, c in (("iEGL", 2), ("ML", 5), ("DWL", 6))},
                 igl_added_mm2=round(float(((out == 4) & (native != 4)).sum() * a2), 4))
        rows.append(r)
        print(r)
    (OUT / "manifest.json").write_text(json.dumps(dict(
        what="EZH2 cKO label maps with the user's round-1 IGL edits (New IGL Step page); MATLAB labels elsewhere",
        sections=rows), indent=1))


if __name__ == "__main__":
    main()
