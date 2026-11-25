"""
Paint the user's IGL Trace answers (https://claude.ai/artifact/3QkoRZC9SYnWdmxMCCJasF, db
collection "igl_trace", saved to outputs/igl_trace/answers_<section>_round<N>.json) onto the
section's anatomy-fixed labels, the same way as the D4 / E4 traces in igl_finalize:

  - "add" strokes (missing IGL, clicked along its middle in order): a smooth line through the
    clicks in click order, painted at the section's ribbon width (median of the ribbon
    description's branch widths) over every tissue label, EGL included (the user traced over
    it); ends within 60 um of IGL or an earlier stroke are joined to it.
  - "area" strokes (the outline of a region of missing IGL, clicked around its edge): the
    closed outline and everything inside it become IGL, whatever the labels say (user, s5_E:
    "the area needs to be filled, right now it is just the perimeter of the missed region").
    outputs/igl_trace/interpretation.json can reread a stroke saved with another kind.
  - "wrong" strokes: listed, not applied (what the IGL really is there comes from the notes).
  - sections marked "nothing to fix" are left as they are.

    python scripts/igl_trace_apply.py s5_E_p27 [...]

Output: outputs/igl_fixed/<section>_igl_fixed_labels.tif, <section>_traces.npz (painted masks
and lines, 4 um grid), the strokes in grid coordinates under the section in
outputs/igl_fixed/user_traces.json, and a log in outputs/igl_fixed/igl_traced.json.
"""
import json
import sys
from pathlib import Path

import numpy as np
import tifffile
from scipy import ndimage as ndi

sys.path.insert(0, str(Path(__file__).parent))
import igl_finalize as F  # noqa: E402
import igl_gapfill_apply as GA  # noqa: E402
import igl_ribbon as RB  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
ANS = ROOT / "outputs/igl_trace"
OUT = ROOT / "outputs/igl_fixed"


def answers(short):
    """the latest saved answers for a section (highest round number)."""
    files = sorted(ANS.glob(f"answers_{short}_round*.json"), key=lambda p: int(p.stem.split("round")[-1]))
    if not files:
        return None, None
    d = json.loads(files[-1].read_text())
    return d.get("data", d), files[-1].name


def outline_fill(g, pts, tear):
    """the closed outline (smooth line through the clicks, back to the first) and its inside,
    within tissue: the region the user outlined."""
    H, W = g["igl"].shape
    P = list(pts) + [list(pts[0])]
    _, _, C, _ = F.trace_band(g, P, tear, ordered=True, egl_full=True, width_um=8.0)
    line = np.zeros((H, W), bool)
    line[np.clip(np.round(C[:, 0]).astype(int), 0, H - 1), np.clip(np.round(C[:, 1]).astype(int), 0, W - 1)] = True
    line = ndi.binary_dilation(line)
    return ndi.binary_fill_holes(line) & (g["tissue"] | g["anom"]) & ~tear, C


def apply(sid):
    short = GA.short(sid)
    D, src = answers(short)
    px = GA.info(sid)[0]
    native = tifffile.imread(GA.RE.FX / f"{sid}_fixed_labels.tif")
    g = RB.grid_labels(native, px)
    H, W = g["igl"].shape
    tear = np.zeros((H, W), bool)
    log = dict(source=src, status=D.get("status") if D else None, note=D.get("note") if D else None, strokes=[])
    if not D or D.get("status") == "nothing":
        log["applied"] = False
        return log
    rib = RB.describe(g)
    w_sec = float(np.median([np.median(b["width"]) for b in rib["branches"]]))
    paint = np.zeros((H, W), bool)
    allow_egl = np.zeros((H, W), bool)
    masks, grid_strokes = {}, []
    I = json.loads((ANS / "interpretation.json").read_text()).get(sid, {}) if (ANS / "interpretation.json").exists() else {}
    for k, s in enumerate(D.get("strokes", [])):
        pts = [[round(fy * H, 2), round(fx * W, 2)] for fx, fy in s["pts"]]
        kind = I.get(str(k), {}).get("kind", s["kind"])
        grid_strokes.append(dict(kind=kind, saved_as=s["kind"], pts=pts))
        if kind == "area" and len(pts) >= 3:
            area, C = outline_fill(g, pts, tear)
            paint |= area
            allow_egl |= area & g["egl"]
            masks[f"stroke{k}_band"], masks[f"stroke{k}_line"] = area, C
            log["strokes"].append(dict(kind="area", points=len(pts), applied=True, area_um2=round(float(area.sum() * RB.GRID ** 2)),
                                       painted_um2={n: round(float((area & g[n]).sum() * RB.GRID ** 2)) for n in ("ml", "wm", "egl", "unass", "anom")}))
            continue
        if kind != "add" or len(pts) < 2:
            log["strokes"].append(dict(kind=kind, points=len(pts), applied=False))
            continue
        band, egl, C, info = F.trace_band(g, pts, tear, connect_um=60.0, ordered=True, egl_full=True, connect_start=True,
                                          also=paint if paint.any() else None, width_um=w_sec)
        paint |= band | egl
        allow_egl |= egl
        masks[f"stroke{k}_band"], masks[f"stroke{k}_line"] = band | egl, C
        log["strokes"].append(dict(kind="add", points=len(pts), applied=True, **info))
    out = F.stamp(native, paint, paint, allow_egl, px=px)
    tifffile.imwrite(OUT / f"{sid}_igl_fixed_labels.tif", out, compression="zlib")
    np.savez_compressed(OUT / f"{sid}_traces.npz", painted=paint, **masks)
    T = json.loads((OUT / "user_traces.json").read_text()) if (OUT / "user_traces.json").exists() else {}
    T[sid] = dict(trace_page=dict(source=f"IGL Trace page, {src}", grid=[H, W], strokes=grid_strokes,
                                  how="fractions of the whole-section image x grid shape (image exactly 2 um per pixel)"))
    (OUT / "user_traces.json").write_text(json.dumps(T, indent=1))
    log.update(applied=True, width_um=round(w_sec), painted_um2=round(float(paint.sum() * RB.GRID ** 2)),
               igl_added_mm2=round(float(((out == 4) & (native != 4)).sum() * px * px / 1e6), 4),
               changed_from={k: round(float(((out == 4) & (native == v)).sum() * px * px)) for k, v in
                             (("unassigned", 1), ("iEGL", 2), ("oEGL", 3), ("ML", 5), ("DWL", 6), ("PCL", 7), ("DCN", 8), ("anomaly", 9))})
    return log


def main(sids):
    p = OUT / "igl_traced.json"
    L = json.loads(p.read_text()) if p.exists() else {}
    for sid in sids:
        L[sid] = apply(sid)
        print(sid, {k: v for k, v in L[sid].items() if k != "strokes"})
        for s in L[sid]["strokes"]:
            print("   ", s)
    p.write_text(json.dumps(L, indent=1))


if __name__ == "__main__":
    main(sys.argv[1:] or ["s5_E_p27", "s5_B_p27"])
