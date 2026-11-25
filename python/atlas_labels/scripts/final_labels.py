"""
Freeze the final label maps of the control sections in one folder, with a manifest of what
was done to each (2026-11-19, after the user's last reviews of the IGL repairs).

    python scripts/final_labels.py           # -> outputs/final_labels/controls/

Per section the best version available:
  s1_2                the user's hand-labelled gold standard (labelled/2018_05_22_s1_2_p27-0021_labelled.tiff),
                      converted to the pipeline's codes (its IGL / ML codes are swapped); the file itself is untouched
  s4_A, s4_C, s5_E    outputs/igl_fixed/<section>_igl_fixed_labels.tif (gap filler with the user's bleach boxes and
                      markers, accepted box / smear fixes, the user's traced and outlined IGL)
  other controls      outputs/igl_gapfill/<section>_gapfilled_labels.tif (gap filler bridges the user judged correct;
                      identical to the anatomy-fixed labels where it made none)
Each is written at native resolution (uint8, zlib) as <short>_labels.tif with its pixel size in the TIFF
resolution tags, and listed in manifest.json / manifest.csv with source, repairs, the user's verdicts and the
IGL area before (anatomy-fixed automatic labels) and after.
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
import tifffile

sys.path.insert(0, str(Path(__file__).parent))
import igl_gapfill_apply as GA  # noqa: E402
import igl_ribbon_eval as RE  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/final_labels/controls"
CODES = {0: "background", 1: "cerebellum, no layer assigned", 2: "iEGL (inner EGL)", 3: "oEGL (outer EGL)", 4: "IGL", 5: "ML",
         6: "DWL (deep white matter)", 7: "PCL", 8: "DCN (deep cerebellar nuclei)", 9: "anomaly (tear / bleach)"}


def repairs(sid):
    """what was done to the section, from the logs, with the user's verdicts."""
    short = GA.short(sid)
    A = json.loads((ROOT / "outputs/igl_gapfill/apply.json").read_text())[sid]
    J = json.loads((ROOT / "outputs/igl_fixed/igl_fixed.json").read_text())
    Tr = json.loads((ROOT / "outputs/igl_fixed/igl_traced.json").read_text()) if (ROOT / "outputs/igl_fixed/igl_traced.json").exists() else {}
    out = []
    if sid in J:
        e = J[sid]
        for b in e["bridges"]:
            out.append(dict(kind="gap bridge", detail=f"{round(b['chord_um'])} um, found by {b['found_by']}, {b['method']}",
                            verdict="correct (IGL Gap Fill, 2nd review)"))
        if "box_fix_B" in e:
            out.append(dict(kind="bleach box fix B", detail=f"box {e['box_fix_B']['box']}, {e['box_fix_B']['changed_um2']} um2 changed",
                            verdict="accepted (IGL Problem Check)"))
        if "smear_fix_folium6_dapi" in e:
            out.append(dict(kind="smear fix, folium 6 (DAPI)", detail=f"{e['smear_fix_folium6_dapi']['changed_um2']} um2 changed",
                            verdict="accepted (IGL Problem Check)"))
        for name, t in e.get("user_traces", {}).items():
            out.append(dict(kind=f"user trace at marker {name}", detail=f"{t['length_um']} um painted {t['width_um']} um wide"
                            + "".join(f"; stroke {s['length_um']} um" for s in t.get("strokes", [])), verdict="right (IGL Gap Fill, round 6)"))
    elif A["joined"] and A["group"] == "control":
        for b in A["candidates"]:
            if b["status"] == "joined":
                out.append(dict(kind="gap bridge", detail=f"{round(b['chord_um'])} um, found by {b['found_by']}, {b['method']}",
                                verdict="correct (IGL Gap Fill, 2nd review)"))
    if sid in Tr and Tr[sid].get("applied"):
        for s in Tr[sid]["strokes"]:
            if s.get("applied"):
                what = f"outlined area filled, {s['area_um2']} um2" if s["kind"] == "area" else f"traced ribbon {s['length_um']} um"
                out.append(dict(kind="user trace (IGL Trace page)", detail=what, verdict="looks good (IGL Trace, round 2)"))
    elif sid in Tr:
        out.append(dict(kind="reviewed on IGL Trace", detail="nothing to fix", verdict=Tr[sid].get("status")))
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for sid in GA.section_list():
        px, group, big = GA.info(sid)
        if group != "control":
            continue
        short = GA.short(sid)
        auto = tifffile.imread(RE.FX / f"{sid}_fixed_labels.tif")
        if short == "s1_2":
            lab, src = RE.gold(), "gold standard (hand-labelled by the user), IGL / ML codes converted"
        elif (ROOT / f"outputs/igl_fixed/{sid}_igl_fixed_labels.tif").exists():
            lab, src = tifffile.imread(ROOT / f"outputs/igl_fixed/{sid}_igl_fixed_labels.tif"), f"outputs/igl_fixed/{sid}_igl_fixed_labels.tif"
        else:
            lab, src = tifffile.imread(ROOT / f"outputs/igl_gapfill/{sid}_gapfilled_labels.tif"), f"outputs/igl_gapfill/{sid}_gapfilled_labels.tif"
        assert lab.shape == auto.shape, sid
        dst = OUT / f"{short}_labels.tif"
        tifffile.imwrite(dst, lab.astype(np.uint8), compression="zlib", resolution=(1 / px, 1 / px),
                         metadata=None, description=json.dumps(dict(section=sid, um_per_px=px, codes=CODES)))
        a2 = px * px / 1e6
        rows.append(dict(section=short, sample=sid, genotype=group, magnification="20x" if big else "10x", um_per_px=px,
                         file=dst.name, source=src, igl_mm2_automatic=round(float((auto == 4).sum() * a2), 4),
                         igl_mm2_final=round(float((lab == 4).sum() * a2), 4),
                         changed_mm2=round(float((lab != auto).sum() * a2), 4) if short != "s1_2" else None,
                         repairs=repairs(sid) if short != "s1_2" else [dict(kind="gold standard", detail="user's manual labels", verdict="reference")]))
        print(f"{short:5s} {rows[-1]['magnification']} IGL {rows[-1]['igl_mm2_automatic']} -> {rows[-1]['igl_mm2_final']} mm2, "
              f"{len(rows[-1]['repairs'])} repairs  <- {src}")
    (OUT / "manifest.json").write_text(json.dumps(dict(codes=CODES, sections=rows), indent=1))
    with open(OUT / "manifest.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["section", "sample", "genotype", "magnification", "um_per_px", "file", "source", "igl_mm2_automatic", "igl_mm2_final",
                    "changed_mm2", "repairs"])
        for r in rows:
            w.writerow([r["section"], r["sample"], r["genotype"], r["magnification"], r["um_per_px"], r["file"], r["source"],
                        r["igl_mm2_automatic"], r["igl_mm2_final"], r["changed_mm2"],
                        " | ".join(f"{x['kind']}: {x['detail']} [{x['verdict']}]" for x in r["repairs"])])
    print(OUT, len(rows), "sections")


if __name__ == "__main__":
    main()
