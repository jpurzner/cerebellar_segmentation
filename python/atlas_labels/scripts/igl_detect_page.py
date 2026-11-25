"""
New IGL Step page: what the image-based IGL step (igl_detect.py) would change on every section.

Controls: held-out predictions (model trained on the other 11 controls); cKO sections: model
trained on all 12 controls. Per section the nuclei image with MATLAB's IGL outlined (cyan), the
IGL the add-only rule would add (green), IGL the model thinks is not IGL (orange, shown only)
and, on controls, the user's earlier corrections (magenta outline). The proposed IGL can be
painted and erased; db collection "igl_detect", one doc per section:
  {strokes: [{mode: "paint" | "erase", r_um, pts: [[fx, fy], ...]}], verdict, note}
fx, fy = fractions of the section image; grid cell (i, j) is inside a stroke when its centre
(i + 0.5, j + 0.5) lies within r_um / 4 grid cells of a segment (igl_detect_edits.replay).

    python scripts/igl_detect_page.py [pass tag, default loso2]   # -> outputs/igl_detect/new_igl_step.html
"""
import base64
import io
import json
import sys
from pathlib import Path

import numpy as np
import tifffile
from PIL import Image
from scipy import ndimage as ndi

sys.path.insert(0, str(Path(__file__).parent))
import igl_detect as D  # noqa: E402
import igl_gapfill_apply as GA  # noqa: E402
import igl_gapfill_page as PG  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/igl_detect"
TEMPLATE = Path(__file__).parent / "igl_detect_template.html"
G = D.G
KK = 2
RULE = dict(t_add=0.5, min_um2=10000.0, t_not=0.2)


def rule(p_igl, auto, tissue, prm=RULE):
    """add-only: confident IGL (p >= t_add) that MATLAB did not call IGL, in pieces of at least
    min_um2; and, shown only, MATLAB IGL the model gives p < t_not, in pieces of at least min_um2."""
    mat = auto == 4
    out = {}
    for name, m in (("add", (p_igl >= prm["t_add"]) & ~mat & tissue), ("not", (p_igl < prm["t_not"]) & mat)):
        lab, n = ndi.label(m, structure=np.ones((3, 3), bool))
        if n:
            size = np.bincount(lab.ravel()) * G * G
            keep = size >= prm["min_um2"]
            keep[0] = False
            m = keep[lab]
        out[name] = m
    return out


def sweep(tag):
    """the add-only rule on the held-out controls at a few settings."""
    good = ["s4_F", "s5_A", "s5_B", "s5_F", "s5_G", "s6_C"]
    rep = ["s4_C", "s5_E", "s4_A", "s5_C", "s6_G"]
    data = {}
    for s in D.controls():
        _, y, auto, tissue = D.load(s)
        data[s] = (y, auto, tissue, np.load(OUT / tag / f"{s}.npz")["p_igl"].astype(np.float32))
    rows = []
    for t in (0.5, 0.7, 0.85):
        for mn in (2000.0, 10000.0):
            prm = dict(RULE, t_add=t, min_um2=mn)
            y, auto, tissue, p = data["s1_2"]
            new = (auto == 4) | rule(p, auto, tissue, prm)["add"]
            gold = D.iou(new & (y > 0), y == 4)
            gadd = float(np.mean([rule(data[s][3], data[s][1], data[s][2], prm)["add"].sum() * G * G for s in good]))
            found, outside = {}, []
            for s in rep:
                y, auto, tissue, p = data[s]
                a = rule(p, auto, tissue, prm)["add"]
                r = (y == 4) & (auto != 4) & tissue
                found[s] = round(float(a[r].mean()), 2) if r.any() else None
                outside.append(float((a & ~r).sum() * G * G))
            rows.append(dict(t_add=t, min_um2=int(mn), gold_iou=round(gold, 3), good_added_um2=int(gadd), found=found,
                             outside_um2=int(np.mean(outside))))
    return rows


def edge(m):
    return m & ~ndi.binary_erosion(m)


def layers(sid, auto, ch_rule, corr):
    """what the editor needs: the nuclei image (2 um per pixel), a transparent overlay with
    MATLAB's IGL outline (cyan), the model's disagreement hint (orange) and the user's earlier
    corrections (magenta), and the proposed IGL on the 4 um grid as a grey PNG
    (0 none, 100 MATLAB's IGL kept, 200 added by the new step)."""
    S, px = D.stains(sid)
    H, W = auto.shape
    img = PG.native_crop(np.stack([S["nuclei"]]), px, (0, H, 0, W), 0)
    b = io.BytesIO()
    Image.fromarray(img.astype(np.uint8)).save(b, "JPEG", quality=82)
    base = "data:image/jpeg;base64," + base64.b64encode(b.getvalue()).decode()
    up = lambda m: np.repeat(np.repeat(m, KK, 0), KK, 1)
    ov = np.zeros((H * KK, W * KK, 4), np.uint8)
    ov[up(ch_rule["not"])] = (255, 150, 30, 70)
    ov[up(ndi.binary_dilation(edge(ch_rule["not"])))] = (255, 150, 30, 230)
    ov[up(edge(auto == 4))] = (59, 224, 255, 255)
    if corr is not None and corr.any():
        ov[up(ndi.binary_dilation(edge(corr)))] = (255, 47, 214, 255)
    b = io.BytesIO()
    Image.fromarray(ov, "RGBA").save(b, "PNG", optimize=True)
    outline = "data:image/png;base64," + base64.b64encode(b.getvalue()).decode()
    prop = np.zeros((H, W), np.uint8)
    prop[auto == 4] = 100
    prop[ch_rule["add"]] = 200
    b = io.BytesIO()
    Image.fromarray(prop, "L").save(b, "PNG", optimize=True)
    proposed = "data:image/png;base64," + base64.b64encode(b.getvalue()).decode()
    return base, outline, proposed


def main(tag="loso2"):
    L = json.loads((OUT / f"{tag}.json").read_text())
    L1 = json.loads((OUT / "loso.json").read_text())
    secs = []
    for sid in GA.section_list():
        s = GA.short(sid)
        group = GA.info(sid)[1]
        _, y, auto, tissue = D.load(s)
        f = OUT / tag / f"{s}.npz" if group == "control" else OUT / "apply2" / f"{s}.npz"
        p = np.load(f)["p_igl"].astype(np.float32)
        ch = rule(p, auto, tissue)
        corr = ((y == 4) & (auto != 4) & tissue) if group == "control" else None
        st = dict(igl_mm2=round(float((auto == 4).sum() * G * G / 1e6), 3), add_um2=int(ch["add"].sum() * G * G),
                  not_um2=int(ch["not"].sum() * G * G))
        if corr is not None and corr.any():
            st.update(corr_um2=int(corr.sum() * G * G), corr_found=round(float(ch["add"][corr].mean()), 2))
        held = next((r for r in L["sections"] if r["section"] == s), None)
        base, outline, proposed = layers(sid, auto, ch, corr)
        secs.append(dict(short=s, group=group, mag="20x" if GA.info(sid)[2] else "10x", w_um=round(auto.shape[1] * G),
                         h_um=round(auto.shape[0] * G), grid=[int(auto.shape[0]), int(auto.shape[1])],
                         base=base, outline=outline, proposed=proposed, stats=st, held=held))
        print(s, st, flush=True)
    sw = sweep(tag)
    pick = next(r for r in sw if r["t_add"] == RULE["t_add"] and r["min_um2"] == int(RULE["min_um2"]))
    gold = next(r for r in L["sections"] if r["section"] == "s1_2")
    gold1 = next(r for r in L1["sections"] if r["section"] == "s1_2")
    cko = [x for x in secs if x["group"] == "cKO"]
    findings = [
        f"<b>On its own the new step is not better than MATLAB at the IGL's edges.</b> Against your hand-labelled s1_2 it scores {gold['igl_iou_new']:.2f} "
        f"(first pass {gold1['igl_iou_new']:.2f}) where MATLAB scores {gold['igl_iou_matlab']:.2f}; on the other held-out controls it agrees with the reviewed IGL at "
        f"{L['summary']['igl_iou_new']:.2f} on average. The second pass, which looks at the layers around each spot, helps a little. So it does not replace MATLAB's IGL.",
        f"<b>It does find what MATLAB loses to bleaching.</b> Without having seen them, the add-only rule (probability at least {RULE['t_add']}, pieces of at least "
        f"{int(RULE['min_um2']):,} µm²) recovers {round(100 * pick['found']['s5_E'])}% of the lobule you outlined in s5_E and {round(100 * pick['found']['s4_C'])}% of your D4 / E4 "
        f"traces and the U in s4_C, leaves your gold standard unchanged ({pick['gold_iou']:.3f}) and adds about {pick['good_added_um2']:,} µm² per good control. "
        "Check those green patches on the controls: they are either MATLAB misses nobody flagged, or the step's errors.",
        "<b>It does not find short bridges</b> (s4_A, s5_C, s6_G: white-matter labels across a wide ribbon): the gap filler stays for those.",
        f"<b>On the EZH2 cKO sections</b> it would add {', '.join(x['short'] + ' ' + format(x['stats']['add_um2'], ',') + ' µm²' for x in cko)}. "
        "Orange marks MATLAB IGL the model thinks is not IGL (for example the inner EGL pulled into the IGL); it is shown only, not changed.",
        "<b>Fix what it got wrong by painting and erasing</b> on each section's image (brush and eraser, any size; erasing works on MATLAB's IGL too). "
        "Your strokes are saved and replayed onto the label grid: they become that section's corrected IGL and the test set for the next round.",
    ]
    data = dict(rule=RULE, sections=secs, loso1=L1["summary"], loso2=L["summary"], sweep=sw, findings=findings,
                per_section=[dict(r, pass1=next(q["igl_iou_new"] for q in L1["sections"] if q["section"] == r["section"])) for r in L["sections"]])
    page = TEMPLATE.read_text().replace("__DATA__", json.dumps(data))
    dst = OUT / "new_igl_step.html"
    dst.write_text(page)
    print(dst, round(dst.stat().st_size / 1e6, 2), "MB")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "loso2")
