"""
IGL Trace page: whole-section views of the controls whose IGL the user graded bad and the
gap filler cannot repair (s5_B, s5_E) and of the 7 EZH2 cKO sections, with tools to trace
missing IGL along its course, outline missing IGL areas and mark IGL that is wrong, in the
image frame of the 4 um label grid (nuclei: ch1 at 10x, ch2 at 20x; NeuN: ch2 at 10x, ch1 at 20x).

    python scripts/igl_trace_page.py [section ...]      # default s5_B, s5_E + the cKO sections
    -> outputs/igl_trace/igl_trace.html

Answers: db collection "igl_trace", one doc per section (short name):
  {strokes: [{kind: "add" | "wrong", pts: [[fx, fy], ...]}], note, status}
fx, fy = fractions of the section image's width / height; grid row = fy * H, column = fx * W
(H, W = the 4 um grid shape, in DATA.sections[].grid).
"""
import base64
import io
import json
import sys
from pathlib import Path

import numpy as np
import tifffile
from PIL import Image, ImageDraw
from scipy import ndimage as ndi

sys.path.insert(0, str(Path(__file__).parent))
import igl_gapfill_apply as GA  # noqa: E402
import igl_gapfill_page as PG  # noqa: E402
import igl_ribbon as RB  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/igl_trace"
TEMPLATE = Path(__file__).parent / "igl_trace_template.html"
G = RB.GRID
KK = 2                      # display pixels per grid pixel (2 um per display pixel)
JPG_Q = 80
IGL_LINE = (59, 224, 255, 255)
END = (255, 255, 255, 255)


def jpg(gray, q=80):
    b = io.BytesIO()
    Image.fromarray(gray.astype(np.uint8)).save(b, "JPEG", quality=q)
    return "data:image/jpeg;base64," + base64.b64encode(b.getvalue()).decode()


def png(img):
    b = io.BytesIO()
    img.save(b, "PNG", optimize=True)
    return "data:image/png;base64," + base64.b64encode(b.getvalue()).decode()


def channel(raw, px, shape, ch):
    H, W = shape
    a = PG.native_crop(raw, px, (0, H, 0, W), ch)
    hh, ww = H * KK, W * KK
    return np.pad(a, ((0, max(hh - a.shape[0], 0)), (0, max(ww - a.shape[1], 0))))[:hh, :ww]


def overlay(igl, tips, shape):
    """transparent layer: outline of the IGL labels (cyan) and the ribbon ends (white rings)."""
    H, W = shape
    edge = igl & ~ndi.binary_erosion(igl)
    rgba = np.zeros((H * KK, W * KK, 4), np.uint8)
    rgba[np.repeat(np.repeat(edge, KK, 0), KK, 1)] = IGL_LINE
    im = Image.fromarray(rgba, "RGBA")
    d = ImageDraw.Draw(im)
    for y, x in tips:
        X, Y = (x + 0.5) * KK, (y + 0.5) * KK
        d.ellipse([X - 7, Y - 7, X + 7, Y + 7], outline=END, width=2)
    return png(im)


def section(sid):
    px, group, big = GA.info(sid)
    A = json.loads((ROOT / "outputs/igl_gapfill/apply.json").read_text())[sid]
    Z = np.load(ROOT / f"outputs/igl_gapfill/{sid}_gapfill.npz")
    code = Z["code"]
    raw = tifffile.imread((GA.EV.SB20 if big else GA.EV.SB10) / f"{sid}_native.tif")
    shape = code.shape
    native = tifffile.imread(GA.RE.FX / f"{sid}_fixed_labels.tif")
    g = RB.grid_labels(native, px)
    rib = RB.describe(g)
    ends = [e["pos"] for e in rib["ends"]]
    reasons = sorted({c["status"].replace("rejected: ", "").split(",")[0] for c in A["candidates"]})
    result, result_text = None, None
    L = json.loads((ROOT / "outputs/igl_fixed/igl_traced.json").read_text()) if (ROOT / "outputs/igl_fixed/igl_traced.json").exists() else {}
    if L.get(sid, {}).get("applied"):
        out = tifffile.imread(ROOT / f"outputs/igl_fixed/{sid}_igl_fixed_labels.tif")
        g2 = RB.grid_labels(out, px)
        rgb = PG.PAL[GA.grid_code(g2)].copy()
        rgb[g2["igl"] & ~g["igl"]] = PG.FILL
        result = png(Image.fromarray(PG.up(rgb, KK)))
        lg = L[sid]
        parts = []
        for x in lg["strokes"]:
            if not x.get("applied"):
                continue
            if x["kind"] == "area":
                parts.append(f"the area inside your outline filled ({x['area_um2']:,} µm²)")
            else:
                parts.append(f"a {x['length_um']:,} µm stroke painted {x['width_um']} µm wide")
        result_text = ("Painted: " + "; ".join(parts) + f". {lg['igl_added_mm2']} mm² of IGL added, over "
                       + ", ".join(f"{v:,} µm² {k}" for k, v in lg["changed_from"].items() if v >= 500) + ".")
    nuc, neun = (2, 1) if big else (1, 2)             # 20x: ch2 nuclei, ch1 NeuN; 10x: ch1 nuclei, ch2 NeuN
    layers = dict(nuc=jpg(channel(raw, px, shape, nuc), JPG_Q), neun=jpg(channel(raw, px, shape, neun), JPG_Q),
                  labels=png(Image.fromarray(PG.up(PG.PAL[code], KK))))
    if result:
        layers["result"] = result
    judged = []
    jf = ROOT / "outputs/igl_gapfill/judged_v2.json"
    if jf.exists():
        for m in json.loads(jf.read_text()):
            if m["sid"] == sid and m["verdict"] and not m["id"].startswith("ref-"):
                judged.append(f"{m['what']}: you said {m['verdict']}")
    return dict(sid=sid, short=GA.short(sid), group=group, grid=[int(shape[0]), int(shape[1])], mag="20x" if big else "10x", judged=judged,
                w_um=round(shape[1] * G), h_um=round(shape[0] * G), layers=layers,
                overlay=overlay(g["igl"], [tuple(map(float, e)) for e in ends], shape),
                n_ends=len(ends), candidates=len(A["candidates"]), reasons=reasons,
                igl_mm2=round(float(g["igl"].sum() * G * G / 1e6), 3), result_text=result_text)


def main(sids):
    OUT.mkdir(parents=True, exist_ok=True)
    secs = [section(s) for s in sids]
    page = TEMPLATE.read_text().replace("__DATA__", json.dumps(dict(sections=secs)))
    dst = OUT / "igl_trace.html"
    dst.write_text(page)
    print(dst, round(dst.stat().st_size / 1e6, 2), "MB", [(s["short"], s["grid"], s["n_ends"]) for s in secs])


if __name__ == "__main__":
    names = sys.argv[1:] or (["s5_B_p27", "s5_E_p27"] + [x for x in GA.section_list() if GA.info(x)[1] == "cKO"])
    main(names)
