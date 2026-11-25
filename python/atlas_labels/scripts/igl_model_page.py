"""
Close-up views of the IGL model's predictions at each user marker, for the IGL Problem
Check page (igl_check_page.py adds them as it["model"]).

Each marker: a 600 x 600 um crop of the IF image at 1 um per pixel (R = p27, G = NeuN,
B = DAPI, rotated to the page frame), centred on the marker, with overlay layers:
  now   - IGL (yellow) and white matter (orange) outlines of the anatomy-fixed labels
  last  - IGL outline after the last fix (igl_fixed), light green
  A/B/C - the model's candidates: centre line + band edge (cyan, magenta, white)
"""
import json
import sys
from pathlib import Path

import numpy as np
import tifffile
from PIL import Image
from scipy import ndimage as ndi
from skimage.draw import line as draw_line

sys.path.insert(0, str(Path(__file__).parent))
import igl_check_page as P  # noqa: E402
import igl_finalize as F  # noqa: E402
import igl_repair_eval as EV  # noqa: E402
import igl_ribbon as RB  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
MOD = ROOT / "outputs/igl_model"
G = RB.GRID
HALF_UM, STEP_UM = 300.0, 1.0
COLS = {"A": (59, 224, 255), "B": (255, 90, 209), "C": (255, 255, 255)}
LABELS = {("s4_A_p27", "m1"): "marker 1 · folium 9 · E6",          # the user's own marker list
          ("s4_C_p27", "C3"): "marker 1 · folium 4 · C3", ("s4_C_p27", "D4"): "marker 2 · folium 4 · D4",
          ("s4_C_p27", "E4"): "marker 3 · folium 7 · E4"}
TARGET_TEXT = {"end": "joins another ribbon end", "ribbon": "joins the side of the ribbon"}


def _frame(tr, centre_grid):
    A, off, crop = np.array(tr["A"]), np.array(tr["offset"]), np.array(tr["crop"], float)
    to_disp = lambda q: np.linalg.solve(A, np.asarray(q, float) * G / (F.PX * P.DS) - off) - crop
    md = to_disp(centre_grid)
    n = int(2 * HALF_UM / STEP_UM)
    s = STEP_UM / (F.PX * P.DS)
    ii, jj = np.mgrid[0:n, 0:n].astype(float)
    disp = np.stack([md[0] + (ii - n / 2) * s, md[1] + (jj - n / 2) * s])
    nat_ds = np.einsum("ij,jhw->ihw", A, disp + crop[:, None, None]) + off[:, None, None]
    native = nat_ds * P.DS
    to_crop = lambda q: (to_disp(q) - md) / s + n / 2
    return native, n, to_crop


def _outline(m, k=1):
    return m & ~ndi.binary_erosion(m, iterations=k)


def _rgba(mask, col, n):
    o = np.zeros((n, n, 4), np.uint8)
    o[mask] = (*col, 255)
    return o


def build(sid, tr):
    z = MOD / f"{sid}_marker_candidates.npz"
    js = MOD / "marker_candidates.json"
    if not z.exists() or not js.exists():
        return []
    arr = np.load(z)
    info = json.loads(js.read_text()).get(sid, {})
    raw = tifffile.imread(P.SB10 / f"{sid}_native.tif").astype(np.float32)
    lab0 = tifffile.imread(EV.RE.FX / f"{sid}_fixed_labels.tif")
    fin_p = F.OUT / f"{sid}_igl_fixed_labels.tif"
    lab1 = tifffile.imread(fin_p) if fin_p.exists() else None
    tis = lab0[::P.DS, ::P.DS] > 0
    norm = []
    for ci in (0, 2, 1):
        lo, hi = np.percentile(raw[ci, ::P.DS, ::P.DS][tis], [1, 99.7])
        norm.append((ci, lo, hi))
    out = []
    for name, m in info.items():
        c = np.array(m["marker_grid"])
        native, n, to_crop = _frame(tr, c)
        chans = []
        for ci, lo, hi in norm:
            v = ndi.map_coordinates(raw[ci], native, order=1, mode="constant", cval=0)
            chans.append(np.clip((v - lo) / max(hi - lo, 1e-6), 0, 1))
        rgb = (np.stack(chans, -1) * 255).astype(np.uint8)
        l0 = ndi.map_coordinates(lab0, native, order=0, mode="constant", cval=0)
        rgb[l0 == 0] = (rgb[l0 == 0] * 0.35).astype(np.uint8)
        now = np.zeros((n, n, 4), np.uint8)
        now[_outline(np.isin(l0, (6, 8)), 2)] = (255, 140, 40, 255)
        now[_outline(l0 == 4, 2)] = (255, 236, 60, 255)
        layers = dict(now=P.uri(Image.fromarray(now, "RGBA"), "PNG", optimize=True))
        if lab1 is not None:
            l1 = ndi.map_coordinates(lab1, native, order=0, mode="constant", cval=0)
            last = np.zeros((n, n, 4), np.uint8)
            last[_outline(np.isin(l1, (6, 8)), 2)] = (255, 140, 40, 255)
            last[_outline(l1 == 4, 2)] = (90, 255, 150, 255)
            layers["last"] = P.uri(Image.fromarray(last, "RGBA"), "PNG", optimize=True)
        gq = native * F.PX / G
        cands = []
        for q in m["candidates"]:
            k = q["id"]
            band = arr[f"{name}_{k}_band"]
            bm = ndi.map_coordinates(band.astype(np.uint8), gq, order=0, mode="constant", cval=0) > 0
            ov = np.zeros((n, n, 4), np.uint8)
            ov[_outline(bm, 2)] = (*COLS[k], 255)
            C = np.array([to_crop(p) for p in arr[f"{name}_{k}_curve"]])
            ln = np.zeros((n, n), bool)
            for p0, p1 in zip(C[:-1], C[1:]):
                rr, cc = draw_line(int(round(p0[0])), int(round(p0[1])), int(round(p1[0])), int(round(p1[1])))
                ok = (rr >= 0) & (rr < n) & (cc >= 0) & (cc < n)
                ln[rr[ok], cc[ok]] = True
            ln = ndi.binary_dilation(ln, iterations=2)
            ov[ln] = (*COLS[k], 255)
            for p in (C[0], C[-1]):
                yy, xx = np.ogrid[:n, :n]
                ring = np.abs(np.hypot(yy - p[0], xx - p[1]) - 7) < 1.6
                ov[ring] = (*COLS[k], 255)
            layers[k] = P.uri(Image.fromarray(ov, "RGBA"), "PNG", optimize=True)
            cands.append(dict(id=k, length_um=round(q["length_um"]), text=f"{TARGET_TEXT[q['target']]}, "
                              f"{round(q['length_um'])} µm bridge, width {q['width_um'][0]}→{q['width_um'][1]} µm"))
        mk = F.MARKERS_ALL[sid][name]
        out.append(dict(name=name, label=LABELS.get((sid, name), name), img=P.uri(Image.fromarray(rgb, "RGB"), "JPEG", quality=88), n=n, layers=layers,
                        cands=cands, x=mk[0], y=mk[1]))
    return out
