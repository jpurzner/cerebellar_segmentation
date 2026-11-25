"""
IGL and reference line (2026-11-21): per section, the IGL being called (current labels: reviewed
controls, the user's corrected cKO sections, MATLAB for s1_5) and the reference line the flattening
measures depth from (outputs/flatten/<short>.npz), coloured by what is under it:
  - EGL labelled under it (within 20 um, past any tear / unassigned tissue);
  - no EGL label, but the column below reaches IGL: the cortex is there, so the EGL is probably
    there too and unlabelled;
  - no EGL label and no IGL below: cut edge, base, peduncle, DCN.
Plus the fissure cuts, and optionally the EGL labels and every 50th um column.

    python scripts/igl_reference_page.py      # -> outputs/igl_reference/page/
"""
import json
import sys
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage as ndi

sys.path.insert(0, str(Path(__file__).parent))
import flatten as F  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
PAGE = ROOT / "outputs/igl_reference/page"
TEMPLATE = Path(__file__).parent / "igl_reference_template.html"


def stretch(v, mask, lo=1.0, hi=99.5):
    a, b = np.percentile(v[mask], [lo, hi])
    return np.clip((v - a) / max(b - a, 1e-6) * 255, 0, 255).astype(np.uint8)


def overlay(mask, rgb, fill_a, edge_a):
    edge = mask & ~ndi.binary_erosion(mask)
    out = np.zeros(mask.shape + (4,), np.uint8)
    out[mask] = (*rgb, fill_a)
    out[edge] = (*rgb, edge_a)
    return out


def egl_raw(L, valid, cover_um=20.0, skip_um=40.0):
    """per column, without any smoothing along the line: EGL label within cover_um of the surface (past
    leading tear / unassigned cells), and the thickness of that first EGL run (um)."""
    N, Dn = L.shape
    skip = np.isin(L, (1, 6)) & valid
    lead = np.minimum(np.argmin(skip, axis=1), int(skip_um / F.STEP))
    k = int(cover_um / F.STEP)
    rows = np.arange(N)[:, None]
    idx = np.minimum(lead[:, None] + np.arange(k)[None, :], Dn - 1)
    isE = (L == 2) & valid
    has = isE[rows, idx].any(1)
    thick = np.zeros(N)
    for i in np.flatnonzero(has):
        col = isE[i, lead[i]:]
        s0 = int(np.argmax(col))
        e0 = s0 + (int(np.argmin(col[s0:])) if not col[s0:].all() else len(col) - s0)
        thick[i] = (e0 - s0) * F.STEP
    return has, thick


def build(short, row):
    Z = np.load(F.OUT / f"{short}.npz")
    auto, best, tissue = F.labels(short)
    gh, gw = tissue.shape
    sid = F.sid_of(short)
    R, k, _ = F.stains_reduced(sid)
    nuc = ndi.zoom(R["nuclei"], (2 * gh / R["nuclei"].shape[0], 2 * gw / R["nuclei"].shape[1]), order=1)
    t2 = np.repeat(np.repeat(tissue, 2, 0), 2, 1)[:nuc.shape[0], :nuc.shape[1]]
    base = (stretch(nuc, t2) * 0.8).astype(np.uint8)
    files = {}
    Image.fromarray(base).save(PAGE / "img" / f"{short}_base.jpg", quality=80, optimize=True)
    files["base"] = f"img/{short}_base.jpg"
    for name, m, rgb, fa, ea in (("igl", best == 4, (255, 214, 10), 120, 255), ("egl", best == 2, (59, 224, 255), 110, 230),
                                 ("cut", Z["mid"], (255, 60, 200), 255, 255)):
        Image.fromarray(overlay(m, rgb, fa, ea), "RGBA").save(PAGE / "img" / f"{short}_{name}.png", optimize=True)
        files[name] = f"img/{short}_{name}.png"
    Q, X, valid, ok, Lb = Z["Q"], Z["X"], Z["valid"], Z["ok"], Z["Lb"]
    N = len(Q)
    igl_below = ((Lb == 4) & valid).any(1)
    raw, thick = egl_raw(Lb, valid)
    cls = np.where(raw, 0, np.where(igl_below, 1, 2)).astype(np.int8)
    piece = Z["piece"]
    um = lambda m: int(m.sum() * F.STEP)
    stats = dict(length_um=int(N * F.STEP), egl_um=um(cls == 0), no_egl_igl_um=um(cls == 1), no_egl_no_igl_um=um(cls == 2),
                 egl_thin_um=um(raw & (thick < 12)), egl_thick_med_um=float(np.median(thick[raw])) if raw.any() else 0.0,
                 egl_gaps=int(sum(1 for a, b in F._runs(~raw & igl_below) if (b - a) * F.STEP >= 20)),
                 igl_um2=int((best == 4).sum() * F.G ** 2), cut_um2=int(Z["mid"].sum() * F.G ** 2),
                 pieces=int(piece.max() + 1))
    step = 2
    line = []
    for i in range(0, N, step):
        line += [int(round(Q[i, 0] * F.G)), int(round(Q[i, 1] * F.G)), int(cls[i]), int(piece[i])]
    cols = []
    for i in range(0, N, 25):
        L = int(valid[i].sum())
        pts = np.vstack([Q[i:i + 1], X[i, :L:5]]) * F.G
        cols.append([int(round(v)) for v in pts.ravel()])
    data = dict(short=short, group=row["group"], base=[int(base.shape[0]), int(base.shape[1])], grid=[int(gh), int(gw)],
                line=line, line_step=step, cols=cols, files=files, stats=stats)
    (PAGE / "data" / f"{short}.json").write_text(json.dumps(data, separators=(",", ":")))
    return dict(short=short, group=row["group"], data=f"data/{short}.json", stats=stats)


def findings(secs):
    def rng(grp, key, frac=True):
        v = [s["stats"][key] / (s["stats"]["length_um"] if frac else 1000) for s in secs if s["group"] == grp]
        return f"{round(min(v) * 100)}–{round(max(v) * 100)}%" if frac else f"{min(v):.1f}–{max(v):.1f} mm"
    return [
        "The line itself follows the tissue edge, which does not need the EGL. The EGL labels come in at three places: "
        "fissure midlines are only found where two labelled EGL sheets meet; the cuts are routed through labelled EGL; "
        "and stretches with no EGL label under them are greyed out and left out of the IGL fit.",
        f"<b>Controls:</b> an EGL label within 20 µm of the line along {rng('control', 'egl_um')} of it, counted point by point; "
        f"{rng('control', 'no_egl_igl_um', False)} per section has no EGL label although the column below reaches IGL, and "
        f"{rng('control', 'egl_thin_um', False)} has an EGL label thinner than 12 µm. <b>cKO:</b> {rng('cKO', 'egl_um')} labelled, "
        f"{rng('cKO', 'no_egl_igl_um', False)} with IGL below and no EGL label, {rng('cKO', 'egl_thin_um', False)} thinner than 12 µm.",
        "Not in those numbers: a closed fissure with no EGL label in it is never cut, so the line jumps across its mouth and "
        "the IGL lining it sits deep inside the tissue with no line near it. Look for IGL bands facing each other with no "
        "magenta line between them.",
        "The repair this sets up, from the IGL outward: where IGL lies under the surface, the EGL band should be there between the "
        "surface and the ML; where two IGL bands face each other with no line between them, there is a closed fissure, and its "
        "midline lies halfway between the two IGL bands, past the ML on each side.",
    ]


def main():
    (PAGE / "img").mkdir(parents=True, exist_ok=True)
    (PAGE / "data").mkdir(parents=True, exist_ok=True)
    summ = json.loads((F.OUT / "summary.json").read_text())
    rows = {r["section"]: r for r in summ["sections"]}
    order = ["s1_2"] + sorted(s for s, r in rows.items() if r["group"] == "control" and s != "s1_2") + \
        sorted(s for s, r in rows.items() if r["group"] == "cKO")
    secs = []
    for s in order:
        secs.append(build(s, rows[s]))
        st = secs[-1]["stats"]
        print(f"{s:5s} {rows[s]['group']:7s} line {st['length_um'] / 1000:5.1f} mm | EGL under it {st['egl_um'] / 1000:5.1f} | "
              f"no EGL, IGL below {st['no_egl_igl_um'] / 1000:4.1f} | no EGL, no IGL {st['no_egl_no_igl_um'] / 1000:4.1f} mm", flush=True)
    (PAGE / "index.html").write_text(TEMPLATE.read_text().replace(
        "__DATA__", json.dumps(dict(collection="igl_reference_check", sections=secs, findings=findings(secs)), separators=(",", ":"))))
    tot = sum(f.stat().st_size for f in PAGE.rglob("*") if f.is_file())
    print(f"page: {len(list(PAGE.rglob('*.*')))} files, {tot / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
