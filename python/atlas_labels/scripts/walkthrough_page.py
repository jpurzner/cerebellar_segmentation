"""
Pipeline walkthrough (2026-11-24): the steps that make the reference line and the corrected labels, and at each
spot the user marked on flatten_check_v8 that is still open, what each step did there: close-ups of the stains,
the input labels, the stain calls, the cuts coloured by the step that drew them, the line with its columns, and
the flattened strip.

    python scripts/walkthrough_page.py      # -> outputs/flat_fix/walkthrough/ (index.html, img/)
"""
import json
import sys
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw
from scipy import ndimage as ndi

sys.path.insert(0, str(Path(__file__).parent))
import flatten as F  # noqa: E402
import egl_repair as ER  # noqa: E402

ROOT = F.ROOT
OUT = ROOT / "outputs/flat_fix/walkthrough"
FIX = ROOT / "outputs/flat_fix"
S = 3                                   # image pixels per 4 um grid cell
LAB = {0: (11, 14, 18), 1: (84, 92, 104), 2: (59, 224, 255), 3: (46, 74, 108), 4: (255, 214, 10), 5: (176, 124, 255), 6: (255, 70, 70)}
CUT = dict(drop=(150, 150, 150), stain=(255, 60, 200), ext=(255, 150, 30), way=(255, 255, 255), peak=(255, 235, 0), trim=(150, 120, 255))
# the open marks of the flatten v8 notes, on the flatten v8 line (page v9); i_v7 = the index the user marked on page v8
MARKS = [
    dict(key="s2_2-11", short="s2_2", i=5617, i_v7=5868, ctr=(175, 445), half=72),
    dict(key="s2_4-19", short="s2_4", i=9595, i_v7=9424, ctr=(358, 535), half=60),
    dict(key="s3_2-8", short="s3_2", i=4236, i_v7=4236, ctr=(76, 295), half=60),
    dict(key="s3_2-9", short="s3_2", i=4653, i_v7=4721, ctr=(88, 333), half=60),
    dict(key="s3_2-11", short="s3_2", i=5513, i_v7=5703, ctr=(155, 452), half=60),
]


def stretch(v, m, lo=1.0, hi=99.5):
    a, b = np.percentile(v[m], [lo, hi]) if m.any() else (0.0, 1.0)
    return np.clip((v - a) / max(b - a, 1e-6), 0, 1)


def win(a, r0, c0, h, fill=0):
    """a[r0-h:r0+h, c0-h:c0+h], padded with fill outside the grid."""
    out = np.full((2 * h, 2 * h) + a.shape[2:], fill, dtype=a.dtype)
    H, W = a.shape[:2]
    ra, rb, ca, cb = max(r0 - h, 0), min(r0 + h, H), max(c0 - h, 0), min(c0 + h, W)
    out[ra - (r0 - h):rb - (r0 - h), ca - (c0 - h):cb - (c0 - h)] = a[ra:rb, ca:cb]
    return out


def up(rgb):
    return np.repeat(np.repeat(rgb, S, 0), S, 1)


def paint(lab):
    rgb = np.zeros(lab.shape + (3,), np.uint8)
    for k, c in LAB.items():
        rgb[lab == k] = c
    return rgb


def ring(img, r, c, h):
    d = ImageDraw.Draw(img)
    y, x = (r - (h[0] - h[2])) * S + S / 2, (c - (h[1] - h[2])) * S + S / 2
    for w, col in ((5, (0, 0, 0)), (2, (255, 214, 10))):
        d.ellipse([x - 11, y - 11, x + 11, y + 11], outline=col, width=w)


def save(img, name):
    p = OUT / "img" / name
    img.save(p, optimize=True)
    return f"img/{name}"


def mark_panels(m, cache):
    s, i = m["short"], m["i"]
    if s not in cache:
        lab, tis = ER.repaired(s)
        _, best, _ = F.labels(s)
        E = np.load(ER.OUT / f"{s}.npz")
        Z = np.load(F.OUT2 / f"{s}.npz")
        R = np.load(FIX / f"{s}.npz")
        St = ER.stain_maps(s, lab, tis)
        Mm, Eg = lab == 3, lab == 2
        nm, pm = np.median(St["neun_a"][Mm]), np.median(St["p27a"][Mm])
        um, ue = np.median(St["nuc_a"][Mm]), np.median(St["nuc_a"][Eg])
        sig = tis & ((St["neun_a"] - nm) + (St["p27a"] - pm) < ER.PRM["arach_thr"]) & (St["nuc_a"] < um + ER.PRM["arach_bright"] * (ue - um))
        tp = ER.pia_tissue(s)
        mid_p, dropped = F.prune_cuts(lab, F.cut_tissue(tp, F.PRM), E["mid"], F.PRM, egl_s=E["egl_stain"])
        Rd, k, _ = F.stains_reduced(F.sid_of(s))
        cache[s] = dict(lab=lab, tis=tis, best=best, E=E, Z=Z, R=R, sig=sig, tp=tp, mid_p=mid_p, dropped=dropped, Rd=Rd)
    C = cache[s]
    lab, tis, E, Z, R = C["lab"], C["tis"], C["E"], C["Z"], C["R"]
    r0, c0 = m["ctr"]
    h = m["half"]
    gh, gw = lab.shape
    hh = (r0, c0, h)
    files = {}
    qr, qc = np.round(Z["Q"][i]).astype(int)
    # 1. stains: DAPI grey, NeuN green, p27 magenta
    Rd = C["Rd"]
    fy, fx = Rd["nuclei"].shape[0] / gh, Rd["nuclei"].shape[1] / gw
    sl = lambda a: a[max(int((r0 - h) * fy), 0):int((r0 + h) * fy), max(int((c0 - h) * fx), 0):int((c0 + h) * fx)]
    tw = win(tis, r0, c0, h)
    tz = np.asarray(Image.fromarray(tw.astype(np.uint8) * 255).resize(sl(Rd["nuclei"]).shape[::-1], Image.NEAREST)) > 0
    d_, n_, p_ = (stretch(sl(Rd[k2]).astype(np.float32), tz) for k2 in ("nuclei", "neun", "p27"))
    rgb = np.stack([0.62 * d_ + 0.75 * p_, 0.62 * d_ + 0.85 * n_, 0.62 * d_ + 0.75 * p_], -1)
    img = Image.fromarray((np.clip(rgb, 0, 1) * 255).astype(np.uint8)).resize((2 * h * S, 2 * h * S), Image.BILINEAR)
    ring(img, qr, qc, hh)
    files["stain"] = save(img, f"{m['key']}_1stain.png")
    # 2. input labels
    img = Image.fromarray(up(paint(win(C["best"], r0, c0, h))))
    ring(img, qr, qc, hh)
    files["input"] = save(img, f"{m['key']}_2input.png")
    # 3. stain calls on dim DAPI
    dg = (np.asarray(Image.fromarray((d_ * 255).astype(np.uint8)).resize((2 * h, 2 * h), Image.BILINEAR)) * 0.35).astype(np.uint8)
    rgb = np.repeat(dg[..., None], 3, -1)
    for key, col in (("egl_stain", (59, 224, 255)), ("igl_stain", (255, 214, 10))):
        rgb[win(E[key], r0, c0, h)] = col
    rgb[win(C["sig"], r0, c0, h)] = CUT["trim"]
    rgb[~tw] = LAB[0]
    img = Image.fromarray(up(rgb))
    ring(img, qr, qc, hh)
    files["calls"] = save(img, f"{m['key']}_3calls.png")
    # 4. cuts by the step that drew them, on the repaired labels (dimmed)
    rgb = (paint(win(lab, r0, c0, h)) * 0.33).astype(np.uint8)
    rgb[win(tis & ~C["tp"], r0, c0, h)] = CUT["trim"]
    fin = Z["mid"]
    dist = lambda a: ndi.distance_transform_edt(~a) if a.any() else np.full(a.shape, 1e9)
    d_er, d_ms, d_ext = dist(E["mid"]), dist(E["mid_stain"]), dist(E["mid_ext"])
    way = C["mid_p"] & (d_er > 2)
    stain_c = fin & ~way & ~Z["peak_add"] & (d_ms <= d_ext) & (d_ms <= 3)
    ext_c = fin & ~way & ~Z["peak_add"] & ~stain_c
    for msk, col in ((C["dropped"], CUT["drop"]), (stain_c, CUT["stain"]), (ext_c, CUT["ext"]), (fin & way, CUT["way"]), (Z["peak_add"] & fin, CUT["peak"])):
        rgb[win(msk, r0, c0, h)] = col
    img = Image.fromarray(up(rgb))
    ring(img, qr, qc, hh)
    files["cuts"] = save(img, f"{m['key']}_4cuts.png")
    near = np.zeros(lab.shape, bool)
    near[max(r0 - h, 0):r0 + h, max(c0 - h, 0):c0 + h] = True
    cut_stats = dict(stain=int((stain_c & near).sum()), ext=int((ext_c & near).sum()), way=int((fin & way & near).sum()),
                     drop=int((C["dropped"] & near).sum()), peak=int((Z["peak_add"] & fin & near).sum()), trim_um2=int((tis & ~C["tp"] & near).sum() * 16))
    # 5. line and columns on the corrected labels
    img = Image.fromarray(up(paint(win(R["lab_new"], r0, c0, h))))
    d = ImageDraw.Draw(img)
    L, valid, X, Q = R["L_new"], Z["valid"], Z["X"], Z["Q"]
    reach = (np.isin(L, (4, 5)) & valid).any(1)
    px = lambda p: ((p[..., 1] - (c0 - h)) * S + S / 2, (p[..., 0] - (r0 - h)) * S + S / 2)
    inside = (np.abs(Q[:, 0] - r0) < h) & (np.abs(Q[:, 1] - c0) < h)
    for j in np.flatnonzero(inside)[::5]:
        n = int(valid[j].sum())
        if n < 2:
            continue
        xs, ys = px(X[j, :n])
        d.line(list(zip(xs.tolist(), ys.tolist())), fill=(77, 255, 136) if reach[j] else (255, 60, 60), width=1)
    xs, ys = px(Q[inside])
    for x_, y_ in zip(xs, ys):
        d.rectangle([x_ - 1, y_ - 1, x_ + 1, y_ + 1], fill=(255, 255, 255))
    ring(img, qr, qc, hh)
    files["cols"] = save(img, f"{m['key']}_5cols.png")
    # 6. the strip: corrected labels and DAPI, +-240 um along the line, 240 um deep
    a, b = max(i - 120, 0), min(i + 120, len(L))
    Dk = 120
    lab_s = paint(L[a:b, :Dk].T)
    dapi = Z["S_nuclei"][a:b, :Dk].astype(np.float32).T
    vm = valid[a:b, :Dk].T
    dapi = (stretch(dapi, vm) * 255).astype(np.uint8)
    dapi[~vm] = 0
    bar = np.zeros((6, b - a, 3), np.uint8)
    bar[:] = (40, 46, 56)
    bar[:, R["geo"][a:b].astype(bool)] = (200, 208, 218)
    bar[:, R["reg"][a:b].astype(bool) & ~R["geo"][a:b].astype(bool)] = (77, 255, 136)
    gap = np.full((2, b - a, 3), 11, np.uint8)
    strip = np.concatenate([bar, gap, lab_s, gap, np.repeat(dapi[..., None], 3, -1)], 0)
    img = Image.fromarray(np.repeat(np.repeat(strip, 2, 0), 2, 1))
    d = ImageDraw.Draw(img)
    x = (i - a) * 2
    d.polygon([(x - 7, 0), (x + 7, 0), (x, 10)], fill=(255, 214, 10), outline=(0, 0, 0))
    files["strip"] = save(img, f"{m['key']}_6strip.png")
    lo, hi = max(i - 100, 0), min(i + 100, len(L))
    stats = dict(no_igl_um=int((~reach[lo:hi]).sum() * 2), geo_um=int(R["geo"][lo:hi].sum() * 2), cuts=cut_stats,
                 mm=round(i * 2 / 1000, 2), mm_v7=round(m["i_v7"] * 2 / 1000, 2), strip_mm=[round(a * 2 / 1000, 2), round(b * 2 / 1000, 2)],
                 igl_on_meninges_um2=int((win(C["best"] == 4, r0, c0, h) & win(C["sig"], r0, c0, h)).sum() * 16),
                 crop_um=2 * h * 4)
    return files, stats


def main():
    (OUT / "img").mkdir(parents=True, exist_ok=True)
    notes = json.loads((FIX / "notes_v8.json").read_text())
    cache, out = {}, []
    for m in MARKS:
        files, stats = mark_panels(m, cache)
        out.append(dict(key=m["key"], short=m["short"], note=notes[m["short"]].get("note", ""), files=files, stats=stats))
        print(m["key"], json.dumps(stats), flush=True)
    (OUT / "marks.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
