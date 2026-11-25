"""
Score a gap-filler setting against the user's frozen verdicts (igl_gapfill_judged.py):
for every judged fill, the share of it the filler paints now (only the part that no fill with the
opposite verdict at the same section shares, since routes from one end overlap). An item counts as filled when
that share is at least 0.3. Correct / join items should be filled, wrong / refuse ones not.
s4_A / s4_C items from the bleach-box run are scored with the user's confirmed regions.

    python scripts/igl_gapfill_score.py ['{"search": "band"}']     # optional prm overrides (JSON)
"""
import json
import sys
from pathlib import Path

import numpy as np
import tifffile

sys.path.insert(0, str(Path(__file__).parent))
import igl_finalize as F  # noqa: E402
import igl_gapfill as GF  # noqa: E402
import igl_gapfill_apply as GA  # noqa: E402
import igl_ribbon as RB  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/igl_gapfill"


def run(sid, prm, with_conf):
    if sid.endswith("_p27") and sid.startswith("s4_"):
        native = tifffile.imread(GA.RE.FX / f"{sid}_fixed_labels.tif")
    else:
        native = tifffile.imread(GA.RE.FX / f"{sid}_fixed_labels.tif")
    px = GA.info(sid)[0]
    g = RB.grid_labels(native, px)
    neun, tear = GA.images(sid, g["igl"].shape)
    conf = F.confirmed_regions(sid, g) if with_conf else None
    new, br, _ = GF.fill_gaps(g, prm=prm, neun=neun, tear=tear, bleach=conf)
    return new & ~g["igl"], br


def score(prm=None, verbose=True):
    tag = "v2" if (OUT / "judged_v2.json").exists() else "v1"
    meta = json.loads((OUT / f"judged_{tag}.json").read_text())
    Mk = np.load(OUT / f"judged_{tag}.npz")
    cache, res = {}, []
    for m in meta:
        conf = m["source"].startswith("s4 run")
        k = (m["sid"], conf)
        if k not in cache:
            cache[k] = run(m["sid"], prm, conf)
        paint, br = cache[k]
        f = Mk[m["key"]].copy()
        if m["igl_there"] is not None:          # only the part no opposite verdict shares
            for o in meta:
                if o["sid"] == m["sid"] and o["igl_there"] is not None and o["igl_there"] != m["igl_there"]:
                    f &= ~Mk[o["key"]]
        share = float((paint & f).sum() / max(f.sum(), 1))
        res.append(dict(m, share=round(share, 2), filled=share >= 0.3))
    # one row per page item (several fills may share a verdict)
    items = {}
    for r in res:
        it = items.setdefault(r["id"], dict(id=r["id"], verdict=r["verdict"], igl_there=r["igl_there"], filled=False, share=0.0))
        it["filled"] |= r["filled"]
        it["share"] = max(it["share"], r["share"])
    rows = list(items.values())
    judged = [r for r in rows if r["igl_there"] is not None]
    agree = sum(r["filled"] == r["igl_there"] for r in judged)
    other = {}
    for (sid, conf), (paint, br) in cache.items():
        jm = np.zeros(paint.shape, bool)
        for m in meta:
            if m["sid"] == sid:
                jm |= Mk[m["key"]]
        other[f"{sid[:30]}{' (conf)' if conf else ''}"] = [round(b["chord_um"]) for b in br if b["status"] == "joined" and "_mask" in b
                                                         and not _overlaps(b, jm)]
    if verbose:
        for r in rows:
            ok = "" if r["igl_there"] is None else ("agrees" if r["filled"] == r["igl_there"] else "DISAGREES")
            print(f"  {r['id']:12s} {str(r['verdict']):8s} painted {r['share']:.2f}  {ok}")
        print(f"agreement {agree}/{len(judged)}; bridges elsewhere (chord um):", {k: v for k, v in other.items() if v})
    return dict(agree=agree, n=len(judged), rows=rows, other=other)


def _overlaps(b, jm):
    y0, y1, x0, x1 = b["_box"]
    return bool((jm[y0:y1, x0:x1] & b["_mask"]).sum() >= 0.3 * b["_mask"].sum())


if __name__ == "__main__":
    prm = json.loads(sys.argv[1]) if len(sys.argv) > 1 else None
    res = score(prm)
    if len(sys.argv) > 2:                     # tag: save as outputs/igl_gapfill/score_<tag>.json
        (OUT / f"score_{sys.argv[2]}.json").write_text(json.dumps(dict(prm=prm, **res), indent=1, default=float))
