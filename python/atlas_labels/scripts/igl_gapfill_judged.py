"""
Freeze the user's verdicts on the IGL Gap Fill page (db collection "gapfill_real", 2026-11-19)
together with the fill each verdict was about, so later filler versions can be scored against
them: outputs/igl_gapfill/judged_v1.npz (one full-grid mask per item) + judged_v1.json.

    python scripts/igl_gapfill_judged.py verdicts_v1.json [verdicts_v2.json]   # {id: {verdict, note}}
    (with the second file: + the round-2 bridges -> judged_v2.*)

Items: the 16 bridges of the first real run (apply_prelayout.json; s1_5-3 and s2_2-1 were shown
as refused ref-s1_5-1 / ref-s2_2-1), and the fills / marker joins refused at s4_A / s4_C
(outputs/igl_fixed/<slide>_layout_refused.npz, grouped as on the page).
Verdict meaning: correct / join = the IGL should be there; wrong / refuse = it should not.
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import igl_gapfill_page as PG  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/igl_gapfill"


def _current_fill(sid, chord_um):
    J = json.loads((PG.FIXED / "igl_fixed.json").read_text())
    if sid in J:
        BZ = np.load(PG.FIXED / f"{sid}_bridges.npz")
        for i, b in enumerate(J[sid]["bridges"]):
            if round(b["chord_um"]) == chord_um:
                return BZ[f"fill_{i}"]
        return None
    A = json.loads((OUT / "apply.json").read_text())
    Z = np.load(OUT / f"{sid}_gapfill.npz")
    for b in A[sid]["candidates"]:
        if b["status"] == "joined" and round(b["chord_um"]) == chord_um:
            y0, y1, x0, x1 = [int(v) for v in Z[f"box_{b['id']}"]]
            f = np.zeros(Z["code"].shape, bool)
            f[y0:y1, x0:x1] = Z[f"fill_{b['id']}"]
            return f
    return None


def main(vfile, round2=None):
    V = json.loads(Path(vfile).read_text())
    A = json.loads((OUT / "apply_prelayout.json").read_text())
    J = json.loads((PG.FIXED / "igl_fixed.json").read_text())
    masks, meta = {}, []

    def add(key, sid, f, src, what):
        v = V.get(key, {})
        good = {"correct": True, "join": True, "wrong": False, "refuse": False}.get(v.get("verdict"))
        masks[f"k{len(meta)}"] = f
        meta.append(dict(key=f"k{len(meta)}", id=key, sid=sid, source=src, what=what, verdict=v.get("verdict"), igl_there=good, note=v.get("note")))

    if round2:                            # round 1 is frozen already: start from it
        meta = json.loads((OUT / "judged_v1.json").read_text())
        M1 = np.load(OUT / "judged_v1.npz")
        masks = {m["key"]: M1[m["key"]] for m in meta}
        A = {}
    for sid, e in A.items():
        Z = np.load(OUT / f"{sid}_gapfill.npz")
        jo = [b for b in e["candidates"] if b["status"] == "joined"]
        for n, b in enumerate(jo, 1):
            k = b["id"]
            y0, y1, x0, x1 = [int(v) for v in Z[f"box_{k}"]]
            f = np.zeros(Z["code"].shape, bool)
            f[y0:y1, x0:x1] = Z[f"fill_{k}"]
            key = f"{e['short']}-{n}"
            key = {"s1_5-3": "ref-s1_5-1", "s2_2-1": "ref-s2_2-1"}.get(key, key)
            add(key, sid, f, "first real run", f"{b['chord_um']:.0f} um bridge ({b['found_by']}, {b['target']})")
    for slide in (() if round2 else ("s4_A_p27", "s4_C_p27")):
        R = np.load(PG.FIXED / f"{slide}_layout_refused.npz")
        code = np.load(OUT / f"{slide}_gapfill.npz")["code"]
        groups = []
        for it in J[slide]["refused_items"]:
            f = R[it["key"]]
            c = PG.cut_mask(code, f)
            for gr in groups:
                a = gr[0][1]
                if a.any() and c.any() and (a & c).sum() / (a | c).sum() >= 0.7:
                    gr.append((f, c, it["what"]))
                    break
            else:
                groups.append([(f, c, it["what"])])
        for n, gr in enumerate(groups, 1):
            for f, c, what in gr:
                add(f"ref-{slide[:4]}-{n}", slide, f, "s4 run with bleach boxes and markers", what)
    if round2:
        # second review (page of round 2, outputs/igl_gapfill/igl_gap_fill_round2.html): bridges
        # p-<section>-<n>; their fills from the current outputs, matched by section and length
        V2 = json.loads(Path(round2).read_text())
        h = (OUT / "igl_gap_fill_round2.html").read_text()
        D = json.loads(h[h.index("const DATA = ") + 13:h.index(";\nconst ALL")])
        for it in D["items"]:
            v = V2.get(it["id"])
            if not v:
                continue
            f = _current_fill(it["sid"], it["chord_um"])
            if f is not None:
                V[it["id"]] = v
                add(it["id"], it["sid"], f, "s4 run, round 2" if it["sid"] in J else "round 2", f"{it['chord_um']} um bridge ({it['found_by']})")
    tag = "v2" if round2 else "v1"
    np.savez_compressed(OUT / f"judged_{tag}.npz", **masks)
    (OUT / f"judged_{tag}.json").write_text(json.dumps(meta, indent=1))
    for m in meta:
        print(f"{m['id']:12s} {str(m['verdict']):8s} {m['what']}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
