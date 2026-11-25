"""
IGL Anomaly Check page: the whole-IGL model's flags on the 10x controls it was not fitted on.

    python scripts/igl_anomaly_page.py

Per section: the IF image in the page frame (anterior left, dorsal up), the flagged
stretches of IGL ribbon (numbered F1, F2, ...; magenta = shape, cyan = signal), and breaks
(K1, K2, ...): ribbon ends whose cheapest continuation (gap model) is
  likely   - cheap without crossing white matter (score <= break threshold)
  possible - cheap only when crossing a short white-matter-labelled patch (score <= 10),
             the signature of bleached IGL labelled white matter, but also seen at 2 of 15
             natural ends of the intact controls
The user marks each Real problem / Fine and can drop markers on anything missed.
-> outputs/igl_whole/igl_anomaly.html
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import igl_check_page as P  # noqa: E402
import igl_whole_model as W  # noqa: E402
import wm_library as WL  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = W.OUT
SIDS = ["s4_A_p27", "s4_C_p27", "s4_F_p27", "s5_A_p27", "s5_B_p27", "s5_C_p27", "s5_E_p27"]
SIGNAL = {"neun", "dapi", "p27", "dneun", "ddapi", "dp27"}
THR_WM = 10.0
PX = 0.65


def to_disp(tr, q):
    A, off, crop = np.array(tr["A"]), np.array(tr["offset"]), np.array(tr["crop"], float)
    q = np.atleast_2d(np.asarray(q, float))
    d = np.linalg.solve(A, (q * W.G / (PX * P.DS) - off).T).T - crop
    h, w = tr["shape"]
    return [[round(float(p[1] / w), 5), round(float(p[0] / h), 5)] for p in d]


def main():
    lib = {e["id"]: e for e in json.loads((WL.OUT / "library.json").read_text())["entries"]}
    res = json.loads((OUT / "flags.json").read_text())
    brk = json.loads((OUT / "breaks.json").read_text())
    thr = json.loads((OUT / "break_threshold.json").read_text())
    items = []
    for sid in SIDS:
        n = sid.replace("_p27", "")
        it, tr, _ = P.build(sid, lib)
        Pr = W.load_profile(n)
        fl = []
        for k, f in enumerate(sorted(res["sections"][n]["flags"], key=lambda f: -f["level"]), 1):
            pts = Pr["pos"][f["i0"]:f["i1"] + 1]
            fl.append(dict(id=f"F{k}", kind="signal" if f["why"][0]["feature"] in SIGNAL else "shape", type=f["type"],
                           length_um=f["length_um"], level=f["level"], pts=to_disp(tr, pts[::2]),
                           also=[w["type"] for w in f["why"][1:] if w["type"] != f["type"]][:2]))
        bks, k = [], 0
        for e in brk[n]:
            b, bw = e["best"], e["best_wm"]
            if b and b["score"] <= thr["threshold"]:
                kind, use = "likely", b
            elif bw and bw["score"] <= THR_WM and bw["wm_um"] > 0:
                kind, use = "possible", bw
            else:
                continue
            k += 1
            bks.append(dict(id=f"K{k}", kind=kind, score=use["score"], length_um=use["length_um"], wm_um=use.get("wm_um", 0),
                            end=to_disp(tr, e["pos"])[0], curve=to_disp(tr, np.array(use["curve"])[::2])))
        s = res["sections"][n]
        items.append(dict(id=sid, w=it["w"], h=it["h"], img=it["img"], outline=it["outline"], labels=it["labels"],
                          folia=it["folia"], cell_px=it["cell_px"], flags=fl, breaks=bks,
                          flagged_pct=s["flagged_per_mm"], ribbon_mm=s["ribbon_mm"]))
        print(sid, len(fl), "flags,", len(bks), "breaks:", [(b["id"], b["kind"], b["score"], b["wm_um"]) for b in bks], flush=True)
    summary = dict(tau_local=dict(zip(res["local_features"], res["tau_local"])), scale=dict(zip(res["local_features"], res["local_scale"])),
                   sections={k: dict(genotype=v["genotype"], training=v["training"], flagged_pct=v["flagged_per_mm"], n_flags=v["n_flags"])
                             for k, v in res["sections"].items()},
                   synthetic=json.loads((OUT / "synthetic_eval.json").read_text())["summary"],
                   breaks=dict(threshold=thr["threshold"], found=thr["gap_ends_found"], natural=thr["natural_ends_flagged"],
                               no_end=thr["gap_edges_without_an_end"]))
    page = (Path(__file__).parent / "igl_anomaly_template.html").read_text()
    page = page.replace("__DATA__", json.dumps(dict(items=items, summary=summary)))
    p = OUT / "igl_anomaly.html"
    p.write_text(page)
    print(p, f"{p.stat().st_size / 1e6:.2f} MB")


if __name__ == "__main__":
    main()
