"""
Apply the IGL ribbon model (igl_model.complete) at the breaks the user marked on s4_A and
s4_C, and keep the three best distinct predictions per marker for the user to choose from.

    python scripts/igl_model_apply.py

Per marker (igl_finalize.MARKERS_ALL), on the anatomy-fixed labels (4 um grid):
  starts  - ribbon ends within 300 um of the marker (normal and light pruning)
  targets - other ribbon ends within 500 um, and points along the ribbon skeleton within
            350 um of the marker, entered from either direction; a target is skipped when
            the IGL already joins it to the start by a route < 3x the straight distance
            (that would be the start's own ribbon looping back)
  forbidden tissue - EGL, background, tears, and white matter outside the confirmed regions
            (bleach box + 100 um, marker + 150 um: bleach is often labelled white matter)
  each start-target pair is completed by the model (elastic curve from the smooth and the
  depth-route starts, lower energy wins); kept if it passes within 100 um of the marker and
  is at most 600 um long; ranked by model energy + length / 300 um (shorter breaks are
  more likely) + (distance from the marker / 50 um)^2 (the user marked the break); the best
  three that differ by > 20 um are kept as A, B, C.
Band painted for each: centre line +/- half the width interpolated between the two ends,
over ML / PCL / unassigned / anomaly (not tears) / white matter in confirmed regions.
-> outputs/igl_model/<sid>_marker_candidates.npz, marker_candidates.json
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import tifffile
from scipy import ndimage as ndi
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).parent))
import igl_finalize as F  # noqa: E402
import igl_check_page as P  # noqa: E402
import igl_model as M  # noqa: E402
import igl_model_eval as ME  # noqa: E402
import igl_repair_eval as EV  # noqa: E402
import igl_ribbon as RB  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/igl_model"
G = RB.GRID
R_START, R_END, R_TGT, MAX_LEN, NEAR_MARK = 300.0, 500.0, 350.0, 600.0, 100.0
LEN_SCALE, MARK_SCALE = 300.0, 50.0
PRM = json.loads((OUT / "benchmark.json").read_text())["summary"]["best_params"]


def route_map(igl, a, c, r_px):
    from skimage.graph import MCP_Geometric
    H, W = igl.shape
    y0, x0 = max(int(c[0]) - r_px, 0), max(int(c[1]) - r_px, 0)
    y1, x1 = min(int(c[0]) + r_px + 1, H), min(int(c[1]) + r_px + 1, W)
    cost = np.where(igl[y0:y1, x0:x1], 1.0, np.inf)
    ia = (int(round(a[0])) - y0, int(round(a[1])) - x0)
    out = np.full(igl.shape, np.inf)
    if not (0 <= ia[0] < cost.shape[0] and 0 <= ia[1] < cost.shape[1]) or not np.isfinite(cost[ia]):
        # start just outside the cleaned IGL: enter it at the nearest IGL pixel
        yy, xx = np.nonzero(np.isfinite(cost))
        if not len(yy):
            return out
        k = np.argmin((yy - ia[0]) ** 2 + (xx - ia[1]) ** 2)
        ia = (int(yy[k]), int(xx[k]))
    cum, _ = MCP_Geometric(cost).find_costs([ia])
    out[y0:y1, x0:x1] = cum
    return out


def skeleton_targets(rib, dE, width):
    pts = []
    for b in rib["branches"]:
        P_ = b["path"].astype(float)
        if len(P_) < 8:
            continue
        S = np.stack([ndi.gaussian_filter1d(P_[:, k], 2.0, mode="nearest") for k in range(2)], 1)
        T = np.gradient(S, axis=0)
        T /= np.maximum(np.hypot(*T.T), 1e-9)[:, None]
        for j in range(3, len(S) - 3, 6):
            lo, hi = max(j - 12, 0), min(j + 13, len(P_))
            seg = b["path"][lo:hi]
            pts.append(dict(pos=S[j], t=T[j], depth=float(np.median(dE[seg[:, 0], seg[:, 1]])),
                            width=float(np.median(width[seg[:, 0], seg[:, 1]]))))
    return pts


def band(C, wa, wb, shape, paintable):
    return ME.band(C, max(wa, 25.0), max(wb, 25.0), shape, paintable)


def candidates_at(sec, c, name):
    rib, aux, dE, width = sec["rib"], sec["aux"], sec["dE"], sec["rib"]["width"]
    ends = list(rib["ends"]) + [e for e in aux["ends"] if all(np.hypot(*(e["pos"] - q["pos"])) * G > 20 for q in rib["ends"])]
    starts = [e for e in ends if np.hypot(*(e["pos"] - c)) * G <= R_START]
    tg_sk = [t for t in sec["sk"] if np.hypot(*(t["pos"] - c)) * G <= R_TGT]
    cands = []
    for e in starts:
        a = np.asarray(e["pos"], float)
        ta = np.asarray(e["dir"], float)
        route = route_map(sec["igl"], a, c, int((R_END + R_TGT) / G))
        targets = []
        for e2 in ends:
            if e2 is e:
                continue
            b = np.asarray(e2["pos"], float)
            targets.append(("end", b, -np.asarray(e2["dir"], float), e2["depth"], e2["width"]))
        for t in tg_sk:
            for sgn in (1, -1):
                targets.append(("ribbon", t["pos"], sgn * t["t"], t["depth"], t["width"]))
        for kind, b, tb_in, db, wb in targets:
            v = b - a
            dist = np.hypot(*v)
            if dist * G < 20 or dist * G > R_END:
                continue
            if v @ ta < -0.17 * dist or v @ tb_in < -0.17 * dist:      # must lie ahead of both ends
                continue
            iy, ix = np.clip(np.round(b).astype(int), 0, np.array(route.shape) - 1)
            if np.isfinite(route[iy, ix]) and route[iy, ix] < 3 * dist:
                continue
            Hc = M.hermite(a, ta, b, tb_in, 40)
            if np.hypot(*np.diff(Hc, axis=0).T).sum() * G > MAX_LEN:
                continue
            avoid = sec["igl"].copy()
            yy, xx = np.ogrid[:avoid.shape[0], :avoid.shape[1]]
            for p in (a, b):
                avoid &= (yy - p[0]) ** 2 + (xx - p[1]) ** 2 > 16
            C, info = M.complete(a, ta, b, tb_in, sec["field"], e["depth"], db, PRM, avoid)
            L = float(np.hypot(*np.diff(C, axis=0).T).sum() * G)
            dm = float(np.min(np.hypot(*(C - c).T)) * G)
            if L > MAX_LEN or dm > NEAR_MARK:
                continue
            cands.append(dict(curve=C, score=info["energy"] + L / LEN_SCALE + (dm / MARK_SCALE) ** 2, energy=info["energy"], length_um=L,
                              marker_dist_um=dm, target=kind, start_init=info["start"], a=a, b=b,
                              width=(float(e["width"]), float(wb)), depth=(float(e["depth"]), float(db))))
    cands.sort(key=lambda q: q["score"])
    keep = []
    for q in cands:
        if all(ME.curve_err(q["curve"], k["curve"])[0] > 20 for k in keep):
            keep.append(q)
        if len(keep) == 3:
            break
    for k, q in zip("ABC", keep):
        q["id"] = k
        q["band"] = band(q["curve"], *q["width"], sec["igl"].shape, sec["paintable"]) & ~sec["g"]["igl"]
    print(f"  {name}: {len(starts)} starts, {len(cands)} completions pass the marker -> "
          + ", ".join(f"{q['id']} {q['target']} {q['length_um']:.0f}um score {q['score']:.1f}" for q in keep), flush=True)
    return keep, len(starts), len(cands)


def section(sid, tr):
    native = tifffile.imread(EV.RE.FX / f"{sid}_fixed_labels.tif")
    g = RB.grid_labels(native, F.PX)
    _, tear = EV.images(sid.replace("_p27", ""), g["igl"].shape)
    box = np.load(ROOT / f"outputs/igl_ribbon/{sid}_ribbon_repair_grid.npz")["bleach"]
    conf = ndi.distance_transform_edt(~box) * G <= 100
    yy, xx = np.indices(g["igl"].shape)
    marks = {}
    for name, (x, y) in F.MARKERS_ALL[sid].items():
        c = F.disp_to_grid(tr, x, y)
        marks[name] = c
        conf |= ((yy - c[0]) ** 2 + (xx - c[1]) ** 2) * G ** 2 <= 150 ** 2
    rib = RB.describe(g)
    aux = RB.describe(g, {**RB.P, "prune_factor": 1.0})
    dE = rib["depth_egl"]
    anom_ok = g["anom"] & ~tear
    forbid = g["egl"] | ~(g["tissue"] | anom_ok) | (g["wm"] & ~conf) | tear
    paintable = (g["ml"] | g["unass"] | anom_ok | (g["wm"] & conf)) & ~g["egl"]
    sec = dict(g=g, rib=rib, aux=aux, dE=dE, igl=rib["igl"], field=M.Field(dE, forbid), paintable=paintable,
               sk=skeleton_targets(rib, dE, rib["width"]), conf=conf)
    return sec, marks


def main():
    trs = {t["sid"]: t for t in json.loads((P.OUT / "transforms.json").read_text())}
    summary = {}
    for sid in F.MARKERS_ALL:
        t0 = time.time()
        sec, marks = section(sid, trs[sid])
        print(sid, flush=True)
        arrays, info = {}, {}
        for name, c in marks.items():
            keep, ns, nc = candidates_at(sec, c, name)
            info[name] = dict(marker_grid=[round(float(v), 1) for v in c], starts=ns, completions=nc,
                              candidates=[{k: (round(v, 2) if isinstance(v, float) else v) for k, v in q.items()
                                           if k in ("id", "score", "energy", "length_um", "marker_dist_um", "target", "start_init")}
                                          | dict(width_um=[round(v) for v in q["width"]], depth_um=[round(v) for v in q["depth"]],
                                                 band_um2=round(float(q["band"].sum() * G ** 2)))
                                          for q in keep])
            for q in keep:
                arrays[f"{name}_{q['id']}_curve"] = q["curve"]
                arrays[f"{name}_{q['id']}_band"] = q["band"]
        arrays["conf"] = sec["conf"]
        np.savez_compressed(OUT / f"{sid}_marker_candidates.npz", **arrays)
        summary[sid] = info
        print(f"  {time.time() - t0:.0f}s", flush=True)
    (OUT / "marker_candidates.json").write_text(json.dumps(summary, indent=1, default=float))


if __name__ == "__main__":
    main()
