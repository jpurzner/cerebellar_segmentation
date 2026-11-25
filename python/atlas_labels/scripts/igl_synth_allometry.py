"""
How the layers scale with folium size in real sections (the rule a generator should follow).

    python scripts/igl_synth_allometry.py     # -> outputs/igl_synth/allometry.json

At every point of the cerebellum's medial axis (the centre line of each folium and lobule,
excluding the base cut), with h = the folium's half-width (distance to the pial surface):
  r_wm   white-matter radius = distance to the nearest non-white-matter pixel (0 if the axis
         point is not white matter)
  r_core distance to the nearest ML (within 150 um of the EGL) or EGL pixel: the radius
         of IGL + white matter
  IGL thickness = r_core - r_wm,   ML thickness = h - r_core
(h here is measured to the EGL: the medial axis is that of the tissue inside the EGL, so
closed fissures, where EGL meets EGL, still separate the folia.)
Binned by h (25 um bins), median and IQR, for the IF controls, IF other controls and Allen
P7 within 600 um of the midline.
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi
from skimage.morphology import medial_axis

sys.path.insert(0, str(Path(__file__).parent))
import igl_synth_stats as S  # noqa: E402

G = S.G
OUT = S.OUT
BINS = np.arange(25, 451, 25)


def measure(g):
    """medial axis of the tissue inside the EGL; half-width h = distance to the EGL (or to
    background where there is no EGL); ML counted only within 150 um of the EGL (the Allen
    'ML' class also takes pale gaps inside the white matter)."""
    inner = g["tissue"] & ~g["egl"]
    ax, dist = medial_axis(inner, return_distance=True)
    h = dist * G
    rows = np.nonzero(g["egl"])[0]
    base_row = rows.max() if len(rows) else g["tissue"].shape[0]
    yy = np.arange(g["tissue"].shape[0])[:, None] * np.ones((1, g["tissue"].shape[1]))
    dE = ndi.distance_transform_edt(~g["egl"]) * G
    ml_layer = g["ml"] & (dE <= 150)
    keep = ax & (yy < base_row - 25) & (dE <= h + 8)        # axis points whose nearest surface is EGL, not the base cut
    r_wm = ndi.distance_transform_edt(g["wm"]) * G
    r_core = ndi.distance_transform_edt(~(ml_layer | g["egl"] | ~g["tissue"])) * G
    sel = keep & (h >= BINS[0])
    hv, wv, cv = h[sel], r_wm[sel], np.minimum(r_core[sel], h[sel])
    return hv, wv, np.maximum(cv - wv, 0), np.maximum(hv - cv, 0)


def main():
    acc = {}
    for group, name, load in S.real_sets():
        if group == "IF EZH2 cKO":
            group = "IF EZH2 cKO"
        g = load()
        hv, wv, iv, ov = measure(g)
        a = acc.setdefault(group, [[], [], [], []])
        for lst, v in zip(a, (hv, wv, iv, ov)):
            lst.append(v)
        print(group, name, len(hv), flush=True)
    res = dict(bins=BINS.tolist(), groups={})
    for grp, (H, Wm, Ig, Om) in acc.items():
        H, Wm, Ig, Om = map(np.concatenate, (H, Wm, Ig, Om))
        out = dict(n=int(len(H)), wm=[], igl=[], outer=[], count=[], wm_present=[])
        for lo, hi in zip(BINS[:-1], BINS[1:]):
            s_ = (H >= lo) & (H < hi)
            out["count"].append(int(s_.sum()))
            for k, v in (("wm", Wm), ("igl", Ig), ("outer", Om)):
                out[k].append(np.percentile(v[s_], [25, 50, 75]).round(1).tolist() if s_.sum() >= 20 else None)
            out["wm_present"].append(round(float((Wm[s_] > 0).mean()), 3) if s_.sum() >= 20 else None)
        # simple fits on the bin medians: outer ~ const, igl ~ a + b h, wm ~ max(0, rho (h - h0))
        mid = (BINS[:-1] + BINS[1:]) / 2
        ok = [i for i in range(len(mid)) if out["igl"][i] is not None]
        if len(ok) >= 3:
            x = mid[ok]
            out["fit_igl"] = np.polyfit(x, [out["igl"][i][1] for i in ok], 1).round(3).tolist()
            out["fit_outer"] = np.polyfit(x, [out["outer"][i][1] for i in ok], 1).round(3).tolist()
            out["fit_wm"] = np.polyfit(x, [out["wm"][i][1] for i in ok], 1).round(3).tolist()
        res["groups"][grp] = out
        print(grp, "ML thickness medians", [v[1] if v else None for v in out["outer"]])
        print(grp, "IGL thickness medians", [v[1] if v else None for v in out["igl"]])
        print(grp, "WM radius medians", [v[1] if v else None for v in out["wm"]], "share with WM", out["wm_present"])
    (OUT / "allometry.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
