"""
Smoothed approximations of real IGLs: the IGL as a layer at a (smooth) depth below the EGL.

    python scripts/igl_synth_laminar.py     # -> outputs/igl_synth/laminar.json, laminar_examples.json

For every real section (IF controls, IF other controls, Allen P7 within 600 um), with dE =
distance to the EGL and C = cerebellum minus EGL:
  wm_const    IGL = C, not white matter, dE >= m            (1 parameter: ML thickness m)
  band_const  IGL = C, a <= dE <= b                         (2 parameters: the IGL band; no white matter used)
  wm_smooth   IGL = C, not white matter, dE >= a_s(x)
  band_smooth IGL = C, a_s(x) <= dE <= b_s(x)
              with mu_s, sd_s = local mean and standard deviation of the IGL's depth below the EGL,
              smoothed over s um (normalised Gaussian convolution over IGL pixels), and
              a_s = mu_s - sqrt(3) sd_s, b_s = mu_s + sqrt(3) sd_s (the exact edges of a band
              filled uniformly in depth); where no IGL is near, the section's constant fit
Smoothing scales s = 50, 100, 200, 400, 800 um. Score = Dice with the real IGL.
The constant parameters are fitted per section by grid search; their distribution across
sections is itself a result (what ML thickness and IGL band the data imply).
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi

sys.path.insert(0, str(Path(__file__).parent))
import igl_synth_stats as S  # noqa: E402

G = S.G
OUT = S.OUT
SCALES = [50, 100, 200, 400, 800]
EXAMPLES = ["s5_G", "Allen P7 +0"]


def dice(a, b):
    return float(2 * (a & b).sum() / max(a.sum() + b.sum(), 1))


def smooth_field(values, mask, sigma_px, fill=np.nan, min_w=1e-3):
    num = ndi.gaussian_filter(np.where(mask, values, 0.0), sigma_px)
    den = ndi.gaussian_filter(mask.astype(float), sigma_px)
    out = np.where(den > min_w, num / np.maximum(den, 1e-9), fill)
    return out, den


def band_fields(dE, igl, sigma_px, res):
    m1, den = smooth_field(dE, igl, sigma_px, fill=np.nan)
    m2, _ = smooth_field(dE ** 2, igl, sigma_px, fill=np.nan)
    sd = np.sqrt(np.maximum(m2 - m1 ** 2, 0))
    a = m1 - np.sqrt(3) * sd
    b = m1 + np.sqrt(3) * sd
    a = np.where(np.isnan(a), res["band_const"]["a"], a)
    b = np.where(np.isnan(b), res["band_const"]["b"], b)
    return a, b


def fit(g):
    dE = ndi.distance_transform_edt(~g["egl"]) * G
    C = g["tissue"] & ~g["egl"]
    igl = g["igl"] & C
    res = {}
    # wm_const: grid search on m
    best = (0, None)
    for m in np.arange(10, 160, 2):
        d = dice(C & ~g["wm"] & (dE >= m), igl)
        if d > best[0]:
            best = (d, m)
    res["wm_const"] = dict(dice=round(best[0], 4), m=float(best[1]))
    # band_const: grid search on a, b
    best = (0, None, None)
    for a in np.arange(10, 140, 4):
        for b in np.arange(a + 20, 400, 8):
            d = dice(C & (dE >= a) & (dE <= b), igl)
            if d > best[0]:
                best = (d, a, b)
    res["band_const"] = dict(dice=round(best[0], 4), a=float(best[1]), b=float(best[2]))
    for sc in SCALES:
        a_s, b_s = band_fields(dE, igl, sc / G, res)
        res[f"wm_smooth_{sc}"] = dict(dice=round(dice(C & ~g["wm"] & (dE >= a_s), igl), 4))
        res[f"band_smooth_{sc}"] = dict(dice=round(dice(C & (dE >= a_s) & (dE <= b_s), igl), 4))
    return res, dict(dE=dE, C=C)


def example_masks(g, res):
    """the real IGL and its approximations, as SVG paths (2 um per unit, like the explainer)."""
    import igl_model_explainer as X
    dE = ndi.distance_transform_edt(~g["egl"]) * G
    C = g["tissue"] & ~g["egl"]
    igl = g["igl"] & C
    up = lambda m: ndi.zoom(m.astype(np.uint8), 2, order=0) > 0
    out = dict(n=[int(g["tissue"].shape[0] * 2), int(g["tissue"].shape[1] * 2)],
               tissue=X.paths(up(g["tissue"])), egl=X.paths(up(g["egl"])), wm=X.paths(up(g["wm"])), igl=X.paths(up(igl)), approx={})
    out["approx"]["band_const"] = X.paths(up(C & (dE >= res["band_const"]["a"]) & (dE <= res["band_const"]["b"])))
    out["approx"]["wm_const"] = X.paths(up(C & ~g["wm"] & (dE >= res["wm_const"]["m"])))
    for sc in SCALES:
        a_s, b_s = band_fields(dE, igl, sc / G, res)
        out["approx"][f"band_smooth_{sc}"] = X.paths(up(C & (dE >= a_s) & (dE <= b_s)))
        out["approx"][f"wm_smooth_{sc}"] = X.paths(up(C & ~g["wm"] & (dE >= a_s)))
    return out


def main():
    rows, ex = [], {}
    for group, name, load in S.real_sets():
        if group == "IF EZH2 cKO":
            continue
        g = load()
        res, _ = fit(g)
        rows.append(dict(group=group, name=name, **{k: v for k, v in res.items()}))
        print(f"{group:10s} {name:16s} wm_const {res['wm_const']['dice']:.3f} (m={res['wm_const']['m']:.0f})  band_const {res['band_const']['dice']:.3f} "
              f"(a={res['band_const']['a']:.0f}, b={res['band_const']['b']:.0f})  band_smooth " +
              " ".join(f"{s}:{res[f'band_smooth_{s}']['dice']:.3f}" for s in SCALES), flush=True)
        if name in EXAMPLES:
            ex[name] = example_masks(g, res)
    (OUT / "laminar.json").write_text(json.dumps(rows, indent=1))
    (OUT / "laminar_examples.json").write_text(json.dumps(ex))
    for grp in ("IF control", "IF other", "Allen P7"):
        R = [r for r in rows if r["group"] == grp]
        keys = ["wm_const", "band_const"] + [f"wm_smooth_{s}" for s in SCALES] + [f"band_smooth_{s}" for s in SCALES]
        print(grp, {k: round(float(np.median([r[k]["dice"] for r in R])), 3) for k in keys},
              "m", np.percentile([r["wm_const"]["m"] for r in R], [10, 50, 90]).round(0).tolist(),
              "a", np.percentile([r["band_const"]["a"] for r in R], [10, 50, 90]).round(0).tolist(),
              "b", np.percentile([r["band_const"]["b"] for r in R], [10, 50, 90]).round(0).tolist())


if __name__ == "__main__":
    main()
