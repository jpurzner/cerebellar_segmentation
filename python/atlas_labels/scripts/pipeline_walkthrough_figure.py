"""
Stage-by-stage figures of the Nissl cerebellum pipeline (allen_nissl_template.py)
on one midline section per method:
  P7  (EGL-based: P4/P7)      -> outputs/walkthrough/walkthrough_P7.png
  P56 (IGL-based: P14/P56)    -> outputs/walkthrough/walkthrough_P56.png
Figures are 1800 x 900 px.
"""
import sys
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage as ndi
from skimage import filters, morphology
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
import allen_nissl_template as T  # noqa: E402

OUT = T.ROOT / "outputs/walkthrough"
OUT.mkdir(parents=True, exist_ok=True)
CLS_CMAP = matplotlib.colors.ListedColormap(
    [(0, 0.75, 0.9), (1, 0.85, 0.1), (1, 0.5, 0.1), (0.55, 0.9, 0.55)])


def seed_window(age, files):
    cfg = T.AGES[age]
    sd = cfg["seed"]
    sec = sd["sec"]
    rgb = np.asarray(Image.open(files[sec]).convert("RGB"))
    g = rgb.astype(float).mean(2)
    if sd["kind"] == "crop":
        cs = Image.open(T.SRC / f"cerebellum_crop_{sd['img_id']}.jpg").size[::-1]
        seed = np.zeros(g.shape, bool)
        r0, c0 = sd["offset"]
        seed[r0 + 40:r0 + cs[0] - 40, c0 + 40:c0 + cs[1] - 40] = True
        w = T.bbox_of(seed, 0, g.shape)
        env0 = T.envelope(ndi.gaussian_filter(g[T.sl(w)], 1.5), None)
    else:
        atlas = T.atlas_cb_mask(sd["atlas_image_id"], sec, g.shape, sd.get("ontology", "admba"))
        w = T.bbox_of(atlas, T.WIN, g.shape)
        env0 = T.envelope(ndi.gaussian_filter(g[T.sl(w)], 1.5),
                          ndi.binary_erosion(atlas[T.sl(w)], iterations=10))
    prior = np.zeros(g.shape, bool)
    prior[T.sl(w)] = env0
    w = T.bbox_of(prior, T.WIN, g.shape)
    return sec, rgb[T.sl(w)], g[T.sl(w)], ndi.binary_dilation(prior[T.sl(w)], iterations=10)


def contour(ax, m, color, lw=0.8, ls="-"):
    if m.any():
        ax.contour(m.astype(float), [0.5], colors=color, linewidths=lw, linestyles=ls)


def figure(age):
    T.configure(age)
    files = T.section_files()
    sec, rgb, gw, prior = seed_window(age, files)
    gs, g1 = ndi.gaussian_filter(gw, 1.5), ndi.gaussian_filter(gw, 1.0)

    if T.MODE == "adult":
        env, core, igl, thin, t = T.envelope_adult(gs, prior)
        pia, lab, fiss = T.segment_adult(g1, env, core, igl, t)
        hull_src = ndi.binary_dilation(core, iterations=int(T.ML_TYP_UM / T.PX_UM)) & env
        D = ndi.distance_transform_edt(~core) * T.PX_UM
        axis_all = morphology.skeletonize(env & ~core)
    else:
        env = T.envelope(gs, prior)
        pia, lab, fiss, t = T.segment_layers(g1, env)
        hull_src = lab == T.LBL["EGL"]
    keep, info = T.base_cut(env, hull_src)
    final = lab.copy()
    if "skipped" not in info:
        final[~keep] = 0
    hull = morphology.convex_hull_image(hull_src) if hull_src.any() else hull_src

    b = T.bbox_of(env, 30, env.shape)
    S = T.sl(b)
    cls = np.digitize(g1, filters.threshold_multiotsu(g1[env], 4)).astype(float)
    cls[~env] = np.nan

    fig, ax = plt.subplots(2, 3, figsize=(18, 9), dpi=100)
    ax = ax.ravel()
    ax[0].imshow(rgb[S])
    ax[0].set_title(f"1. Allen {age} Nissl, section {sec} (midline)")

    ax[1].imshow(cls[S], cmap=CLS_CMAP, vmin=0, vmax=3, interpolation="nearest")
    names = ("darkest", "dark", "pale", "palest")
    ax[1].set_title("2. 4 grey classes inside the cerebellum\n(" + " / ".join(names) + ")")

    ax[2].imshow(rgb[S])
    contour(ax[2], env[S], "red", 1.0)
    if T.MODE == "adult":
        contour(ax[2], igl[S], "yellow", 0.6)
        contour(ax[2], (core & ~igl)[S], "orange", 0.6)
        ax[2].set_title("3. envelope (red) = tissue within ML depth of the\n"
                        "dark IGL ribbon (yellow); WM inside IGL (orange)")
    else:
        contour(ax[2], (lab == T.LBL["EGL"])[S], "cyan", 0.6)
        ax[2].set_title("3. envelope (red) = dark-cell mass, closed + filled;\n"
                        "EGL network touching the surface (cyan)")

    ax[3].imshow(gw[S], cmap="gray")
    if T.MODE == "adult":
        yy, xx = np.nonzero(axis_all[S])
        ax[3].plot(xx, yy, ",", color="deepskyblue")
        ax[3].set_title("4. fissures: ML medial axis (blue); kept where both\n"
                        "sides are IGL, not the surface (red)")
    else:
        ax[3].set_title("4. fissures: midline of EGL lying > 50 um below the\n"
                        "envelope (fused EGL pairs), joined to the surface (red)")
    yy, xx = np.nonzero(fiss[S])
    ax[3].plot(xx, yy, ",", color="red")

    ax[4].imshow(rgb[S])
    contour(ax[4], hull[S], "white", 0.8, "--")
    contour(ax[4], env[S], "red", 0.8)
    if "skipped" not in info:
        A = np.array(info["A_px"]) - [b[0], b[2]]
        B = np.array(info["B_px"]) - [b[0], b[2]]
        ax[4].plot([A[1], B[1]], [A[0], B[0]], "-", color="magenta", lw=2)
        rem = (env & ~keep)[S]
        ax[4].imshow(np.ma.masked_where(~rem, rem), cmap="cool", alpha=0.5)
        ax[4].set_title(f"5. brainstem cut: convex hull (white) of {'IGL+ML' if T.MODE == 'adult' else 'EGL'};\n"
                        f"base edge = lobule I-X line (magenta), {info['removed_mm2']:.2f} mm2 removed")
    else:
        ax[4].set_title("5. brainstem cut: " + info["skipped"])

    ax[5].imshow(final[S], cmap=T.CMAP, vmin=0, vmax=6, interpolation="nearest")
    ax[5].set_title("6. labels: cyan EGL, green ML, yellow IGL, orange WM")
    for a in ax:
        a.set_xticks([]); a.set_yticks([])
    fig.suptitle(f"{age}: how the outline and layers are extracted (anterior left, dorsal up)")
    fig.tight_layout()
    p = OUT / f"walkthrough_{age}.png"
    fig.savefig(p)
    plt.close(fig)
    print(p, Image.open(p).size)


if __name__ == "__main__":
    for age in (sys.argv[1:] or ["P7", "P56"]):
        figure(age)
