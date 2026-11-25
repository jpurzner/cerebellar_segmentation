"""
Figure for fit_p7_template.py:  outputs/p7_template_fit/<slide>/<slide>_fit.png (1800 x 1200)

  a IF section (p27 red, NeuN green, DAPI blue) + fitted template outline
  b automatic MATLAB labels (anomaly grey)     c gold standard
  d template after similarity fit              e template after B-spline
  f fit score vs template medial-lateral offset
"""
import json
import sys
from pathlib import Path

import numpy as np
import tifffile
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
import fit_p7_template as F  # noqa: E402

CM = matplotlib.colors.ListedColormap(
    [(1, 1, 1), (0.1, 0.1, 0.4), (0.4, 0.7, 0.95), (0, 0.75, 0.9), (1, 0.85, 0.1),
     (0.55, 0.9, 0.55), (1, 0.5, 0.1), (1, 1, 1), (1, 1, 1), (0.55, 0.55, 0.55)])


def composite(slide):
    raw = tifffile.imread(F.SANDBOX / f"{slide}_native.tif").astype(np.float32)
    C, H, W = raw.shape
    b = F.BLOCK
    h, w = H // b, W // b
    sm = raw[:, :h * b, :w * b].reshape(C, h, b, w, b).mean((2, 4))
    rgb = np.stack([np.clip(sm[i] / np.percentile(sm[i], 99.5), 0, 1) for i in range(3)], -1)
    return rgb


def main(slide="2018_05_22_s1_2_p27-0021"):
    d = F.OUT_ROOT / slide
    z = np.load(d / f"{slide}_p7template_grid.npz")
    rep = json.loads((d / f"{slide}_fit.json").read_text())
    um = rep["saved_variant_um"]
    bs = z[f"bspline_{um:.0f}um"]
    sim = z["similarity"]
    rgb = composite(slide)

    fig, ax = plt.subplots(2, 3, figsize=(18, 12), dpi=100)
    ax = ax.ravel()
    ax[0].imshow(rgb)
    ax[0].contour((sim > 0).astype(float), [0.5], colors="white", linewidths=0.8, linestyles="--")
    ax[0].contour((bs > 0).astype(float), [0.5], colors="red", linewidths=0.9)
    ax[0].contour((bs == 4).astype(float), [0.5], colors="yellow", linewidths=0.5)
    ax[0].set_title("a  IF section (p27 R, NeuN G, DAPI B)\n"
                    "template outline: similarity (white dashed), B-spline (red); IGL edge (yellow)")
    gd = rep.get("dice_vs_gold", {})

    def dtxt(k):
        x = gd.get(k)
        return "" if not x else "\nDice vs gold: " + "  ".join(f"{n} {v:.2f}" for n, v in x.items())

    ax[1].imshow(z["auto"], cmap=CM, vmin=0, vmax=9, interpolation="nearest")
    ax[1].set_title("b  automatic labels (fit input; grey = anomaly)" + dtxt("auto_segmentation"))
    if z["gold"].size:
        ax[2].imshow(z["gold"], cmap=CM, vmin=0, vmax=9, interpolation="nearest")
        ax[2].set_title("c  gold standard (IGL/ML swap corrected)")
    b = rep["best"]
    ax[3].imshow(sim, cmap=CM, vmin=0, vmax=9, interpolation="nearest")
    ax[3].set_title(f"d  P7 template sec {b['sec']} ({b['offset_um']:+.0f} um), similarity fit\n"
                    f"rot {np.rad2deg(b['params'][0]):.0f} deg, scale {np.exp(b['params'][1]):.2f}, "
                    f"mirror {b['mirror']}" + dtxt("template_similarity"))
    ax[4].imshow(bs, cmap=CM, vmin=0, vmax=9, interpolation="nearest")
    ax[4].set_title(f"e  template after B-spline warp ({um:.0f} um control spacing)"
                    + dtxt(f"template_bspline_{um:.0f}um"))
    for a in ax[:5]:
        a.set_xticks([]); a.set_yticks([])

    secs = rep["sections"]
    off = np.array([s["offset_um"] for s in secs])
    sc = np.array([s["score"] for s in secs])
    ax[5].plot(off, sc, ".", color="0.3")
    ax[5].plot([b["offset_um"]], [b["score"]], "o", color="red", ms=9, mfc="none")
    ax[5].set_xlabel("template section offset from midline (um)")
    ax[5].set_ylabel("similarity-fit score (mean soft Dice)")
    ax[5].set_title("f  which template section fits best")
    ho = rep.get("holdout")
    if ho:
        txt = "hold-out: accuracy inside hidden\n%d disks of r = %.0f um (3 seeds)\n" % (
            ho["n_disks"], ho["radius_um"]) + "\n".join(f"{k}: {v:.3f}" for k, v in ho["mean_acc"].items())
        ax[5].text(0.02, 0.02, txt, transform=ax[5].transAxes, fontsize=9, va="bottom",
                   bbox=dict(fc="white", ec="0.7"))
    fig.suptitle(f"{slide}: Allen P7 Nissl template fitted to the IF section "
                 "(cyan EGL, green ML, yellow IGL, orange WM)")
    fig.tight_layout()
    p = d / f"{slide}_fit.png"
    fig.savefig(p)
    plt.close(fig)
    from PIL import Image
    print(p, Image.open(p).size)


if __name__ == "__main__":
    main(*sys.argv[1:2])
