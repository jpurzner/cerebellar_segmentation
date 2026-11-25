"""
Apply ribbon_join.py to the eight 20x slides at native resolution and draw a
before/after figure.

    python scripts/ribbon_join_apply.py

Decisions are made at 4x downsampling and stamped onto native pixels that had the
source label (same as depth_model/anatomy_fixes.py apply_native).
Input:  depth_model/out/anatomy_fixes/<slide>_fixed_labels.tif
Output: outputs/ribbon_join/native/<slide>_joined.tif       rules P + M + W + F + S (default)
        outputs/ribbon_join/native/<slide>_joined_PMW.tif   rules P + M + W (no fissure rules)
        outputs/ribbon_join/native_scores.json              s1_2 gold accuracy at native res
        outputs/ribbon_join/joined_s2_5_s2_2.png            figure (1800 x 1200)
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
import ribbon_join as RJ  # noqa: E402
import ribbon_join_eval as E  # noqa: E402

NAT = E.OUT / "native"
LABEL_COLORS = [(0, 0, 0), (0.1, 0.1, 0.45), (0.4, 0.7, 0.95), (0, 0.85, 0.9), (1, 0.85, 0.1),
                (0.55, 0.9, 0.55), (1, 0.5, 0.1), (0.9, 0.1, 0.1), (0.85, 0.2, 0.85), (0.45, 0.45, 0.45)]
CM = matplotlib.colors.ListedColormap(LABEL_COLORS)


def apply_native(native, ds, decided_ds):
    H, W = native.shape
    ys = np.minimum(np.arange(H) // ds, decided_ds.shape[0] - 1)
    xs = np.minimum(np.arange(W) // ds, decided_ds.shape[1] - 1)
    orig_up = native[::ds, ::ds][np.ix_(ys, xs)]
    new_up = decided_ds[np.ix_(ys, xs)]
    change = (new_up != orig_up) & (native == orig_up)
    out = native.copy()
    out[change] = new_up[change]
    return out


def gold_native():
    g = tifffile.imread(E.SEG / "labelled" / f"{E.GOLD_SLIDE}_labelled.tiff").astype(np.int16)
    gs = g.copy()
    gs[g == 4], gs[g == 5] = 5, 4
    return gs


def main():
    NAT.mkdir(parents=True, exist_ok=True)
    scores = {}
    keep = {}
    for s in E.SLIDES:
        native = tifffile.imread(E.FIXED / f"{s}_fixed_labels.tif").astype(np.uint8)
        ds = native[::E.DS, ::E.DS]
        for tag, rules in (("", ("P", "M", "W", "F", "S")), ("_PMW", ("P", "M", "W"))):
            j, _ = RJ.join(ds, E.PX, {"rules": rules})
            out = apply_native(native, E.DS, j)
            tifffile.imwrite(NAT / f"{s}_joined{tag}.tif", out, compression="zlib")
            if s == E.GOLD_SLIDE:
                g = gold_native()
                m = np.isin(g, E.SCORED)
                scores["+".join(rules)] = dict(before=float((native[m] == g[m]).mean()),
                                               after=float((out[m] == g[m]).mean()),
                                               iou_before=E.ious(native, g), iou_after=E.ious(out, g))
            if s in ("2018_05_22_s2_5_p27-0026", "2018_05_22_s2_2_p27") and tag == "":
                keep[s] = (ds, j)
        print(f"{s}: written", flush=True)
    (E.OUT / "native_scores.json").write_text(json.dumps(scores, indent=1))
    for k, v in scores.items():
        print(f"s1_2 native, {k}: acc {v['before']:.4f} -> {v['after']:.4f}; IoU {v['iou_before']} -> {v['iou_after']}")

    fig, ax = plt.subplots(2, 3, figsize=(18, 12), dpi=100)
    for r, (s, (before, after)) in enumerate(keep.items()):
        ch = before != after
        ax[r, 0].imshow(before, cmap=CM, vmin=0, vmax=9, interpolation="nearest")
        ax[r, 0].set_title(f"{s[-14:]}  before (anatomy-fixed labels)")
        ax[r, 1].imshow(after, cmap=CM, vmin=0, vmax=9, interpolation="nearest")
        ax[r, 1].set_title("after rules P + M + W + F + S")
        base = np.where(before > 0, 0.25, 1.0)
        rgb = np.dstack([base] * 3)
        cols = np.array(LABEL_COLORS)
        rgb[ch] = cols[after[ch]]
        ax[r, 2].imshow(rgb, interpolation="nearest")
        ax[r, 2].set_title(f"changed pixels in their new colour ({100 * ch.sum() / (before > 0).sum():.1f}% of tissue)")
    for a in ax.ravel():
        a.set_xticks([]); a.set_yticks([])
    fig.suptitle("Ribbon joining by layer order (P, M, W, F, S): EGL sky/cyan, ML green, IGL yellow, DWL orange, PCL red")
    fig.tight_layout()
    p = E.OUT / "joined_s2_5_s2_2.png"
    fig.savefig(p)
    plt.close(fig)
    from PIL import Image
    print(p, Image.open(p).size)


if __name__ == "__main__":
    main()
