"""Sweep input-preprocessing variants on s1_2 with-DCN MATLAB at 20x.

The 20x problem: individual nuclei are resolved with gaps that the MATLAB
pipeline's internal Gaussian (sigma=5 px = 2.56 um) doesn't close. Nuclear
gaps cause layer masks to fragment, and bwareafilt(>=200000 px) then drops
the small pieces.

Hypothesis: preprocessing the input channels (Gaussian blur, max-filter
dilation, median) before MATLAB sees them will merge nuclei into solid
layer signals.

Variants tested:
  - native               (baseline, no preprocessing)
  - gauss_sig1um         (light Gaussian)
  - gauss_sig2um         (medium Gaussian)
  - maxfilt_r2px         (1 um disk dilation - preserves intensity peaks)
  - maxfilt_r3px         (1.5 um disk dilation)
  - median_5px           (median 5x5 smoothing - edge preserving)

Each variant runs cerebellum_threshold_segment20x_with_dcn on s1_2 native
resolution. We compare IoU vs gold for each layer.

Output:
  /python/matlab_preproc_sweep/<variant>_panel.png
  /python/matlab_preproc_sweep/_INDEX.png
  /python/matlab_preproc_sweep/<variant>.log
"""
from __future__ import annotations
import subprocess
import time
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage as ndi
import tifffile

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from slide_manifest import manifest

ROOT = Path(__file__).resolve().parent
LABELLED_DIR = ROOT.parent / "labelled"
OUT = ROOT.parent / "matlab_preproc_sweep"
OUT.mkdir(exist_ok=True)
PX_UM = 0.5119049
MATLAB = "/Applications/MATLAB_R2022b.app/bin/matlab"
REPO_ROOT = "/Users/jpurzner/Dropbox/imaging_analysis/cerebellar_segmentation"
KOVESI = "/Users/jpurzner/Dropbox/imaging_analysis/MatlabFns"

LABEL_COLORS = np.array([
    [0,   0,   0], [0,   0, 128], [0, 128, 255], [0, 255, 255],
    [255, 255, 0], [128, 255, 128],[255, 128, 0], [255,   0, 0],
    [255, 0, 255],   # 8 DCN — magenta
], dtype=np.uint8)

LABEL_NAMES = {0:"bg",1:"cereb",2:"iEGL",3:"oEGL",4:"IGL",5:"ML",6:"DWL",7:"PCL",8:"DCN"}


# Variants: (name, function-on-stack)
def variant_native(stk):     return stk.copy()

def _gauss(stk, sig_um):
    sig_px = sig_um / PX_UM
    out = np.zeros_like(stk)
    for c in range(3):
        out[c] = ndi.gaussian_filter(stk[c].astype(np.float32),
                                      sigma=sig_px).astype(stk.dtype)
    return out

def variant_gauss1(stk):     return _gauss(stk, 1.0)
def variant_gauss2(stk):     return _gauss(stk, 2.0)

def _maxfilt(stk, r_px):
    """Disk-shaped maximum filter — dilates bright nuclei in intensity space."""
    yy, xx = np.ogrid[-r_px:r_px+1, -r_px:r_px+1]
    foot = (xx*xx + yy*yy) <= r_px*r_px
    out = np.zeros_like(stk)
    for c in range(3):
        out[c] = ndi.maximum_filter(stk[c], footprint=foot)
    return out

def variant_maxr2(stk):      return _maxfilt(stk, 2)
def variant_maxr3(stk):      return _maxfilt(stk, 3)

def variant_median5(stk):
    out = np.zeros_like(stk)
    for c in range(3):
        out[c] = ndi.median_filter(stk[c], size=5)
    return out


VARIANTS = [
    ("native",       variant_native),
    ("gauss_sig1um", variant_gauss1),
    ("gauss_sig2um", variant_gauss2),
    ("maxfilt_r2px", variant_maxr2),
    ("maxfilt_r3px", variant_maxr3),
    ("median_5px",   variant_median5),
]


def to_unit(im):
    im = im.astype(np.float32); lo, hi = im.min(), im.max()
    return (im - lo) / max(hi - lo, 1e-9)


def make_input(stack, sandbox, name):
    out = sandbox / f"s1_2_{name}.tif"
    tifffile.imwrite(out, stack, photometric="minisblack")
    return out


def run_matlab(input_tif, sandbox, tag):
    safe = tag.replace("-", "_").replace(".", "_")
    script_name = f"driver_{safe}.m"
    script_path = sandbox / script_name
    code = f"""addpath('{REPO_ROOT}');
addpath(genpath('{KOVESI}'));
addpath(genpath('{REPO_ROOT}/BaiSkeletonPruningDCE'));
addpath(genpath('{REPO_ROOT}/Skeleton'));
addpath(genpath('{REPO_ROOT}/WCHT'));
addpath(genpath('{REPO_ROOT}/CurevUtils1.1'));
set(0, 'DefaultFigureVisible', 'off');
img = '{input_tif}';
fprintf('input: %s\\n', img);
t0 = tic;
try
    layer = cerebellum_threshold_segment20x_with_dcn(img);
    fprintf('OK in %.1fs\\n', toc(t0));
catch ME
    fprintf('FAIL: %s\\n', ME.message);
    for k = 1:length(ME.stack)
        fprintf('  at %s line %d\\n', ME.stack(k).name, ME.stack(k).line);
    end
end
"""
    script_path.write_text(code)
    log = OUT / f"{tag}.log"
    with open(log, "w") as f:
        proc = subprocess.run([MATLAB, "-batch",
                                f"cd('{sandbox}'); run('{script_name}')"],
                               stdout=f, stderr=subprocess.STDOUT, text=True)
    return proc.returncode, log


def per_class_iou(pred, gt):
    """Per-class IoU. Returns dict {class_name: iou} for non-trivial classes."""
    out = {}
    for c, name in LABEL_NAMES.items():
        if c == 0: continue   # skip background
        p = (pred == c); g = (gt == c)
        if g.sum() < 100 and p.sum() < 100:
            continue
        inter = (p & g).sum(); union = (p | g).sum()
        out[name] = inter / max(union, 1)
    return out


def egl_iou(pred, gt):
    """Combined EGL (iEGL ∪ oEGL) IoU — most diagnostic for layer quality."""
    pe = (pred == 2) | (pred == 3); ge = (gt == 2) | (gt == 3)
    inter = (pe & ge).sum(); union = (pe | ge).sum()
    return inter / max(union, 1)


def main():
    slides = manifest()
    s = next(s for s in slides if s.slide_id == "2018_05_22_s1_2_p27-0021")
    stack = tifffile.imread(s.input_path)
    print(f"loaded s1_2: shape {stack.shape}")

    # Find gold
    cands = sorted(
        list(LABELLED_DIR.glob(f"{s.slide_id}_corrected.tif*")) +
        list(LABELLED_DIR.glob(f"{s.slide_id}_labelled.tif*")) +
        list(LABELLED_DIR.glob(f"{s.slide_id}_labelld.tif*")),
        key=lambda p: p.stat().st_mtime, reverse=True)
    cands = [c for c in cands if c.suffix in (".tif", ".tiff")]
    if not cands:
        print("no gold standard found, aborting"); return
    gold = tifffile.imread(cands[0]).astype(np.uint8)
    print(f"gold: {cands[0].name}")

    # Display RGB (R=p27, G=NeuN, B=DAPI) — channel order is (p27, neun, dapi) at 20x
    rgb = np.stack([to_unit(stack[0]), to_unit(stack[1]), to_unit(stack[2])], axis=-1)
    rgb /= max(rgb.max(), 1e-9)

    sandbox = OUT / "sandbox"
    sandbox.mkdir(exist_ok=True)

    rows = []
    for name, fn in VARIANTS:
        print(f"\n=== {name} ===")
        stk_v = fn(stack)
        input_path = make_input(stk_v, sandbox, name)
        t0 = time.time()
        rc, log = run_matlab(str(input_path), sandbox, name)
        elapsed = time.time() - t0
        print(f"  MATLAB done in {elapsed:.0f}s, rc={rc}")
        labels_path = input_path.with_name(input_path.stem + "_labels.tif")
        if not labels_path.exists():
            print(f"  labels not found at {labels_path}")
            print(f"  see {log}")
            continue
        labels = tifffile.imread(labels_path).astype(np.uint8)

        H = min(labels.shape[0], gold.shape[0], rgb.shape[0])
        W = min(labels.shape[1], gold.shape[1], rgb.shape[1])
        labels = labels[:H, :W]
        gld = gold[:H, :W]
        rgb_c = rgb[:H, :W]

        iou_egl = egl_iou(labels, gld)
        ious = per_class_iou(labels, gld)
        dcn_pct = 100 * (labels == 8).sum() / labels.size
        print(f"  EGL IoU vs gold: {iou_egl:.3f}")
        for cname, iou in ious.items():
            pp = 100*(labels == [k for k,v in LABEL_NAMES.items() if v==cname][0]).sum()/labels.size
            gp = 100*(gld == [k for k,v in LABEL_NAMES.items() if v==cname][0]).sum()/gld.size
            print(f"    {cname:6s}  IoU={iou:.3f}  pred={pp:5.2f}%  gold={gp:5.2f}%")

        rows.append({"name": name, "labels": labels, "gold": gld, "rgb": rgb_c,
                     "elapsed": elapsed, "iou_egl": iou_egl,
                     "iou_per": ious, "dcn_pct": dcn_pct})

        # Per-variant panel
        fig, axes = plt.subplots(1, 3, figsize=(21, 7))
        axes[0].imshow(rgb_c); axes[0].axis("off")
        axes[0].set_title(f"s1_2 RGB", fontsize=11)
        axes[1].imshow(LABEL_COLORS[labels]); axes[1].axis("off")
        axes[1].set_title(
            f"{name} ({elapsed:.0f}s)\nEGL IoU={iou_egl:.3f}  DCN={dcn_pct:.2f}%",
            fontsize=11)
        axes[2].imshow(LABEL_COLORS[gld]); axes[2].axis("off")
        axes[2].set_title("GOLD", fontsize=11)
        fig.tight_layout()
        fig.savefig(OUT / f"{name}_panel.png", dpi=85)
        plt.close(fig)

    if not rows:
        print("\nno successful runs"); return

    # Summary panel: rows × 2 (variant output | gold)
    print("\nbuilding summary panel...")
    fig, axes = plt.subplots(len(rows), 2, figsize=(14, len(rows)*7))
    if len(rows) == 1: axes = axes[None, :]
    for i, r in enumerate(rows):
        axes[i, 0].imshow(LABEL_COLORS[r["labels"]]); axes[i, 0].axis("off")
        axes[i, 0].set_title(
            f"{r['name']} — EGL IoU={r['iou_egl']:.3f}  DCN={r['dcn_pct']:.2f}%",
            fontsize=11)
        axes[i, 1].imshow(LABEL_COLORS[r["gold"]]); axes[i, 1].axis("off")
        axes[i, 1].set_title("GOLD", fontsize=11)
    fig.suptitle("s1_2 input-preprocessing sweep — with-DCN MATLAB",
                 fontsize=14)
    fig.tight_layout()
    fig.savefig(OUT / "_INDEX.png", dpi=80)
    plt.close(fig)
    print(f"summary: {OUT / '_INDEX.png'}")

    # Final summary table
    print("\n=== summary ===")
    cls_keys = ["iEGL", "oEGL", "IGL", "ML", "DWL", "PCL", "DCN"]
    header = f"{'variant':14s} {'time':>5} {'EGL':>6} " + " ".join(f"{c:>6}" for c in cls_keys)
    print(header)
    for r in rows:
        line = f"  {r['name']:12s} {r['elapsed']:>4.0f}s {r['iou_egl']:>6.3f} "
        for c in cls_keys:
            line += f"{r['iou_per'].get(c, 0.0):>6.3f} "
        print(line)


if __name__ == "__main__":
    main()
