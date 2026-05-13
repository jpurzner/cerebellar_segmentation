"""Test preprocessing variants on all 8 20x slides.

s1_2 already showed preprocessing doesn't improve a working slide. But the
user says other 20x slides have global failures (layers collapse, DWL eats
everything, etc.). Maybe preprocessing rescues those.

Variants to test (skip the catastrophic gauss_sig2 and maxfilt_r3):
  - gauss_sig1um   (light Gaussian)
  - maxfilt_r2px   (1 um disk dilation — best EGL on s1_2)
  - median_5px     (edge-preserving smoothing)

Native baseline comes from /python/matlab_dcn_multi/ (already computed).

For each (slide, variant) we:
  - run with-DCN MATLAB
  - read the raw labels
  - report per-class %s
  - flag global failures (any class outside expected range)

Output:
  /python/matlab_preproc_8slides/<slide>_compare.png   (4-col: RGB|native|variants)
  /python/matlab_preproc_8slides/_INDEX.png
  /python/matlab_preproc_8slides/<slide>_<variant>.log
  /python/matlab_preproc_8slides/SUMMARY.txt          (per-slide failure flags)
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
OUT = ROOT.parent / "matlab_preproc_8slides"
OUT.mkdir(exist_ok=True)
NATIVE_DIR = ROOT.parent / "matlab_dcn_multi" / "sandbox"   # reuse native runs
PX_UM = 0.5119049
MATLAB = "/Applications/MATLAB_R2022b.app/bin/matlab"
REPO_ROOT = "/Users/jpurzner/Dropbox/imaging_analysis/cerebellar_segmentation"
KOVESI = "/Users/jpurzner/Dropbox/imaging_analysis/MatlabFns"

LABEL_COLORS = np.array([
    [0,   0,   0], [0,   0, 128], [0, 128, 255], [0, 255, 255],
    [255, 255, 0], [128, 255, 128],[255, 128, 0], [255,   0, 0],
    [255, 0, 255],
], dtype=np.uint8)
LABEL_NAMES = {0:"bg",1:"cereb",2:"iEGL",3:"oEGL",4:"IGL",5:"ML",6:"DWL",7:"PCL",8:"DCN"}


# Variants
def variant_native(stk):     return stk.copy()

def _gauss(stk, sig_um):
    sig_px = sig_um / PX_UM
    out = np.zeros_like(stk)
    for c in range(3):
        out[c] = ndi.gaussian_filter(stk[c].astype(np.float32),
                                      sigma=sig_px).astype(stk.dtype)
    return out

def variant_gauss1(stk):     return _gauss(stk, 1.0)

def _maxfilt(stk, r_px):
    yy, xx = np.ogrid[-r_px:r_px+1, -r_px:r_px+1]
    foot = (xx*xx + yy*yy) <= r_px*r_px
    out = np.zeros_like(stk)
    for c in range(3):
        out[c] = ndi.maximum_filter(stk[c], footprint=foot)
    return out

def variant_maxr2(stk):      return _maxfilt(stk, 2)

def variant_median5(stk):
    out = np.zeros_like(stk)
    for c in range(3):
        out[c] = ndi.median_filter(stk[c], size=5)
    return out


# Don't include native — we'll load it from matlab_dcn_multi/
NEW_VARIANTS = [
    ("gauss_sig1um", variant_gauss1),
    ("maxfilt_r2px", variant_maxr2),
    ("median_5px",   variant_median5),
]


def to_unit(im):
    im = im.astype(np.float32); lo, hi = im.min(), im.max()
    return (im - lo) / max(hi - lo, 1e-9)


def make_input(stack, sandbox, name):
    out = sandbox / f"{name}.tif"
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


def class_pcts(labels):
    """Return dict {class_name: pct}."""
    return {LABEL_NAMES[c]: 100*(labels == c).sum()/labels.size
            for c in LABEL_NAMES}


def flag_failures(pcts):
    """Return list of failure-mode strings, empty if no global failure."""
    flags = []
    if pcts["DWL"] > 25:   flags.append(f"DWL_HIGH({pcts['DWL']:.0f}%)")
    if pcts["DWL"] < 4:    flags.append(f"DWL_LOW({pcts['DWL']:.1f}%)")
    if pcts["iEGL"] > 12:  flags.append(f"iEGL_HIGH({pcts['iEGL']:.0f}%)")
    if pcts["iEGL"] < 3:   flags.append(f"iEGL_LOW({pcts['iEGL']:.1f}%)")
    if pcts["oEGL"] > 8:   flags.append(f"oEGL_HIGH({pcts['oEGL']:.1f}%)")
    if pcts["oEGL"] < 1:   flags.append(f"oEGL_LOW({pcts['oEGL']:.2f}%)")
    if pcts["IGL"] < 5:    flags.append(f"IGL_LOW({pcts['IGL']:.1f}%)")
    if pcts["DCN"] > 5:    flags.append(f"DCN_HIGH({pcts['DCN']:.0f}%)")
    return flags


def main():
    slides = manifest()
    twenty = [s for s in slides if s.mag == "20x"]
    print(f"running {len(NEW_VARIANTS)} variants on {len(twenty)} 20x slides")

    sandbox = OUT / "sandbox"
    sandbox.mkdir(exist_ok=True)

    summary_rows = []
    for s in twenty:
        print(f"\n=== {s.slide_id} ===")
        stack = tifffile.imread(s.input_path)
        rgb = np.stack([to_unit(stack[0]), to_unit(stack[1]),
                        to_unit(stack[2])], axis=-1)
        rgb /= max(rgb.max(), 1e-9)

        # Load native baseline (already computed in matlab_dcn_multi)
        native_path = NATIVE_DIR / f"{s.slide_id}_native_labels.tif"
        if not native_path.exists():
            print(f"  WARNING: no native baseline at {native_path}")
            native = None
        else:
            native = tifffile.imread(native_path).astype(np.uint8)

        slide_results = {"native": native}
        for vname, vfn in NEW_VARIANTS:
            tag = f"{s.slide_id}_{vname}"
            stk_v = vfn(stack)
            input_path = make_input(stk_v, sandbox, tag)
            t0 = time.time()
            rc, log = run_matlab(str(input_path), sandbox, tag)
            elapsed = time.time() - t0
            labels_path = input_path.with_name(input_path.stem + "_labels.tif")
            if not labels_path.exists():
                print(f"  {vname}: failed ({elapsed:.0f}s, rc={rc})")
                slide_results[vname] = None; continue
            labels = tifffile.imread(labels_path).astype(np.uint8)
            slide_results[vname] = labels
            print(f"  {vname}: {elapsed:.0f}s OK")

        # Compute %s + failure flags per variant
        print(f"  {'variant':14s} {'iEGL':>5} {'oEGL':>5} {'IGL':>5} {'ML':>5} "
              f"{'DWL':>5} {'PCL':>5} {'DCN':>5}  flags")
        per_variant = {}
        for vname in ["native"] + [v for v,_ in NEW_VARIANTS]:
            labs = slide_results[vname]
            if labs is None:
                print(f"  {vname:12s}  MISSING"); continue
            pcts = class_pcts(labs)
            flags = flag_failures(pcts)
            per_variant[vname] = {"pcts": pcts, "flags": flags, "labels": labs}
            print(f"  {vname:12s}  {pcts['iEGL']:>5.2f} {pcts['oEGL']:>5.2f} "
                  f"{pcts['IGL']:>5.2f} {pcts['ML']:>5.2f} {pcts['DWL']:>5.2f} "
                  f"{pcts['PCL']:>5.2f} {pcts['DCN']:>5.2f}  "
                  f"{', '.join(flags) if flags else 'OK'}")

        summary_rows.append({"slide": s.slide_id, "rgb": rgb,
                             "per_variant": per_variant})

        # Per-slide comparison panel
        n_panels = 1 + len(per_variant)   # RGB + variants
        fig, axes = plt.subplots(1, n_panels, figsize=(n_panels*7, 7))
        axes[0].imshow(rgb); axes[0].axis("off")
        axes[0].set_title(f"{s.slide_id}\nRGB", fontsize=10)
        for i, vname in enumerate(per_variant.keys()):
            d = per_variant[vname]
            ax = axes[i+1]
            H = min(d["labels"].shape[0], rgb.shape[0])
            W = min(d["labels"].shape[1], rgb.shape[1])
            ax.imshow(LABEL_COLORS[d["labels"][:H,:W]]); ax.axis("off")
            flag_str = "FAIL: " + ", ".join(d["flags"]) if d["flags"] else "ok"
            ax.set_title(
                f"{vname}\nDCN={d['pcts']['DCN']:.2f}%  iEGL={d['pcts']['iEGL']:.1f}%  "
                f"oEGL={d['pcts']['oEGL']:.1f}%\n{flag_str}",
                fontsize=9, color=("red" if d["flags"] else "black"))
        fig.tight_layout()
        fig.savefig(OUT / f"{s.slide_id}_compare.png", dpi=80)
        plt.close(fig)

    # Build summary text
    print("\nbuilding SUMMARY.txt...")
    with open(OUT / "SUMMARY.txt", "w") as f:
        f.write("Per-slide preprocessing comparison — 20x with-DCN MATLAB\n")
        f.write("=" * 80 + "\n\n")
        f.write("Failure indicators (none = OK):\n")
        f.write("  DWL>25%, DWL<4%  : DWL collapsed/exploded\n")
        f.write("  iEGL>12%, iEGL<3%: iEGL detection broken\n")
        f.write("  oEGL>8%, oEGL<1% : oEGL detection broken\n")
        f.write("  IGL<5%           : IGL collapsed\n")
        f.write("  DCN>5%           : DCN exploded\n\n")
        for r in summary_rows:
            f.write(f"--- {r['slide']} ---\n")
            for vname, d in r["per_variant"].items():
                p = d["pcts"]
                flag = ", ".join(d["flags"]) if d["flags"] else "OK"
                f.write(f"  {vname:14s} iEGL={p['iEGL']:5.2f} oEGL={p['oEGL']:5.2f} "
                        f"IGL={p['IGL']:5.2f} ML={p['ML']:5.2f} DWL={p['DWL']:5.2f} "
                        f"DCN={p['DCN']:5.2f}  {flag}\n")
            f.write("\n")
    print(f"summary: {OUT / 'SUMMARY.txt'}")

    # Build contact sheet INDEX
    n = len(summary_rows)
    fig, axes = plt.subplots(n, 1, figsize=(28, n*7))
    if n == 1: axes = [axes]
    for i, r in enumerate(summary_rows):
        img = plt.imread(OUT / f"{r['slide']}_compare.png")
        axes[i].imshow(img); axes[i].axis("off")
        axes[i].set_title(r['slide'], fontsize=12)
    fig.suptitle("Preprocessing variants on 8 20x slides — global failure check",
                 fontsize=16)
    fig.tight_layout()
    fig.savefig(OUT / "_INDEX.png", dpi=70)
    plt.close(fig)
    print(f"contact sheet: {OUT / '_INDEX.png'}")

    # Final summary at console
    print("\n=== which variant FAILS LEAST per slide? ===")
    for r in summary_rows:
        fails = {v: len(d["flags"]) for v, d in r["per_variant"].items()}
        best = min(fails, key=fails.get)
        msg = f"  {r['slide']:30s}  best={best:14s}  ({fails[best]} flags)  "
        msg += "  ".join(f"{v}:{n}" for v,n in fails.items())
        print(msg)


if __name__ == "__main__":
    main()
