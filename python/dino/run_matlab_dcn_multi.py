"""Run with-DCN MATLAB variant on all 8 good 20x slides — native resolution.

Native gives DCN; ds2x loses it (per earlier finding in matlab_dcn_test).
Goal: confirm DCN segmentation is anatomically correct across slides AND
verify the EGL-preservation fix (DCN restricted to deep DWL interior).

Output:
  /python/matlab_dcn_multi/<slide>_with_dcn.png   — per-slide RGB | DCN | gold
  /python/matlab_dcn_multi/_INDEX.png              — contact sheet
  /python/matlab_dcn_multi/<slide>_with_dcn.log    — MATLAB log
"""
from __future__ import annotations
import subprocess
import time
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import tifffile

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from slide_manifest import manifest

ROOT = Path(__file__).resolve().parent
LABELLED_DIR = ROOT.parent / "labelled"
OUT = ROOT.parent / "matlab_dcn_multi"
OUT.mkdir(exist_ok=True)
PX_UM = 0.5119049
MATLAB = "/Applications/MATLAB_R2022b.app/bin/matlab"
REPO_ROOT = "/Users/jpurzner/Dropbox/imaging_analysis/cerebellar_segmentation"
KOVESI = "/Users/jpurzner/Dropbox/imaging_analysis/MatlabFns"

# 9-class colormap including DCN (label 8 = magenta)
LABEL_COLORS = np.array([
    [0,   0,   0],   # 0 background
    [0,   0, 128],   # 1 all_cerebellum
    [0, 128, 255],   # 2 iEGL
    [0, 255, 255],   # 3 oEGL
    [255, 255, 0],   # 4 IGL (parser-name; biologically ML — see PROJECT_HISTORY §4)
    [128, 255, 128], # 5 ML  (parser-name; biologically IGL)
    [255, 128, 0],   # 6 DWL
    [255,   0, 0],   # 7 PCL
    [255, 0, 255],   # 8 DCN — magenta (NEW)
], dtype=np.uint8)

LABEL_NAMES = {0:"bg",1:"cereb",2:"iEGL",3:"oEGL",4:"IGL",5:"ML",6:"DWL",7:"PCL",8:"DCN"}


def to_unit(im):
    im = im.astype(np.float32); lo, hi = im.min(), im.max()
    return (im - lo) / max(hi - lo, 1e-9)


def make_native_input(stack, sandbox, base):
    """Save the native-resolution stack as TIFF for MATLAB to consume."""
    out = sandbox / f"{base}_native.tif"
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


def find_gold(slide_id):
    """Find a hand-corrected gold-standard label TIFF, or return None."""
    cands = sorted(
        list(LABELLED_DIR.glob(f"{slide_id}_corrected.tif*")) +
        list(LABELLED_DIR.glob(f"{slide_id}_labelled.tif*")) +
        list(LABELLED_DIR.glob(f"{slide_id}_labelld.tif*")),
        key=lambda p: p.stat().st_mtime, reverse=True)
    cands = [c for c in cands if c.suffix in (".tif", ".tiff")]
    return cands[0] if cands else None


def class_breakdown(labels, gold):
    """Return list of (class_name, pred_pct, gold_pct) for non-trivial classes."""
    rows = []
    for c in range(9):
        pp = 100*(labels == c).sum()/labels.size
        gp = 100*(gold == c).sum()/gold.size if gold is not None else 0.0
        if pp > 0.05 or gp > 0.05:
            rows.append((LABEL_NAMES.get(c, str(c)), pp, gp))
    return rows


def main():
    slides = manifest()
    # 20x only — DCN is a 20x feature, and the with-DCN script is
    # cerebellum_threshold_segment20x_with_dcn (20x specific).
    twenty = [s for s in slides if s.mag == "20x"]
    print(f"running with-DCN on {len(twenty)} 20x slides")

    sandbox = OUT / "sandbox"
    sandbox.mkdir(exist_ok=True)

    rows = []
    for s in twenty:
        tag = f"{s.slide_id}_dcn"
        print(f"\n=== {s.slide_id} ===")
        try:
            stack = tifffile.imread(s.input_path)
        except Exception as e:
            print(f"  cannot read input: {e}")
            continue

        rgb = np.stack([to_unit(stack[0]), to_unit(stack[1]),
                        to_unit(stack[2])], axis=-1)
        rgb /= max(rgb.max(), 1e-9)

        input_path = make_native_input(stack, sandbox, s.slide_id)
        t0 = time.time()
        rc, log = run_matlab(str(input_path), sandbox, tag)
        elapsed = time.time() - t0
        print(f"  MATLAB done in {elapsed:.0f}s, rc={rc}")

        labels_path = input_path.with_name(input_path.stem + "_labels.tif")
        if not labels_path.exists():
            print(f"  labels not found at {labels_path}")
            print(f"  see {log}")
            continue
        labels = tifffile.imread(labels_path).astype(np.uint8)

        gold_path = find_gold(s.slide_id)
        gold = tifffile.imread(gold_path).astype(np.uint8) if gold_path else None

        # Crop everything to common size
        H = labels.shape[0]; W = labels.shape[1]
        if gold is not None:
            H = min(H, gold.shape[0]); W = min(W, gold.shape[1])
            gold = gold[:H, :W]
        H = min(H, rgb.shape[0]); W = min(W, rgb.shape[1])
        labels = labels[:H, :W]
        rgb_c = rgb[:H, :W]

        # Class breakdown
        print("  class %  with_dcn  gold")
        for name, pp, gp in class_breakdown(labels, gold):
            print(f"    {name:6s}  {pp:6.2f}    {gp:6.2f}")

        dcn_pct = 100 * (labels == 8).sum() / labels.size
        rows.append({
            "slide": s.slide_id, "labels": labels, "rgb": rgb_c,
            "gold": gold, "gold_name": gold_path.name if gold_path else None,
            "elapsed": elapsed, "dcn_pct": dcn_pct,
        })

        # Per-slide panel
        n_cols = 3 if gold is not None else 2
        fig, axes = plt.subplots(1, n_cols, figsize=(7*n_cols, 7))
        axes[0].imshow(rgb_c); axes[0].axis("off")
        axes[0].set_title(f"{s.slide_id}\nRGB", fontsize=10)
        axes[1].imshow(LABEL_COLORS[labels]); axes[1].axis("off")
        axes[1].set_title(
            f"with-DCN MATLAB ({elapsed:.0f}s)\nDCN (magenta) = {dcn_pct:.2f}%",
            fontsize=10)
        if gold is not None:
            axes[2].imshow(LABEL_COLORS[gold]); axes[2].axis("off")
            axes[2].set_title(f"GOLD\n({gold_path.name})", fontsize=10)
        fig.tight_layout()
        out_path = OUT / f"{s.slide_id}_with_dcn.png"
        fig.savefig(out_path, dpi=85)
        plt.close(fig)
        print(f"  saved {out_path.name}")

    # Contact sheet of all per-slide panels
    if rows:
        print("\nbuilding contact sheet...")
        n = len(rows)
        ncol = 2
        nrow = (n + ncol - 1) // ncol
        fig, axes = plt.subplots(nrow, ncol, figsize=(20, nrow*8))
        axes = np.atleast_2d(axes).ravel()
        for i, r in enumerate(rows):
            img = plt.imread(OUT / f"{r['slide']}_with_dcn.png")
            axes[i].imshow(img); axes[i].axis("off")
            flag = "★" if r["gold"] is not None else ""
            axes[i].set_title(
                f"{r['slide']}{flag}  DCN={r['dcn_pct']:.2f}%",
                fontsize=10)
        for i in range(n, len(axes)):
            axes[i].axis("off")
        fig.suptitle(
            "with-DCN MATLAB on 20x slides — magenta = DCN, "
            "★ = gold-corrected available",
            fontsize=14)
        fig.tight_layout()
        out_index = OUT / "_INDEX.png"
        fig.savefig(out_index, dpi=80)
        plt.close(fig)
        print(f"contact sheet: {out_index}")

    # Summary
    print("\n=== summary ===")
    print(f"{'slide':32s} {'time':>5} {'DCN%':>7}")
    for r in rows:
        print(f"  {r['slide']:30s} {r['elapsed']:>4.0f}s {r['dcn_pct']:>6.2f}")


if __name__ == "__main__":
    main()
