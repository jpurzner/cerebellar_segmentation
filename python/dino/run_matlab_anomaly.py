"""Run with-DCN-anomaly V2 MATLAB on all 8 20x slides.

V2 detects ONLY signal-dropout anomalies (bleaching + tearing), classified
by shape. Fold + seam detectors removed (V1 fold was firing on EGL).

The anomaly variant (cerebellum_threshold_segment20x_with_dcn_anomaly):
  - detects DAPI < 0.05 inside cerebellum (signal dropout)
  - classifies dropout components by shape:
      long thin (eccentricity > 0.95 OR aspect > 5) → tearing
      solid roundish/rectangular                    → bleaching
  - subtracts anomaly from all_cerebellum before layer detection
  - labels detected regions as class 9 in the output

Per-slide 4-panel review:
  RGB | RGB+anomaly overlay (R=bleach, G=tear) | with-DCN+anomaly output | gold

Output:
  /python/matlab_anomaly/<slide>_anomaly_review.png
  /python/matlab_anomaly/_INDEX.png
  /python/matlab_anomaly/SUMMARY.txt
  /python/matlab_anomaly/<slide>_dcna.log
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
OUT = ROOT.parent / "matlab_anomaly"
OUT.mkdir(exist_ok=True)
PX_UM = 0.5119049
MATLAB = "/Applications/MATLAB_R2022b.app/bin/matlab"
REPO_ROOT = "/Users/jpurzner/Dropbox/imaging_analysis/cerebellar_segmentation"
KOVESI = "/Users/jpurzner/Dropbox/imaging_analysis/MatlabFns"

# 10-class colormap including DCN (8) + anomaly (9)
LABEL_COLORS = np.array([
    [0,   0,   0],   # 0 background
    [0,   0, 128],   # 1 all_cerebellum
    [0, 128, 255],   # 2 iEGL
    [0, 255, 255],   # 3 oEGL
    [255, 255, 0],   # 4 IGL
    [128, 255, 128], # 5 ML
    [255, 128, 0],   # 6 DWL
    [255,   0, 0],   # 7 PCL
    [255, 0, 255],   # 8 DCN — magenta
    [80,  80,  80],  # 9 anomaly — dark gray
], dtype=np.uint8)
LABEL_NAMES = {0:"bg",1:"cereb",2:"iEGL",3:"oEGL",4:"IGL",5:"ML",6:"DWL",7:"PCL",8:"DCN",9:"anomaly"}


def to_unit(im):
    im = im.astype(np.float32); lo, hi = im.min(), im.max()
    return (im - lo) / max(hi - lo, 1e-9)


def make_input(stack, sandbox, slide_id):
    out = sandbox / f"{slide_id}_native.tif"
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
    layer = cerebellum_threshold_segment20x_with_dcn_anomaly(img);
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
    cands = sorted(
        list(LABELLED_DIR.glob(f"{slide_id}_corrected.tif*")) +
        list(LABELLED_DIR.glob(f"{slide_id}_labelled.tif*")) +
        list(LABELLED_DIR.glob(f"{slide_id}_labelld.tif*")),
        key=lambda p: p.stat().st_mtime, reverse=True)
    cands = [c for c in cands if c.suffix in (".tif", ".tiff")]
    return cands[0] if cands else None


def overlay_anomaly_on_rgb(rgb, bleach, tear, alpha=0.6):
    """Overlay anomaly masks onto RGB: R=bleach, magenta=tear."""
    out = rgb.copy()
    H, W = rgb.shape[:2]
    bleach = bleach[:H, :W] > 0
    tear = tear[:H, :W] > 0
    out[bleach] = (1 - alpha) * out[bleach] + alpha * np.array([1.0, 0.2, 0.0])
    out[tear]   = (1 - alpha) * out[tear]   + alpha * np.array([1.0, 0.0, 1.0])
    return np.clip(out, 0, 1)


def main():
    slides = manifest()
    twenty = [s for s in slides if s.mag == "20x"]
    print(f"running anomaly variant on {len(twenty)} 20x slides")

    sandbox = OUT / "sandbox"
    sandbox.mkdir(exist_ok=True)

    rows = []
    for s in twenty:
        tag = f"{s.slide_id}_dcna"
        print(f"\n=== {s.slide_id} ===")
        try:
            stack = tifffile.imread(s.input_path)
        except Exception as e:
            print(f"  cannot read input: {e}"); continue

        rgb = np.stack([to_unit(stack[0]), to_unit(stack[1]),
                        to_unit(stack[2])], axis=-1)
        rgb /= max(rgb.max(), 1e-9)

        input_path = make_input(stack, sandbox, s.slide_id)
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

        # Per-anomaly-type masks (V2: bleach + tear; fold/seam dropped)
        stem = input_path.stem
        bleach = tifffile.imread(sandbox / f"{stem}_anomaly_bleach.tif")
        tear   = tifffile.imread(sandbox / f"{stem}_anomaly_tear.tif")

        gold_path = find_gold(s.slide_id)
        gold = tifffile.imread(gold_path).astype(np.uint8) if gold_path else None

        H = labels.shape[0]; W = labels.shape[1]
        if gold is not None:
            H = min(H, gold.shape[0]); W = min(W, gold.shape[1])
            gold = gold[:H, :W]
        H = min(H, rgb.shape[0]); W = min(W, rgb.shape[1])
        labels = labels[:H, :W]
        rgb_c = rgb[:H, :W]

        cereb_px = max((labels != 0).sum(), 1)
        bleach_pct = 100 * (bleach > 0).sum() / cereb_px
        tear_pct   = 100 * (tear > 0).sum() / cereb_px
        anomaly_pct = 100 * (labels == 9).sum() / cereb_px

        print(f"  anomaly: bleach={bleach_pct:.2f}% tear={tear_pct:.2f}%  "
              f"total={anomaly_pct:.2f}% (of cereb)")
        print(f"  class %  iEGL={100*(labels==2).sum()/labels.size:.2f}  "
              f"oEGL={100*(labels==3).sum()/labels.size:.2f}  "
              f"IGL={100*(labels==4).sum()/labels.size:.2f}  "
              f"ML={100*(labels==5).sum()/labels.size:.2f}  "
              f"DWL={100*(labels==6).sum()/labels.size:.2f}  "
              f"DCN={100*(labels==8).sum()/labels.size:.2f}  "
              f"anomaly={100*(labels==9).sum()/labels.size:.2f}")

        rows.append({
            "slide": s.slide_id, "labels": labels, "rgb": rgb_c,
            "bleach_pct": bleach_pct, "tear_pct": tear_pct,
            "anomaly_pct": anomaly_pct, "elapsed": elapsed,
            "bleach": bleach, "tear": tear,
            "gold": gold, "gold_name": gold_path.name if gold_path else None,
        })

        # 4-panel review per slide
        n_cols = 4 if gold is not None else 3
        fig, axes = plt.subplots(1, n_cols, figsize=(7*n_cols, 7))
        axes[0].imshow(rgb_c); axes[0].axis("off")
        axes[0].set_title(f"{s.slide_id}\nRGB", fontsize=10)
        overlay = overlay_anomaly_on_rgb(rgb_c, bleach, tear)
        axes[1].imshow(overlay); axes[1].axis("off")
        axes[1].set_title(
            f"Anomaly overlay (R=bleach magenta=tear)\n"
            f"bleach={bleach_pct:.2f}% tear={tear_pct:.2f}%",
            fontsize=10)
        axes[2].imshow(LABEL_COLORS[labels]); axes[2].axis("off")
        axes[2].set_title(
            f"with-DCN+anomaly ({elapsed:.0f}s)\n"
            f"DCN={100*(labels==8).sum()/labels.size:.2f}% "
            f"anomaly={100*(labels==9).sum()/labels.size:.2f}%",
            fontsize=10)
        if gold is not None:
            axes[3].imshow(LABEL_COLORS[gold]); axes[3].axis("off")
            axes[3].set_title(f"GOLD\n({gold_path.name})", fontsize=10)
        fig.tight_layout()
        out_path = OUT / f"{s.slide_id}_anomaly_review.png"
        fig.savefig(out_path, dpi=85)
        plt.close(fig)
        print(f"  saved {out_path.name}")

    # Contact sheet
    if rows:
        print("\nbuilding contact sheet...")
        n = len(rows)
        fig, axes = plt.subplots(n, 1, figsize=(28, n*7))
        if n == 1: axes = [axes]
        for i, r in enumerate(rows):
            img = plt.imread(OUT / f"{r['slide']}_anomaly_review.png")
            axes[i].imshow(img); axes[i].axis("off")
            axes[i].set_title(
                f"{r['slide']}  bleach={r['bleach_pct']:.2f}% "
                f"tear={r['tear_pct']:.2f}%  total={r['anomaly_pct']:.2f}%",
                fontsize=11)
        fig.suptitle(
            "with-DCN + anomaly V2 (DAPI dropout) — 8 20x slides "
            "(R=bleach magenta=tear, gray=anomaly class in segmentation)",
            fontsize=14)
        fig.tight_layout()
        out_index = OUT / "_INDEX.png"
        fig.savefig(out_index, dpi=70)
        plt.close(fig)
        print(f"contact sheet: {out_index}")

    # Summary file
    with open(OUT / "SUMMARY.txt", "w") as f:
        f.write("Anomaly detection V2 summary (% of cerebellum tissue)\n")
        f.write("Detector: DAPI < 0.05 inside cerebellum, classified by shape\n")
        f.write("=" * 70 + "\n\n")
        f.write(f"{'slide':32s} {'bleach':>7} {'tear':>6} {'total':>7}\n")
        for r in rows:
            f.write(f"  {r['slide']:30s} {r['bleach_pct']:>6.2f}% "
                    f"{r['tear_pct']:>5.2f}% {r['anomaly_pct']:>6.2f}%\n")
    print(f"summary: {OUT / 'SUMMARY.txt'}")


if __name__ == "__main__":
    main()
