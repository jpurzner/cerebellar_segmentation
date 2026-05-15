"""Run with-DCN 10x MATLAB variant on all 11 good 10x slides.

Uses cerebellum_threshold_segment2_with_dcn_anomaly (the 10x sibling of
cerebellum_threshold_segment20x_with_dcn). Channel order in 10x input TIFFs
is already (p27, DAPI, NeuN) — that's what the MATLAB script expects, so we
write the stack as-is.

Output:
  /python/matlab_anomaly_10x/<slide>_with_dcn.png
  /python/matlab_anomaly_10x/_INDEX.png
  /python/matlab_anomaly_10x/<slide>_dcn.log
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
OUT = ROOT.parent / "matlab_anomaly_10x"
OUT.mkdir(exist_ok=True)
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
    [255, 0, 255],   # 8 DCN — magenta
    [80,  80,  80],  # 9 anomaly — dark gray (NEW)
], dtype=np.uint8)

LABEL_NAMES = {0:"bg",1:"cereb",2:"iEGL",3:"oEGL",4:"IGL",5:"ML",6:"DWL",7:"PCL",8:"DCN"}


def to_unit(im):
    im = im.astype(np.float32); lo, hi = im.min(), im.max()
    return (im - lo) / max(hi - lo, 1e-9)


def make_input(stack, sandbox, slide_id):
    """Write stack as-is (10x channels are already p27/DAPI/NeuN in TIFF order)."""
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
    layer = cerebellum_threshold_segment2_with_dcn_anomaly(img);
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


def make_rgb_for_display(stack, channel_order):
    """Build an RGB display image. channel_order tells us what each TIFF page is.
    For visual review we always show: R=p27, G=NeuN, B=DAPI (bio convention)."""
    role_to_idx = {role: i for i, role in enumerate(channel_order)}
    p27  = stack[role_to_idx["p27"]]
    neun = stack[role_to_idx["neun"]]
    dapi = stack[role_to_idx["dapi"]]
    rgb = np.stack([to_unit(p27), to_unit(neun), to_unit(dapi)], axis=-1)
    return rgb / max(rgb.max(), 1e-9)


def main():
    slides = manifest()
    ten = [s for s in slides if s.mag == "10x"]
    print(f"running with-DCN on {len(ten)} 10x slides")

    sandbox = OUT / "sandbox"
    sandbox.mkdir(exist_ok=True)

    rows = []
    for s in ten:
        tag = f"{s.slide_id}_dcn"
        print(f"\n=== {s.slide_id} ===")
        try:
            stack = tifffile.imread(s.input_path)
        except Exception as e:
            print(f"  cannot read input: {e}")
            continue

        # Some 10x TIFFs are stored as (H, W, 3) instead of (3, H, W) — handle both
        if stack.ndim == 3 and stack.shape[-1] == 3 and stack.shape[0] != 3:
            stack = np.moveaxis(stack, -1, 0)
        if stack.ndim != 3 or stack.shape[0] < 3:
            print(f"  unexpected stack shape: {stack.shape}")
            continue

        rgb = make_rgb_for_display(stack, s.channel_order)

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

        gold_path = find_gold(s.slide_id)
        gold = tifffile.imread(gold_path).astype(np.uint8) if gold_path else None

        H = labels.shape[0]; W = labels.shape[1]
        if gold is not None:
            H = min(H, gold.shape[0]); W = min(W, gold.shape[1])
            gold = gold[:H, :W]
        H = min(H, rgb.shape[0]); W = min(W, rgb.shape[1])
        labels = labels[:H, :W]
        rgb_c = rgb[:H, :W]

        print("  class %  with_dcn  gold")
        for c in range(9):
            pp = 100*(labels == c).sum()/labels.size
            gp = 100*(gold == c).sum()/gold.size if gold is not None else 0.0
            if pp > 0.05 or gp > 0.05:
                print(f"    {LABEL_NAMES.get(c, str(c)):6s}  {pp:6.2f}    {gp:6.2f}")

        dcn_pct = 100 * (labels == 8).sum() / labels.size
        rows.append({
            "slide": s.slide_id, "labels": labels, "rgb": rgb_c,
            "gold": gold, "gold_name": gold_path.name if gold_path else None,
            "elapsed": elapsed, "dcn_pct": dcn_pct,
        })

        n_cols = 3 if gold is not None else 2
        fig, axes = plt.subplots(1, n_cols, figsize=(7*n_cols, 7))
        axes[0].imshow(rgb_c); axes[0].axis("off")
        axes[0].set_title(f"{s.slide_id}\nRGB (R=p27 G=NeuN B=DAPI)", fontsize=10)
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

    # Contact sheet
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
            "with-DCN MATLAB on 10x slides — magenta = DCN, "
            "★ = gold-corrected available",
            fontsize=14)
        fig.tight_layout()
        out_index = OUT / "_INDEX.png"
        fig.savefig(out_index, dpi=80)
        plt.close(fig)
        print(f"contact sheet: {out_index}")

    print("\n=== summary ===")
    print(f"{'slide':20s} {'time':>5} {'DCN%':>7}")
    for r in rows:
        print(f"  {r['slide']:18s} {r['elapsed']:>4.0f}s {r['dcn_pct']:>6.2f}")


if __name__ == "__main__":
    main()
