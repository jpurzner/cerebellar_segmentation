"""Test the with-DCN MATLAB variant on s1_2 native + ds2x. Visualize results
showing the new DCN class.
"""
from __future__ import annotations
import subprocess
import time
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage as ndi
import tifffile
from scipy.ndimage import zoom

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from dataset import matlab_segments_to_labels
from slide_manifest import manifest

ROOT = Path(__file__).resolve().parent
LABELLED_DIR = ROOT.parent / "labelled"
OUT = ROOT.parent / "matlab_dcn_test"
OUT.mkdir(exist_ok=True)
PX_UM = 0.5119049
MATLAB = "/Applications/MATLAB_R2022b.app/bin/matlab"
REPO_ROOT = "/Users/jpurzner/Dropbox/imaging_analysis/cerebellar_segmentation"
KOVESI = "/Users/jpurzner/Dropbox/imaging_analysis/MatlabFns"

# 9-class colormap including DCN (label 8)
LABEL_COLORS = np.array([
    [0,   0,   0],   # 0 background
    [0,   0, 128],   # 1 all_cerebellum
    [0, 128, 255],   # 2 iEGL
    [0, 255, 255],   # 3 oEGL
    [255, 255, 0],   # 4 IGL
    [128, 255, 128], # 5 ML
    [255, 128, 0],   # 6 DWL
    [255,   0, 0],   # 7 PCL
    [255, 0, 255],   # 8 DCN — magenta (new)
], dtype=np.uint8)

LABEL_NAMES = {0:"bg",1:"cereb",2:"iEGL",3:"oEGL",4:"IGL",5:"ML",6:"DWL",7:"PCL",8:"DCN"}


def to_unit(im):
    im = im.astype(np.float32); lo, hi = im.min(), im.max()
    return (im - lo) / max(hi - lo, 1e-9)


def make_input(stack, kind, sandbox, base):
    if kind == "native":
        out = sandbox / f"{base}_native.tif"
        tifffile.imwrite(out, stack, photometric="minisblack")
        return out
    if kind == "ds2x":
        out_stack = np.zeros((3, stack.shape[1] // 2, stack.shape[2] // 2),
                             dtype=stack.dtype)
        for c in range(3):
            out_stack[c] = zoom(stack[c], 0.5, order=1)
        out = sandbox / f"{base}_ds2x.tif"
        tifffile.imwrite(out, out_stack, photometric="minisblack")
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


def main():
    slides = manifest()
    s1_2 = next(s for s in slides if s.slide_id == "2018_05_22_s1_2_p27-0021")
    sandbox = OUT / "sandbox"; sandbox.mkdir(exist_ok=True)

    stack = tifffile.imread(s1_2.input_path)
    cands = sorted(
        list(LABELLED_DIR.glob(f"{s1_2.slide_id}_corrected.tif*")) +
        list(LABELLED_DIR.glob(f"{s1_2.slide_id}_labelled.tif*")) +
        list(LABELLED_DIR.glob(f"{s1_2.slide_id}_labelld.tif*")),
        key=lambda p: p.stat().st_mtime, reverse=True)
    cands = [c for c in cands if c.suffix in (".tif", ".tiff")]
    gold = tifffile.imread(cands[0]).astype(np.uint8)

    rgb = np.stack([to_unit(stack[0]), to_unit(stack[1]), to_unit(stack[2])], axis=-1)
    rgb /= max(rgb.max(), 1e-9)

    rows = []
    base = s1_2.slide_id
    for kind in ["native", "ds2x"]:
        tag = f"{kind}_with_dcn"
        print(f"\n=== {tag} ===")
        input_path = make_input(stack, kind, sandbox, f"{base}_{tag}")
        t0 = time.time()
        rc, log = run_matlab(str(input_path), sandbox, tag)
        elapsed = time.time() - t0
        print(f"  MATLAB done in {elapsed:.0f}s, rc={rc}")
        # raw labels file
        labels_path = input_path.with_name(input_path.stem + "_labels.tif")
        if not labels_path.exists():
            print(f"  labels not found at {labels_path}")
            print(f"  see {log}")
            continue
        labels = tifffile.imread(labels_path).astype(np.uint8)
        print(f"  labels shape: {labels.shape}, classes: {sorted(set(labels.ravel()[::1000]))}")
        # If downsampled, upsample
        if labels.shape != stack.shape[1:]:
            from skimage.transform import resize
            labels = resize(labels, stack.shape[1:], order=0, preserve_range=True,
                             anti_aliasing=False).astype(np.uint8)
        H = min(labels.shape[0], gold.shape[0]); W = min(labels.shape[1], gold.shape[1])
        labels = labels[:H, :W]; gld = gold[:H, :W]; rgb_c = rgb[:H, :W]

        # Class breakdown
        print("  class %  with_dcn  gold")
        for c in range(9):
            pp = 100*(labels==c).sum()/labels.size
            gp = 100*(gld==c).sum()/gld.size
            if pp > 0.05 or gp > 0.05:
                print(f"    {LABEL_NAMES.get(c,c):6s}  {pp:6.2f}    {gp:6.2f}")

        rows.append({"tag": tag, "labels": labels, "elapsed": elapsed})

    if not rows:
        print("\nno successful runs"); return

    # Render: rows × (RGB | with_dcn output | gold)
    print("\nrendering panel...")
    fig, axes = plt.subplots(len(rows), 3, figsize=(21, len(rows)*7))
    if len(rows) == 1: axes = axes[None, :]
    for r_i, r in enumerate(rows):
        H = r["labels"].shape[0]; W = r["labels"].shape[1]
        axes[r_i, 0].imshow(rgb[:H, :W]); axes[r_i, 0].axis("off")
        axes[r_i, 0].set_title(f"{r['tag']}\nRGB", fontsize=11)
        axes[r_i, 1].imshow(LABEL_COLORS[r["labels"]]); axes[r_i, 1].axis("off")
        dcn_pct = 100*(r["labels"]==8).sum()/r["labels"].size
        axes[r_i, 1].set_title(f"with-DCN MATLAB ({r['elapsed']:.0f}s)\n"
                                f"DCN (magenta) = {dcn_pct:.2f}%", fontsize=10)
        axes[r_i, 2].imshow(LABEL_COLORS[gold[:H, :W]]); axes[r_i, 2].axis("off")
        axes[r_i, 2].set_title("GOLD (your hand-corrected, no DCN class)", fontsize=10)
    fig.suptitle("MATLAB with new DCN class — magenta is the deep cerebellar nuclei",
                  fontsize=14)
    fig.tight_layout()
    out_path = OUT / "s1_2_with_dcn.png"
    fig.savefig(out_path, dpi=85)
    plt.close(fig)
    print(f"saved {out_path}")


if __name__ == "__main__":
    main()
