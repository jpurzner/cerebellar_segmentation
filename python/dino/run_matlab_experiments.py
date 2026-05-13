"""Test MATLAB on input variants and tightened iEGL parameter for s1_2.

Generates input variants (downsampled, blurred), then runs actual MATLAB on
each, then renders comparison panels vs gold standard.

Variants tested:
  A) native + default MATLAB params
  B) downsampled 2× + default MATLAB params
  C) blurred σ=4 µm + default MATLAB params
  D) native + tight-iEGL MATLAB variant (sens 0.3→0.15, thresh 0.2→0.35)
"""
from __future__ import annotations
import shutil
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
OUT = ROOT.parent / "matlab_experiments"
OUT.mkdir(exist_ok=True)
PX_UM = 0.5119049
MATLAB = "/Applications/MATLAB_R2022b.app/bin/matlab"
REPO_ROOT = "/Users/jpurzner/Dropbox/imaging_analysis/cerebellar_segmentation"
KOVESI = "/Users/jpurzner/Dropbox/imaging_analysis/MatlabFns"

LABEL_COLORS = np.array([
    [0,   0,   0],   [0,   0, 128], [0, 128, 255], [0, 255, 255],
    [255, 255, 0], [128, 255, 128],[255, 128, 0], [255,   0, 0],
], dtype=np.uint8)


def to_unit(im):
    im = im.astype(np.float32); lo, hi = im.min(), im.max()
    return (im - lo) / max(hi - lo, 1e-9)


def egl_iou(pred, gt):
    pe = (pred == 2) | (pred == 3); ge = (gt == 2) | (gt == 3)
    inter = (pe & ge).sum(); union = (pe | ge).sum()
    return inter / max(pe.sum(), 1), inter / max(ge.sum(), 1), inter / max(union, 1)


def make_input_variant(stack, kind, sandbox_dir, base_name):
    """Save a TIFF variant + return its path."""
    if kind == "native":
        out = sandbox_dir / f"{base_name}_native.tif"
        tifffile.imwrite(out, stack, photometric="minisblack")
        return out
    if kind == "ds2x":
        out_stack = np.zeros((3, stack.shape[1] // 2, stack.shape[2] // 2),
                             dtype=stack.dtype)
        for c in range(3):
            out_stack[c] = zoom(stack[c], 0.5, order=1)
        out = sandbox_dir / f"{base_name}_ds2x.tif"
        tifffile.imwrite(out, out_stack, photometric="minisblack")
        return out
    if kind.startswith("blur"):
        sig_um = float(kind.replace("blur", ""))
        sig_px = sig_um / PX_UM
        out_stack = np.zeros_like(stack)
        for c in range(3):
            out_stack[c] = ndi.gaussian_filter(stack[c].astype(np.float32),
                                               sigma=sig_px).astype(stack.dtype)
        out = sandbox_dir / f"{base_name}_{kind}.tif"
        tifffile.imwrite(out, out_stack, photometric="minisblack")
        return out
    raise ValueError(kind)


def run_matlab(input_tif, fn_name, log_path, sandbox_dir, tag):
    """Run MATLAB on input_tif using fn_name segment function.
    Driver scripts must have valid MATLAB identifier filenames (no hyphens,
    underscores OK but no leading underscore)."""
    # Make sandbox subdirectory exist + use a CLEAN filename
    safe_tag = tag.replace("-", "_").replace(".", "_")
    script_name = f"driver_{safe_tag}.m"
    script_path = sandbox_dir / script_name
    matlab_code = f"""addpath('{REPO_ROOT}');
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
    layer = {fn_name}(img);
    fprintf('OK in %.1fs\\n', toc(t0));
catch ME
    fprintf('FAIL: %s\\n', ME.message);
    for k = 1:length(ME.stack)
        fprintf('  at %s line %d\\n', ME.stack(k).name, ME.stack(k).line);
    end
end
"""
    script_path.write_text(matlab_code)
    # `run` likes to be in the same dir as the script. cd first.
    cmd = [MATLAB, "-batch",
           f"cd('{sandbox_dir}'); run('{script_name}')"]
    with open(log_path, "w") as f:
        proc = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, text=True)
    return proc.returncode


def main():
    slides = manifest()
    s1_2 = next(s for s in slides if s.slide_id == "2018_05_22_s1_2_p27-0021")

    sandbox = OUT / "sandbox"
    sandbox.mkdir(exist_ok=True)

    print(f"loading s1_2: {s1_2.input_path}")
    stack = tifffile.imread(s1_2.input_path)
    print(f"  shape {stack.shape}")

    matlab_gt = matlab_segments_to_labels(tifffile.imread(s1_2.gt_path))

    cands = sorted(
        list(LABELLED_DIR.glob(f"{s1_2.slide_id}_corrected.tif*")) +
        list(LABELLED_DIR.glob(f"{s1_2.slide_id}_labelled.tif*")) +
        list(LABELLED_DIR.glob(f"{s1_2.slide_id}_labelld.tif*")),
        key=lambda p: p.stat().st_mtime, reverse=True)
    cands = [c for c in cands if c.suffix in (".tif", ".tiff")]
    gold = tifffile.imread(cands[0]).astype(np.uint8)
    print(f"gold: {cands[0].name}")

    # Define experiments
    experiments = [
        ("A_native_default",   "native", "cerebellum_threshold_segment20x"),
        ("B_ds2x_default",     "ds2x",   "cerebellum_threshold_segment20x"),
        ("C_blur4um_default",  "blur4",  "cerebellum_threshold_segment20x"),
        ("D_native_tightiegl", "native", "cerebellum_threshold_segment20x_tight_iegl"),
        ("E_blur4um_tightiegl","blur4",  "cerebellum_threshold_segment20x_tight_iegl"),
    ]

    rows = []
    base = s1_2.slide_id
    for tag, input_kind, fn_name in experiments:
        print(f"\n=== {tag} ===")
        # generate input variant
        input_path = make_input_variant(stack, input_kind, sandbox, f"{base}_{tag}")
        print(f"  input variant: {input_path.name}")

        log_path = OUT / f"{tag}.log"
        t0 = time.time()
        rc = run_matlab(str(input_path), fn_name, log_path, sandbox, tag)
        elapsed = time.time() - t0
        print(f"  MATLAB finished in {elapsed:.0f}s, rc={rc}")

        # MATLAB writes outputs to the same directory as the input
        seg_path = input_path.with_name(input_path.stem + "_segments.tif")
        if not seg_path.exists():
            print(f"  segments not found at {seg_path}")
            with open(log_path) as f:
                print("\n--- log tail ---")
                print(f.read()[-1500:])
            continue
        seg_rgb = tifffile.imread(seg_path)
        seg = matlab_segments_to_labels(seg_rgb)

        # If downsampled, upsample back to full res for comparison
        if seg.shape != matlab_gt.shape:
            from skimage.transform import resize
            seg = resize(seg, matlab_gt.shape, order=0, preserve_range=True,
                          anti_aliasing=False).astype(np.uint8)
        H = min(seg.shape[0], matlab_gt.shape[0], gold.shape[0])
        W = min(seg.shape[1], matlab_gt.shape[1], gold.shape[1])
        seg = seg[:H, :W]; mat = matlab_gt[:H, :W]; gld = gold[:H, :W]

        Pm, Rm, IoUm = egl_iou(seg, mat)
        Pg, Rg, IoUg = egl_iou(seg, gld)
        print(f"  EGL IoU vs MATLAB GT: {IoUm:.3f}")
        print(f"  EGL IoU vs GOLD     : {IoUg:.3f}")

        # per-class
        names = {2:"iEGL",3:"oEGL",4:"IGL",5:"ML",6:"DWL",7:"PCL"}
        for c in [2,3,4,5,6,7]:
            pp = 100*(seg==c).sum()/seg.size
            mp = 100*(mat==c).sum()/mat.size
            gp = 100*(gld==c).sum()/gld.size
            print(f"    {names[c]:5s}  pred={pp:5.2f}  matlab={mp:5.2f}  gold={gp:5.2f}")

        rows.append({"tag": tag, "input": input_kind, "fn": fn_name,
                     "pred": seg, "matlab": mat, "gold": gld,
                     "iou_matlab": IoUm, "iou_gold": IoUg,
                     "elapsed": elapsed})

    if not rows:
        print("no successful runs"); return

    # Render combined panel
    print("\nrendering combined panel...")
    rgb = np.stack([to_unit(stack[0]), to_unit(stack[1]), to_unit(stack[2])], axis=-1)
    H = rows[0]["pred"].shape[0]; W = rows[0]["pred"].shape[1]
    rgb = rgb[:H, :W]
    fig, axes = plt.subplots(len(rows) + 1, 4, figsize=(28, (len(rows)+1)*7))
    # top row: RGB | empty | MATLAB GT | GOLD
    axes[0, 0].imshow(rgb); axes[0, 0].axis("off"); axes[0, 0].set_title("RGB", fontsize=12)
    axes[0, 1].axis("off"); axes[0, 1].text(0.5, 0.5, "(reference row)", ha="center")
    axes[0, 2].imshow(LABEL_COLORS[rows[0]["matlab"]]); axes[0, 2].axis("off")
    axes[0, 2].set_title("MATLAB GT (reference)", fontsize=12)
    axes[0, 3].imshow(LABEL_COLORS[rows[0]["gold"]]); axes[0, 3].axis("off")
    axes[0, 3].set_title("GOLD (your hand-corrected)", fontsize=12)
    for r_i, r in enumerate(rows):
        ax = axes[r_i + 1]
        ax[0].imshow(rgb); ax[0].axis("off"); ax[0].set_title(r["tag"], fontsize=11)
        ax[1].imshow(LABEL_COLORS[r["pred"]]); ax[1].axis("off")
        ax[1].set_title(f"MATLAB out\nIoU(M)={r['iou_matlab']:.2f}  IoU(G)={r['iou_gold']:.2f}",
                         fontsize=10)
        ax[2].imshow(LABEL_COLORS[r["matlab"]]); ax[2].axis("off")
        ax[2].set_title("MATLAB GT", fontsize=10)
        ax[3].imshow(LABEL_COLORS[r["gold"]]); ax[3].axis("off")
        ax[3].set_title("GOLD", fontsize=10)
    fig.suptitle("MATLAB experiments on s1_2: input variants × iEGL params",
                  fontsize=14)
    fig.tight_layout()
    panel_path = OUT / "s1_2_matlab_experiments.png"
    fig.savefig(panel_path, dpi=80)
    plt.close(fig)
    print(f"\nsaved {panel_path}")

    print("\n=== summary ===")
    print(f"{'tag':22s} {'time':>6} {'IoU(MATLAB)':>12} {'IoU(GOLD)':>11}")
    for r in rows:
        print(f"  {r['tag']:20s} {r['elapsed']:>5.0f}s "
              f"{r['iou_matlab']:>12.3f} {r['iou_gold']:>11.3f}")


if __name__ == "__main__":
    main()
