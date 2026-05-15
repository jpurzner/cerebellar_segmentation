"""Export binary EGL boundary masks for each slide.

For each slide outputs:
  <slide>_egl_outer.tif   — EGL pixels touching background (PIA-side line)
  <slide>_egl_inner.tif   — EGL pixels touching ML/PCL/IGL/DWL (deep-side line)
  <slide>_egl_skel.tif    — 1-px-wide skeleton through the middle of the EGL ribbon

All masks are uint8 binary (0 or 255). Designed to feed into a model
that's learning the boundary of the cerebellum.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import tifffile
from scipy.ndimage import binary_dilation
from skimage.morphology import skeletonize

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from slide_manifest import manifest

ROOT = Path(__file__).resolve().parent
SAND = ROOT.parent / "matlab_anomaly" / "sandbox"
SAND_10X = ROOT.parent / "matlab_anomaly_10x" / "sandbox"
OUT = ROOT.parent / "egl_boundary_masks"
OUT.mkdir(exist_ok=True)


def find_latest_labels(slide_id, sandbox_dirs):
    """Find the most recently produced labels file across given sandbox dirs."""
    candidates = []
    for sb in sandbox_dirs:
        if not sb.exists(): continue
        candidates.extend(sb.glob(f"{slide_id}*_labels.tif"))
    if not candidates: return None
    return max(candidates, key=lambda p: p.stat().st_mtime)


def export_boundary_masks(slide_id, labels):
    """Compute and save outer/inner/skel masks for one slide."""
    egl_mask = (labels == 2) | (labels == 3)
    if egl_mask.sum() == 0:
        print(f"  {slide_id}: no EGL pixels, skipping"); return None

    bg_mask = (labels == 0)
    cereb_non_egl = (labels >= 4)   # IGL, ML, DWL, PCL, DCN, anomaly

    # Outer boundary: EGL pixels touching background (pia interface)
    outer = egl_mask & binary_dilation(bg_mask, iterations=1)
    # Inner boundary: EGL pixels touching deep tissue (ML/PCL/IGL/DWL/DCN)
    inner = egl_mask & binary_dilation(cereb_non_egl, iterations=1)
    # Skeleton: 1-px-wide centerline through the EGL ribbon
    skel = skeletonize(egl_mask)

    # Save as binary uint8 (0 / 255)
    tifffile.imwrite(OUT / f"{slide_id}_egl_outer.tif",
                      (outer * 255).astype(np.uint8), compression="zlib")
    tifffile.imwrite(OUT / f"{slide_id}_egl_inner.tif",
                      (inner * 255).astype(np.uint8), compression="zlib")
    tifffile.imwrite(OUT / f"{slide_id}_egl_skel.tif",
                      (skel * 255).astype(np.uint8), compression="zlib")
    return {"outer": int(outer.sum()), "inner": int(inner.sum()),
            "skel": int(skel.sum()), "shape": labels.shape}


def main():
    slides = manifest()
    print(f"exporting EGL boundary masks for {len(slides)} slides")

    for s in slides:
        # Try both 20x and 10x sandboxes
        labels_path = find_latest_labels(s.slide_id, [SAND, SAND_10X])
        if labels_path is None:
            print(f"  {s.slide_id}: no labels found, skip"); continue

        labels = tifffile.imread(labels_path).astype(np.uint8)
        info = export_boundary_masks(s.slide_id, labels)
        if info is None: continue
        print(f"  {s.slide_id} ({s.mag}, {info['shape']}): "
              f"outer={info['outer']} inner={info['inner']} skel={info['skel']} "
              f"(from {labels_path.name})")

    print(f"\noutput: {OUT}/")
    print("  *_egl_outer.tif  pia-side EGL boundary (1-px-wide)")
    print("  *_egl_inner.tif  deep-side EGL boundary")
    print("  *_egl_skel.tif   EGL centerline skeleton")


if __name__ == "__main__":
    main()
