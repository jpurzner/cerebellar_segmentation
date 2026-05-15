"""Export full EGL binary masks (entire ribbon, not boundary or skeleton).

For each slide outputs:
  <slide>_egl_full.tif   — uint8 binary, 255 wherever set_bin == 2 (iEGL) or 3 (oEGL)

Single-channel, same dimensions as the source image. Binary 0/255.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import tifffile

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from slide_manifest import manifest

ROOT = Path(__file__).resolve().parent
SAND = ROOT.parent / "matlab_anomaly" / "sandbox"
SAND_10X = ROOT.parent / "matlab_anomaly_10x" / "sandbox"
OUT = ROOT.parent / "egl_full_masks"
OUT.mkdir(exist_ok=True)


def find_latest_labels(slide_id, sandbox_dirs):
    """Find the most recently produced labels file across given sandbox dirs."""
    candidates = []
    for sb in sandbox_dirs:
        if not sb.exists(): continue
        candidates.extend(sb.glob(f"{slide_id}*_labels.tif"))
    if not candidates: return None
    return max(candidates, key=lambda p: p.stat().st_mtime)


def main():
    slides = manifest()
    print(f"exporting full EGL masks for {len(slides)} slides")

    for s in slides:
        labels_path = find_latest_labels(s.slide_id, [SAND, SAND_10X])
        if labels_path is None:
            print(f"  {s.slide_id}: no labels found, skip"); continue

        labels = tifffile.imread(labels_path).astype(np.uint8)
        egl_full = (labels == 2) | (labels == 3)
        tifffile.imwrite(OUT / f"{s.slide_id}_egl_full.tif",
                          (egl_full * 255).astype(np.uint8), compression="zlib")
        n_egl = int(egl_full.sum())
        print(f"  {s.slide_id} ({s.mag}, {labels.shape}): "
              f"egl={n_egl} px ({100*n_egl/labels.size:.2f}% of image) "
              f"-> {s.slide_id}_egl_full.tif")

    print(f"\noutput: {OUT}/")


if __name__ == "__main__":
    main()
