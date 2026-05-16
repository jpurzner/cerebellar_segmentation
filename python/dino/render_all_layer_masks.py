"""Export per-layer binary masks for every slide.

For each slide, outputs one uint8 binary TIFF per anatomical layer:
  <slide>_iEGL.tif      label 2
  <slide>_oEGL.tif      label 3
  <slide>_EGL.tif       labels 2 | 3 (combined EGL ribbon)
  <slide>_IGL.tif       label 4
  <slide>_ML.tif        label 5
  <slide>_DWL.tif       label 6
  <slide>_PCL.tif       label 7
  <slide>_DCN.tif       label 8
  <slide>_anomaly.tif   label 9

All single-channel uint8 binary (0 / 255), same dimensions as source.
Designed to feed into a shape-based layer-completion model.
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
OUT = ROOT.parent / "layer_masks"
OUT.mkdir(exist_ok=True)

# Layer name -> list of integer labels that comprise that layer
# (EGL is the union of iEGL + oEGL)
LAYERS = {
    "iEGL":    [2],
    "oEGL":    [3],
    "EGL":     [2, 3],
    "IGL":     [4],
    "ML":      [5],
    "DWL":     [6],
    "PCL":     [7],
    "DCN":     [8],
    "anomaly": [9],
}


def find_latest_labels(slide_id, sandbox_dirs):
    """Find the most recently produced labels file across given sandbox dirs."""
    candidates = []
    for sb in sandbox_dirs:
        if not sb.exists(): continue
        candidates.extend(sb.glob(f"{slide_id}*_labels.tif"))
    if not candidates: return None
    return max(candidates, key=lambda p: p.stat().st_mtime)


def export_one_slide(slide_id, labels):
    """Export all layer binary masks for one slide. Returns dict of pixel counts."""
    counts = {}
    for layer_name, label_ids in LAYERS.items():
        mask = np.zeros(labels.shape, dtype=bool)
        for L in label_ids:
            mask |= (labels == L)
        out_path = OUT / f"{slide_id}_{layer_name}.tif"
        tifffile.imwrite(out_path, (mask * 255).astype(np.uint8),
                          compression="zlib")
        counts[layer_name] = int(mask.sum())
    return counts


def main():
    slides = manifest()
    print(f"exporting {len(LAYERS)} layer masks for {len(slides)} slides "
          f"({len(LAYERS) * len(slides)} TIFFs total)")

    summary_rows = []
    for s in slides:
        labels_path = find_latest_labels(s.slide_id, [SAND, SAND_10X])
        if labels_path is None:
            print(f"  {s.slide_id}: no labels, skip"); continue

        labels = tifffile.imread(labels_path).astype(np.uint8)
        counts = export_one_slide(s.slide_id, labels)
        total = labels.size

        # Brief per-slide summary
        line = f"  {s.slide_id} ({s.mag}, {labels.shape})"
        for name in ["EGL", "IGL", "ML", "DWL", "PCL", "DCN", "anomaly"]:
            pct = 100 * counts[name] / total
            line += f"  {name}={pct:.1f}%"
        print(line)
        summary_rows.append((s.slide_id, s.mag, total, counts))

    # Per-layer + per-slide summary
    print(f"\noutput: {OUT}/")
    print(f"   {len(summary_rows)} slides x {len(LAYERS)} layers = "
          f"{len(summary_rows) * len(LAYERS)} TIFFs")

    summary_path = OUT / "_SUMMARY.txt"
    with open(summary_path, "w") as f:
        f.write("Per-slide per-layer pixel counts (and % of image)\n")
        f.write("=" * 90 + "\n\n")
        for sid, mag, total, counts in summary_rows:
            f.write(f"=== {sid} ({mag}) total={total} ===\n")
            for name in ["iEGL", "oEGL", "EGL", "IGL", "ML",
                          "DWL", "PCL", "DCN", "anomaly"]:
                pct = 100 * counts[name] / total
                f.write(f"  {name:>8s}: {counts[name]:>10d} ({pct:>5.2f}%)\n")
            f.write("\n")
    print(f"summary: {summary_path}")


if __name__ == "__main__":
    main()
