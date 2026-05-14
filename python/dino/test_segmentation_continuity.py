"""Diagnostic continuity tests for segmentation outputs.

Reports per-slide topology metrics and flags problematic segmentations.
Doesn't modify any labels — purely a post-hoc validation tool.

Tests per layer (EGL = iEGL+oEGL, ML, IGL, DWL):
  - Number of connected components
  - Area of largest component
  - Fraction of total layer area in largest component (continuity score)
  - Number of large isolated components (>500 px²)

Failure flags:
  - EGL: >5 components OR <80% in largest (should be near-continuous ribbon)
  - ML: >5 components OR <80% in largest
  - IGL: >5 components OR <80% in largest
  - DWL: >40 isolated branches (some branching is normal — it's a tree)

Outputs:
  /python/test_continuity/SUMMARY.txt        per-slide flag report
  /python/test_continuity/<slide>_metrics.txt  per-slide detailed metrics
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import tifffile
from skimage import measure

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from slide_manifest import manifest

ROOT = Path(__file__).resolve().parent
SAND = ROOT.parent / "matlab_anomaly" / "sandbox"
OUT = ROOT.parent / "test_continuity"
OUT.mkdir(exist_ok=True)

# Layer label mapping
LAYERS = {
    "EGL":  [2, 3],   # iEGL + oEGL together
    "ML":   [5],
    "IGL":  [4],
    "DWL":  [6],
}


def layer_mask(labels, layer_label_list):
    """Build binary mask from one or more label values."""
    m = np.zeros(labels.shape, dtype=bool)
    for L in layer_label_list:
        m |= (labels == L)
    return m


def layer_metrics(labels, layer_name, label_list):
    """Compute connectivity metrics for one layer."""
    mask = layer_mask(labels, label_list)
    total_area = mask.sum()
    if total_area == 0:
        return {"n_components": 0, "total_area": 0, "largest_area": 0,
                "largest_pct": 0.0, "n_large_components": 0,
                "components_by_size": []}
    cc = measure.label(mask, connectivity=2)
    props = measure.regionprops(cc)
    areas = sorted([p.area for p in props], reverse=True)
    largest = areas[0] if areas else 0
    largest_pct = 100 * largest / total_area
    n_large = sum(1 for a in areas if a >= 500)
    return {
        "n_components": len(areas),
        "total_area": total_area,
        "largest_area": largest,
        "largest_pct": largest_pct,
        "n_large_components": n_large,
        "components_by_size": areas[:10],   # top 10 for inspection
    }


def flag_failures(layer_name, m):
    """Return list of failure flag strings for a layer."""
    flags = []
    if layer_name in ("EGL", "ML", "IGL"):
        if m["n_large_components"] > 5:
            flags.append(f"{layer_name}_FRAGMENTED({m['n_large_components']} components)")
        if m["largest_pct"] < 80 and m["total_area"] > 10000:
            flags.append(f"{layer_name}_DISCONTINUOUS(largest={m['largest_pct']:.0f}%)")
    elif layer_name == "DWL":
        # DWL is naturally a branching tree — many components are normal
        # but >40 large isolated branches suggests over-detection
        if m["n_large_components"] > 40:
            flags.append(f"DWL_OVER-FRAGMENTED({m['n_large_components']} large components)")
    return flags


def find_latest_labels(slide_id):
    """Find the most recently produced labels file for a slide."""
    candidates = sorted(SAND.glob(f"{slide_id}*_labels.tif"),
                        key=lambda p: p.stat().st_mtime, reverse=True)
    return candidates[0] if candidates else None


def main():
    slides = manifest()
    twenty = [s for s in slides if s.mag == "20x"]
    print(f"running continuity tests on {len(twenty)} 20x slides")

    summary_lines = ["Per-slide layer continuity report\n",
                     "=" * 80 + "\n\n"]

    for s in twenty:
        labels_path = find_latest_labels(s.slide_id)
        if labels_path is None:
            print(f"  {s.slide_id}: no labels file found, skip")
            continue

        labels = tifffile.imread(labels_path).astype(np.uint8)
        per_layer = {}
        all_flags = []

        for layer_name, label_list in LAYERS.items():
            m = layer_metrics(labels, layer_name, label_list)
            per_layer[layer_name] = m
            all_flags.extend(flag_failures(layer_name, m))

        # Print to console
        flag_str = "OK" if not all_flags else " | ".join(all_flags)
        print(f"\n=== {s.slide_id} ({labels_path.name}) ===")
        print(f"  status: {flag_str}")
        print(f"  {'layer':>5} {'n_comp':>7} {'large':>6} {'largest':>8} {'largest%':>9}")
        for layer_name, m in per_layer.items():
            print(f"  {layer_name:>5} {m['n_components']:>7} "
                  f"{m['n_large_components']:>6} {m['largest_area']:>8} "
                  f"{m['largest_pct']:>8.1f}%")

        # Per-slide detailed report
        with open(OUT / f"{s.slide_id}_metrics.txt", "w") as f:
            f.write(f"Continuity metrics for {s.slide_id}\n")
            f.write(f"Source: {labels_path.name}\n")
            f.write(f"Status: {flag_str}\n\n")
            for layer_name, m in per_layer.items():
                f.write(f"--- {layer_name} ---\n")
                f.write(f"  total area: {m['total_area']} px\n")
                f.write(f"  components: {m['n_components']} ({m['n_large_components']} >= 500 px)\n")
                f.write(f"  largest:    {m['largest_area']} px ({m['largest_pct']:.1f}% of total)\n")
                f.write(f"  top 10 sizes: {m['components_by_size']}\n\n")

        # Summary
        summary_lines.append(f"{s.slide_id}: {flag_str}\n")
        for layer_name, m in per_layer.items():
            summary_lines.append(
                f"  {layer_name}: n={m['n_components']} large={m['n_large_components']} "
                f"largest_pct={m['largest_pct']:.1f}%\n")
        summary_lines.append("\n")

    # Write summary
    with open(OUT / "SUMMARY.txt", "w") as f:
        f.writelines(summary_lines)
    print(f"\nsummary: {OUT / 'SUMMARY.txt'}")
    print(f"per-slide details: {OUT}/<slide>_metrics.txt")


if __name__ == "__main__":
    main()
