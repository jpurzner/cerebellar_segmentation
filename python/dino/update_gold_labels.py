"""For each slide with a hand-corrected labels TIFF in /python/labelled/,
re-compute the patch-grid labels in the chosen cache directory.

Backs up the original cache to <cache>_pre_gold/ before overwriting.
"""
from __future__ import annotations
import argparse, shutil
from pathlib import Path
import numpy as np
import tifffile

from slide_manifest import manifest

PATCH = 14
LABELLED = Path(__file__).resolve().parent.parent / "labelled"


def find_corrected(slide_id: str) -> Path | None:
    """Match any of: _corrected.tif, _labelled.tif*, _labeld*.tif*, etc."""
    cands = sorted(
        list(LABELLED.glob(f"{slide_id}_corrected.tif*")) +
        list(LABELLED.glob(f"{slide_id}_labelled.tif*")) +
        list(LABELLED.glob(f"{slide_id}_labelld.tif*")) +
        list(LABELLED.glob(f"{slide_id}_labeld.tif*")),
        key=lambda p: p.stat().st_mtime, reverse=True)
    # ignore SVG
    cands = [c for c in cands if c.suffix in (".tif", ".tiff")]
    return cands[0] if cands else None


def update_one_cache(cache_dir: Path, slide_id: str, corrected_path: Path):
    cache_path = cache_dir / f"{slide_id}_features.npz"
    if not cache_path.exists():
        print(f"  cache missing: {cache_path}")
        return False
    d = dict(np.load(cache_path, allow_pickle=True))
    H, W = int(d["H"]), int(d["W"])
    gy_total, gx_total = d["patch_grid"]
    gy_total, gx_total = int(gy_total), int(gx_total)

    corrected = tifffile.imread(corrected_path).astype(np.uint8)
    if corrected.shape != (H, W):
        # crop to common
        Hc = min(corrected.shape[0], H); Wc = min(corrected.shape[1], W)
        new_label = np.zeros((H, W), dtype=np.uint8)
        new_label[:Hc, :Wc] = corrected[:Hc, :Wc]
        corrected = new_label

    # Re-derive patch grid via majority vote, accounting for scale factor
    # if cache used target um/px different from native
    target_um = float(d.get("target_um_per_px", d.get("pixel_size_um", 0)))
    pixel_um = float(d.get("pixel_size_um", 0))
    if target_um == 0 or pixel_um == 0:
        scale = 1.0
    else:
        scale = pixel_um / target_um   # how much we resampled to get DINO input

    # patch i,j at DINO res maps to original-res pixels
    # patch_y in DINO = i*PATCH .. (i+1)*PATCH; in original: /scale
    # i.e. each patch covers PATCH/scale original pixels
    label_grid = np.zeros((gy_total, gx_total), dtype=np.uint8)
    n_classes = max(8, corrected.max() + 1)
    ylo = (np.arange(gy_total) * PATCH / scale).astype(int)
    yhi = ((np.arange(gy_total) + 1) * PATCH / scale).astype(int)
    xlo = (np.arange(gx_total) * PATCH / scale).astype(int)
    xhi = ((np.arange(gx_total) + 1) * PATCH / scale).astype(int)
    for i in range(gy_total):
        ya, yb = ylo[i], min(yhi[i], H)
        if ya >= yb: continue
        for j in range(gx_total):
            xa, xb = xlo[j], min(xhi[j], W)
            if xa >= xb: continue
            cell = corrected[ya:yb, xa:xb]
            if cell.size == 0: continue
            counts = np.bincount(cell.ravel(), minlength=n_classes)
            label_grid[i, j] = int(counts.argmax())

    d["labels"] = label_grid
    d["has_labels"] = True
    d["gold_label_source"] = corrected_path.name
    np.savez_compressed(cache_path, **d)
    print(f"  updated {cache_path.name}  scale={scale:.3f}  "
          f"labels classes: {sorted(set(label_grid.ravel()))}")
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", type=Path,
                    default=Path(__file__).resolve().parent / "cache_um0.5",
                    help="cache directory to update")
    ap.add_argument("--backup", action="store_true",
                    help="backup cache to <cache>_pre_gold/ before updating")
    args = ap.parse_args()

    if not args.cache.exists():
        print(f"cache dir not found: {args.cache}")
        return
    if args.backup:
        bk = args.cache.parent / f"{args.cache.name}_pre_gold"
        if not bk.exists():
            print(f"backing up {args.cache} -> {bk}")
            shutil.copytree(args.cache, bk)
        else:
            print(f"backup already exists: {bk}")

    slides = manifest()
    print(f"\nscanning for corrected labels in {LABELLED}/...")
    updated = []
    for s in slides:
        c = find_corrected(s.slide_id)
        if c is None:
            continue
        print(f"\n{s.slide_id}: corrected label = {c.name}")
        ok = update_one_cache(args.cache, s.slide_id, c)
        if ok: updated.append(s.slide_id)

    print(f"\nupdated {len(updated)} slides in {args.cache.name}: {updated}")


if __name__ == "__main__":
    main()
