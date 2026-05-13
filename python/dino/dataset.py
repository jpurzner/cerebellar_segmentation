"""Slide loading + label preprocessing.

Loads:
  - Input TIF (3 channels: p27, NeuN, DAPI)
  - MATLAB segments TIF (RGB, decoded into 8-class label image)

Provides tile cropping for fitting the giant ~5000x5000 px images through
DINOv2 in memory-bounded chunks.
"""
from __future__ import annotations
from pathlib import Path
from typing import Iterator, Tuple
import numpy as np
import tifffile
from skimage import exposure

# Layer label IDs (matching cerebellum_pipeline.py set_bin)
# 0 = background, 1 = all_cerebellum, 2 = iEGL, 3 = oEGL, 4 = IGL,
# 5 = ML, 6 = DWL, 7 = PCL
N_CLASSES = 8

# MATLAB segments.tif RGB color → label ID, empirically calibrated jet(8) palette
_MATLAB_RGB_TO_LABEL = [
    ((0, 0, 0),       0),  # background
    ((0, 0, 128),     1),  # all_cerebellum (dark blue, jet level 1)
    ((0, 0, 255),     1),  # also seen as label 1 in some renders
    ((0, 128, 255),   2),  # iEGL (azure)
    ((0, 255, 255),   3),  # oEGL (cyan)
    ((255, 255, 0),   4),  # IGL (yellow)
    ((128, 255, 128), 5),  # ML (light green)
    ((255, 128, 0),   6),  # DWL (orange)
    ((255, 0, 0),     7),  # PCL (red)
]

PX_UM_20X = 0.5119049
PX_UM_10X = 1.0239


def matlab_segments_to_labels(seg_rgb: np.ndarray) -> np.ndarray:
    """Decode MATLAB jet-rendered segmentation back to label IDs."""
    s = seg_rgb.astype(int)
    H, W, _ = s.shape
    labels = np.zeros((H, W), dtype=np.uint8)
    for rgb, lbl in _MATLAB_RGB_TO_LABEL:
        match = np.all(np.abs(s - np.array(rgb)) < 16, axis=-1)
        labels[match] = lbl
    return labels


def load_slide(input_tif: Path,
                channel_order=("p27", "neun", "dapi")) -> np.ndarray:
    """Load 3-channel TIF, return (3, H, W) float32 in [0,1] after CLAHE.
    Output channels are reordered to (p27, neun, dapi) order regardless of input.
    """
    stack = tifffile.imread(input_tif)
    if stack.shape[0] not in (3, 4) and stack.shape[-1] in (3, 4):
        stack = np.moveaxis(stack, -1, 0)
    name_to_idx = {n: i for i, n in enumerate(channel_order)}

    def to_unit(im):
        im = im.astype(np.float32)
        lo, hi = im.min(), im.max()
        return (im - lo) / max(hi - lo, 1e-9)

    p27  = exposure.equalize_adapthist(to_unit(stack[name_to_idx["p27"]]),  clip_limit=0.01)
    neun = exposure.equalize_adapthist(to_unit(stack[name_to_idx["neun"]]), clip_limit=0.01)
    dapi = exposure.equalize_adapthist(to_unit(stack[name_to_idx["dapi"]]), clip_limit=0.01)
    return np.stack([p27, neun, dapi], axis=0).astype(np.float32)


def slide_to_dino_input(stack: np.ndarray) -> np.ndarray:
    """Convert (3, H, W) [p27, NeuN, DAPI] in [0,1] to a DINO-ready RGB float32
    tensor normalised with ImageNet stats."""
    rgb = stack.transpose(1, 2, 0)   # H, W, 3
    # ImageNet normalization
    mean = np.array([0.485, 0.456, 0.406], dtype=np.float32)
    std  = np.array([0.229, 0.224, 0.225], dtype=np.float32)
    rgb = (rgb - mean) / std
    return rgb.transpose(2, 0, 1).astype(np.float32)   # back to C, H, W


def tile_iter(H: int, W: int,
              tile: int = 1024, overlap: int = 128,
              snap_to: int = 14) -> Iterator[Tuple[int, int, int, int]]:
    """Yield (y0, y1, x0, x1) for tiles of size `tile` with given overlap.
    Coordinates are aligned to multiples of `snap_to` (DINO patch size = 14).
    """
    step = tile - overlap
    for y0 in range(0, H, step):
        y1 = min(H, y0 + tile)
        # snap y1 to multiple of snap_to (rounding down)
        y1_snap = max(y0 + snap_to, (y1 // snap_to) * snap_to)
        for x0 in range(0, W, step):
            x1 = min(W, x0 + tile)
            x1_snap = max(x0 + snap_to, (x1 // snap_to) * snap_to)
            yield (y0, y1_snap, x0, x1_snap)


def load_label_for_slide(seg_tif: Path) -> np.ndarray:
    """Load MATLAB segments.tif and decode into label ids."""
    seg = tifffile.imread(seg_tif)
    return matlab_segments_to_labels(seg)
