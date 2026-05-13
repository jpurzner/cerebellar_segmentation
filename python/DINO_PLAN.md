# DINOv2-based cerebellar layer segmentation — implementation plan

## Goal
Replace the threshold-based EGL/iEGL/oEGL/IGL/ML/DWL detectors with a small
segmentation head trained on top of frozen **DINOv2** features. Aim for IoU
≥0.6 across all slides (vs current classical port mean 0.33 with σ ≈ 0.15).

## Why DINOv2
- Self-supervised on 142M natural images; feature representations transfer
  surprisingly well to histology / IF microscopy.
- "ViT-S/14" backbone gives **dense patch features at 14×14 px resolution**
  (small enough for thin EGL ribbons that are ~30-50 µm wide ≈ 60-100 px at 20×).
- Frozen-backbone + small linear head is a well-established, low-risk recipe
  (DINOv2 paper, Sec 7). No backbone fine-tuning needed for our scale.
- Runs CPU-only at inference (slow but feasible). Apple-Silicon MPS for speed.

## Architecture

```
Input (3, H, W)  uint16 microscopy stack [p27, NeuN, DAPI]
  │
  ├── Channel normalisation (CLAHE per channel) → float32 in [0,1]
  ├── DINO-ready RGB stack: clamp + scale to ImageNet stats
  │   We map (p27, NeuN, DAPI) → (R, G, B) and normalise the way ImageNet wants
  │
  ▼
DINOv2 ViT-S/14 (frozen) → patch tokens (N_patches, 384)
  │
  ▼
Per-patch features upsampled to original resolution (bilinear)
  │
  ▼
Small linear / MLP head (384 → 64 → 7 layer-classes)
  │
  ▼
Per-pixel logits → softmax → label map (8 classes inc. background)
```

Tile-and-stitch for memory: process the 5000×5000 px image in 1024×1024 px
tiles with 128 px overlap, blend logits across tiles.

## Training data

We have **14 20× slides** with MATLAB segmentation outputs (segments.tif). Plus
**4 10× slides**. Total **18 weak labels**.

- Use MATLAB segments as **weak supervision**. They're imperfect (s1_3, s1_4
  have catastrophic failures) but mostly correct.
- Filter labels: drop slides where MATLAB area fractions are pathological
  (e.g. <2% EGL or >50% DWL). Keep ~10 slides for training.
- Held-out: 2-3 slides for validation (NEVER in training).
- Optionally: augment with the **port v9 outputs** as a second weak label
  source. Train on regions where port and MATLAB AGREE (high-confidence) and
  weight UNCERTAIN regions less.

## Training loop

1. Load 1024×1024 tile + corresponding label tile.
2. Forward through DINOv2 frozen (cached for speed).
3. Linear head → per-pixel logits (8 classes).
4. **Loss = weighted cross-entropy + Dice**, with class weights inversely
   proportional to class frequency (EGL is rare → high weight).
5. Adam, lr 1e-3, batch_size 8, 50-100 epochs.
6. Validate every 5 epochs on held-out slides; stop on plateau.

Total compute estimate: ~2-4 hours on M1/M2 Mac (CPU/MPS).

## Evaluation

Same metric as we've been using:
- IoU vs MATLAB GT (with caveat that MATLAB itself has failures)
- Per-class IoU (iEGL, oEGL, IGL, ML, DWL, PCL)
- Visual 4-panel comparison plots like figs/v_all_14/

Cross-validation over slides (leave-one-out on held-out slides).

## Layer-quant compatibility

DINO output is per-pixel labels. Same as port output. So existing layer_quant
pipeline (1D row-wise profiles per linear cut) works unchanged.

## Risks & mitigations

| risk | mitigation |
|---|---|
| MATLAB labels are noisy on some slides | Filter slides; weight regions by port-vs-MATLAB agreement |
| 14 slides too few | Heavy augmentation (rotation, scaling, intensity jitter); use random crops as separate examples |
| Inherits MATLAB blind spots (deep folds) | Add radial-spoke-derived "EGL hint" as extra input channel; OR augment training labels with phasesym-extended EGL |
| MPS / DINOv2 setup issues | Fall back to CPU-only inference; pre-compute features once and cache |

## Concrete file layout

```
python/dino/
├── README.md                       (this plan + run instructions)
├── extract_features.py             (load image, run DINOv2 frozen, save .npz)
├── train_head.py                   (load cached features + labels, train linear head)
├── predict.py                      (run trained head on new image, save segmentation)
├── evaluate.py                     (compute IoU vs MATLAB GT, compare to port)
└── dataset.py                      (image loading, tile cropping, label preprocessing)
```

## Step-by-step execution

1. `pip install torch torchvision dinov2` (or use HuggingFace `dinov2-small`)
2. Run `extract_features.py` for all 14 slides → cached features (~50 GB).
   Re-run only if model changes.
3. Run `train_head.py` on 11 slides, holding out 3 → trained head (~10 MB).
4. Run `predict.py` on held-out 3 slides → segmentation outputs.
5. Run `evaluate.py` → IoU table + diff plots.
6. If IoU > current port mean (0.33), iterate on architecture/loss.
   If much higher (>0.6), declare success and integrate.

## Time / cost
- Setup + scaffolding: 2-3 hours
- Feature extraction (14 slides × 5000×5000 / 14×14 patches): ~30 min
- Training: 1-2 hours
- Evaluation: 30 min
- **Total first cut: ~half a working day**

## Open decisions
- ViT-S/14 (smaller, fast) vs ViT-B/14 (larger, more accurate)? Start with ViT-S.
- Linear head vs 2-layer MLP? Linear first — interpretable, fast, baseline.
- Use radial-spoke output as extra input channel? Probably yes (helps deep folds).
- Backbone fine-tuning? No — keeps things simple, avoids over-fitting to 14 slides.
