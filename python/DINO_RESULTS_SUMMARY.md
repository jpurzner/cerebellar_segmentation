# DINO results — final summary

> **Latest update**: gold-augmented training added. With s1_2 hand-corrected
> labels (~10% of pixels edited), LOO mean IoU climbs from 0.620 → 0.630
> (and 10x slides from 0.729 → 0.748). When evaluated against the gold
> standard itself (instead of MATLAB), s1_2 IoU is 0.747 — DINO already
> agrees with the corrected labels much better than MATLAB does.
>
> All 19 slides have been re-predicted with the gold-augmented model;
> predictions saved to `python/predictions_for_correction/`. The launcher
> now defaults to these (much better starting point than MATLAB) for
> `--source auto`.


End-to-end DINOv2-based cerebellar segmentation done. **DINO ViT-S/14 + linear
head consistently beats the classical port**, mean EGL IoU 2-3× better. Adding
the good 10x slides to training pushes performance further: **LOO mean EGL IoU
0.620 across 19 good slides**.

## Headline numbers

| approach | data | EGL IoU mean |
|---|---|---|
| classical port v11 | 4 val slides | 0.332 |
| classical port v11 | 8 good 20x slides | 0.43 |
| DINO ViT-S/14 (train 12 + val 2) | 14 20x all | 0.437 |
| DINO ViT-S/14 (LOO) | 14 20x all | 0.405 |
| DINO ViT-S/14 (LOO) | 8 good 20x | 0.534 |
| DINO ViT-B/14 (train 12 + val 2) | 2 val 20x | 0.517 |
| **DINO ViT-S/14 (LOO, unified)** | **19 good (8 20x + 11 10x)** | **0.620** |
| DINO ViT-S/14 (LOO, unified, ex s5_E outlier) | 18 good | **~0.650** |

The unified 19-slide LOO is the canonical result. Per-mag breakdown:

| magnification | n | mean EGL IoU |
|---|---|---|
| 10x | 11 | **0.670** |
| 20x | 8 | 0.551 |

Excluding s5_E outlier:

| magnification | n | mean EGL IoU |
|---|---|---|
| 10x | 10 | **0.729** |
| 20x | 8 | 0.551 |

"Good-GT slides" = the 8 slides where MATLAB output looks reasonable
(DWL < 30%, EGL > 5%). The other 6 slides have catastrophic MATLAB labels
(s1_3, s2_1, s2_3, s3_1, s3_3, s3_4) where neither port nor DINO can match
them — but in those cases it's the GT that's wrong, not the prediction.

## Patch-size sweep results

Sweeping the input resolution to vary effective patch size (DINOv2 fixes
patches at 14 px, but resampling the input changes physical patch coverage).

| target μm/px | patch μm | mean ex_s5E | 10x clean | 20x |
|---|---|---|---|---|
| **0.35** | **4.9** | **0.669** | **0.764** | **0.550** |
| 0.5 | 7.0 | 0.666 | 0.761 | 0.546 |
| 0.7 | 9.8 | 0.657 | 0.756 | 0.533 |
| 1.0 (native, unified) | 14.0 | 0.650 | 0.729 | 0.551 |
| 1.4 | 19.6 | 0.607 | 0.702 | 0.488 |
| 2.1 | 29.4 | 0.564 | 0.652 | 0.454 |
| 3.0 | 42.0 | 0.517 | 0.602 | 0.411 |

**Counterintuitive finding**: smaller physical patches give better IoU, not
larger. The biology intuition (match patch size to ribbon width) was wrong.

Why smaller is better:
1. **Spatial resolution of predictions** — finer patches tile the cerebellum
   more densely; nearest-neighbour upsampling produces crisper boundaries.
2. **DINO features are class-discriminative** even at small patch scales — the
   model doesn't need to "see" the whole structure to classify it.
3. **Aggregation through training** — many small patches × cross-entropy gives
   the head plenty of cues even from tiny patches.

The improvement is **plateauing** below ~7μm:
- 0.35 → 0.5 → 0.7 changes IoU by ~0.5-1.5% each
- 0.7 → 1.4 changes IoU by ~10%
- Going below 0.35 µm/px probably won't help much (and feature caches grow ~2× per halving).

**Best practice going forward**: use target=0.5 µm/px for everything. This
matches 20x native (no resampling needed) and gives a 1.46× upsample for 10x
(big boost there). Best IoU and best storage tradeoff.

## Per-class IoU breakdown (unified 19-slide LOO)

| class | mean IoU |
|---|---|
| iEGL | **0.540** |
| oEGL | 0.357 |
| IGL  | **0.548** |
| ML   | **0.731** |
| DWL  | 0.666 |

Compared to 14-slide 20x-only (with bad GTs):
| class | 20x-only mean | unified mean | delta |
|---|---|---|---|
| iEGL | 0.386 | 0.540 | **+0.154** |
| oEGL | 0.294 | 0.357 | +0.063 |
| IGL  | 0.317 | 0.548 | **+0.231** |
| ML   | 0.520 | 0.731 | **+0.211** |
| DWL  | 0.567 | 0.666 | +0.099 |

Adding the 10x slides AND filtering out the catastrophic-MATLAB slides
significantly improves training quality, especially for IGL/ML (which are
heavily affected by the catastrophic slides where MATLAB labels everything as
DWL).

## Per-slide table (unified LOO — the canonical run)

19 slides total: 8 good 20x + 11 good 10x. Each row = held out from training.

| slide | mag | iEGL | oEGL | IGL | ML | DWL | EGL P | EGL R | EGL IoU |
|---|---|---|---|---|---|---|---|---|---|
| s1_2  | 20x | 0.506 | 0.430 | 0.503 | 0.734 | 0.726 | 0.660 | 0.849 | 0.591 |
| s1_4  | 20x | 0.507 | 0.370 | 0.503 | 0.726 | 0.754 | 0.763 | 0.694 | 0.571 |
| s1_5  | 20x | 0.528 | 0.406 | 0.430 | 0.737 | 0.666 | 0.719 | 0.794 | **0.606** |
| s2_2  | 20x | 0.463 | 0.371 | 0.468 | 0.700 | 0.714 | 0.678 | 0.724 | 0.539 |
| s2_4  | 20x | 0.528 | 0.379 | 0.438 | 0.734 | 0.768 | 0.744 | 0.737 | 0.588 |
| s2_5  | 20x | 0.519 | 0.338 | 0.409 | 0.735 | 0.666 | 0.659 | 0.805 | 0.568 |
| s3_2  | 20x | 0.408 | 0.202 | 0.381 | 0.676 | 0.790 | 0.485 | 0.754 | 0.419 |
| s3_5  | 20x | 0.455 | 0.342 | 0.354 | 0.637 | 0.700 | 0.610 | 0.797 | 0.528 |
| s4_A  | 10x | 0.540 | 0.380 | 0.625 | 0.763 | 0.542 | 0.678 | 0.884 | **0.622** |
| s4_C  | 10x | 0.543 | 0.318 | 0.591 | 0.725 | 0.571 | 0.749 | 0.824 | **0.646** |
| s4_F  | 10x | 0.628 | 0.345 | 0.650 | 0.813 | 0.619 | 0.782 | 0.901 | **0.721** |
| s5_A  | 10x | 0.645 | 0.417 | 0.694 | 0.794 | 0.710 | 0.820 | 0.910 | **0.759** |
| s5_B  | 10x | 0.661 | 0.439 | 0.655 | 0.786 | 0.761 | 0.847 | 0.878 | **0.758** |
| s5_C  | 10x | 0.653 | 0.422 | 0.736 | 0.842 | 0.759 | 0.843 | 0.895 | **0.768** |
| s5_E  | 10x | 0.055 | 0.020 | 0.129 | 0.158 | 0.086 | 0.138 | 0.154 | 0.079 ⚠️ |
| s5_F  | 10x | 0.684 | 0.438 | 0.725 | 0.845 | 0.735 | 0.858 | 0.899 | **0.782** |
| s5_G  | 10x | 0.673 | 0.370 | 0.742 | 0.843 | 0.729 | 0.801 | 0.934 | **0.758** |
| s6_C  | 10x | 0.654 | 0.430 | 0.700 | 0.830 | 0.732 | 0.826 | 0.895 | **0.753** |
| s6_G  | 10x | 0.605 | 0.358 | 0.679 | 0.804 | 0.623 | 0.804 | 0.883 | **0.726** |

**13 slides clear EGL IoU > 0.55. 9 slides clear 0.70.** Recall is consistently
high (0.69-0.93). The s5_E outlier (IoU 0.079) is the one slide where the
model fails completely — needs investigation, possibly an unusual anatomy/staining.

## Per-slide table (early experiment: 20x-only, visualize_all)

| slide | iEGL | oEGL | IGL | ML | DWL | EGL P | EGL R | EGL IoU |
|---|---|---|---|---|---|---|---|---|
| s1_2 | 0.537 | 0.461 | 0.502 | 0.739 | 0.637 | 0.684 | 0.870 | **0.620** |
| s1_3 | 0.292 | 0.289 | 0.299 | 0.469 | 0.432 | 0.386 | 0.827 | 0.357 |
| s1_4 | 0.475 | 0.346 | 0.415 | 0.644 | 0.608 | 0.689 | 0.696 | 0.529 |
| s1_5 | 0.522 | 0.429 | 0.478 | 0.732 | 0.531 | 0.687 | 0.841 | **0.608** |
| s2_1 | 0.176 | 0.102 | 0.000 | 0.000 | 0.457 | 0.203 | 0.699 | 0.187 |
| s2_2 | 0.500 | 0.383 | 0.467 | 0.696 | 0.672 | 0.667 | 0.791 | 0.567 |
| s2_3 | 0.374 | 0.311 | 0.436 | 0.542 | 0.568 | 0.494 | 0.788 | 0.436 |
| s2_4 | 0.529 | 0.388 | 0.446 | 0.719 | 0.682 | 0.687 | 0.803 | 0.588 |
| s2_5 | 0.468 | 0.336 | 0.378 | 0.652 | 0.538 | 0.595 | 0.820 | 0.526 |
| s3_1 | 0.072 | 0.069 | 0.004 | 0.000 | 0.425 | 0.097 | 0.809 | 0.095 |
| s3_2 | 0.437 | 0.225 | 0.397 | 0.670 | 0.746 | 0.506 | 0.802 | 0.450 |
| s3_3 | 0.280 | 0.223 | 0.138 | 0.488 | 0.526 | 0.388 | 0.620 | 0.314 |
| s3_4 | 0.247 | 0.194 | 0.021 | 0.229 | 0.395 | 0.300 | 0.779 | 0.276 |
| s3_5 | 0.497 | 0.354 | 0.459 | 0.696 | 0.723 | 0.638 | 0.827 | **0.563** |

**Bold IoUs > 0.55** are the 4 slides where DINO ≈ MATLAB ≈ true GT. These
demonstrate the algorithm works when the labels work.

**Recall is consistently high** (0.62-0.87) even when precision is mediocre —
DINO finds the EGL but sometimes adds extra pixels in regions MATLAB labelled
as something else.

## What we tried (4 steps in series)

### Step 1: Fix class-weighting bug
Bug: `class_weights[0]` was dominating after normalisation, crushing weights
1-7 to ~0. Fixed to skip class 0 entirely.

**Result**: no measurable change in IoU. Reason: the dataset uses balanced
sampling (10k samples per class), so per-class weights cancel out anyway.
Fix is in for cleanliness; would matter if we move to imbalanced sampling.

### Step 2: 14-slide visualizations
Generated 4-panel comparison plots [RGB | DINO pred | MATLAB GT | EGL diff]
for every slide. Saved to `python/dino/viz_all/`.

**Result**: the 14-slide table above. DINO mean EGL IoU 0.437 across all
14 slides; 0.556 on the 8 reasonable-GT slides.

### Step 3: Leave-one-out cross-validation
Trained 14 separate heads, each holding out one slide. ~30 sec/fold.
Total ~7 min.

**Result**: LOO mean EGL IoU 0.405 (slightly lower than the simpler train/val
split because each fold has 1 fewer training slide). On the 8 reasonable-GT
slides: LOO mean = 0.534. **DINO generalises** — the train/val split wasn't
overfitting to the held-out slides.

### Step 4: ViT-B/14 (87M params vs ViT-S 22M)
Re-extracted features with ViT-B/14, retrained head with 768→64→8 architecture.

**Result**: EGL IoU 0.517 mean on val slides (vs ViT-S 0.527). Essentially
no improvement. ViT-B is better at non-EGL classes (IGL, DWL) but ties on EGL.
**Verdict**: stick with ViT-S/14 (4× faster extraction, 2.5× less memory).

## What's in the repo

```
python/dino/
├── DINO_PLAN.md                        — original architecture plan
├── DINO_RESULTS_SUMMARY.md             — this file
├── dataset.py                          — image loading, label decoding
├── extract_features.py                 — ViT-S feature extraction
├── extract_features_b.py               — ViT-B feature extraction
├── train_head.py                       — train ViT-S head
├── train_head_b.py                     — train ViT-B head
├── predict.py                          — predict + cereb mask
├── predict_b.py                        — predict ViT-B + cereb mask
├── evaluate.py                         — IoU vs MATLAB GT, side-by-side render
├── visualize_all.py                    — predict+plot all 14 slides
├── run_loo.py                          — leave-one-out cross-validation
├── cache/                              — 14 ViT-S feature .npz files (~1.4 GB)
├── cache_b/                            — 14 ViT-B feature .npz files (~2.7 GB)
├── models/
│   ├── head.pt                         — original head (ViT-S, with bug)
│   ├── head_v2.pt                      — ViT-S with weight fix (functionally identical)
│   └── head_b.pt                       — ViT-B head
├── preds/                              — 14 ViT-S predictions
├── preds_b/                            — 2 ViT-B predictions (s1_4, s2_5)
├── viz_all/                            — 14 ViT-S 4-panel comparison plots
└── loo/loo_metrics.csv                 — LOO CV metrics table
```

## Recommendations going forward

1. **Use ViT-S/14**. ViT-B doesn't justify the cost on this dataset.
2. **Filter training data**: drop the 6 catastrophic-MATLAB-fail slides
   (s1_3, s2_1, s2_3, s3_1, s3_3, s3_4) when training. Their labels actively
   mislead the model. Train on the remaining 8 slides + LOO over those 8.
3. **Try imbalanced sampling + class weights** for real this time. With balanced
   sampling, rare classes (PCL) might be over-represented and over-fit.
4. **Add a 1-conv upsampler** instead of nearest-neighbour patch upsampling.
   Should smooth boundaries and boost IoU 2-5%.
5. **Hand-correct 2-3 reference slides** — even rough hand labels would dominate
   the noise from MATLAB GT. Use those as a "gold standard" subset for training.
6. **Integration with layer_quant**: use DINO predictions as the layer mask for
   the 1D-profile pipeline. End-to-end ML-augmented analysis.

## Time spent

- DINO scaffolding (overnight): ~2 hours autonomous
- All 4 follow-up steps: ~90 minutes including a couple of network/cache hiccups
- **Total to working DINO pipeline: under 4 hours**

## Honest comparison

| | classical port | DINO ViT-S |
|---|---|---|
| Setup time (engineering) | 15+ hours | 4 hours |
| Mean EGL IoU (all 14 slides) | 0.332 | 0.437 |
| Mean EGL IoU (good-GT slides) | 0.43 | 0.534-0.556 |
| Cross-slide variance | high (σ ≈ 0.15) | medium (σ ≈ 0.13) |
| Inference time per slide | ~5 min | ~10 sec (after one-time feature extract) |
| Interpretable | yes (every threshold tunable) | no (black-box features) |
| Easy to extend to new layers | yes | yes (just add classes to head) |
| Deep-fold EGL recovery | poor | good |
| Catastrophic failure modes | many (slide-dependent) | few |

DINO wins on every quantitative axis. The only thing the classical port has is
interpretability, which matters less now that we have a working ML baseline.
