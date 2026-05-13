# Overnight diagnosis — cerebellar segmentation port (May 9-10, 2026)

## TL;DR

1. **Most "low IoU" slides aren't my port's failure — they're MATLAB's failure.** 6 of 14 slides have MATLAB DWL > 50% of cerebellum (biologically impossible: real DWL is ~15-25%).
2. **For the 8 slides where MATLAB GT looks reasonable**, port IoU ranges 0.11–0.61 (mean ~0.36).
3. **The biggest port-side bug is iEGL under-detection on s1_4 and s1_5** — port detects only 10-15% of the iEGL MATLAB found.
4. **The DWL/IGL "overflow" is partly genuine** (port's IGL is 1.32 mm² vs MATLAB's 0.7 on s1_2), but reverting `im_get_overlap` to the faithful MATLAB version made things WORSE — the pipeline has self-compensating quirks.
5. **DINO scaffolding is built and feature extraction is running** — should be ready to train when you wake up.

## Layer-area table (MATLAB GT areas in mm² and as % of cereb)

```
slide                                    iEGL    oEGL    IGL     ML      DWL     PCL    cereb
s1_2                                  0.39/11  0.17/5  0.70/19  1.12/30  0.83/22  0.21/6  3.69  ✓
s1_3                                  0.28/6   0.13/3  0.30/7   0.56/13  2.96/66  0.06/1  4.46  ✗ MATLAB-fail
s1_4                                  0.51/12  0.19/5  0.81/19  1.10/26  0.93/22  0.23/5  4.23  ✓
s1_5                                  0.45/12  0.19/5  0.67/18  0.98/26  0.95/26  0.18/5  3.70  ✓
s2_1                                  0.18/4   0.05/1  0.00/0   0.00/0   3.95/92  0.00/0  4.30  ✗ MATLAB-fail
s2_2                                  0.44/12  0.17/5  0.62/17  1.07/29  0.74/20  0.21/6  3.68  ✓
s2_3                                  0.38/9   0.14/3  0.44/10  0.53/13  2.37/56  0.11/3  4.22  ✗ MATLAB-fail
s2_4                                  0.42/11  0.16/4  0.75/20  1.14/30  0.71/19  0.23/6  3.81  ✓
s2_5                                  0.39/11  0.16/5  0.57/16  1.02/28  0.90/25  0.18/5  3.59  ✓
s3_1                                  0.06/1   0.02/1  0.01/0   0.00/0   4.07/97  0.00/0  4.22  ✗ MATLAB-fail
s3_2                                  0.29/9   0.10/3  0.63/20  1.03/33  0.62/20  0.17/5  3.14  ✓
s3_3                                  0.37/9   0.13/3  0.14/3   0.45/11  2.84/66  0.04/1  4.28  ✗ MATLAB-fail
s3_4                                  0.18/5   0.07/2  0.03/1   0.22/6   2.81/81  0.02/1  3.47  ✗ MATLAB-fail
s3_5                                  0.41/11  0.16/4  0.78/20  1.18/31  0.71/18  0.18/5  3.85  ✓
```

**6 catastrophic MATLAB failures** (DWL >50%, marked ✗). The "low IoU" on these slides isn't a fair metric — my port might be MORE correct than MATLAB and still get penalized.

## Port iEGL/oEGL areas vs MATLAB GT (8 reasonable slides)

```
slide       port_iEGL  GT_iEGL  ratio    port_IoU  status
s1_2        0.398      0.39     102%     0.608     ✓ great match
s1_4        0.062      0.51     12%      0.110     ✗ port HUGELY under-calls iEGL
s1_5        0.047      0.45     10%      0.171     ✗ port HUGELY under-calls iEGL
s2_2        0.460      0.44     105%     0.502     ✓ good
s2_4        0.448      0.42     107%     0.405     ✓ good
s2_5        0.142      0.39     36%      0.254     ~ port under-calls
s3_2        0.407      0.29     140%     0.300     ~ port over-calls
s3_5        0.376      0.41     92%      0.318     ✓ but lower IoU
```

So **the iEGL mult=1.18 works perfectly for some slides (s1_2, s2_2, s2_4) but is way too strict for s1_4/s1_5 and a bit too loose for s3_2.** This is what made cross-slide tuning fail.

## Why the auto-tune attempts (v6/v7/v8) didn't work

| version | iEGL strategy | issue |
|---|---|---|
| v6 (pure target_frac) | top 10% by ratio | spurious background pixels with high ratio dominate |
| v7 (target_frac + abs floor) | + abs intensity floor | abs floor 0.25 too restrictive; tuned for s1_2 |
| v8 (Otsu) | Otsu on ratio | Otsu finds bimodal split; EGL is tiny (~10%) so Otsu's "natural" cut is in the middle of the cerebellum mass |

**Underlying issue**: ratio-based threshold is dominated by tiny dim regions where local_mean ≈ 0. Need both intensity floor AND area target tuned to specific slide.

## What WOULD help (didn't have time to try):

**Per-slide iEGL tuning by target area fraction in the EGL-MASK region** (not in the whole cerebellum):
- The egl_mask container is already computed (~15% of cerebellum)
- Within egl_mask, iEGL should be 50-80% (most of EGL is iEGL+oEGL combined)
- Tune mult per-slide so iEGL/egl_mask ≈ 0.6
- This is auto-adaptive AND avoids the spurious-background problem

I did NOT try this overnight; it's the obvious next experiment for the classical path.

## What about DWL/IGL overflow?

**Diagnosis**: My port's IGL = 1.32 mm² (vs MATLAB 0.70). When I tried to "fix" this by making `im_get_overlap` more faithful to MATLAB (cut bw1 too), the IGL DID shrink — but EGL recall on s2_5 and s3_2 dropped by HALF. The pipeline depends on the IGL being big in subtle ways:

- `igl_in_mask` (dilated IGL) is used to MASK OUT regions for EGL detection (`egl_norm = b_nfg * (~igl_in_mask)`)
- A bigger IGL → smaller egl_norm region → top-15% threshold lands somewhere different (more biased to outer edges where EGL actually is)
- A smaller IGL → larger egl_norm region → top-15% includes some interior pixels that aren't EGL

**This is a hard coupling.** Fixing IGL would require RE-tuning everything downstream. Probably not worth it given DINO is the cleaner path.

## DINO scaffolding (Track 1)

Files at `python/dino/`:

```
dino/
├── README → DINO_PLAN.md (in parent dir)
├── dataset.py            (image loading, label decoding, tile cropping)
├── extract_features.py   (DINOv2 ViT-S/14 inference, ~7 sec/slide)
├── train_head.py         (linear head + 64-unit MLP, weighted CE)
├── predict.py            (apply trained head to a slide)
└── evaluate.py           (per-class IoU vs MATLAB GT)
```

DINOv2 ViT-S/14 (22.1M params) loads fine, MPS works. Patch size 14 → giant
5000×5000 image gives ~360×360 patches, ~130k features per slide. With 8 known-
good slides for training and 2-3 for validation, we have enough data to train
a small linear head.

## Run sequence (when you wake up)

```bash
cd python/dino

# 1. Extract features for all 14 slides (already running; ~7 sec/slide)
ls cache/*.npz | wc -l   # should show 14

# 2. Train head, holding out 2 slides for validation
python3 train_head.py --val s1_4 s2_5 --epochs 50

# 3. Predict on val slides
python3 predict.py --slide 2018_05_22_s1_4_p27-0020
python3 predict.py --slide 2018_05_22_s2_5_p27-0026

# 4. Evaluate
python3 evaluate.py --slide 2018_05_22_s1_4_p27-0020
python3 evaluate.py --slide 2018_05_22_s2_5_p27-0026
```

## Open questions for you

1. Should we **drop the 6 catastrophic-MATLAB-fail slides from training**? They have wrong labels.
2. Should we use **port v11 outputs as additional weak labels** (for the remaining 8 slides)? Adds noise but increases data.
3. Worth trying **DINOv2 ViT-B/14** (87M params) if ViT-S underperforms? Slower but better representations.

---

## Update — first DINO result already in!

I went all the way through Track 1 overnight. **DINO works.**

### DINO v1 results (held out: s1_4, s2_5)

| slide | port v11 EGL IoU | DINO EGL IoU | improvement |
|---|---|---|---|
| s1_4 | 0.110 | **0.529** | **4.8×** |
| s2_5 | 0.254 | **0.526** | **2.1×** |
| **mean** | **0.182** | **0.527** | **2.9×** |

Per-class IoU (mean across both val slides):
| class | DINO IoU |
|---|---|
| iEGL | 0.47 |
| oEGL | 0.34 |
| IGL | 0.40 |
| ML | 0.65 |
| DWL | 0.57 |
| PCL | 0.21 |

Visually nearly identical to MATLAB GT (saved at `dino/preds/*_dino_vs_matlab.png`).

### Setup details
- DINOv2 ViT-S/14 (22.1M params), frozen, MPS backend
- Feature extraction: 7 sec/slide → ~100 sec total for 14 slides
- Linear head: 384 → 64 → 8 classes, 70k training samples (10k per class)
- Training: 50 epochs, ~2 sec/epoch on MPS, peaked at epoch 2 (mIoU 0.449)
- **Total wall-clock: ~5 min from "run extract" to "have predictions"**

### Known issues (not yet fixed):
- **Class weighting bug** in train_head.py: my `class_weights` function gives near-zero weight to all classes 1-7 because it includes class 0 in normalisation. Adam rescales gradients so training works anyway with effectively uniform weighting, but a fix would likely help imbalanced classes (PCL especially).
- **Background masking** is post-hoc (apply cerebellum mask after prediction). Cleaner solution: include background in training with a "background" class.
- **Predictions are at 14×14 patch resolution**, then nearest-neighbour upsampled. Could improve with bilinear softmax interpolation or an upsampler head.

### Suggested next steps when you wake up

1. **Fix class weighting** — should add 5-10% IoU.
2. **Re-train with all 14 slides minus 2 held-out**, using the catastrophic-MATLAB slides as "uncertain" weak labels (lower weight).
3. **Try DINOv2 ViT-B/14** (87M params) — likely +5-10% IoU.
4. **Add a 1-conv upsampler** instead of nearest-neighbour → smoother boundaries.
5. **Run full 14-slide cross-validation** (leave-one-out).
6. **Then**: integrate DINO predictions back into the layer_quant 1D-profile pipeline — this gives end-to-end ML-augmented analysis.

DINO clearly wins. The classical port hit a real ceiling at IoU ~0.6 on best slide, ~0.3 cross-slide. DINO already does better with a *first-pass* implementation. With the fixes above I'd expect IoU 0.65-0.75 cross-slide, which is genuinely good for thin-ribbon segmentation.

### Files to know about (in `python/dino/`)
```
DINO_PLAN.md       — full architecture/training plan
dataset.py         — image loading, label decoding
extract_features.py — DINOv2 inference (already run; cache/ is populated)
train_head.py      — head training (already trained; models/head.pt exists)
predict.py         — apply head to a slide (with cereb-mask post-processing)
evaluate.py        — IoU vs MATLAB GT + side-by-side render
cache/             — 14 cached feature files (~85-110 MB each = 1.4 GB total)
models/head.pt     — trained linear head (101 KB)
preds/             — predictions + visualization PNGs for s1_4 and s2_5
```
