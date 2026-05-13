# Cerebellar Segmentation — Project History

A complete record of approaches attempted, results, failures, and outstanding
decisions. Written 2026-05-12 after the conversation hit context limits.

> **CRITICAL — DO NOT MODIFY THE GOLD STANDARD WITHOUT EXPLICIT USER PERMISSION**
>
> File: `python/labelled/2018_05_22_s1_2_p27-0021_labelled.tiff`
> This is the user's hand-corrected ground truth (~10 hours of work in napari,
> 9.96% of pixels edited). It is the only true ground truth in the project and
> is irreplaceable.

---

## 1. Project Goal

Build a working segmentation pipeline for P7 mouse cerebellar
immunofluorescence images that classifies cerebellar layers:

| label | abbr  | meaning                                       | typical marker            |
|------:|-------|-----------------------------------------------|---------------------------|
| 0     | bg    | background                                    | —                         |
| 1     | cereb | all cerebellum (mask)                         | tissue mask               |
| 2     | iEGL  | inner external granule layer (postmitotic)    | p27+ DAPI dense           |
| 3     | oEGL  | outer external granule layer (proliferating)  | DAPI dense, p27−          |
| 4     | IGL   | internal granule layer (mature granule cells) | NeuN+ DAPI dense          |
| 5     | ML    | molecular layer                                | sparse                    |
| 6     | DWL   | deep white layer (white matter)                | DAPI sparse, p27 sparse   |
| 7     | PCL   | Purkinje cell layer                            | thin between IGL/ML       |
| 8     | DCN   | deep cerebellar nuclei (NEW)                   | NeuN+ clusters inside DWL |

**Channel orders** (caused multiple bugs early on):
- 10x slides: `(p27, DAPI, NeuN)`
- 20x slides: `(p27, NeuN, DAPI)`

**Pixel sizes:**
- 20x: 0.5119049 µm/px
- 10x: ~1.024 µm/px

---

## 2. Slide Inventory

19 "good" slides used for training/eval (defined in `python/dino/slide_manifest.py`):
- 8 good 20x slides
- 11 good 10x slides

Excluded as bad: 20x: s1_3, s2_1, s2_3, s3_1, s3_3, s3_4. 10x: s4_B, s4_E, s4_G, s6_E.

The single hand-corrected slide is `2018_05_22_s1_2_p27-0021` (20x).

---

## 3. Approaches Attempted (chronological)

### 3.1 MATLAB baseline (existing, starting point)

- `cerebellum_threshold_segment2.m` (10x)
- `cerebellum_threshold_segment20x.m` (20x)

Threshold-based segmentation using `adaptthresh`, `bwareafilt`, `imfill`,
custom skeleton/curvature utilities. Produces RGB segments via `label2rgb`.

**Status:** generally good on 10x, weaker on 20x.

**Known weaknesses (per user observation):**
- iEGL/oEGL boundary biased toward iEGL on 20x
- NeuN+ clusters inside DWL (Deep Cerebellar Nuclei) get lumped into DWL

### 3.2 Classical Python port (`cerebellum_pipeline.py`)

A full Python re-implementation following MATLAB's logic (DAPI threshold →
EGL ring → IGL fill → ML/DWL/PCL completion).

**Result: BAD.** Mean IoU vs MATLAB on 10x slides where MATLAB is good:
**0.357**. Some slides had DWL=58% / IGL=0% (DWL over-fill bug never fully
fixed). User's verdict: *"the python port was pretty bad, it never even came
close to the matlab script in terms of segmentation."*

**Status:** ABANDONED.

Files: `cerebellum_pipeline.py`, `port_*.log` files. Multiple variant scripts
(`run_port_*.py`, `validate_*.log`) are vestiges of debugging.

### 3.3 Radial spokes / kernel approaches

Attempted to fix the Python port's EGL detection by drawing radial rays from
pia and IGL boundary, then voting along those rays.

Files: `radial_spokes.py`, `igl_spokes.py`, `spoke_profiles.py`,
`run_radial_v[1-5].py`, `diag_*` scripts.

**Result: did not converge to MATLAB-quality output.**
**Status:** ABANDONED.

### 3.4 DINOv2 ViT-S/14 + linear head (THE WINNING APPROACH)

Self-supervised vision transformer features + small head trained on MATLAB
labels as weak supervision.

- Feature extractor: `dinov2_vits14` (384-dim)
- Patch tokens via `forward_features()['x_norm_patchtokens']`
- Cache stored at `python/dino/cache_um0.5/` (target 0.5 µm/px input)
- Head: `LinearHead` in `python/dino/train_head.py`
- Output classes: 8 (does not include DCN)

**Patch-size sweep:** tested target µm/px ∈ {0.35, 0.5, 0.7, 1.0, 1.5, 3.0}.
Best: **0.5 µm/px** (≈ 7 µm patch).

**Trained models** (in `python/dino/models/`):
- `head_full_gold.pt` — DINOv2 head, all 19 slides, gold-boost ×4
- `head_v3_full_gold.pt` — DINOv3 head (worse, see below)

**Predictions:** `python/predictions_for_correction/<slide>_pred.tif`,
review panels in `python/review_panels/`.

**Status:** WORKING. The DINOv2 + gold-boost head is the current best
end-to-end model and what should be used for new slides.

### 3.5 DINOv3 ViT-S/16

License-gated weights — user obtained pre-signed URLs from Meta and
downloaded with curl (wget broken due to missing libunistring).

**Result:** visually WORSE than DINOv2. User: *"DinoV3 seems visually worse
than the V2."*

Files: `python/dino/cache_v3_um0.5/`, `python/dino/predict_and_render_v3.py`,
`python/predictions_v3/`, `python/review_panels_v3/`.

**Status:** retained for reference, not used.

### 3.6 Hand correction with napari (gold standard creation)

`python/dino/launch_napari.py` — opens DINOv2 prediction in napari for editing.

- Auto-save to `python/labelled/<slide_id>_corrected.tif`
- Ctrl-S keyboard shortcut + 💾 button (added after save bug)

**Bug:** initially napari saved an SVG instead of TIFF (wrong layer selected).
First edit attempt was lost. Fix: explicit Save button + Ctrl-S; also: do a
small test-save before doing major edits.

**Result for s1_2:** 263k → 2.3M pixels edited (9.96% of image).
Top transitions:
- iEGL → oEGL: 320k pixels (this is the "iEGL biased" problem the user saw)
- cereb → DWL: 196k
- PCL → IGL: 181k

Files in `python/labelled/`:
- `2018_05_22_s1_2_p27-0021_labelled.tiff` ← **the gold standard**
- `*_labelld.tiff`, `*_labeld.tiff` ← earlier saves, superseded
- `*_labeld.svg` ← the broken save (60 MB, leave for forensics, do not load)

### 3.7 Gold-boost training

Re-trained DINOv2 head with the gold-corrected s1_2 weighted ×4 in the loss.
Cache at `python/dino/cache_um0.5/` (gold-updated); backup of pre-gold cache
at `python/dino/cache_um0.5_pre_gold/`.

**Result:** measurable improvement on s1_2; small or neutral on other slides
(only one slide is gold-corrected).

### 3.8 MATLAB input-variant experiments (`run_matlab_experiments.py`)

Wrote a Python wrapper that runs real MATLAB via subprocess on s1_2 with
input variants:

| tag                  | input              | function                                    | result                       |
|----------------------|--------------------|---------------------------------------------|------------------------------|
| A_native_default     | native 0.51 µm/px  | `cerebellum_threshold_segment20x`           | OK; IoU(GOLD) = 0.735        |
| B_ds2x_default       | downsampled 2×     | `cerebellum_threshold_segment20x`           | OK; **IoU(GOLD) = 0.805** ✓  |
| C_blur4um_default    | Gaussian σ=4 µm    | `cerebellum_threshold_segment20x`           | too aggressive               |
| D_native_tightiegl   | native             | `cerebellum_threshold_segment20x_tight_iegl`| **CRASH (rc=-9, 1 sec)**     |
| E_blur4um_tightiegl  | blur σ=4           | `cerebellum_threshold_segment20x_tight_iegl`| crashed too                  |

**Headline result: ds2x downsampling improved IoU vs gold from 0.735 → 0.805
on s1_2.**

Subprocess pattern that finally worked:
- Write `.m` driver file with valid MATLAB identifier (no hyphens, no
  leading underscore).
- Invoke MATLAB with `-batch "cd('sandbox'); run('driver_xxx.m')"`.
- Earlier attempts with multi-line `-batch "..."` failed with "No MATLAB
  command specified".

Output: `python/matlab_experiments/s1_2_matlab_experiments.png`.

### 3.9 Tight-iEGL MATLAB variant (FAILED)

`cerebellum_threshold_segment20x_tight_iegl.m`: tightened iEGL params:
- line 342: adaptthresh sensitivity `0.3 → 0.15`
- line 350: threshold `>= 0.2 → >= 0.35`

**Result: rc=-9 within 1 second.** Likely a script-modification issue (could
be that the variable referenced by the threshold is renamed elsewhere). Not
yet investigated.

### 3.10 DCN class added (`cerebellum_threshold_segment20x_with_dcn.m`)

NEW MATLAB variant adding DCN as label 8.

Key changes:
- DCN detection block before set_bin assignment, using NeuN within DWL:
  ```matlab
  dcn_final = bwareafilt(dcn_seed > 0, [round(2000*scaling^2) Inf])
  set_bin(dcn_final == 1) = 8;   % between DWL and IGL priority
  ```
- `region_names = {..., 'DCN'};`
- Raw uint8 label save:
  ```matlab
  imwrite(uint8(set_bin), [file(1:(end-4)) '_labels.tif'])
  ```
- `imwrite(label2rgb(set_bin, jet(9), 'k'), ...)` — explicit `jet(9)` to fit
  9 classes (default `jet` is 7-row → "Index in position 1 exceeds array
  bounds" crash).
- `col_map` extended to 8 rows, magenta `[1 0 1]` for DCN.

Tested via `run_matlab_dcn.py` on s1_2 native + ds2x:
- **native:** 35 DCN clusters, ~2.3% of image, magenta correctly placed in
  central tract (anatomically plausible) ✓
- **ds2x:** only 0.31% DCN — downsampling smooths small clusters away. ✗

Output: `python/matlab_dcn_test/s1_2_with_dcn.png`.

**Trade-off (not yet resolved):**
- native gives DCN but worse layer boundaries (IoU = 0.735)
- ds2x gives better boundaries (IoU = 0.805) but loses DCN

---

## 4. Critical Bug: IGL/ML label swap

**Status:** discovered, NOT YET FIXED.

`python/dino/dataset.py::matlab_segments_to_labels()` parses the MATLAB
RGB segments TIFF by color. The mapping has IGL and ML SWAPPED:

| color    | parser says | MATLAB actually says |
|----------|-------------|----------------------|
| yellow   | IGL (4)     | **ML (5)**           |
| lt green | ML (5)      | **IGL (4)**          |

iEGL/oEGL/DWL/PCL are correctly mapped — only labels 4 and 5 are swapped.

**Verification:** running the with-DCN variant produces a raw uint8 label
file (not just RGB). Direct comparison:
- raw label 4 (MATLAB defines this as `c_final = IGL`) = 4,291,164 pixels =
  18.44%, displays as lt green
- raw label 5 (`ml_final = ML`) = 2,666,004 pixels = 11.46%, displays as yellow

The parser inverted these.

**Implications:**
- The gold standard (s1_2) was painted by the user looking at napari, which
  showed the (wrongly-named) labels. So the gold standard has the swap baked
  in: pixels labeled "IGL" are biologically ML and vice versa.
- The DINOv2 head learned the swap. Its predictions are internally
  consistent but biologically mis-named for IGL/ML.
- All IoU numbers are still meaningful (consistent within the project).
- Only the biological interpretation of "IGL" vs "ML" is wrong.

**Open decision (a/b):**
- (a) Fix `matlab_segments_to_labels` and re-train. Requires either
  re-painting the gold standard or remapping its labels 4↔5. Re-training
  DINO heads. Updating any scripts that hard-code label numbers.
- (b) Document and ignore. Treat label 4 = "ML" and label 5 = "IGL"
  internally; only translate at final reporting time.

---

## 5. Bugs Encountered & Fixes

| bug                                              | fix                                            |
|--------------------------------------------------|------------------------------------------------|
| MATLAB `-batch` with multi-line script failed    | Write `.m` driver, `cd + run('script.m')`      |
| MATLAB driver names with `-` rejected            | sanitize tag: replace `-` and `.` with `_`     |
| `label2rgb(set_bin, 'jet', 'k')` index OOB at 8  | use `label2rgb(set_bin, jet(9), 'k')`          |
| napari saved SVG instead of TIFF                 | added explicit Save button + Ctrl-S shortcut   |
| DINOv3 weights HTTP 403                          | user obtained pre-signed URLs from Meta        |
| `wget` broken (libunistring.2.dylib missing)     | use `curl` instead                             |
| Channel-order mismatch 10x vs 20x                | added `channel_order` to slide_manifest        |
| Python port DWL over-fill                        | NEVER FIXED → port abandoned                   |
| IGL/ML label swap                                | **NOT FIXED — open decision (§4)**             |
| Tight-iEGL MATLAB variant rc=-9                  | NOT INVESTIGATED                               |

---

## 6. Current Pipeline State

**Best end-to-end model for new slides:**
1. Run DINOv2 ViT-S/14 → `python/dino/cache_um0.5/<slide>_features.npz`
2. Apply head: `python/dino/models/head_full_gold.pt`
3. Outputs: `python/predictions_for_correction/<slide>_pred.tif`
4. Optional: hand-correct in `python/dino/launch_napari.py <short_id>`

**Best MATLAB variant per s1_2:**
- ds2x + default → IoU vs gold = 0.805 (no DCN)
- native + with_dcn → IoU = ~0.735 + magenta DCN class
- (no single config has both)

**MATLAB GT used for DINO training:**
- Generated from native + default (the original output).
- Has the IGL/ML swap baked into the parser, so DINO learned the swap.

---

## 7. Outstanding Decisions (waiting on user)

1. **IGL/ML label swap** — option (a) fix and re-train, or (b) document &
   ignore? (See §4.)
2. **DCN strategy:**
   - native + DCN (loses ML/DWL precision)
   - ds2x (better boundaries, no DCN)
   - merge: run both, take ds2x for layer assignments + native DCN mask as
     overlay (suggested but not yet implemented)
3. **Whether to regenerate "improved MATLAB GT"** from ds2x runs and
   re-train DINOv2 head on those.
4. **Whether to hand-correct more slides** — only s1_2 is gold; gold-boost
   only really helps when there's gold to boost.

---

## 8. Outstanding Tasks Not Yet Done

- Run ds2x on more 20x slides to confirm cross-slide generalization.
- Sweep downsample factors: 1.25×, 1.5×, 2×, 2.5×.
- Test gentler blur (σ=1, σ=2 µm); σ=4 too aggressive.
- Investigate why `cerebellum_threshold_segment20x_tight_iegl` rc=-9.
- Implement merge strategy for DCN (if user picks that option).

---

## 9. Key File Locations

| purpose                              | path                                                                       |
|--------------------------------------|----------------------------------------------------------------------------|
| Original MATLAB (10x)                | `cerebellum_threshold_segment2.m`                                          |
| Original MATLAB (20x)                | `cerebellum_threshold_segment20x.m`                                        |
| MATLAB tight-iEGL variant (BROKEN)   | `cerebellum_threshold_segment20x_tight_iegl.m`                             |
| MATLAB with DCN variant              | `cerebellum_threshold_segment20x_with_dcn.m`                               |
| Slide manifest                       | `python/dino/slide_manifest.py`                                            |
| MATLAB→labels parser (HAS BUG §4)    | `python/dino/dataset.py::matlab_segments_to_labels`                        |
| DINOv2 head training                 | `python/dino/train_head.py`                                                |
| DINOv2 cache (best model input)      | `python/dino/cache_um0.5/`                                                 |
| Best DINOv2 head                     | `python/dino/models/head_full_gold.pt`                                     |
| DINOv2 predictions                   | `python/predictions_for_correction/`                                       |
| napari editor                        | `python/dino/launch_napari.py`                                             |
| **Gold standard (DO NOT MODIFY)**    | `python/labelled/2018_05_22_s1_2_p27-0021_labelled.tiff`                   |
| MATLAB experiments runner            | `python/dino/run_matlab_experiments.py`                                    |
| MATLAB experiments output            | `python/matlab_experiments/s1_2_matlab_experiments.png`                    |
| MATLAB DCN runner                    | `python/dino/run_matlab_dcn.py`                                            |
| MATLAB DCN output                    | `python/matlab_dcn_test/s1_2_with_dcn.png`                                 |
| Earlier docs                         | `python/DINO_PLAN.md`, `python/DINO_RESULTS_SUMMARY.md`, `python/FINDINGS_OVERNIGHT.md` |

---

## 10. Useful Commands

```bash
# launch napari for hand correction (s1_2 already done — DO NOT REOPEN)
cd python/dino && python launch_napari.py <short_slide_id>

# re-run MATLAB input-variant experiments
cd python/dino && python run_matlab_experiments.py

# run MATLAB with DCN on s1_2
cd python/dino && python run_matlab_dcn.py

# re-train DINOv2 head with current gold-boost
cd python/dino && python train_head.py --gold-boost 4

# predict + render review panels
cd python/dino && python predict_and_render_v3.py  # (DINOv3 — worse, DINOv2 has its own runner)
```

---

*Written 2026-05-12 to preserve project state across the conversation
context window. If you read this and something is unclear or stale, ask
before guessing — and never edit the gold standard without permission.*
