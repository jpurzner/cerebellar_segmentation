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

---

# Addendum 2026-05-16: Iterations after May 12

This addendum captures the major iterations between May 12 and the
agent transition. The pipeline grew substantially since the original
doc was written. The current pipeline head is
`cerebellum_threshold_segment20x_with_dcn_anomaly.m` (1362 lines) and
a 10x sibling `cerebellum_threshold_segment2_with_dcn_anomaly.m`.

For a focused orientation see **`HANDOFF.md`** — this addendum is the
chronological record.

## 11. IGL membership + PCL size fixes (commit `5ccc0ca`)

User reported "non-NeuN positive cells gathering near the PC cell
layer" — DAPI-bright PC and basket cells were being absorbed into IGL
because the original `c_final = ca_whole | c_whole | cb_whole` was a
pure OR. The `cb_whole` (DAPI-confirmed IGL membership) pulled in any
DAPI-bright region in the smoothed-NeuN-top-40% area.

**Fix 1 — IGL membership requires raw NeuN.** Dropped pure-DAPI
confirmation entirely. New rule:
`c_final = c_whole | (ca_whole & c_nf > c_level)`. p27-confirmed
regions now also need raw NeuN above Otsu.

**Fix 2 — PCL gated on PC soma size.** Original `pc_layer_bin = (pc_layer_mask .* b_nfg) > 0.4`
caught every DAPI-bright pixel in the IGL/EGL ring. New:
`bwareafilt(seed, [500*scaling² 8000*scaling²])` — only blobs of
PC-soma size (132-2100 µm² at 20x = 13-52 µm diameter) survive.
Excludes basket/stellate interneurons (~8-10 µm) and large clumped
artifacts.

Validated on 8 slides. s2_2 IGL pulled back by 1.23%, PCL down 1.2%,
DWL grew correspondingly. Other layers stable.

## 12. V3 EGL/IGL adjudication (commit `c17f592`)

Addresses the EZH2 cKO failure mode where ML shrinks and EGL/IGL
physically touch. Dual-evidence rule:
- **EGL → IGL** if currently EGL AND raw NeuN > 0.65 AND > 50 µm from pia AND not in EGL territory (< 50 µm from pia)
- **IGL → EGL** if currently IGL AND raw NeuN < 0.30 AND < 50 µm from pia AND not in IGL territory (> 30 µm from pia)

Critical: uses **raw NeuN (`mat2gray(double(c))`)**, NOT `c_nf` (which
goes through CLAHE and gets inflated in EGL regions where signal is
uniformly low — CLAHE stretches the noise to span [0,1]).

Calibrated empirically on s2_2:
- iEGL NeuN q50 = 0.45, q95 = 0.87
- IGL NeuN q50 = 0.79, q5 = 0.38
- They overlap 0.40-0.50 → per-pixel NeuN cannot cleanly discriminate

Validated on 8 slides:
- Clean slides barely change (~0% reassignment)
- s2_2: 148K pixels reassigned, biggest EGL→IGL shift (-0.61% iEGL, +0.62% IGL)

## 13. V4 ribbon discontinuity (commits `87213d5` → `8a9c859`)

**V4 V1 (`87213d5`):** disc-closing-based ribbon discontinuity detector
for EGL and IGL ribbons. Bridges gaps up to 50 µm. Reassigns gap-fill
pixels currently labeled as the "other" layer with matching signal.

**V4 V2 (`8a9c859`):** Three major changes per user feedback:
1. **Directional close instead of disc** — multi-angle line strel
   union (6 angles, 0/30/60/90/120/150°, length 50 µm). Bridges
   along ribbon direction without expanding perpendicular (which a
   disc strel does, pulling in too much adjacent tissue).
2. **EGL reclaims from IGL OR DWL OR DCN** — not just IGL. User
   reported "DWL being in the EGL" and "DCN are expanding into the
   adjacent EGL".
3. **ML added as third ribbon** — same directional close logic.

Plus a sibling iteration for DWL fingers (DWL is structurally a
branching tree, but same multi-angle directional close bridges along
branches without inflating). DWL fingers reclaim from ML where
combined raw intensity < 0.25 AND > 50 µm from pia.

## 14. DCN-at-surface anatomical reclaim (commit `8177097`)

DCN is by definition deep cerebellar nuclei. A DCN label within 50 µm
of pia is anatomically impossible — DCN-detector over-extension where
imclose/imfill pushed the NeuN+ seed past the dwl_interior boundary
into adjacent EGL.

Reassigns surface-DCN pixels to iEGL (if p27+) or oEGL (if p27-).
Clean mass conservation: DCN -0.08 to -0.40% per slide, iEGL+oEGL gain
exact same.

## 15. DWL fingers + DWL-at-surface (commit `f96e2af`)

**DWL fingers** — same multi-angle directional close as ribbons,
applied to DWL tree. Reclaims ML pixels in DWL gap-fill regions where
combined raw intensity < 0.25 AND > 50 µm from pia.

**DWL-at-surface** — parallel to DCN-at-surface. DWL is deep white
matter, should never touch pia. Threshold tighter than DCN (30 µm vs
50 µm) since DWL legitimately sits closer to surface in folium tips.

Validated on 8 slides:
- s2_5 (user-reported "persistent DWL on outer surface"): DWL drop
  -0.81%, iEGL+oEGL gain +0.80
- All slides showed DWL-at-surface mis-detections (0.32-0.81%)

## 16. ML/DWL refinement attempts — ABANDONED

User: "still a lot of DWL getting called ML". Tried 5 variants in
sequence, all either too aggressive (5+% on clean s1_2) or too
conservative (~0% reassigned). Detailed in HANDOFF.md §"Iter-1".

**Lesson:** signal distributions DWL/ML are too similar (q50 0.28 vs
0.30), and ML is anatomically adjacent to IGL/DWL across most of its
boundary. Cannot be solved with signal + simple topology rules.
Needs atlas or learned discriminator.

## 17. Label-1 (missing-value) reassignment

**V1 (commit `552db4f`):** Naive nearest-neighbor propagation within
50 µm. Worked but had a critical bug: deep "cereb-only" pixels could
get EGL labels propagated INTO them (because a stray EGL pixel was
within 50 µm). User: "we are now putting an EGL across the base of
the cerebellum where there is DWL".

**V3 (commit `1d28675`):** Depth-gated propagation:
- pia_dist < 30 µm (surface) → {iEGL, oEGL, IGL, PCL}
- pia_dist 30-80 µm (mid)    → all biological labels
- pia_dist > 80 µm (deep)    → {IGL, ML, DWL, PCL, DCN}  (NO EGL)

Validated on 8 slides: every slide shows oEGL DOWN 0.83-1.97% with
ML/IGL UP correspondingly. s1_2 max move: oEGL -0.93%. Clean.

## 18. Continuity test suite (commit `402acad`)

`python/dino/test_segmentation_continuity.py`. Diagnostic-only.
Reports per-layer: n_components, n_large_components (>=500 px),
largest_area, largest_pct.

**Critical finding:** ML is HEAVILY fragmented across all 8 slides:
- ML largest component holds only 15-50% of total ML
- IGL largest 42-54%
- EGL: typically 56-96% (better)
- DWL: 2-14 large components (naturally branching)

Drove the design of iter-4 orphan cleanup.

## 19. Orphan fragment cleanup (commit `07a90be`)

For each ribbon (EGL/IGL/ML), find connected components ≤ 500 µm²
and reassign to dominant boundary label. Implemented as a local
function `orphan_cleanup` (sibling of `thresh_by_area` at end of file).

**Bug to be aware of:** original used `sum(set_bin == 5)` which on
2D arrays returns a 1xN vector (per-column sum), causing fprintf to
fire 4860 times. Fixed by using `nnz(...)`.

Validated on 8 slides:
- ML n_components cut ~half (e.g., s2_5: 101 → 55)
- IGL n_components also cut (s2_5: 23 → 12)
- PCL gains modestly (+0.15 to +0.49%) — many small ML islands
  surrounded by PCL territory get reclassified as PCL (likely real
  misclassified PC bodies)
- Largest_pct still doesn't change much — orphan cleanup addresses
  small fragments only, not the few big splits

## 20. iter-5 EGL-at-depth — TRIED AND REVERTED

Attempted "EGL labeled at pia_dist > 100 µm → reassign". Caught 5-9%
of image — over-corrected. User reported regressions: EGL lost from
deep folium invaginations + IGL bleeding to pia. Reverted (commits
`f5ae524` → `71c3d36`).

**Lesson:** pia_dist alone is unreliable because the cereb mask fills
sulci, making deep folium-tip EGL pixels look topologically "deep"
even though they're locally at the surface.

## 21. iter-6 EGL-at-cut-surface reclaim (commit `05d0116`)

User feedback after iter-5 revert: "we still need to fix the EGL over
the cut surface of the DWL / DCN".

**Anatomical insight:** EGL has at least 50 µm of ML/PCL/IGL between
it and DWL/DCN — even in thin folium tips. So EGL DIRECTLY ADJACENT
to DWL/DCN (within 30 µm, with no ML buffer) is anatomically
impossible — it's on a microtome cut surface.

**Dual-evidence rule:**
1. EGL pixel within 30 µm of DWL or DCN
2. Signal supports DWL: combined raw < 0.30 (dim)
   OR signal supports DCN: raw NeuN > 0.4 AND combined > 0.20

Real EGL has high DAPI AND low NeuN — adjacent EGL pixels with
EGL-like signal stay (could be folium-tip legitimate). Only clearly-
misclassified cut-surface EGL gets reassigned.

Validated on 8 slides — clean conservation: EGL pulls back 0.5-1.2%
per slide, DWL+DCN absorb. s1_2 max move = DWL +0.78%. User: "this
is the best so far".

## 22. 10x port (commit `d091f1f`)

Created `cerebellum_threshold_segment2_with_dcn_anomaly.m` as a port
of the full 20x anomaly pipeline. Adjustments:
- **Channel order swap:** at 10x page 2 = DAPI, page 3 = NeuN
  (opposite of 20x). Local variables a/b/c map the same (a=p27,
  b=DAPI, c=NeuN), just different imread page indices.
- **Pixel size:** 0.5119049 → 1.0239 throughout. All pia_dist_um,
  line_len_um, max_dist_px, orphan_max_area calculations updated.

Also created `python/dino/run_matlab_anomaly_10x.py` — batch runner
adapted from `run_matlab_dcn_10x.py` with the new function name and
output dir (`python/matlab_anomaly_10x/`).

Validated on 11 10x slides. All rc=0, layer percentages healthy:
- iEGL: 6.4-8.6%, oEGL: 2.3-3.1%, IGL: 17-23%, ML: 11-17%,
  DWL: 5-7%, PCL: 0.8-2%, DCN: 0.05-1.27%
- Anomaly: 1.8-2.5% (higher than 20x because smaller image area)

## 23. EGL boundary tooling (commits `d091f1f`, `67124c1`)

Three Python scripts for exporting EGL geometry, intended to feed a
downstream model that's learning the cerebellum boundary:

| script | output | use |
|--------|--------|-----|
| `render_egl_lines.py` | 3-panel PNG per slide | review/QC |
| `render_egl_boundary_mask.py` | binary outer/inner/skel TIFFs | boundary modeling |
| `render_egl_full_mask.py` | binary full-ribbon TIFFs | training |

Output dirs (all gitignored): `python/egl_lines/`,
`python/egl_boundary_masks/`, `python/egl_full_masks/`.

For 19 slides total (8 20x + 11 10x). Single-channel uint8, 0/255.

---

*Addendum written 2026-05-16 during agent handoff. Continue documenting
new iterations here. The companion `HANDOFF.md` is the focused
orientation doc for new agents.*
