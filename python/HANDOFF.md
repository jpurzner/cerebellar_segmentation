# Cerebellar Segmentation — Agent Handoff (2026-05-16)

Quick orientation doc for the next agent. The full timeline lives in
`PROJECT_HISTORY.md`; the user-facing repo overview is in `../README.md`.
This file is a focused summary of **what works, what doesn't, and what
the layers actually are**, so the next agent can be productive
immediately without re-deriving things.

---

## TL;DR — current state

- **MATLAB pipeline is the canonical output target.** Two scripts:
  - `cerebellum_threshold_segment20x_with_dcn_anomaly.m` (20x, 1362 lines)
  - `cerebellum_threshold_segment2_with_dcn_anomaly.m` (10x, ported)
- **DINOv2 weak-supervision pipeline** exists in `python/dino/` and works
  but is mostly stale relative to the MATLAB improvements.
- **9-class output** (background, cereb-only, 6 anatomical layers, DCN,
  anomaly). Plus a label-9 "anomaly" class for detected tears/bleach.
- **Last commit on `main`:** `67124c1` (full EGL binary mask exporter).
  Repo: <https://github.com/jpurzner/cerebellar_segmentation>.
- **Critical: DO NOT MODIFY the gold standard**
  `python/labelled/2018_05_22_s1_2_p27-0021_labelled.tiff` without
  explicit permission. Single source of truth, ~10 hours of manual paint.

---

## Cerebellar anatomy (what the layers are)

P7 mouse cerebellum has a layered cortex over a deep white matter core.
**Cross-section, going from the pial surface inward:**

```
   PIA (cerebellum surface)
   ┌──────────────────────────────┐
   │  oEGL  ── DAPI dense, p27−  │  ← proliferating granule precursors
   │  iEGL  ── DAPI dense, p27+  │  ← postmitotic granule precursors
   ├──────────────────────────────┤
   │  ML    ── sparse cells,      │  ← molecular layer (parallel fibres
   │           moderate NeuN      │    + Purkinje dendrites)
   ├──────────────────────────────┤
   │  PCL   ── DAPI-bright soma   │  ← Purkinje cell layer (one cell
   │           between IGL/ML     │    thick)
   ├──────────────────────────────┤
   │  IGL   ── NeuN+ DAPI dense   │  ← internal granule layer (mature
   │                              │    granule cells)
   ├──────────────────────────────┤
   │  DWL   ── DAPI sparse, low   │  ← deep white matter (axon tracts)
   │           p27 and NeuN       │
   │  DCN   ── NeuN+ clusters     │  ← deep cerebellar nuclei (large
   │           WITHIN DWL         │    NeuN+ neurons embedded in DWL)
   └──────────────────────────────┘
   DEEP (cerebellar core)
```

The cortex (EGL → ML → PCL → IGL) wraps each folium — so in a slide
you see this ribbon-like pattern repeated around every folium, with
DWL/DCN in the central tract feeding the folia.

### Layer label scheme (used everywhere)

| label | abbr     | colormap (LABEL_COLORS) | bio                                |
|------:|----------|-------------------------|-------------------------------------|
| 0     | bg       | black                   | outside the tissue                 |
| 1     | cereb    | dark blue               | tissue mask, no specific layer     |
| 2     | iEGL     | sky blue                | inner external granule, p27+       |
| 3     | oEGL     | cyan                    | outer external granule, p27−       |
| 4     | IGL      | yellow                  | internal granule (NeuN+)           |
| 5     | ML       | light green             | molecular layer                    |
| 6     | DWL      | orange                  | deep white matter                  |
| 7     | PCL      | red                     | Purkinje cell soma                 |
| 8     | DCN      | magenta                 | NeuN+ deep nuclei in DWL           |
| 9     | anomaly  | dark gray               | tear / bleach (auto-detected)      |

> **Known bug carried throughout:** `python/dino/dataset.py::matlab_segments_to_labels`
> historically swapped labels 4 (IGL) and 5 (ML) when parsing the MATLAB
> RGB segments TIFF. The DINOv2 head was trained on the swapped labels;
> the gold standard was painted in napari with the swapped labels visible.
> All MATLAB raw `_labels.tif` files use the CORRECT convention
> (label 4 = IGL, label 5 = ML). See PROJECT_HISTORY §4 for the empirical
> verification. The label codes in the LABEL_COLORS array above are the
> MATLAB-raw correct convention.

---

## Ribbon vs finger taxonomy

This taxonomy is **central** to several of the recent algorithmic fixes.
Understanding it is essential before touching the V4 / V6 / orphan-cleanup
logic.

### Ribbons (continuous bands)

These layers form unbroken bands that wrap each folium:

- **EGL** (iEGL + oEGL) — outer band hugging the pia. Topologically a
  closed-curve ribbon around the entire cerebellum periphery,
  branching into each folium and around each sulcus. Typical thickness
  10-30 µm.
- **ML** — middle band, between EGL and IGL/PCL. Same topology as EGL
  (follows folium structure). Thickness ~20-50 µm.
- **IGL** — band of mature granule cells just inside PCL. Thickness
  50-150 µm. Slightly less ribbon-like than EGL/ML because it can be
  thicker / less uniform.
- **PCL** — single-cell-thick line of Purkinje somata between IGL and
  ML. Effectively a discontinuous string of cell bodies more than a
  continuous ribbon — handled separately.

### Fingers / branching structures

- **DWL** — the cerebellar arbor vitae (white-matter tree). Branches
  into each folium from the central peduncle. NOT a smooth ribbon;
  has finger-like extensions of varying width.
- **DCN** — three discrete clusters of large NeuN+ neurons embedded
  within the DWL. Round/oval blobs, not connected.

### Why this taxonomy matters

The **V4 ribbon discontinuity detector** (lines ~895–985 in
`cerebellum_threshold_segment20x_with_dcn_anomaly.m`) uses
multi-angle line strel closing (6 angles, every 30°) to bridge gaps
along the local ribbon/finger direction without expanding perpendicular
(which a disc strel would do). It runs on EGL, IGL, ML, AND DWL — each
with the same directional-close logic.

---

## Per-layer signal characteristics (empirical, from s2_2)

For pipeline tuning, these distributions matter. All are raw `mat2gray`
normalized values (no CLAHE) for the NeuN channel:

| layer | NeuN q5 | NeuN q50 | NeuN q95 | what it means |
|-------|---------|----------|----------|----------------|
| iEGL  | 0.26    | 0.45     | 0.87     | EGL has SIGNIFICANT NeuN at P7 (not zero like adult) |
| oEGL  | 0.20    | 0.36     | 0.51     | moderate, lower than iEGL |
| **IGL** | **0.38** | **0.79** | **1.00** | mature granule cells, very high NeuN |
| ML    | 0.20    | 0.30     | 0.52     | sparse cells, low-moderate NeuN |
| DWL   | 0.20    | 0.28     | 0.47     | very sparse cells |
| PCL   | 0.22    | 0.35     | 0.68     | DAPI-bright but NeuN range varies |
| DCN   | 0.31    | 0.43     | 0.70     | NeuN+ neurons in DWL |

**Key implications:**
- IGL and EGL **overlap** in the 0.40–0.50 NeuN range — per-pixel NeuN
  alone CANNOT cleanly discriminate EGL vs IGL at P7. Adjudication
  needs dual evidence (signal + topology).
- DWL and ML have near-identical NeuN distributions (q50 ≈ 0.28 vs
  0.30). **Signal alone cannot discriminate ML vs DWL** — boundary
  must come from topology.
- This is why several attempted fixes failed catastrophically (iter-1,
  iter-5 V4): pure per-pixel-signal rules can't separate overlapping
  classes.

---

## Pipeline architecture (current MATLAB)

`cerebellum_threshold_segment20x_with_dcn_anomaly.m` runs in this order:

1. **Read channels** (page 1=p27, page 2=NeuN, page 3=DAPI at 20x —
   reversed at 10x). Assign to local `a`, `b`, `c` (where `b`=DAPI,
   `c`=NeuN always).
2. **Preprocess** each channel: `adapthisteq` (CLAHE) →
   `background_subtract` → `imgaussfilt(5)` → `anisodiff` → `phasesym`.
3. **Detect PC bodies** — round bright p27+ spots (40-300 px circular).
4. **Build `all_cerebellum` mask** — `imclose` + `imfill('holes')` on
   thresholded combined-channel signal.
5. **Anomaly detection (V4)** — tear/bleach via tissue-aware flat-field
   + per-pixel `combined_ff < 0.13` + shape filter (tear = high
   eccentricity, bleach = round + filled). Outputs `*_anomaly_*.tif`.
6. **Subtract anomaly from `all_cerebellum`** — pipeline runs on
   cleaned tissue.
7. **Layer detection** (largely unchanged from original): IGL via
   smoothed NeuN top 40%, DWL via low p27 minus low NeuN, EGL via
   DAPI ring outside IGL, iEGL/oEGL split via adaptive p27 thresholds,
   ML by exclusion.
8. **DCN detection** — NeuN+ clusters within `dwl_interior` (DWL
   eroded by 10 px to avoid spillage).
9. **Initial `set_bin` assignment** — priority order (bottom-up
   overwrites): cereb → ml → pcl_layer → dwl → dcn → IGL → oEGL →
   iEGL → pc_bodies. (Final overrides happen later.)
10. **V3 EGL/IGL adjudication** (lines ~836–894) — dual-evidence rule
    (NeuN + pia_dist) for swapping pixels at the EGL/IGL contact
    zone. Targeted at EZH2-cKO mice where ML shrinks and EGL/IGL touch.
11. **V4 ribbon discontinuity (lines ~895–990)** — directional close on
    EGL, IGL, ML, AND DWL fingers. Reclaims gap-fill pixels with
    matching signal.
12. **DCN-at-surface reclaim** — anatomically impossible (DCN is
    deep), reassigned to iEGL/oEGL by local p27.
13. **DWL-at-surface reclaim** — same idea, tighter threshold (30 µm
    vs DCN's 50 µm).
14. **EGL-at-cut-surface reclaim (iter-6, line ~1020)** — EGL within
    30 µm of DWL/DCN with mismatched signal → reassigned. Targets
    the cut-surface failure mode where microtome exposed deep tissue
    is mis-labeled as EGL.
15. **Label-1 reclaim (depth-gated)** — propagates unassigned
    "cereb-only" pixels to nearest specifically-labeled neighbor,
    gated by depth zone (no EGL deep, no IGL at surface).
16. **Orphan fragment cleanup** — finds connected components < 500 µm²
    for each ribbon (EGL/IGL/ML) and reassigns each to its dominant
    boundary neighbor label.
17. **Final overrides:** `pc_bin_filt` → label 7, `anomaly_mask` →
    label 9. These are the most reliable detections and trump
    everything above.

`cerebellum_threshold_segment2_with_dcn_anomaly.m` is the 10x sibling
— same logic, swapped page 2/3 for channel read, 0.5119 → 1.0239 for
µm conversions.

---

## What we tried (chronological summary of recent iterations)

The full history is in PROJECT_HISTORY.md. Recent post-foundation
iterations:

### Iter-1: ML/DWL boundary refinement — **ABANDONED**

User reported "lots of DWL getting called ML". Tried 5 variants:
- V1: IGL-proximity gate (within 30 µm of IGL) — too aggressive,
  caught 5-6% of image because ML is BY DEFINITION adjacent to IGL.
- V2: Loosened V4 DWL fingers threshold + raw NeuN gate — too
  conservative, found ~0 new pixels.
- V3: "ML too deep" (pia_dist > 100 µm) — too aggressive, caught
  8-10% because most ML legitimately has pia_dist > 100 (folium tips
  + branching tract).
- V4: ML pixels at IGL boundary specifically — still caught ~5% (ML
  is adjacent to IGL across most of its inner boundary).
- V5: Sandwich rule (ML adjacent to BOTH IGL and DWL) — fired on ~0%
  (too rare).

**Lesson learned:** signal overlap is too narrow (DWL median NeuN
0.28, ML median 0.30) AND ML is anatomically adjacent to IGL by
design. Pure signal + topology rules can't reliably refine the
ML/DWL boundary. Likely requires either atlas-based reasoning or
better gold-standard data to train a discriminator.

### Iter-2: Missing-value (label 1) reassignment — **WORKING** (`552db4f`)

**V1 (committed):** Propagate "cereb-only" pixels to nearest labeled
neighbor within 50 µm. Worked but had a problem: deep cereb-only
pixels could get EGL labels propagated INTO them (because a stray
EGL pixel was within 50 µm).

**V3 fix (committed `1d28675`):** Depth-gated propagation:
- pia_dist < 30 µm (surface) → {iEGL, oEGL, IGL, PCL}
- pia_dist 30-80 µm (mid)    → any biological label
- pia_dist > 80 µm (deep)    → {IGL, ML, DWL, PCL, DCN}  (NO EGL)

Validated clean on 8 slides. Each slide had cereb 4-6% → 0-1%, with
most going to oEGL (most cereb-only pixels are near surface).

### Iter-3: Continuity test suite — **WORKING** (`402acad`)

`python/dino/test_segmentation_continuity.py`. Diagnostic only.
Reports per-layer: n_components, largest component %.

**Key finding:** ML is HEAVILY fragmented across all slides (largest
component holds only 15-50% of total ML). IGL also fragmented
(largest 42-54%). EGL relatively continuous (56-96%). DWL within
naturally-branching range.

This data drove iter-4.

### Iter-4: Orphan fragment cleanup — **WORKING** (`07a90be`)

For each ribbon (EGL/IGL/ML), find components ≤ 500 µm² and reassign
each to dominant boundary label. Local helper function
`orphan_cleanup` (sibling of `thresh_by_area` at end of file). Cut
ML n_components roughly in half on most slides. Mass conservation
clean.

**One bug** to be aware of: I originally used `sum(set_bin == 5)`
which returns a 1xN vector for a 2D array — fprintf then fires
once per element (~4860 times). Fixed by using `nnz(...)`.

### Iter-5 V4: EGL-at-depth reclaim — **REVERTED** (`71c3d36`)

Tried to fix "EGL crossing base of cerebellum" using pia_dist > 100 µm
+ signal mismatch. Caught 5-9% of image — OVER-corrected. User
reported NEW regressions: EGL lost from deep invaginations + IGL
bleeding to pia. Reverted.

**Lesson:** depth-based EGL reassignment is brittle because deep
folium tips legitimately have EGL near the surface but at large
`pia_dist` from the WHOLE-cerebellum boundary (the local pia
distance is small, but bwdist measures global).

### Iter-6: EGL-at-cut-surface reclaim — **WORKING** (`05d0116`)

Replaced iter-5 with an **anatomical adjacency** rule: EGL pixels
within 30 µm of DWL/DCN with signal supporting DWL or DCN → reassign.
Real EGL has high DAPI + low NeuN AND has at least 50 µm of
ML/PCL/IGL between it and DWL. So EGL DIRECTLY ADJACENT to DWL/DCN
must be on a microtome cut surface.

Validated on 8 slides — clean conservation: EGL pulls back 0.5-1.2%
per slide, DWL+DCN absorb. Other layers untouched. User said "this
is the best so far".

### 10x port — **WORKING** (`d091f1f`)

Created `cerebellum_threshold_segment2_with_dcn_anomaly.m` —
identical to 20x version with channel-order swap (page 2=DAPI,
page 3=NeuN at 10x) and `0.5119049 → 1.0239` for µm conversions.
Validated on 11 10x slides — clean layer percentages, all rc=0.

### EGL boundary tooling — **WORKING** (`d091f1f`, `67124c1`)

Three Python scripts for exporting EGL geometry for downstream
modeling (CNN training on cerebellum boundaries):

- `python/dino/render_egl_lines.py` — 3-panel PNG per slide
  (RGB | boundary+skel overlay | segmentation).
- `python/dino/render_egl_boundary_mask.py` — binary uint8 TIFFs:
  `<slide>_egl_outer.tif`, `_egl_inner.tif`, `_egl_skel.tif`.
- `python/dino/render_egl_full_mask.py` — full ribbon mask
  (`<slide>_egl_full.tif`, single channel uint8, 255 where EGL).

All outputs gitignored. Output dirs: `python/egl_lines/`,
`python/egl_boundary_masks/`, `python/egl_full_masks/`.

---

## What hasn't worked (don't repeat these mistakes)

1. **Per-pixel NeuN to discriminate EGL/IGL at P7.** Distributions
   overlap badly (0.40-0.50). Needs dual evidence.
2. **Per-pixel ML/DWL discrimination.** Signal distributions
   essentially identical. Iter-1 wasted 5 variants on this.
3. **Pure pia_dist for EGL-deep detection.** Folium-tip EGL has
   large pia_dist when the global cereb mask fills sulci. Caught
   real EGL as false positives. iter-5 V4.
4. **Aggressive disc dilation on ribbon detection.** Expands
   perpendicular AND along, pulling in adjacent tissue. Fixed by
   multi-angle line strel union (V4 V2).
5. **Pure-DAPI gating in IGL membership.** Original `c_final =
   ca_whole | c_whole | cb_whole` was OR of three signals.
   `cb_whole` (DAPI-confirmed) pulled in PC bodies and basket
   cells. Fixed in `5ccc0ca`.
6. **DINOv3 vs DINOv2.** DINOv3 was tried; user said "visually
   worse than V2". The DINOv2 ViT-S/14 pipeline at target 0.5 µm/px
   patch size is the better config.
7. **Python port of the MATLAB pipeline.** Abandoned. Mean IoU 0.357
   vs MATLAB on 10x; DWL over-fill never fixed.
8. **Tight-iEGL MATLAB variant** (`cerebellum_threshold_segment20x_tight_iegl.m`)
   — crashes with rc=-9. Never investigated.

---

## Known limitations / open TODOs

(From PROJECT_HISTORY plus recent user feedback)

1. **ML/DWL boundary** — heavy mis-classification persists. Needs
   atlas-based or learned approach.
2. **Big-piece splits in ML/IGL** — orphan cleanup only addresses
   small fragments. ML largest component is still 15-50% of total
   on most slides. Suggested fix: merge across PCL/DWL "bridges"
   (pixels of layer L that touch each other only through PCL or
   DCN should be considered connected).
3. **s2_2 small invagination lost to DWL** — specific failure mode
   on EZH2-cKO slides where ML is shrunk.
4. **s2_2 superior right DWL mislabel** — separate slide-specific
   issue.
5. **s3_2 fold disrupting EGL** — anomaly detector doesn't catch
   tissue folds. V1 had a fold detector but was dropped (false-
   positive on EGL). Could revisit with smarter logic.
6. **IGL/ML swap bug** (PROJECT_HISTORY §4) — Python `dataset.py`
   parser swaps labels 4 and 5. Should be fixed before any new
   DINO training.

---

## File map

```
cerebellum_threshold_segment20x.m                    ★ 20x original (10-year-old)
cerebellum_threshold_segment20x_with_dcn.m           20x + DCN class only
cerebellum_threshold_segment20x_with_dcn_anomaly.m   ★★ 20x CURRENT HEAD (1362 lines)
cerebellum_threshold_segment20x_tight_iegl.m         broken iEGL-tightening variant
cerebellum_threshold_segment2.m                      ★ 10x original
cerebellum_threshold_segment2_with_dcn.m             10x + DCN
cerebellum_threshold_segment2_with_dcn_anomaly.m     ★★ 10x CURRENT HEAD (port of 20x)

python/PROJECT_HISTORY.md                             full chronological log
python/HANDOFF.md                                     ★ this file
python/labelled/2018_05_22_s1_2_p27-0021_labelled.tiff
                                                      ★★★ GOLD — DO NOT MODIFY

python/dino/
├── slide_manifest.py                                 8 good 20x + 11 good 10x
├── dataset.py                                        ⚠ has IGL/ML swap bug
├── launch_napari.py                                  hand-correction editor
├── extract_features*.py                              DINOv2 feature cache
├── train_*.py                                        head training
├── predict*.py                                       head inference + render
├── run_matlab_anomaly.py                             20x batch runner (8 slides)
├── run_matlab_anomaly_10x.py                         10x batch runner (11 slides)
├── run_matlab_dcn.py / multi / 10x                   older with-DCN-only runners
├── render_egl_lines.py                               line-diagram visualization
├── render_egl_boundary_mask.py                       binary boundary TIFFs
├── render_egl_full_mask.py                           binary ribbon TIFFs
├── test_segmentation_continuity.py                   diagnostic continuity tests
└── explore_tears_v[1-4].py                           anomaly-detector design

python/matlab_anomaly/                                20x outputs (auto-generated)
python/matlab_anomaly_10x/                            10x outputs (auto-generated)
python/egl_lines/                                     line-diagram PNGs (gitignored)
python/egl_boundary_masks/                            binary boundary TIFFs
python/egl_full_masks/                                binary ribbon TIFFs
```

---

## Pixel sizes + channel orders (gets confused often)

| magnification | pixel size | TIFF page 1 | page 2 | page 3 |
|---------------|------------|-------------|--------|--------|
| 20x           | 0.5119 µm  | p27         | NeuN   | DAPI   |
| 10x           | 1.0239 µm  | p27         | DAPI   | NeuN   |

In the MATLAB scripts, local variables always map as: `a = p27`,
`b = DAPI`, `c = NeuN` (regardless of magnification). Only the
`imread(file, N)` page indices differ.

---

## Running the pipeline

### MATLAB

```matlab
addpath('/path/to/cerebellar_segmentation');
addpath(genpath('/path/to/cerebellar_segmentation/BaiSkeletonPruningDCE'));
addpath(genpath('/path/to/cerebellar_segmentation/Skeleton'));
addpath(genpath('/path/to/cerebellar_segmentation/WCHT'));
addpath(genpath('/path/to/cerebellar_segmentation/CurevUtils1.1'));
addpath(genpath('/path/to/MatlabFns'));

% 20x:
layer = cerebellum_threshold_segment20x_with_dcn_anomaly('myslide.tif');
% 10x:
layer = cerebellum_threshold_segment2_with_dcn_anomaly('myslide.tif');
```

Writes alongside the input: `_labels.tif` (uint8 0..9),
`_segments.tif` (RGB), `_anomaly_*.tif`, `_montage.tif`, etc.

### Python (batch on all good slides)

```bash
cd python/dino
# 20x: 8 slides ~25 min
/path/to/conda/python run_matlab_anomaly.py
# 10x: 11 slides ~15 min
/path/to/conda/python run_matlab_anomaly_10x.py
```

### Continuity tests

```bash
/path/to/conda/python test_segmentation_continuity.py
# reports per-layer n_components, largest_pct for each slide
# output: python/test_continuity/SUMMARY.txt + per-slide reports
```

### Export EGL masks (for downstream modeling)

```bash
/path/to/conda/python render_egl_full_mask.py      # full ribbon
/path/to/conda/python render_egl_boundary_mask.py  # outer/inner/skel
/path/to/conda/python render_egl_lines.py          # 3-panel review PNGs
```

---

## Repo status (16 commits ahead of original)

Each commit is anchored to a specific user concern. Recent ones in
order:

| commit  | what |
|---------|------|
| `67124c1`| Full EGL binary mask exporter |
| `d091f1f`| 10x pipeline port + EGL boundary tooling |
| `05d0116`| **EGL-at-cut-surface reclaim (iter-6) — user "best so far"** |
| `71c3d36`| Revert iter-5 V4 (over-correction) |
| `f5ae524`| iter-5 V4 EGL-at-depth (reverted) |
| `07a90be`| Orphan fragment cleanup (iter-4) |
| `1d28675`| Label-1 reclaim depth-gated (iter-2 V3) |
| `552db4f`| Label-1 reclaim V1 (had EGL-at-base bug) |
| `402acad`| Continuity test suite (iter-3) |
| `8a9c859`| Ribbon adjudication V4 V2 (directional close, 3 ribbons) |
| `87213d5`| Ribbon V4 V1 (disc closing) |
| `c17f592`| V3 EGL/IGL adjudication |
| `8177097`| DCN-at-surface |
| `f96e2af`| DWL fingers + DWL-at-surface |
| `5ccc0ca`| IGL membership + PCL size gate |
| `4cded9f`| README |
| `8fddebb`| Foundation (129 files) |

`git log --oneline` for current.
