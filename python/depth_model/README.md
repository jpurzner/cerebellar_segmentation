# Depth-band layer model and anatomy fixes

First attempt at the shape-based layer completion described in
`python/HANDOFF.md` ("Next phase"). Two parts:

- **Depth-band relabelling** (`ribbon_unfold_laplace.py` … `batch_apply_model.py`).
  Maps cortex pixels into an unrolled `(s, ψ)` frame and relabels pixels whose
  depth disagrees with their label. Evaluated out of sample, it does **not**
  improve agreement with the gold slide (see Results).
- **Anatomy fixes** (`anatomy_fixes.py`). Four rules based on where the IGL
  is. They raise gold-slide accuracy from 83.1 % to 86.8 % and are the part
  worth using.

## Evaluation

`tune_relabel.py` scores any relabelling on:

- **Gold accuracy:** pixel accuracy and per-layer IoU on s1_2, the only
  hand-painted slide (`python/labelled/…_labelled.tiff`, read only). That file
  has IGL and ML painted swapped, so labels 4 and 5 are swapped back on load.
  Unassigned (label 1) and background pixels in the gold are not scored.
- **Continuity:** share of each layer's area in its largest connected piece,
  as in `python/dino/test_segmentation_continuity.py`. The gold slide itself
  scores only 56 % for ML and 52 % for IGL, below that test's 80 % threshold,
  so continuity is used as a sanity check, not as an objective.

## Anatomy fixes

White matter and DCN sit under the IGL; ML sits between EGL and IGL. R1 and
R2 ask whether IGL lies on the straight line from a pixel to its nearest EGL.

- **R1:** DWL with no IGL between it and the EGL, within 125 µm of the EGL →
  ML. A chunk is skipped if converting it would split off a significant piece
  of the white-matter tree.
- **R2:** DCN with no IGL between it and the EGL, within 100 µm → nearest
  cortex label.
- **R3:** DWL islands cut off from the main tree → the cortex label around
  them. Islands that reconnect to the tree through tear (anomaly) pixels are
  left alone: on the 10x slides these are real white matter beyond a tear.
- **R4:** tissue enclosed by IGL, after closing gaps of up to 25 µm in the
  IGL + white-matter ring → DWL. This recovers folium white-matter cores the
  pipeline missed: when their IGL ring is broken, the MATLAB step's
  `imfill(igl_in,'holes')` never fills them, and
  `bwareafilt(dwl_final,[200000 Inf])` drops any detached white matter.

Decisions are made at 4× downsampling and applied at native resolution only
to pixels carrying the rule's source label.

```bash
python3 python/depth_model/anatomy_fixes.py
```

Writes `out/anatomy_fixes/<slide>_fixed_labels.tif` for the 8 good 20x and
11 10x slides, plus `summary.csv`.

Results on s1_2 against gold, at native resolution:

| | before | R1–R3 | all four |
|---|---|---|---|
| pixel accuracy | 0.831 | 0.848 | **0.868** |
| ML IoU | 0.659 | 0.693 | 0.740 |
| DWL IoU | 0.736 | 0.796 | 0.871 |
| EGL IoU (mean of iEGL, oEGL) | 0.557 | 0.561 | 0.569 |

Of the changed area on s1_2, 64 % now matches gold and 7 % was right before
and is now wrong. R4 recovers 73 % of a white-matter branch in the upper
right that the pipeline had labelled ML, IGL, EGL and PCL.

On the other slides the rules change 4–8.4 % of cortex + white-matter pixels
on 20x and 0–2.4 % on 10x. R4 adds nothing on 10x. On the 20x slides each R4
fill is one or two large patches with the same signal as the s1_2 branch
(low NeuN and p27, higher DAPI than ML), although R4 picks them by position
only. Worth checking by eye: s3_5 (largest fill, 0.25 mm²) and s1_5 (small
fill, half from EGL, with EGL-like DAPI).

Caveat: all thresholds were chosen on the one gold slide. They sit on
plateaus (R1 gives +1.3 points anywhere from 100 to 150 µm; R4 gives the same
gain with 16 or 25 µm gap closing, and fails at 49 µm, where it seals across
fissures), but nothing outside s1_2 checks that the relabelled pixels are
correct.

An intensity-only detector for missed white matter was also tried and
dropped: the best setting recovered about 20 % of the s1_2 branch, and small
threshold changes swung it to flagging only wrong pixels.

## Depth-band relabelling

### Method

1. **Unfold** (`ribbon_unfold_laplace.py`). Solve Laplace's equation on the
   tissue with ψ = 0 on surface oEGL pixels and ψ = 1 on DWL. Each cortex
   pixel is traced along −∇ψ to the oEGL and takes the arclength of its
   landing point as `s`. Tissue fragments under 200 px are dropped.
2. **Score** (`depth_bands_diagnostic.py`, `batch_score_slides.py`).
3. **Build** (`build_reference_model.py`). Pool `(s/s_max, ψ, layer)` over
   slides into per-bin ψ quantiles per layer.
4. **Apply** (`apply_reference_model.py`, `batch_apply_model.py`). Per-layer
   Gaussian likelihood with a prior boost on the current label; relabel if the
   top posterior is ≥ 0.7, its margin is ≥ 0.3, and the change isn't
   iEGL ↔ oEGL.

### Results

`tune_relabel.py` swept 270 settings: three ways of estimating the bands
(pooled per arclength bin, pooled ignoring arclength, and per slide from its
own labels along its own arclength) × confidence 0.5–0.9 × margin 0–0.3 ×
prior boost 0–1.5. Pooled models for a 10x target leave that slide out.
Against the gold slide (baseline accuracy 0.831 at 4× downsampling):

| band estimate | best | worst |
|---|---|---|
| pooled per arclength bin | −2.9 points | −12.8 points |
| pooled, ignoring arclength | ±0 (changes nothing) | −3.4 points |
| per slide, own labels | ±0 (changes < 0.2 %) | −2.6 points |

No setting improves on the unmodified labels. The earlier claim that the
`10x_all` model "relabels 0.65 % of the 10x slides and catches more errors
on bad slides" was an in-sample effect (see below).

### Why it doesn't work

- **Arclength isn't comparable across slides.** Where `s = 0` falls and the
  order of the oEGL components are arbitrary per slide, so a pooled
  per-`s`-bin band mixes different anatomy, with as few as 6 PCL pixels per
  bin. In sample, the bands partly fit the target's own labels, which is why
  the stored `10x_all` model changes only about 1 % of its own slides; left
  out, the same model changes 8–13 %, much of it noise (e.g. IGL → iEGL).
- **ψ overlaps too much between layers.** ML's spread covers most of 0–1.
- **The largest errors are on the DWL boundary**, and ψ pins every
  DWL-labelled pixel at 1, so a chunk of ML wrongly labelled DWL becomes
  "deep" by construction. The anatomy fixes avoid ψ for exactly this reason.

### Usage

Run from the repo root; outputs go to `python/depth_model/out/` (gitignored).

```bash
python3 python/depth_model/batch_score_slides.py --sandbox python/matlab_anomaly_10x/sandbox --output_root python/depth_model/out/ribbon_unfold_10x --ds 4
```

```bash
python3 python/depth_model/build_reference_model.py --strategy all --tag 10x_all --ribbon_root python/depth_model/out/ribbon_unfold_10x
```

```bash
python3 python/depth_model/tune_relabel.py --quick
```

`tune_relabel.py` and `anatomy_fixes.py` read unfold outputs from
`out/ribbon_unfold_20x/` (8 good 20x slides) and `out/ribbon_unfold_10x/`
(11 10x slides).

## Remaining errors on s1_2 (not addressed)

- **oEGL labelled iEGL** (35 % of gold oEGL): a consistent shift of the
  iEGL/oEGL boundary toward the surface, within about 10 px. Looks like a
  convention difference between the p27 threshold and the painting; needs a
  decision on which is right, not a geometry rule.
- **Gold DWL labelled ML at the cut surface:** a median 380 µm from any
  pipeline DWL, and indistinguishable from ML by depth or IGL position. Needs
  image signal.
- **ML ↔ PCL:** boundary placement within about 10 px.

## Known limitations

- No fill-in yet for unassigned (label 1) or anomaly pixels.
- `s2_2_p27_adjud2` and `adjud3` fail in the unfold: no oEGL component reaches
  200 px.
