# Cerebellar Segmentation

Automated segmentation of postnatal day 7 (P7) mouse cerebellar
immunofluorescence images into anatomical layers, using either a
classical MATLAB threshold pipeline or a DINOv2 self-supervised
linear-probe pipeline.

The project started as a MATLAB-only tool and now has two coexisting
pipelines: the original threshold-based MATLAB code and a Python /
DINOv2 pipeline that uses MATLAB output as weak supervision plus
hand-corrected gold-standard labels.

---

## What the pipelines segment

Input: a 3-channel widefield fluorescence image of a P7 mouse cerebellar
section. Channels are `p27`, `NeuN`, `DAPI` (channel order differs
between 10x and 20x scans — see *Image conventions* below).

Output: per-pixel layer labels.

| label | layer | abbr  | biology                                                |
|-------|-------|-------|--------------------------------------------------------|
| 0     | background     | bg    | outside the tissue                                              |
| 1     | all_cerebellum | cereb | tissue mask, no specific layer assigned                          |
| 2     | inner external granule layer  | iEGL  | postmitotic granule precursors, p27+ DAPI-dense                 |
| 3     | outer external granule layer  | oEGL  | proliferating granule precursors, p27− DAPI-dense                |
| 4     | internal granule layer        | IGL   | mature granule cells, NeuN+ DAPI-dense                          |
| 5     | molecular layer               | ML    | parallel fibres + Purkinje dendrites, sparse cells               |
| 6     | deep white layer              | DWL   | white matter, sparse DAPI                                        |
| 7     | Purkinje cell layer           | PCL   | Purkinje somata, DAPI-bright between IGL and ML                  |
| 8     | deep cerebellar nuclei (NEW)  | DCN   | NeuN+ clusters within DWL — added in `_with_dcn` variant         |
| 9     | anomaly (NEW)                 | —     | bleached/torn regions — added in `_with_dcn_anomaly` variant     |

> **Note on label 4 vs 5 in DINOv2 outputs.** A historical bug in
> `python/dino/dataset.py::matlab_segments_to_labels` swapped IGL and ML
> when parsing the MATLAB RGB segments TIFF. The DINOv2 head was
> trained on the swapped labels, the gold standard was painted in
> napari with the swapped labels visible, so the swap is internally
> consistent in the Python world. The raw MATLAB `set_bin` files
> produced by the `with_dcn` variants use the correct
> (`4 = IGL`, `5 = ML`) convention. See `python/PROJECT_HISTORY.md` §4
> for the full story.

---

## Image conventions

| magnification | pixel size (µm/px) | TIFF page order               |
|---------------|--------------------|-------------------------------|
| 10x           | 1.0239             | page 1 = p27, page 2 = DAPI, page 3 = NeuN |
| 20x           | 0.5119             | page 1 = p27, page 2 = NeuN, page 3 = DAPI |

The MATLAB scripts read the channels into local variables `a` (p27),
`b` (DAPI), `c` (NeuN) regardless of magnification — only the page
indices in `imread(file, page)` differ between the 10x and 20x scripts.

---

## Repository layout

```
cerebellum_threshold_segment.m              shared helpers (older)
cerebellum_threshold_segment2.m             ★ 10x classical pipeline (original)
cerebellum_threshold_segment2_with_dcn.m    10x + DCN class
cerebellum_threshold_segment20x.m           ★ 20x classical pipeline (original)
cerebellum_threshold_segment20x_with_dcn.m  20x + DCN class
cerebellum_threshold_segment20x_with_dcn_anomaly.m
                                            20x + DCN + tear/bleach anomaly
                                            detection + raw-NeuN IGL gate +
                                            PC-soma-size PCL gate (current head)
cerebellum_threshold_segment20x_tight_iegl.m
                                            tighter iEGL threshold variant
                                            (currently crashes — kept for
                                            reference)

background_subtract.m, smallfill.m, im_get_overlap.m,
thresh_by_area.m, get_thresh_by_area_adapt.m, ...
                                            shared MATLAB helpers used by
                                            the segmentation scripts

python/
├── PROJECT_HISTORY.md                      full chronological log of
│                                           approaches tried and what failed
├── DINO_PLAN.md, DINO_RESULTS_SUMMARY.md   intermediate planning notes
├── FINDINGS_OVERNIGHT.md                   one-night sweep notes
│
├── cerebellum_pipeline.py                  abandoned classical Python port
├── radial_spokes.py, igl_spokes.py, ...    abandoned spoke-based attempts
│
└── dino/
    ├── slide_manifest.py                   list of 8 good 20x + 11 good 10x
    │                                       slides with paths and channel order
    ├── dataset.py                          MATLAB labels → integer-class arrays
    ├── extract_features*.py                cache DINOv2 features per slide
    ├── train_head*.py                      train linear/MLP head on cached
    │                                       features (gold-boost weighting)
    ├── predict*.py                         apply head, post-process, save
    │                                       per-pixel label TIFFs
    ├── evaluate.py, run_loo*.py            leave-one-out cross-validation
    │
    ├── launch_napari.py                    open a slide in napari for hand
    │                                       correction; auto-saves to
    │                                       python/labelled/<slide>_corrected.tif
    │
    ├── run_matlab_dcn.py                   run with-DCN MATLAB on s1_2
    ├── run_matlab_dcn_multi.py             run with-DCN on all 8 20x slides
    ├── run_matlab_dcn_10x.py               run with-DCN on all 11 10x slides
    ├── run_matlab_anomaly.py               run with-DCN+anomaly on 8 20x
    ├── run_matlab_experiments.py           input-variant sweep (downsample, blur)
    ├── run_matlab_preproc_*.py             preprocessing variant tests
    │
    └── explore_tears_v[1-4].py             anomaly-detector design iterations
                                            (Python-only ground truthing
                                             before MATLAB port)
```

Outputs (gitignored):

```
python/dino/cache_um*/        cached DINOv2 patch features per slide
python/dino/models/           trained head weights (.pt)
python/predictions_*/         per-slide DINOv2 predictions (TIFF)
python/labelled/              hand-corrected gold standard (currently 1 slide)
python/matlab_*/              outputs of MATLAB runs (per-slide PNGs + masks)
python/anomaly_explore_*/     diagnostic panels from explore_tears_v[1-4].py
```

---

## Pipeline 1: Classical MATLAB

Open MATLAB, add the project root to the path, then call the
appropriate function on a 3-channel TIFF:

```matlab
addpath('/path/to/cerebellar_segmentation');
addpath(genpath('/path/to/cerebellar_segmentation/BaiSkeletonPruningDCE'));
addpath(genpath('/path/to/cerebellar_segmentation/Skeleton'));
addpath(genpath('/path/to/cerebellar_segmentation/WCHT'));
addpath(genpath('/path/to/cerebellar_segmentation/CurevUtils1.1'));
addpath(genpath('/path/to/MatlabFns'));   % Peter Kovesi's image library

% 20x with the latest detection improvements:
layer = cerebellum_threshold_segment20x_with_dcn_anomaly('myslide.tif');
% 10x with DCN:
layer = cerebellum_threshold_segment2_with_dcn('myslide.tif');
```

Each script writes alongside the input:

| filename suffix             | content                                             |
|-----------------------------|-----------------------------------------------------|
| `_labels.tif`               | uint8 label image, values 0..9 (per the table above) |
| `_segments.tif`             | RGB visualisation via `label2rgb(set_bin, jet(N), 'k')` |
| `_montage.tif`              | montage of intermediate channel maps                |
| `_anf_overlay.tif`          | p27 channel with `a_final` (iEGL) outline           |
| `_anomaly_*.tif`            | per-type anomaly masks (only `_with_dcn_anomaly`)    |

`layer` is a MATLAB table of per-region area + mean channel intensity.

### How the classical pipeline works

A high-level walkthrough of the algorithm — see
`cerebellum_threshold_segment20x.m` for the canonical implementation:

1. **Preprocess** each channel: `mat2gray` → `adapthisteq` (CLAHE) →
   `background_subtract` → `imgaussfilt(5)` → `anisodiff` →
   `phasesym` (edge detection).
2. **Detect Purkinje cell bodies** as small round blobs where p27 is
   moderately high AND p27 > DAPI AND p27 > NeuN (lines 91–117).
3. **Build the whole-cerebellum mask** by thresholding the sum of
   normalised DAPI + p27, then large morphological closing + hole
   filling (lines 125–148).
4. **Detect IGL** by smoothing NeuN heavily (3× `ordfilt2 + imgaussfilt`)
   then keeping the top 40% of intensity inside cerebellum, then
   confirming each candidate region with raw NeuN intensity (and
   optionally p27 / DAPI in the older OR-style logic).
5. **Detect DWL** as deep regions of moderate-low p27 with very low NeuN.
6. **EGL detection**: threshold DAPI in regions outside the IGL mask.
7. **iEGL/oEGL split**: adaptive thresholding on p27 (`a_bin`)
   intersected with the EGL ring → iEGL; complement DAPI thresholded
   inside the EGL ring → oEGL.
8. **PCL (layer)** = the ring between IGL and EGL where DAPI is bright.
9. **Final assignment** (priority is bottom-up — later overwrites
   earlier):

   ```
   1 = all_cerebellum
   5 = ml_final
   7 = pc_layer_bin       (PCL, the wider DAPI-bright ring)
   6 = dwl_final
   8 = dcn_final          (only in _with_dcn variants)
   4 = c_final            (IGL — overwrites PCL/DWL/DCN)
   3 = b_final            (oEGL)
   2 = a_final            (iEGL)
   7 = pc_bin_filt        (small round PC bodies — overrides EVERYTHING)
   9 = anomaly_mask       (only in _with_dcn_anomaly — overrides everything)
   ```

### Recent additions (`_with_dcn_anomaly` variant)

The `_with_dcn_anomaly.m` file extends the standard 20x pipeline with:

- **DCN class (label 8).** NeuN+ clusters within the DWL interior
  (eroded 5 px to keep DCN strictly inside white matter). Catches the
  large NeuN+ neurons of the deep cerebellar nuclei that would
  otherwise be lumped into DWL or IGL.
- **Tear / bleach detection (label 9).** Tissue-aware flat-field
  correction on each channel, then per-pixel `combined_ff < 0.13`
  inside an eroded morphological tissue mask, then shape filter:
  eccentricity ≥ 0.95 OR aspect ratio ≥ 4 → tear; ecc < 0.85 AND
  extent > 0.5 → bleach. Validated on s1_2's hand-labelled tear
  ground truth (79% recall on tear #1, 51% on tear #2).
- **IGL membership now requires raw NeuN.** Original code accepted IGL
  via OR of (NeuN, p27, DAPI). DAPI-confirmation incorrectly absorbed
  PCs and basket interneurons. New rule:
  `c_final = c_whole | (ca_whole & c_nf > c_level)`. Drops pure-DAPI
  confirmation entirely.
- **PCL gated on PC soma size.** Original `pc_layer_bin` was every
  DAPI-bright pixel in the IGL/EGL ring, which absorbed small
  basket/stellate interneurons (8–10 µm). New filter:
  `bwareafilt(seed, [500*scaling² 8000*scaling²])` → keeps only blobs
  of 13–52 µm diameter, the Purkinje soma range.

---

## Pipeline 2: DINOv2 weak supervision

This pipeline uses the MATLAB output as weak supervision to train a
linear / MLP head on top of frozen DINOv2 ViT-S/14 patch features.
After the head is trained, predictions are typically smoother and more
robust to anomalies than the raw MATLAB output, but still inherit any
systematic mis-classifications in the training labels (e.g. the IGL/ML
swap).

```bash
cd python/dino

# Cache features at target pixel size 0.5 um/px (best for P7 cerebellum)
python extract_features_unified.py --target-um 0.5

# Train head on all slides with 4x weight on hand-corrected ones
python train_full_gold.py --gold-boost 4

# Predict + render review panels for all slides
python predict_all_for_correction.py
```

`launch_napari.py <slide_short_id>` opens napari with the DINO
prediction (or hand-corrected save if it exists) loaded as an
editable label layer. Use the brush / fill / erase tools, then
Ctrl-S to save to `python/labelled/<slide_id>_corrected.tif`.

> **Critical: do NOT modify the gold standard
> `python/labelled/2018_05_22_s1_2_p27-0021_labelled.tiff` without
> express user permission.** It is the only ground truth available
> and represents ~10 hours of manual painting.

---

## Slide manifest

`python/dino/slide_manifest.py` defines 19 "good" slides (8 good 20x +
11 good 10x). Slides excluded as catastrophic-MATLAB on initial review:

- 20x: `s1_3, s2_1, s2_3, s3_1, s3_3, s3_4`
- 10x: `s4_B, s4_E, s4_G, s6_E`

The hand-corrected gold standard is `2018_05_22_s1_2_p27-0021` (20x).

---

## Outstanding issues

See `python/PROJECT_HISTORY.md` for full chronology and rationales.

- **EGL/IGL contact in EZH2 cKO.** EZH2 conditional knockouts have a
  shrunken ML, so EGL and IGL physically touch with no separation.
  The current pipeline can't always tell which is which when this
  happens. Best example: `2018_05_22_s2_2_p27` — in one folium EGL
  pixels get assigned IGL, elsewhere the opposite. Possible
  fix: per-pixel NeuN re-classification of any pixel currently in
  IGL ∪ EGL. Atlas-based deformation is a fallback.
- **IGL/ML label swap** in `dataset.py::matlab_segments_to_labels`
  (PROJECT_HISTORY §4). The DINOv2 head learned the swap; the gold
  standard was painted with the swap. Internally consistent but
  biologically misnamed for label 4 vs 5.
- **Tight-iEGL MATLAB variant crashes** (rc=-9), not yet investigated.
- **Tear detector recall ~80% on the largest tears, ~50% on smaller
  ones**, with ~50% of detections occurring in unlabelled (but
  visually plausible) small tears. May want to either tighten the
  detector or expand the gold standard.

---

## Dependencies

### MATLAB
- Image Processing Toolbox
- Statistics and Machine Learning Toolbox
- [Peter Kovesi's MATLAB image processing functions](https://www.peterkovesi.com/matlabfns/)
  — provides `phasesym`, `nonmaxsup`, `anisodiff`, `edgelink`,
  `filledgegaps`. Expected at `~/Dropbox/imaging_analysis/MatlabFns`
  per the runner scripts; adjust path as needed.
- Bundled (in repo, gitignored): `BaiSkeletonPruningDCE`, `Skeleton`,
  `WCHT`, `CurevUtils1.1` — used by some helpers.

### Python (conda recommended)
- Python 3.10+
- `torch`, `torchvision` — for DINOv2
- `numpy`, `scipy`, `scikit-image`, `scikit-learn`
- `tifffile`, `napari`, `matplotlib`
- `gh` (GitHub CLI) — only for `git push` if you need a high
  `http.postBuffer` workaround for the foundation commit (`524288000`)

### External tools
- MATLAB R2022b at `/Applications/MATLAB_R2022b.app/bin/matlab`
  (path is hard-coded in `run_matlab_*.py`; edit for other
  installations).

---

## Project status

Active development. The classical MATLAB pipeline is the canonical
output target for now; the DINOv2 pipeline is a working alternative
that benefits from hand correction but inherits any training-label
errors. Recent work has focused on reducing systematic
mis-classifications in the MATLAB pipeline (DCN class, anomaly
detection, IGL/PCL membership rules) rather than expanding the
DINOv2 side.
