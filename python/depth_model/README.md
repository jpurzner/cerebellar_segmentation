# Depth-band layer model

First attempt at the shape-based layer completion described in
`python/HANDOFF.md` ("Next phase"). Each cortex pixel is mapped into an
unrolled `(s, ψ)` frame, where `s` is arclength along the oEGL and `ψ` is
depth from oEGL (0) to DWL (1). In that frame the cortical layers should be
horizontal bands. A reference model of those bands, pooled over well-segmented
slides, is then used to relabel pixels whose depth disagrees with their label.

## Method

1. **Unfold** (`ribbon_unfold_laplace.py`). Solve Laplace's equation on the
   tissue with ψ = 0 on surface oEGL pixels (oEGL touching background,
   unclassified tissue or anomaly) and ψ = 1 on DWL. The pial edge itself is
   not used as the boundary because the segmentation often loses it. Each
   cortex pixel is traced along −∇ψ to the oEGL and takes the arclength of its
   landing point as `s`. Tissue fragments under 200 px are dropped (they made
   the solve singular).
2. **Score** (`depth_bands_diagnostic.py`, `batch_score_slides.py`). Per-slide
   ψ quantiles per layer and s-bin, plus a scorecard ranking slides on layer
   coverage, depth ordering, interface separation, PCL/IGL collapse into DWL,
   and streamline reach.
3. **Build** (`build_reference_model.py`). Pool `(s/s_max, ψ, layer)` over the
   chosen slides into per-bin ψ quantiles per layer, with smoothed interface
   curves.
4. **Apply** (`apply_reference_model.py`, `batch_apply_model.py`). Per-layer
   Gaussian likelihood in `(s/s_max, ψ)`, with a prior boost on the current
   label. A pixel is relabelled only if the top posterior is ≥ 0.7, its margin
   over the runner-up is ≥ 0.3, and the change is not an iEGL ↔ oEGL swap
   (those two are not separable by depth).

## Usage

Run from the repo root. Outputs go to `python/depth_model/out/` (gitignored).

```bash
python3 python/depth_model/batch_score_slides.py --sandbox python/matlab_anomaly_10x/sandbox --output_root python/depth_model/out/ribbon_unfold_10x --ds 4
```

```bash
python3 python/depth_model/build_reference_model.py --strategy all --tag 10x_all --ribbon_root python/depth_model/out/ribbon_unfold_10x
```

```bash
python3 python/depth_model/batch_apply_model.py --model 10x_all --ribbon_root python/depth_model/out/ribbon_unfold_10x
```

Without `--sandbox` / `--output_root`, `batch_score_slides.py` runs on the 20x
sandbox (`python/matlab_anomaly/sandbox`) into `out/ribbon_unfold`.
Single-slide entry points: `ribbon_unfold_laplace.py <labels.tif> [ds]`,
`depth_bands_diagnostic.py <ribbon_coords.npz>`,
`apply_reference_model.py --slide <slide_id>`, and
`compare_reference_models.py --a <tag> --b <tag>`.

## Results so far

- **Reference model:** `10x_all`, pooled from all 11 10x slides (1.7 M cortex
  pixels). The 10x slides scored lower than the 20x iteration variants on the
  coverage-weighted scorecard (mean 0.76 vs 0.80), mostly because PCL is hard
  to trace at 10x, but their layer depths are anatomically cleaner: IGL is not
  pushed into DWL. Interface medians: oEGL|iEGL 0.025, iEGL|ML 0.116,
  ML|PCL 0.243, PCL|IGL 0.450, all shallower than the 20x-derived model.
- **Applied to the 11 10x slides:** 0.23–1.19 % of cortex pixels relabelled
  per slide (0.65 % pooled). Largest corrections: ML → IGL (3,503 px),
  iEGL → IGL (1,879), ML → oEGL (1,501), IGL → ML (902), iEGL → ML (846).
- **20x sandbox:** 49 of 51 files scored. Most are iteration variants of about
  8 physical sections, so they are not independent samples.

## Known limitations

- **No true fill-in yet.** The model only relabels pixels that already carry a
  cortex label. Predicting labels for unassigned (label 1) or anomaly pixels
  would need the dense ψ field saved for those pixels too.
- `s2_2_p27_adjud2` and `adjud3` fail: no oEGL component reaches 200 px, so
  `oegl_arclength()` has nothing to parameterise.
- On the 10x slides about 10 % of oEGL-labelled pixels sit at ψ ≈ 1, deep in
  the cortex. Not yet checked whether that is an annotation convention or
  noise; it does not move the oEGL band median.
- 20x apply numbers recorded before the overwrite-mask fix (which had let
  iEGL ↔ oEGL swaps and low-margin flips through) are inflated and should be
  rerun.
