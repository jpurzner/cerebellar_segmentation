# Atlas-guided layer labels and layer quantification (P7 Ezh2 cKO)

Downstream of the MATLAB / Python segmentation in this repository. It takes the anatomy-fixed label maps of the
P7 p27 / NeuN / DAPI sections, repairs them, unrolls the cortex into flattened strips, and measures layer
thickness, stain intensity and the shape and orientation of the nuclei in the molecular layer, wild type against
Ezh2 cKO.

This folder holds the code only. Its working copy, with the inputs and outputs (label maps, flattened strips,
Allen reference images, figures), is `~/Dropbox/Q_research/py_projects/cerebellar_atlas_labels`; each script's
paths point there and to `python/` of this repository (`SEG`). Every script starts with a docstring giving its
purpose, inputs, outputs and the decisions behind it.

## Sections, pixel sizes, genotypes

| set | sections | camera | µm / px |
|---|---|---|---|
| 2017 litter, 10x panoramas | s4_A, s4_C, s4_F, s5_A, s5_B, s5_C, s5_E, s5_F, s5_G, s6_C, s6_G | AxioCam MR3, 0.63× adapter | 1.0239 |
| 2018 litter, 20x panoramas | s1_2, s1_4, s1_5, s2_2, s2_4, s2_5, s3_2, s3_5 | AxioCam MR3, 0.63× adapter | 0.5119 |

Both sets come from the same Axio Imager.M2 (6.45 µm camera pixels / (objective × 0.63)); `scripts/zvi_read.py`
reads the `.zvi` metadata and `scripts/px_scale_check.py` checks the 10x scale against nuclear texture and the
Allen P7 atlas. The 10x labels and strips were built at 0.65 µm/px; `quant_layers.um_scale` converts every length
and threshold to the true scale.

Genotypes follow the original analysis notebook (`R_projects/cerebellar_layer_quant/analysis/p27_NeuN_Ezh2_cko_2024_1_09.Rmd`):
10x mice A, B, E cKO and C, D, F, G WT (the section letter is the mouse); 20x slide position 1–5 cKO, cKO, WT,
cKO, cKO. `quant_layers.genotype()` encodes this.

## Stages and main scripts

1. **Allen reference** — `download_allen_nissl.py`, `allen_nissl_template.py` (cerebellum outline and layers from the
   Developing Mouse Brain Nissl series), `fit_p7_template.py` (similarity + B-spline fit of the P7 template to a section).
2. **IGL repair** — `igl_detect.py` (IGL from all three stains), `igl_gapfill.py` / `igl_gapfill_apply.py` (gap filler),
   `igl_trace_apply.py`, `igl_finalize.py`, `folium_repair.py`; review pages `*_page.py` with their `*_template.html`.
   `final_labels.py` freezes the reviewed control label maps with a manifest.
3. **EGL repair and flattening** — `egl_repair.py`, `flatten.py` (columns at right angles to the outer EGL face),
   `walk.py` (the clockwise reference line), `shape_check.py`, `flat_fix.py` (layer Viterbi + stain model on the
   strips; the final 4 µm label grids).
4. **Quantification** — `quant_layers.py` (per-column thickness, per-layer stain, p27+ inner EGL, depth and
   layer-coordinate profiles), `ml_nuclei.py` (DAPI nuclei in the ML, eccentricity and angle to the line
   perpendicular to the IGL), `ml_nuclei_crops.py` (the same on the old 20x folium crops),
   `quant_figures.py`, `ml_nuclei_figures.py`, `figstyle.py` (style, statistics, layout check), `quant_report.py`.

Statistics: each mouse is the mean of its sections; cKO − WT by least squares on mouse means with litter as a
factor (`figstyle.tests`), with the pooled and 2017-only Welch tests beside it.

```bash
python3 scripts/quant_layers.py
python3 scripts/ml_nuclei.py && python3 scripts/ml_nuclei.py --down s1_2 s1_4 s1_5 s2_2 s2_4 s2_5 s3_2 s3_5
python3 scripts/quant_figures.py && python3 scripts/ml_nuclei_figures.py && python3 scripts/quant_report.py
```

Python 3 with numpy, scipy, pandas, scikit-image, tifffile, SimpleITK and matplotlib.
