# Shape-based layer correction

Post-processing for the MATLAB label maps (`matlab_anomaly/sandbox`, `matlab_anomaly_10x/sandbox`)
aimed at the main remaining failure: pieces of EGL, ML and IGL lost to neighbouring layers,
mostly **ML labelled DWL** and **DWL fingers labelled ML**.

Each slide calibrates its own model; nothing is pooled across slides. That avoids the
cross-slide comparability problem that sank the pooled depth-band model in `../depth_model/`.

```bash
python3 python/shape_model/eval_gold.py          # score against the s1_2 gold (read only)
python3 python/shape_model/run_shape_correct.py  # all good 20x + 10x slides -> out/
```

`out/<slide>_shape_labels.tif` uses the canonical label index at native resolution.
`out/<slide>_overview.png` and `out/<slide>_zoom.png` are the review panels.

## Method (`shape_correct.correct`)

Work at about 1 µm/px (20x downsampled ×2, 10x native).

1. **EGL rescue from DCN.** The iter-6 cut-surface step turned real EGL into DCN when the
   EGL touched DWL (on s1_2 about half of all DCN pixels are gold EGL). DCN pixels are
   returned to EGL when they are EGL-dense (DAPI above 0.6 on a scale where slide-median
   DWL = 0 and EGL = 1), NeuN-low, within 100 µm of IGL, and connected to the existing EGL
   across gaps under 20 µm. iEGL vs oEGL is set by p27. 92 % of these pixels are gold EGL.
   The gap tolerance catches the broken DCN strips along s1_2's left folium and the s2_2
   fissure: +0.005 ML, +0.006 DWL.
2. **Anchors.** The trusted EGL is:
   - components touching the outer surface;
   - ribbons of 8000 µm² or more (EGL lining a closed fissure);
   - fragments chained to either of those by gaps under 20 µm (ribbon continuity).

   Pixels with IGL-level NeuN are then dropped from the anchor. On s4_C, iter-6 labelled an
   IGL band as EGL; as an anchor it turned the neighbouring white-matter trunk into ML.

   The IGL is the second anchor.
3. **Anchor-distance relabel, 3 iterations.**
   1. Every ML, PCL and DWL/DCN pixel gets two coordinates: distance to the surface
      (the EGL anchor, or the tissue edge where EGL is missing) and distance to the IGL.
   2. A smoothed 2-D histogram of the slide's own labels gives P(ML, PCL or DWL | coordinates).
   3. A pixel is relabelled when the new class has posterior ≥ 0.7 and its current class
      ≤ 0.15. ML↔DWL uses a looser 0.6 / 0.2, because that is where the errors are.

   Guards:
   - PCL never becomes DWL.
   - DWL at a cut surface (within 40 µm of the edge, more than 40 µm from any EGL) is left alone.
   - **Betweenness:** DWL can only become ML or PCL if its nearest EGL and nearest IGL lie in
     opposite directions (cos < −0.3), because ML sits between EGL and IGL. On gold, 86 % of
     ML pixels pass this test and only 15 % of DWL pixels do.
4. **DWL sanity.** Newly created DWL is undone if:
   - the resulting DWL is thinner than 15 µm half-width (slivers between IGL and PCL);
   - it doesn't connect to existing DWL or DCN (white matter is one tree).
5. **Never cut the white-matter tree.** DWL→ML changes are undone when they sever a piece of
   the original tree that is at least 10,000 µm² and has at least 20 % of its border on IGL,
   i.e. the IGL-wrapped core of a folium.
   - On s4_C (10x), the pieces that would have been cut off were 77k and 202k µm² with
     48–64 % IGL borders.
   - The piece the s1_2 fix legitimately separates has no IGL border and isn't gold DWL.
   - The 1–2k µm² crumbs left on s2_2 after restoring ML don't count as cuts.

Unrecognised EGL fragments are **not** deleted by default (`relabel_orphan_egl=False`); see
the open issues below.

## Results on the gold slide (s1_2, ds2, DCN merged into DWL)

| IoU | ML | PCL | DWL | EGL | iEGL | oEGL | IGL |
|---|---|---|---|---|---|---|---|
| MATLAB iter-6 | 0.631 | 0.476 | 0.722 | 0.731 | 0.558 | 0.446 | 0.906 |
| shape-corrected | **0.719** | **0.508** | **0.819** | **0.749** | 0.566 | 0.451 | 0.906 |
| same, with orphan EGL deletion on (measured before the gap-tolerant rescue) | 0.730 | 0.509 | 0.840 | 0.759 | 0.577 | 0.451 | 0.906 |

- **Changed pixels:** 140k, 4.0 % of tissue. 80 % are now correct; 3 % were correct before.
- **Gold ML labelled DWL:** 77.7k → 27.9k. **Gold DWL labelled ML:** 65.1k → 32.3k.

### What is still wrong on s1_2

13 % of gold-labelled tissue (ds2), in order:

| Error | Pixels | Notes |
|---|---|---|
| oEGL called iEGL | 113k | Not fixable with the rules below. |
| PCL↔ML | 62k | Thin boundary noise, not chunks. |
| ML↔DWL | 60k | Mostly 3 blobs, see below. |
| IGL↔DWL | 45k | |
| ML↔IGL | 35k | |

The ML↔DWL blobs:
- A folium tip where gold ML borders white matter with no EGL. The betweenness gate blocks
  this by design.
- The white-matter finger, which needs orphan EGL deletion.
- A NeuN-dark finger inside the IGL.

## What was tried and dropped (s1_2 gold unless noted)

- **Depth from tissue edge plus local column pooling** made every class worse. Depth from the
  edge mixes folia, and sulci and tips distort it.
- **Distance to background as a third histogram axis** was too sparse: the misassigned chunk
  just confirmed itself.
- **Signed "which side of the IGL" feature** (two variants) gave no gain, because it is
  noisy right next to the IGL.
- **Betweenness as a histogram feature** instead of a gate was slightly worse.
- **Blocking ML→DWL unless EGL and IGL are on the same side** removed a third of the correct
  finger fixes and didn't improve precision.
- **Blocking DWL creation within 10–30 µm of PCL** also killed the finger fix, because iter-6
  puts stray PCL inside DWL fingers.
- **EGL rescue from DWL** (not just DCN) was 76–85 % precise on s1_2. On s1_4 it put EGL
  inside cell-dense P7 white matter, and ML then grew from it. It is off by default.
- **Protecting EGL fragments by DAPI alone** failed: stray EGL patches in s1_2 white matter
  and fissure EGL fragments in s2_2 have the same DAPI and NeuN. Gap-linking (20 µm)
  separates them; DAPI doesn't.
- **Blocking DWL→ML by original white-matter thickness** didn't separate the s4_C trunk
  (49–56 µm half-width) from the correct s1_2 fixes (30–72 µm).
- **Re-splitting iEGL/oEGL by the p27/DAPI ratio.**
  - The misassigned oEGL is DAPI-saturated with middling p27, and the ratio does
    distinguish gold oEGL (~0.35) from iEGL (~0.8).
  - But oEGL IoU rises only as much as iEGL falls (sum 1.01 → 1.01 at best).
  - Splitting by depth within the EGL ribbon was much worse (sum 0.62–0.78).
  - The MATLAB p27 split is the best of everything tried.
- **IGL rescue (EGL/DCN with IGL-level NeuN → IGL), aimed at the s4_C band.**
  - From EGL: wrong. Differentiating iEGL already has some NeuN, so it moved surface and
    fissure iEGL into IGL on s2_4 and s5_E.
  - From DCN only: never reached the s4_C band, because it isn't connected to IGL.
  - Removed.
- **Refining IGL boundaries by local NeuN** (5–15 µm bands, 5–15 µm smoothing) made every
  class worse (IGL 0.906 → 0.80–0.85, PCL 0.51 → 0.35–0.41). The painted gold boundaries
  are smoother than the signal.
- **Smoothing label boundaries** (per-class Gaussian vote, σ 1–5 px) made things slightly to
  clearly worse.
- **IGL → DWL where NeuN is DWL-like and connected to DWL** was 0 % precise on gold: every
  converted pixel was gold IGL or ML.

## Open issues

- **s4_C (10x) trunk.** About 10k px of a white-matter trunk still becomes ML. It is beside a
  band that iter-6 labelled DCN/EGL but whose signal is IGL. The fix belongs upstream:
  IGL-signal pixels shouldn't be labelled DCN or EGL. A post-hoc IGL rescue was tried and
  dropped (see above).

- **Orphan EGL policy.**
  - What it is: EGL fragments that are not surface-attached, not ribbon-sized and more than
    20 µm from the ribbon.
  - Deleting them fixes white-matter fingers on s1_2, worth +0.016 ML and +0.027 DWL IoU.
  - In the centre of s2_5 the same rule removes what looks like fissure EGL and fills it
    with DWL.
  - Off by default until we decide which way to go.
- **s2_5 centre and upper right, s3_2 fold.** These change the most, and the anatomy is
  ambiguous there. Review the zoom panels.
- **oEGL labelled iEGL** (gold oEGL 36 % called iEGL). Neither shape nor the p27/DAPI ratio
  fixes it. It may need a decision on where the painted oEGL/iEGL boundary should sit.
- **The continuity test is misleading as a target.** The gold slide's own ML and IGL are
  each in two large halves (largest piece 56 % and 52 %).

## Relation to `../depth_model/`

A parallel effort. Its pooled depth-band relabelling didn't beat the unmodified labels.
Its `anatomy_fixes.py` R1 rule (DWL with no IGL between it and the EGL → ML) uses the same
geometric idea as the betweenness gate here. Its reported numbers (ML 0.66 → 0.69, DWL
0.74 → 0.79) are at native resolution with a slightly different gold preparation, so they
aren't directly comparable with the table above.
