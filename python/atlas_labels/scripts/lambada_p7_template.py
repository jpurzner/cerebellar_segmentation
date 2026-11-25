"""
Extract a P7 cerebellum shape template (outline + internal features) from the
LAMBADA P7 developmental atlas (25 um isotropic, light-sheet).

LAMBADA citation (CC-BY 4.0): De Launoit E, Skriabine S, Doumazane E,
Rousseau C, Renier N. LAMBADA, a Lightsheet-Aligned Mouse Brain Annotated
Developmental Atlas.

What LAMBADA gives at P7 (verified numerically, see stdout):
  * annotation label 512 (CB) = smooth cerebellar envelope with shallow
    label-0 notches at fissure mouths (no deep fissures, no layers)
  * protein_mbp             = myelin -> white-matter core / arbor vitae
  * reference               = autofluorescence; fissures appear as dark valleys
  * the volume is left/right mirrored, so one hemisphere is enough

Per sagittal slice (midline -> lateral) this script writes, at UP x the atlas
resolution (25/UP um per px), arrays in DISPLAY orientation:
    rows = superior -> inferior,  cols = anterior -> posterior (Allen-plate style)

Label codes follow the user's segmentation scheme (ribbon_unfold_laplace.LBL):
    0 = background, 1 = CEREB (cortex band between pia and WM, not split),
    6 = DWL (white matter, MBP-positive)
Boundary codes: 1 = pial outline, 6 = white-matter outline.
Fissure candidates (dark-valley ridges touching the surface) are saved
separately and NOT burned into the labels until confirmed.
"""
import json
from pathlib import Path

import numpy as np
import nibabel as nib
from scipy import ndimage as ndi
from skimage import filters, measure, morphology
from PIL import Image
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/Users/jpurzner/Dropbox/Q_research/py_projects/cerebellar_atlas_labels")
SRC = ROOT / "LAMBADA/LAMBADA-P7/LAMBADA_25um_{}_P7_v1.0.nii.gz"
OUT = ROOT / "outputs/lambada_p7"
OUT.mkdir(parents=True, exist_ok=True)

VOX_UM = 25.0
UP = 4                      # upsampling factor -> 6.25 um / px
PX_UM = VOX_UM / UP
CB_ID = 512
X_MID = 154                 # midline slice (atlas mirrored about 154.5)
MBP_LOG_THR = 9.6           # Otsu on log1p(MBP) for CB voxels with log1p(MBP) > 6
MARGIN = 6                  # voxels of padding around the CB bbox
FIG_OFFSETS = [0, 8, 16, 24, 32, 40]   # voxels from midline (x 25 um)


def load(name, dtype=float):
    return np.asarray(nib.load(str(SRC).format(name)).dataobj).astype(dtype)


def to_display(a2d):
    """LAMBADA slice [y(A-P, anterior +), z(S-I, superior +)] ->
    rows superior->inferior, cols anterior->posterior."""
    return a2d.T[::-1, ::-1]


def largest_cc(mask):
    lab = measure.label(mask)
    if lab.max() == 0:
        return mask
    return lab == (np.argmax(np.bincount(lab.ravel())[1:]) + 1)


def process_slice(ann, ref, mbp, x, bbox):
    y0, y1, z0, z1 = bbox
    cb = to_display(ann[x, y0:y1, z0:z1] == CB_ID).astype(float)
    r = to_display(ref[x, y0:y1, z0:z1])
    b = to_display(np.log1p(mbp[x, y0:y1, z0:z1]))

    # Sub-voxel envelope: smooth the binary label, upsample, threshold at 0.5.
    env = ndi.zoom(ndi.gaussian_filter(cb, 0.6), UP, order=1) > 0.5
    env = largest_cc(env)
    ru = ndi.zoom(r, UP, order=3)
    bu = ndi.zoom(ndi.gaussian_filter(b, 0.5), UP, order=3)

    # White matter: MBP above threshold, inside envelope, drop specks.
    wm = env & (bu > MBP_LOG_THR)
    wm = morphology.remove_small_objects(wm, max_size=(2 * UP) ** 2)
    wm = ndi.binary_fill_holes(wm)

    # Fissure candidates: dark ridges in autofluorescence, inside the envelope,
    # outside the WM, touching the pial surface.
    sr = filters.sato(ru, sigmas=range(UP // 2, UP + 1), black_ridges=True)
    sr[~env] = 0
    fiss = np.zeros_like(env)
    if env.any():
        thr = np.percentile(sr[env], 92)
        cand = morphology.skeletonize(env & (sr > thr) & ~wm)
        din = ndi.distance_transform_edt(env)
        lab = measure.label(cand, connectivity=2)
        for p in measure.regionprops(lab):
            cc = tuple(p.coords.T)
            if din[cc].min() <= 1.5 * UP and p.area >= 3 * UP:
                fiss[cc] = True

    label = np.zeros(env.shape, np.uint8)
    label[env] = 1
    label[wm] = 6
    bnd = np.zeros(env.shape, np.uint8)
    bnd[env & ~ndi.binary_erosion(env)] = 1
    bnd[wm & ~ndi.binary_erosion(wm)] = 6

    hull = morphology.convex_hull_image(env) if env.any() else env
    stats = dict(
        x=int(x), offset_um=float((x - X_MID) * VOX_UM),
        cb_area_mm2=float(env.sum() * PX_UM ** 2 / 1e6),
        wm_area_mm2=float(wm.sum() * PX_UM ** 2 / 1e6),
        wm_frac=float(wm.sum() / max(env.sum(), 1)),
        solidity=float(env.sum() / max(hull.sum(), 1)),
        n_wm_comp=int(measure.label(wm).max()),
        fissure_px=int(fiss.sum()),
        n_fissure=int(measure.label(fiss, connectivity=2).max()),
    )
    contours = dict(
        pia=[c * PX_UM for c in measure.find_contours(env.astype(float), 0.5)],
        wm=[c * PX_UM for c in measure.find_contours(wm.astype(float), 0.5)],
    )
    return dict(ref=ru, mbp=bu, env=env, wm=wm, fiss=fiss, label=label,
                bnd=bnd, stats=stats, contours=contours)


def main():
    ann = load("annotation", np.uint16)
    ref = load("reference")
    mbp = load("protein_mbp")
    cb3 = ann == CB_ID

    # Fixed bbox over the hemisphere so every slice shares one frame.
    idx = np.argwhere(cb3[X_MID:])
    _, ymin, zmin = idx.min(0)
    _, ymax, zmax = idx.max(0)
    bbox = (ymin - MARGIN, ymax + MARGIN + 1, zmin - MARGIN, zmax + MARGIN + 1)
    x_last = X_MID + int(idx[:, 0].max())
    xs = list(range(X_MID, x_last + 1))
    print(f"bbox y[{bbox[0]}:{bbox[1]}] z[{bbox[2]}:{bbox[3]}], slices x={xs[0]}..{xs[-1]}")

    results = {}
    labels, bnds, fisss = [], [], []
    for x in xs:
        res = process_slice(ann, ref, mbp, x, bbox)
        results[x] = res
        labels.append(res["label"])
        bnds.append(res["bnd"])
        fisss.append(res["fiss"])

    labels = np.stack(labels)
    bnds = np.stack(bnds)
    fisss = np.stack(fisss)
    H, W = labels.shape[1:]
    print(f"slice frame {H}x{W} px at {PX_UM} um/px; {len(xs)} slices")

    # ---- multi-page TIFF stacks (page k = offset k*25 um from midline)
    def save_stack(arr, path):
        pages = [Image.fromarray(a) for a in arr]
        pages[0].save(path, save_all=True, append_images=pages[1:],
                      compression="tiff_lzw")
    save_stack(labels, OUT / "p7_cb_template_labels_stack.tif")
    save_stack(bnds, OUT / "p7_cb_template_boundaries_stack.tif")
    save_stack(fisss.astype(np.uint8), OUT / "p7_cb_fissure_candidates_stack.tif")

    # ---- single midline slice for quick use
    mid = results[X_MID]
    Image.fromarray(mid["label"]).save(OUT / "p7_cb_midline_labels.tif", compression="tiff_lzw")
    Image.fromarray(mid["bnd"]).save(OUT / "p7_cb_midline_boundaries.tif", compression="tiff_lzw")

    # ---- contours (um, display frame) + stats
    np.savez_compressed(
        OUT / "p7_cb_contours.npz",
        x=np.array(xs), offset_um=np.array([(x - X_MID) * VOX_UM for x in xs]),
        px_um=PX_UM,
        **{f"pia_{x}": max(results[x]["contours"]["pia"], key=len, default=np.zeros((0, 2)))
           for x in xs},
        **{f"wm_{x}": (np.concatenate(results[x]["contours"]["wm"])
                       if results[x]["contours"]["wm"] else np.zeros((0, 2)))
           for x in xs},
    )
    stats = [results[x]["stats"] for x in xs]
    (OUT / "p7_cb_slice_stats.json").write_text(json.dumps(stats, indent=1))

    # ---- 3D label volume (25 um, atlas space, both hemispheres) for 3D use
    vol = np.zeros(ann.shape, np.uint8)
    vol[cb3] = 1
    vol[cb3 & (ndi.gaussian_filter(np.log1p(mbp), 0.5) > MBP_LOG_THR)] = 6
    nib.save(nib.Nifti1Image(vol, nib.load(str(SRC).format("annotation")).affine),
             str(OUT / "p7_cb_template_labels_25um.nii.gz"))

    # ---- overview figure (<= 1800 x 900 px). NOT to be opened by the agent.
    sel = [X_MID + o for o in FIG_OFFSETS if X_MID + o in results]
    fig, axes = plt.subplots(3, len(sel), figsize=(18, 9), dpi=100)
    for j, x in enumerate(sel):
        r = results[x]
        for i in range(3):
            ax = axes[i, j]
            ax.set_xticks([]); ax.set_yticks([])
        axes[0, j].imshow(r["ref"], cmap="gray", vmin=0, vmax=np.percentile(r["ref"], 99.5))
        axes[1, j].imshow(r["mbp"], cmap="magma", vmin=6, vmax=12.5)
        for i in (0, 1):
            axes[i, j].contour(r["env"], [0.5], colors="red", linewidths=0.8)
            axes[i, j].contour(r["wm"], [0.5], colors="yellow", linewidths=0.8)
        fy, fx = np.nonzero(r["fiss"])
        axes[0, j].plot(fx, fy, ",", color="cyan")
        cmap = matplotlib.colors.ListedColormap(
            [(0, 0, 0)] + [(0.1, 0.1, 0.4)] * 1 + [(0, 0, 0)] * 4 + [(1, 0.55, 0.1)])
        axes[2, j].imshow(r["label"], cmap=cmap, vmin=0, vmax=6, interpolation="nearest")
        axes[0, j].set_title(f"+{(x - X_MID) * VOX_UM:.0f} um from midline", fontsize=10)
    axes[0, 0].set_ylabel("autofluor.\nred=pia yellow=WM\ncyan=fissure cand.", fontsize=9)
    axes[1, 0].set_ylabel("MBP (log)", fontsize=9)
    axes[2, 0].set_ylabel("labels\n1=CEREB 6=DWL", fontsize=9)
    fig.suptitle("LAMBADA P7 cerebellum template - sagittal, anterior LEFT, dorsal UP", fontsize=12)
    fig.tight_layout()
    fig.savefig(OUT / "p7_cb_template_overview.png")
    plt.close(fig)
    w, h = Image.open(OUT / "p7_cb_template_overview.png").size
    print(f"overview figure {w}x{h}px")

    print("\nper-slice stats (every 4th):")
    print(" off_um  CB_mm2  WM_mm2  WMfrac  solid  nWM  nFiss fissPx")
    for s in stats[::4]:
        print(f" {s['offset_um']:6.0f}  {s['cb_area_mm2']:6.3f}  {s['wm_area_mm2']:6.3f}  "
              f"{s['wm_frac']:6.3f}  {s['solidity']:.3f}  {s['n_wm_comp']:3d}  "
              f"{s['n_fissure']:5d} {s['fissure_px']:6d}")


if __name__ == "__main__":
    main()
