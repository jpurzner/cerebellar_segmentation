#!/Users/jpurzner/miniconda3/bin/python3
"""Launch napari pre-loaded with image + MATLAB segmentation as starting labels.

Usage:
  python launch_napari.py s1_2          # any substring of slide_id
  python launch_napari.py s5_A
  python launch_napari.py s4_F --dino   # also overlay DINO prediction as separate layer

Edit the "labels" layer with brush/fill/erase. Save with napari File menu, OR
the script auto-saves the labels layer to corrected/<slide_id>_corrected.tif
when you close the viewer.

Layer label values (matches MATLAB set_bin):
  1 = all_cerebellum (default)
  2 = iEGL  (azure on MATLAB jet(8) palette)
  3 = oEGL  (cyan)
  4 = IGL   (yellow)
  5 = ML    (light green)
  6 = DWL   (orange)
  7 = PCL   (red)
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path
import numpy as np
import tifffile

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from slide_manifest import manifest, Slide
from dataset import load_slide, matlab_segments_to_labels

CORRECTED_DIR = ROOT.parent / "labelled"   # /python/labelled/ (sibling of /python/dino/)
CORRECTED_DIR.mkdir(exist_ok=True)
DINO_PRED_DIR = ROOT / "preds_unified"
# Newer DINO predictions for editing-as-start-point (gold-augmented model)
DINO_FOR_CORRECTION_DIR = ROOT.parent / "predictions_for_correction"

LABEL_NAMES = {
    0: "background", 1: "all_cerebellum", 2: "iEGL", 3: "oEGL",
    4: "IGL", 5: "ML", 6: "DWL", 7: "PCL",
}


def find_slide(query: str) -> Slide:
    slides = manifest()
    matches = [s for s in slides if query in s.slide_id]
    if not matches:
        print(f"no slide matches '{query}'. available:")
        for s in slides:
            print(f"  {s.mag}  {s.slide_id}")
        sys.exit(1)
    if len(matches) > 1:
        print(f"multiple matches for '{query}':")
        for s in matches: print(f"  {s.slide_id}")
        sys.exit(1)
    return matches[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("slide", help="substring of slide_id, e.g. 's1_2' or 's5_A'")
    ap.add_argument("--dino", action="store_true",
                    help="also overlay current DINO prediction as a separate layer")
    ap.add_argument("--source", choices=["matlab", "corrected", "dino", "auto"], default="auto",
                    help="where to load initial labels from. 'auto' picks: corrected if exists, else dino, else matlab")
    args = ap.parse_args()

    slide = find_slide(args.slide)
    print(f"slide: {slide.slide_id}  ({slide.mag}, pixel={slide.pixel_size_um}um)")
    print(f"  input: {slide.input_path}")
    print(f"  GT:    {slide.gt_path}")

    print("loading channels...")
    stack = load_slide(slide.input_path, channel_order=slide.channel_order)
    H, W = stack.shape[1], stack.shape[2]
    print(f"  shape ({H},{W})  channel_order={slide.channel_order}")

    # Initial labels — find newest existing corrected/labelled file for this slide
    corrected_path = CORRECTED_DIR / f"{slide.slide_id}_corrected.tif"
    candidates = sorted(
        list(CORRECTED_DIR.glob(f"{slide.slide_id}_corrected.tif*")) +
        list(CORRECTED_DIR.glob(f"{slide.slide_id}_labe*.tif*")) +
        list(CORRECTED_DIR.glob(f"{slide.slide_id}_labelled*.tif*")),
        key=lambda p: p.stat().st_mtime, reverse=True)
    dino_correction_path = DINO_FOR_CORRECTION_DIR / f"{slide.slide_id}_dino_pred.tif"
    # 'auto' resolution: prior corrections > DINO > MATLAB
    src = args.source
    if src == "auto":
        if candidates: src = "corrected"
        elif dino_correction_path.exists(): src = "dino"
        else: src = "matlab"
        print(f"auto-selected source: {src}")

    if src == "corrected" and candidates:
        resume_path = candidates[0]
        print(f"loading PRIOR labels (most recent): {resume_path.name}")
        init_labels = tifffile.imread(resume_path)
    elif src == "dino":
        # prefer the new gold-augmented predictions if they exist
        if dino_correction_path.exists():
            print(f"loading DINO (gold-augmented) prediction: {dino_correction_path.name}")
            init_labels = tifffile.imread(dino_correction_path)
        else:
            dino_path = DINO_PRED_DIR / f"{slide.slide_id}_pred.npz"
            if not dino_path.exists():
                print(f"DINO pred not found, falling back to MATLAB")
                init_labels = matlab_segments_to_labels(tifffile.imread(slide.gt_path))
            else:
                d = np.load(dino_path)
                init_labels = d["label_full"]
    else:
        print(f"loading MATLAB segments as initial labels: {slide.gt_path.name}")
        gt_rgb = tifffile.imread(slide.gt_path)
        init_labels = matlab_segments_to_labels(gt_rgb)
    if init_labels.shape != (H, W):
        print(f"  resizing labels {init_labels.shape} -> ({H},{W})")
        Hc = min(init_labels.shape[0], H); Wc = min(init_labels.shape[1], W)
        init_labels_full = np.zeros((H, W), dtype=init_labels.dtype)
        init_labels_full[:Hc, :Wc] = init_labels[:Hc, :Wc]
        init_labels = init_labels_full

    # Optional: DINO overlay as a second labels layer
    dino_labels = None
    if args.dino:
        dino_path = DINO_PRED_DIR / f"{slide.slide_id}_pred.npz"
        if dino_path.exists():
            d = np.load(dino_path)
            dino_labels = d["label_full"]
            if dino_labels.shape != (H, W):
                Hc = min(dino_labels.shape[0], H); Wc = min(dino_labels.shape[1], W)
                tmp = np.zeros((H, W), dtype=dino_labels.dtype)
                tmp[:Hc, :Wc] = dino_labels[:Hc, :Wc]
                dino_labels = tmp
            print(f"loaded DINO pred from {dino_path.name}")
        else:
            print(f"--dino requested but no pred at {dino_path}")

    # ===== Open napari =====
    import napari
    viewer = napari.Viewer(title=f"{slide.slide_id}  ({slide.mag})")

    # Add 3 single-channel image layers with proper colormaps so user can blend
    # them as a multi-channel display
    # channel_order tells us the role of each channel in stack[0..2]
    role_to_color = {"p27": "red", "neun": "green", "dapi": "blue"}
    for i, role in enumerate(slide.channel_order):
        viewer.add_image(stack[i], name=role, colormap=role_to_color.get(role, "gray"),
                          blending="additive",
                          contrast_limits=[stack[i].min(), np.percentile(stack[i], 99.5)])

    # Add MATLAB labels as the editable layer
    labels_layer = viewer.add_labels(init_labels.astype(np.int32),
                                      name=f"labels ({args.source})",
                                      opacity=0.45)

    # Add DINO pred as reference (read-only-ish: low opacity)
    if dino_labels is not None:
        viewer.add_labels(dino_labels.astype(np.int32), name="DINO pred (ref)",
                           opacity=0.25, visible=False)

    # Print legend
    print("\n=== label legend ===")
    for v, name in LABEL_NAMES.items():
        if v in [0, 1]: continue
        print(f"  {v} = {name}")
    print("\n=== keyboard shortcuts in napari Labels layer ===")
    print("  2 = paint mode | 3 = fill | 4 = pick | 5 = erase")
    print("  [ / ] = decrease/increase brush size")
    print("  ctrl-Z = undo | shift-A = toggle preserve labels")
    print("  press 'L' to bring the labels layer to front")

    # Add a Save button via magicgui — simpler and more reliable than close hook
    from magicgui.widgets import PushButton
    save_btn = PushButton(text="💾 Save corrected labels")
    info_msg = [""]
    def do_save(_=None):
        try:
            corrected = labels_layer.data.astype(np.uint8)
            tifffile.imwrite(corrected_path, corrected, compression="zlib")
            uniq = sorted(set(int(v) for v in np.unique(corrected)))
            msg = f"saved {corrected_path.name} | classes present: {uniq}"
            info_msg[0] = msg
            print(msg)
        except Exception as e:
            print(f"save failed: {e}")
    save_btn.changed.connect(do_save)
    viewer.window.add_dock_widget(save_btn, area="left", name="save")

    # Also bind Ctrl+S → save
    @viewer.bind_key("Ctrl-S", overwrite=True)
    def save_keybind(viewer):
        do_save()

    # Best-effort on-close save (might fail silently on macOS)
    def on_close(_=None):
        try:
            corrected = labels_layer.data.astype(np.uint8)
            tifffile.imwrite(corrected_path, corrected, compression="zlib")
            print(f"on-close auto-save: {corrected_path}")
        except Exception:
            pass
    try:
        viewer.window._qt_window.destroyed.connect(on_close)
    except Exception:
        pass

    print("\n💾  to save: click the 'Save corrected labels' button in the left panel")
    print("    OR press Ctrl-S")
    print(f"    saves to: {corrected_path}")

    napari.run()


if __name__ == "__main__":
    main()
