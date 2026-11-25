"""Compare the P56 Nissl segmentation (outputs/allen_p56) with the Allen P56
developing-atlas annotation (CCF ids: CB subtree, *gr = granular, *mo = molecular)
on every annotated plate inside the processed range.  Dice per structure."""
import json, re, sys, urllib.request
from pathlib import Path
import numpy as np
from PIL import Image, ImageDraw
import xml.etree.ElementTree as ET
from svgpathtools import parse_path
sys.path.insert(0, str(Path(__file__).parent))
import allen_nissl_template as T

T.configure("P56")
stats = {s["sec"]: s for s in json.loads((T.OUT / "p56_nissl_stats.json").read_text())}
G = json.loads((T.ROOT / "allen_reference/CCF_structure_graph_1.json").read_text())
acr = {}
def walk(n):
    acr[str(n["id"])] = n["acronym"]; [walk(c) for c in n.get("children", []) or []]
[walk(x) for x in G["msg"]]
cbids = T.ccf_cb_ids()
q = json.loads(urllib.request.urlopen(
    "http://api.brain-map.org/api/v2/data/AtlasImage/query.json?criteria=atlas_data_set(atlases[id$eq181276165])"
    "&num_rows=500", timeout=60).read())["msg"]
plates = sorted([(i["section_number"], i["id"]) for i in q if i.get("annotated") and i["section_number"] in stats])
files = T.section_files()
dice = lambda a, b: 2 * (a & b).sum() / max(a.sum() + b.sum(), 1)
rows = []
for sec, aid in plates:
    cache = T.ROOT / "allen_reference/atlas_svgs_d2" / f"P56_sec{sec}_{aid}.svg"
    if not cache.exists():
        cache.write_bytes(urllib.request.urlopen(
            f"http://api.brain-map.org/api/v2/svg_download/{aid}?downsample=2", timeout=60).read())
    shape = Image.open(files[sec]).size[::-1]
    ims = {k: Image.new("L", (shape[1], shape[0]), 0) for k in ("cb", "gr", "mo")}
    def rast(el, sc):
        for m in re.finditer(r"scale\(([-\d.e]+)", el.attrib.get("transform", "")):
            sc *= float(m.group(1))
        sid = el.attrib.get("structure_id")
        if sid in cbids and "d" in el.attrib:
            p = parse_path(el.attrib["d"])
            pts = [(z.real * sc, z.imag * sc) for z in (p.point(t) for t in np.linspace(0, 1, 1500))]
            ImageDraw.Draw(ims["cb"]).polygon(pts, fill=1)
            if acr[sid].endswith("gr"): ImageDraw.Draw(ims["gr"]).polygon(pts, fill=1)
            if acr[sid].endswith("mo"): ImageDraw.Draw(ims["mo"]).polygon(pts, fill=1)
        for c in el: rast(c, sc)
    rast(ET.parse(cache).getroot(), 1.0)
    st = stats[sec]
    r0, c0 = st["crop_origin_px"]; h, w = st["crop_shape"]
    A = {k: np.asarray(v).astype(bool)[r0:r0 + h, c0:c0 + w] for k, v in ims.items()}
    A["wm"] = A["cb"] & ~A["gr"] & ~A["mo"]
    lab = np.asarray(Image.open(T.OUT / f"sections/sec{sec:03d}_labels.tif"))
    N = dict(cb=lab > 0, gr=lab == 4, mo=lab == 5, wm=lab == 6)
    px = (T.PX_UM / 1000) ** 2
    rows.append((sec, st["offset_um"], A["cb"].sum() * px, N["cb"].sum() * px,
                 dice(N["cb"], A["cb"]), dice(N["gr"], A["gr"]), dice(N["mo"], A["mo"]),
                 dice(N["wm"], A["wm"]), A["wm"].sum() * px, N["wm"].sum() * px))
print(" sec  off_um  atlasCB nisslCB | Dice CB   IGL    ML    WM | WM atlas nissl (mm2)")
for r in rows:
    print(f" {r[0]:3d} {r[1]:+6.0f}   {r[2]:5.2f}   {r[3]:5.2f}  |  {r[4]:.2f}  {r[5]:.2f}  {r[6]:.2f}  {r[7]:.2f} |  {r[8]:.2f}  {r[9]:.2f}")
arr = np.array([r[4:8] for r in rows])
print("median Dice  CB %.2f  IGL %.2f  ML %.2f  WM %.2f  (n=%d plates)" % (*np.median(arr, 0), len(rows)))
json.dump([dict(zip(["sec", "offset_um", "atlas_cb_mm2", "nissl_cb_mm2", "dice_cb", "dice_igl",
                     "dice_ml", "dice_wm", "atlas_wm_mm2", "nissl_wm_mm2"], map(float, r))) for r in rows],
          open(T.OUT / "p56_validation_vs_atlas.json", "w"), indent=1)
