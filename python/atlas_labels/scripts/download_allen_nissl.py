"""
Download the cerebellum range of an Allen Developing Mouse Brain reference
sagittal Nissl series.

    python scripts/download_allen_nissl.py <SectionDataSet id> <out_dir> [half_range_um]

1. section list -> out_dir/section_images.json
2. tiny thumbnails (downsample=6) of every section -> midline by tissue-area symmetry
3. d2 images (downsample=2, ~4 um/px) for sections within +-half_range_um of the
   midline -> out_dir/secNNNN_id<ID>_d2.jpg  (skips files already present)
"""
import io
import json
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from PIL import Image

API = "http://api.brain-map.org/api/v2"


def fetch(url, tries=4):
    for k in range(tries):
        try:
            return urllib.request.urlopen(url, timeout=120).read()
        except Exception as e:  # noqa: BLE001
            if k == tries - 1:
                raise
            time.sleep(2 * (k + 1))


def main(ds_id, out_dir, half_range_um=2500.0):
    out = Path(out_dir)
    (out / "thumbs").mkdir(parents=True, exist_ok=True)
    meta = json.loads(fetch(f"{API}/data/SectionDataSet/{ds_id}.json?include=section_images"))["msg"][0]
    imgs = sorted(meta["section_images"], key=lambda i: i["section_number"])
    (out / "section_images.json").write_text(json.dumps(imgs, indent=0))
    sec_um = float(meta["section_thickness"])
    print(f"dataset {ds_id}: {len(imgs)} sections, {sec_um} um thick", flush=True)

    def thumb(i):
        p = out / "thumbs" / f"sec{i['section_number']:04d}_id{i['id']}_d6.jpg"
        if not p.exists():
            p.write_bytes(fetch(f"{API}/section_image_download/{i['id']}?downsample=6"))
        a = np.asarray(Image.open(p).convert("L")).astype(float)
        return i["section_number"], float((a < 200).sum())

    with ThreadPoolExecutor(4) as ex:
        area = dict(ex.map(thumb, imgs))
    secs = np.array(sorted(area))
    ar = np.array([area[s] for s in secs])
    best = (-1, None)
    for c in np.arange(np.percentile(secs, 25), np.percentile(secs, 75), 0.5):
        m = np.array([(a1, np.interp(2 * c - s0, secs, ar)) for s0, a1 in zip(secs, ar)
                      if secs.min() <= 2 * c - s0 <= secs.max()])
        if len(m) > 0.6 * len(secs):
            r = np.corrcoef(m[:, 0], m[:, 1])[0, 1]
            if r > best[0]:
                best = (r, c)
    mid = best[1]
    half = half_range_um / sec_um
    want = [i for i in imgs if abs(i["section_number"] - mid) <= half]
    print(f"midline ~ section {mid} (r={best[0]:.3f}); downloading {len(want)} d2 images "
          f"(sections {want[0]['section_number']}..{want[-1]['section_number']})", flush=True)
    (out / "midline.json").write_text(json.dumps(dict(mid_section=mid, r=best[0])))

    done = [0]

    def d2(i):
        p = out / f"sec{i['section_number']:04d}_id{i['id']}_d2.jpg"
        if not p.exists():
            data = fetch(f"{API}/section_image_download/{i['id']}?downsample=2")
            Image.open(io.BytesIO(data)).verify()
            p.write_bytes(data)
        done[0] += 1
        if done[0] % 25 == 0:
            print(f"  {done[0]}/{len(want)}", flush=True)
        return p.stat().st_size

    with ThreadPoolExecutor(4) as ex:
        sizes = list(ex.map(d2, want))
    print(f"done: {len(sizes)} files, {sum(sizes) / 1e6:.0f} MB in {out}", flush=True)


if __name__ == "__main__":
    main(int(sys.argv[1]), sys.argv[2], float(sys.argv[3]) if len(sys.argv) > 3 else 2500.0)
