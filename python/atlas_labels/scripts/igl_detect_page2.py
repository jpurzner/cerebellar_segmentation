"""
New IGL Step, round 2 page: every section starts from its best labels so far (controls: the
repaired labels; the 6 cKO sections the user corrected in round 1: their corrected labels;
s1_5: MATLAB's labels) and shows only the round-2 step's changes on top (green added, orange
removed: igl_detect_round.compose), editable with paint / erase. Round 1 vs round 2 on held-out
sections from outputs/igl_detect/round2_eval.json.

Proposal PNG on the 4 um grid: 0 none, 60 removed by the step, 100 IGL as it stands, 200 added.
db collection "igl_detect_r2", one doc per section: {strokes, verdict, note}; replay with
igl_detect_edits.corrected_igl(proposal >= 80, strokes).

    python scripts/igl_detect_page2.py      # -> outputs/igl_detect/new_igl_step_round2.html
"""
import base64
import io
import json
import sys
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage as ndi

sys.path.insert(0, str(Path(__file__).parent))
import igl_detect as D  # noqa: E402
import igl_detect_round as R  # noqa: E402
import igl_gapfill_apply as GA  # noqa: E402
import igl_gapfill_page as PG  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/igl_detect"
TEMPLATE = Path(__file__).parent / "igl_detect_template2.html"
G = D.G
KK = 2


def b64(img, fmt, **kw):
    b = io.BytesIO()
    img.save(b, fmt, **kw)
    return f"data:image/{fmt.lower()};base64," + base64.b64encode(b.getvalue()).decode()


def edge(m):
    return m & ~ndi.binary_erosion(m)


def main():
    E = json.loads((OUT / "round2_eval.json").read_text())
    edited = D.edited_cko()
    secs = []
    for sid in GA.section_list():
        s = GA.short(sid)
        group = GA.info(sid)[1]
        _, y, auto, tissue = D.load(s)
        if group == "control":
            base, base_name = y == 4, "repaired labels"
        elif s in edited:
            base, base_name = y == 4, "your round-1 edits"
        else:
            base, base_name = auto == 4, "MATLAB's labels"
        # the better model per group on held-out sections (round2_eval.json): round 1 (trained on the
        # 12 controls) for the cKO sections, round 2 (+ the user's 6 corrected cKO sections) for controls
        rnd = 1 if group == "cKO" else 2
        P = R.probs(s, rnd)
        new, add, rem = R.compose(P, base, tissue)
        prop = np.zeros(base.shape, np.uint8)
        prop[rem] = 60
        prop[base & ~rem] = 100
        prop[add] = 200
        S, px = D.stains(sid)
        H, W = base.shape
        img = PG.native_crop(np.stack([S["nuclei"]]), px, (0, H, 0, W), 0)
        up = lambda m: np.repeat(np.repeat(m, KK, 0), KK, 1)
        ov = np.zeros((H * KK, W * KK, 4), np.uint8)
        ov[up(edge(base))] = (59, 224, 255, 255)
        before = (base != (auto == 4)) & tissue
        if before.any():
            ov[up(ndi.binary_dilation(edge(before)))] = (255, 47, 214, 230)
        a = lambda m: int(m.sum() * G * G)
        needs = "check the changes" if (a(add) + a(rem)) >= 20000 else ("small changes only" if (a(add) + a(rem)) else "no changes")
        secs.append(dict(short=s, group=group, mag="20x" if GA.info(sid)[2] else "10x", w_um=round(W * G), h_um=round(H * G), grid=[H, W],
                         base=b64(Image.fromarray(img.astype(np.uint8)), "JPEG", quality=82),
                         outline=b64(Image.fromarray(ov, "RGBA"), "PNG", optimize=True),
                         proposed=b64(Image.fromarray(prop, "L"), "PNG", optimize=True), needs=needs,
                         stats=dict(base=base_name, model=f"round {rnd}", igl_mm2=round(a(new) / 1e6, 3), add_um2=a(add), rem_um2=a(rem))))
        print(s, secs[-1]["stats"], needs, flush=True)
    sm = E["summary"]
    ck, ct = sm["cKO"], sm["control"]
    gold = next(r for r in E["sections"] if r["section"] == "s1_2")
    findings = [
        f"<b>Learning from your cKO edits did not help the cKO sections.</b> On your 6 corrected cKO sections, each held out of its own model, round 2 reproduces "
        f"{round(100 * ck['r2_reproduced'])}% of the cells you changed and overlaps your corrected IGL at {ck['r2_iou']:.3f}; round 1, which learned from the controls only, "
        f"reproduces {round(100 * ck['r1_reproduced'])}% at {ck['r1_iou']:.3f} (MATLAB {ck['iou_matlab']:.3f}), and changes less you left alone "
        f"({ck['r1_extra_um2']:,} against {ck['r2_extra_um2']:,} µm² per section). Likely reason: outside the areas you painted, those sections still carry MATLAB's cKO labels, "
        "so round 2 also learned from IGL that is really EGL wherever you did not erase it. (Round 1's removal thresholds were tuned on these same edits, which flatters it a little.)",
        f"<b>On the controls round 2 is slightly better</b>: IGL overlap {ct['r2_iou']:.3f} against {ct['r1_iou']:.3f}, about {ct['r2_extra_um2']:,} µm² of changes you did not make "
        f"per section against {ct['r1_extra_um2']:,} (s6_G's 17,700 µm² of unwanted additions are gone); s1_2 against your gold standard {gold['r2']['iou']:.3f} (MATLAB {gold['iou_matlab']:.3f}).",
        "<b>So below, each group uses its better model</b>: round 1 on the cKO sections, round 2 on the controls. Every section starts from its best labels so far "
        "(controls: the repaired labels; your 6 cKO sections: your round-1 edits; s1_5: MATLAB) and shows only the step's changes on top: green added, orange removed "
        "(IGL it calls EGL). &ldquo;Check the changes&rdquo; marks sections with at least 20,000 µm² of them.",
        "<b>Fix what is wrong with the brush and eraser</b>, or mark the section fine: painting over orange keeps that IGL, erasing green drops the addition. "
        "To learn from cKO edits properly, the next round should train only on the cells you actually edited, not on the rest of those sections.",
    ]
    data = dict(collection="igl_detect_r2", sections=secs, eval=E, findings=findings)
    page = TEMPLATE.read_text().replace("__DATA__", json.dumps(data))
    dst = OUT / "new_igl_step_round2.html"
    dst.write_text(page)
    print(dst, round(dst.stat().st_size / 1e6, 2), "MB")


if __name__ == "__main__":
    main()
