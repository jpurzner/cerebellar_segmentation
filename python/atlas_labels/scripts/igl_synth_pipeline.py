"""
Overnight driver for the synthetic-IGL study.

    python scripts/igl_synth_pipeline.py

1. prior sampling of the later families (foldma, finger, fold2; resumes where stopped)
2. ABC over all families
3. local search around the closest parameter sets of the best family for each target:
   one round from the prior samples, then two more rounds restarting from the best found
4. report data and page
"""
import json
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/igl_synth"
PY = sys.executable


def sh(*args):
    t0 = time.time()
    print(">>", " ".join(args), flush=True)
    r = subprocess.run([PY, *args], cwd=ROOT)
    print(f"<< exit {r.returncode} after {time.time() - t0:.0f}s", flush=True)
    return r.returncode


def best_family(abc, target):
    fams = abc["families"]
    return min(fams, key=lambda f: (-fams[f]["targets"][target]["sensible"], fams[f]["targets"][target]["best_distance"]))


def main():
    sh("scripts/igl_synth_gen.py", "run", "foldma", "8000", "9")
    sh("scripts/igl_synth_gen.py", "run", "finger", "12000", "9")
    sh("scripts/igl_synth_gen.py", "run", "fold2", "8000", "9")
    sh("scripts/igl_synth_abc.py")
    abc = json.loads((OUT / "abc.json").read_text())
    for t in ("IF", "Allen"):
        fam = best_family(abc, t)
        print("refine target", t, "family", fam, abc["families"][fam]["targets"][t]["best_distance"], flush=True)
        sh("scripts/igl_synth_gen.py", "refine", fam, t, "3000", "9")
        sh("scripts/igl_synth_gen.py", "refine_iter", fam, t, "3000", "prior,refine", "refine2")
        sh("scripts/igl_synth_gen.py", "refine_iter", fam, t, "3000", "prior,refine,refine2", "refine3")
    sh("scripts/igl_synth_report.py")


if __name__ == "__main__":
    main()
