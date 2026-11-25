"""
Third overnight step: after igl_synth_pipeline2.py, refine the medial-white-matter model
(foldma, the model that reached your controls) against the EZH2 cKO sections too, so the
control and cKO parameter ranges can be compared within one model; then rebuild the report.

    python scripts/igl_synth_pipeline3.py
"""
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PY = sys.executable


def sh(*args):
    t0 = time.time()
    print(">>", " ".join(args), flush=True)
    r = subprocess.run([PY, *args], cwd=ROOT)
    print(f"<< exit {r.returncode} after {time.time() - t0:.0f}s", flush=True)


def main():
    while subprocess.run(["pgrep", "-f", r"igl_synth_pipeline2\.py"], capture_output=True, text=True).stdout.strip():
        time.sleep(30)
    print("second round finished", time.ctime(), flush=True)
    fam, t = "foldma", "cKO"
    sh("scripts/igl_synth_gen.py", "refine", fam, t, "3000", "9")
    sh("scripts/igl_synth_gen.py", "refine_iter", fam, t, "3000", "prior,refine", "refine2")
    sh("scripts/igl_synth_gen.py", "refine_iter", fam, t, "3000", "prior,refine,refine2", "refine3")
    sh("scripts/igl_synth_gen.py", "refine_iter", fam, t, "3000", "prior,refine,refine2,refine3", "refine4")
    sh("scripts/igl_synth_report.py")
    print("all done", time.ctime(), flush=True)


if __name__ == "__main__":
    main()
