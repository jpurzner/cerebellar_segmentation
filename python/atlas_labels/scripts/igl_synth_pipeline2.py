"""
Second overnight round for the synthetic-IGL study (after igl_synth_pipeline.py).

    python scripts/igl_synth_pipeline2.py

1. waits for the first round to finish
2. prior sampling of fold3 (fold2 with the EGL-squeeze and jitter corrections)
3. ABC over all families
4. local search (4 rounds) for the two closest families against your controls and against
   Allen P7, and the closest family against the EZH2 cKO sections, so the parameter shifts
   between control-like and cKO-like IGLs can be compared. The first round's finger
   refinement is redone here: n_sub was perturbed to a non-integer, which silently discarded
   every refined sample with a split lobule.
5. report data and page
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


def first_round_running():
    r = subprocess.run(["pgrep", "-f", r"igl_synth_pipeline\.py"], capture_output=True, text=True)
    return bool(r.stdout.strip())


def ranked(abc, target):
    fams = abc["families"]
    return sorted(fams, key=lambda f: (-fams[f]["targets"][target]["sensible"], fams[f]["targets"][target]["best_distance"]))


def main():
    while first_round_running():
        time.sleep(30)
    print("first round finished", time.ctime(), flush=True)
    sh("scripts/igl_synth_gen.py", "run", "fold3", "6000", "9")
    sh("scripts/igl_synth_abc.py")
    abc = json.loads((OUT / "abc.json").read_text())
    plan = []
    for t, k in (("IF", 2), ("Allen", 2), ("cKO", 1)):
        for fam in ranked(abc, t)[:k]:
            plan.append((fam, t))
    print("plan", plan, flush=True)
    (OUT / "refine_plan.json").write_text(json.dumps(plan))
    for fam, t in plan:
        print("refine target", t, "family", fam, abc["families"][fam]["targets"][t]["best_distance"], flush=True)
        sh("scripts/igl_synth_gen.py", "refine", fam, t, "3000", "9")
        sh("scripts/igl_synth_gen.py", "refine_iter", fam, t, "3000", "prior,refine", "refine2")
        sh("scripts/igl_synth_gen.py", "refine_iter", fam, t, "3000", "prior,refine,refine2", "refine3")
        sh("scripts/igl_synth_gen.py", "refine_iter", fam, t, "3000", "prior,refine,refine2,refine3", "refine4")
    sh("scripts/igl_synth_report.py")
    print("all done", time.ctime(), flush=True)


if __name__ == "__main__":
    main()
