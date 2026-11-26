"""Replay the exact reduction checks; opt into the longer complete periods.

    python3 explore/circuit-general/series-fibres/run_checks.py
    python3 explore/circuit-general/series-fibres/run_checks.py --examples --periods
"""
import argparse
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[3]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--examples", action="store_true",
                        help="replay every five/six-loop fibre cubature ordinate")
    parser.add_argument("--periods", action="store_true",
                        help="replay complete prism, K3,3 and five-loop fan periods")
    parser.add_argument("--reciprocal", action="store_true",
                        help="check certified reciprocal reconstruction and the K5 bridge (requires g++)")
    parser.add_argument("--k5-family", action="store_true",
                        help="recompute all three full-source six-loop K5 examples")
    arguments = parser.parse_args()
    scripts = ["test_treeflow.py", "test_series.py", "test_dirichlet.py",
               "test_support.py", "test_fibres.py"]
    for script in scripts:
        print("Checking " + script, flush=True)
        subprocess.run([sys.executable, str(ROOT / "code" / "tests" / script)],
                       cwd=ROOT, check=True)
    environment = dict(os.environ, HTK_FAST="1")
    subprocess.run([sys.executable, str(ROOT / "code/tests/test_engine.py"),
                    "test_tree_independence", "test_exterior_algebra_vs_determinant",
                    "test_sign_rules", "test_extraction_sign_formula",
                    "test_catalog_counts", "test_plan_independence"],
                   cwd=ROOT, env=environment, check=True)
    if arguments.examples:
        subprocess.run([sys.executable, str(Path(__file__).with_name("run_examples.py"))],
                       cwd=ROOT, check=True)
    if arguments.periods:
        subprocess.run([sys.executable,
                        str(ROOT / "explore/circuit-general/fiveloop/check_fk5.py")],
                       cwd=ROOT, check=True)
    if arguments.reciprocal:
        for script in ("series-fibres/test_reciprocal.py", "series-fibres/test_core_sample.py",
                       "series-fibres/test_core_numerator.py", "period-dispatch/test_k5_bridge.py"):
            print("Checking " + script, flush=True)
            subprocess.run([sys.executable, str(ROOT / "explore/circuit-general" / script)],
                           cwd=ROOT, check=True)
    if arguments.k5_family:
        subprocess.run([sys.executable, str(Path(__file__).with_name("run_k5_family.py")), "--recompute"],
                       cwd=ROOT, check=True)
    print("All requested exact checks passed.", flush=True)


if __name__ == "__main__":
    main()
