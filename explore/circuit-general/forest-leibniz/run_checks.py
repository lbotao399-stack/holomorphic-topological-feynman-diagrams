"""Replay exact forest Ward controls and optional complete certificates."""
import argparse
from pathlib import Path
import subprocess
import sys
from time import perf_counter

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--examples", action="store_true",
                        help="replay all current words, higher-loop lifts, and the full period combination")
    arguments = parser.parse_args()
    scripts = [
        HERE / "test_ward_local.py",
        HERE / "test_lifted_ward.py",
        HERE / "ward_three.py",
        HERE / "replay_existing_ward.py",
        ROOT / "explore/circuit-general/symmetry/check_conformal_combination.py",
        ROOT / "code/tests/test_forest_periods.py",
        ROOT / "code/tests/test_product_forests.py",
    ]
    if arguments.examples:
        scripts.extend([
            HERE / "certify_words.py",
            HERE / "lifted_ward.py",
            ROOT / "explore/circuit-general/symmetry/certify_forest_combination.py",
        ])
    start = perf_counter()
    for script in scripts:
        print("Running " + str(script.relative_to(ROOT)), flush=True)
        subprocess.run([sys.executable, str(script)], cwd=ROOT, check=True)
    print("All exact forest Ward checks passed (%.3f seconds)." %
          (perf_counter() - start), flush=True)


if __name__ == "__main__":
    main()
