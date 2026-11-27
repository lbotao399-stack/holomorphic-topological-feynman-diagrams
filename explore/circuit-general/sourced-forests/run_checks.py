"""Run exact sourced-forest controls and optionally replay all finite certificates."""
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
                        help="recompute all exact values and certificate implementation hashes")
    arguments = parser.parse_args()
    commands = [
        ROOT / "code/tests/test_forest_contractions.py",
        ROOT / "code/tests/test_forest_periods.py",
        ROOT / "code/tests/test_product_forests.py",
        HERE / "tests_exact.py",
        HERE / "independent_six.py",
        HERE / "comb_reference.py",
    ]
    if arguments.examples:
        commands.append(HERE / "certify.py")
        commands.append(HERE / "certify_missing_period.py")
    start = perf_counter()
    for script in commands:
        print("Running " + str(script.relative_to(ROOT)), flush=True)
        subprocess.run([sys.executable, str(script)], cwd=ROOT, check=True)
    print("All exact sourced-forest checks passed (%.3f seconds)." %
          (perf_counter() - start), flush=True)


if __name__ == "__main__":
    main()
