"""Replay sourced-ear local/original-graph checks and complete certificates."""
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
                        help="also replay all complete kernels and spinors through eleven loops")
    arguments = parser.parse_args()
    scripts = [ROOT / "code/tests/test_forest_ears.py", HERE / "sourced_ward.py"]
    if arguments.examples:
        scripts.append(HERE / "certify.py")
    start = perf_counter()
    for script in scripts:
        print("Running " + str(script.relative_to(ROOT)), flush=True)
        subprocess.run([sys.executable, str(script)], cwd=ROOT, check=True)
    print("All exact sourced-ear checks passed (%.3f seconds)." % (perf_counter() - start), flush=True)


if __name__ == "__main__":
    main()
