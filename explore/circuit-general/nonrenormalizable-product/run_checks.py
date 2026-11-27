"""Replay complete higher-superpotential Wick calculations and independent checks."""
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
                        help="regenerate and compare both complete all-source certificates")
    arguments = parser.parse_args()
    scripts = [HERE / "test_quartic_wick.py", HERE / "test_nonabelian_wick.py"]
    if arguments.examples:
        scripts.extend([HERE / "quartic_wick.py", HERE / "nonabelian_wick.py"])
    start = perf_counter()
    for script in scripts:
        print("Running " + str(script.relative_to(ROOT)), flush=True)
        subprocess.run([sys.executable, str(script)], cwd=ROOT, check=True)
    print("All complete higher-superpotential product checks passed (%.3f seconds)." %
          (perf_counter() - start), flush=True)


if __name__ == "__main__":
    main()
