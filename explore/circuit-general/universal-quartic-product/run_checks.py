"""Replay the complete universal tensor, higher-W classification, and printed formulas."""
from pathlib import Path
import subprocess
import sys
from time import perf_counter

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def main():
    start = perf_counter()
    for name in ("quartic_tensor.py", "higher_W.py", "test_quartic_tensor.py", "render_tensor.py"):
        print("Running " + str((HERE / name).relative_to(ROOT)), flush=True)
        subprocess.run([sys.executable, str(HERE / name)], cwd=ROOT, check=True)
    print("All universal tensor and higher-W controls passed (%.3f seconds)." % (perf_counter() - start), flush=True)


if __name__ == "__main__":
    main()
