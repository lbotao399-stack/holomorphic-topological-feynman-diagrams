"""Run focused independent tests and replay the complete geometric certificate."""

from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]

if __name__ == "__main__":
    for script in (ROOT / "code/tests/test_forest_cube_expansion.py",
                   ROOT / "code/tests/test_solve_product_forest_cubes.py", HERE / "certify.py"):
        print("Running " + str(script.relative_to(ROOT)), flush=True)
        subprocess.run([sys.executable, str(script)], cwd=ROOT, check=True)
