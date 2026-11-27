"""Replay the complete universal-forest implementation and exact certificates."""
import json
from pathlib import Path
import subprocess
import sys
from time import perf_counter

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def main():
    import flint
    import networkx
    import sympy
    print("Python %s; python-flint %s; networkx %s; sympy %s." %
          (sys.version.split()[0], flint.__version__, networkx.__version__, sympy.__version__), flush=True)
    start = perf_counter()
    scripts = [
        ROOT / "code/tests/test_forest_kron.py",
        ROOT / "code/tests/test_forest_reduction.py",
        ROOT / "code/tests/test_forest_cubes.py",
        HERE / "certify_kron.py",
        HERE / "certify_integrals.py",
        HERE / "certify_cubes.py",
        HERE / "composed_five_loop.py",
    ]
    for script in scripts:
        print("\nRunning %s" % script.relative_to(ROOT), flush=True)
        subprocess.run([sys.executable, str(script)], cwd=ROOT, check=True)
    for name in ("five_loop", "sourced_fork"):
        completed = subprocess.run([sys.executable, str(ROOT / "code/solve_product_forest.py"),
                                    str(HERE / (name + ".json")), "--finite-only"],
                                   cwd=ROOT, check=True, text=True, capture_output=True)
        if json.loads(completed.stdout) != json.loads((HERE / (name + "_solution.json")).read_text()):
            raise ArithmeticError("a complete command-line result changed")
        print("Verified the complete %s command-line result." % name, flush=True)
    from render_results import render
    if render() != (HERE / "FULL_RESULTS.md").read_text():
        raise ArithmeticError("the printed full-spinor results differ from the exact coefficient archive")
    print("Every printed bare kernel and full first-current spinor replays.", flush=True)
    print("Complete universal-forest replay passed in %.3f seconds." % (perf_counter() - start), flush=True)


if __name__ == "__main__":
    main()
