"""Run the exact forest/fan controls and optional full high-loop certificates."""
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
                        help="recompute every stored higher-loop forest and fan coefficient")
    args = parser.parse_args()
    commands = [
        [ROOT / "code/tests/test_forest_model.py"],
        [ROOT / "code/tests/test_fan_moments.py"],
        [HERE / "tests_exact.py"],
        [ROOT / "code/tests/test_product_forests.py"],
        [ROOT / "code/tests/test_known.py", "test_S1_fans", "test_S1_fan_general_sources",
         "test_S1_triangle_displacements"],
    ]
    if args.examples:
        commands += [[HERE / "certify_forests.py", "--check"], [HERE / "certify_fans.py"]]
    start = perf_counter()
    for command in commands:
        print("Running " + " ".join(str(item.relative_to(ROOT)) if isinstance(item, Path)
                                    else item for item in command), flush=True)
        subprocess.run([sys.executable] + list(map(str, command)), cwd=ROOT, check=True)
    print("All exact forest/fan checks passed (%.3f s)." % (perf_counter() - start), flush=True)


if __name__ == "__main__":
    main()
