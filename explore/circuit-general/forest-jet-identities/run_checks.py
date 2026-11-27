"""Replay the complete current identities and rational boundary algorithm.

The default compares fresh boundary integrals against the complete saved
rank-two, rank-three, and mixed-input original-graph tensors.  --recompute
also regenerates the independent sixteen-entry mixed-input period oracle
from the original graphs; it does not replace or shorten that tensor.
"""

import argparse
from pathlib import Path
import subprocess
import sys
from time import perf_counter

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
JETS = HERE.parent / "universal-product-jets"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--recompute", action="store_true",
                        help="also recompute every original graph entry in the full mixed-input tensor")
    args = parser.parse_args()
    commands = [HERE / "check_identities.py", JETS / "routing.py", JETS / "test_boundary_product.py"]
    if args.recompute:
        commands.append(HERE / "mixed_input_control.py")
    start = perf_counter()
    for path in commands:
        print("\nRunning", path.relative_to(ROOT), flush=True)
        subprocess.run([sys.executable, str(path)], cwd=ROOT, check=True)
    print("\nAll complete source-simplex, boundary, and original-current checks passed.", flush=True)
    print("Elapsed %.3f seconds." % (perf_counter()-start), flush=True)


if __name__ == "__main__":
    main()
