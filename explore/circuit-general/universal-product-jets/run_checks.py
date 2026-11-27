"""Replay complete higher-jet and rational-boundary product checks."""
from pathlib import Path
import argparse
import subprocess
import sys
from time import perf_counter

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--recompute", action="store_true", help="also redo all original six-graph integrations")
    args = parser.parse_args()
    start = perf_counter()
    commands = [("routing.py",), ("verify_certificates.py",), ("test_covariants.py",), ("render_jets.py",), ("test_boundary_product.py",)]
    if args.recompute:
        commands.append(("jet_tensor.py", "--orders", "2", "3"))
    for command in commands:
        print("Running " + " ".join(command), flush=True)
        subprocess.run([sys.executable, str(HERE / command[0]), *command[1:]], check=True)
    print("All complete tensor, original-current, and rational-boundary checks passed.", flush=True)
    print("Elapsed %.3f seconds." % (perf_counter() - start), flush=True)


if __name__ == "__main__":
    main()
