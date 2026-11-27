"""Validate every stored finite coefficient; optionally replay the heavy audits.

The quick mode checks complete census inputs, all stored source/current
identities, all printed tensors, exact archive digests, and independently
integrated controls. It does not claim to recompute all 200 integrals.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from time import perf_counter

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def dependencies(name):
    data = json.loads((HERE / name).read_text())
    if not data.get("implementation_sha256"):
        raise ArithmeticError("an independent audit has no exact implementation hashes")
    for relative, expected in data["implementation_sha256"].items():
        path = (ROOT / relative).resolve()
        if not path.is_relative_to(ROOT) or hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise ArithmeticError("an independent audit dependency changed: " + relative)
    print("Every exact dependency hash verified for %s." % name, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--deep", action="store_true")
    options = parser.parse_args()
    started = perf_counter()
    subprocess.run([sys.executable, str(HERE / "check_reconstruction.py")], cwd=ROOT, check=True)
    dependencies("independent_core_certificate.json")
    dependencies("chain_first_audit.json")
    subprocess.run([sys.executable, str(HERE / "independent_core.py"), "--compare-only"], cwd=ROOT, check=True)
    from render_results import render
    if render() != (HERE / "FULL_RESULTS.md").read_text():
        raise ArithmeticError("a printed complete spinor differs from its exact stored entries")
    print("All 200 printed bare kernels and all 2200 printed original-current spinors agree exactly.", flush=True)
    if options.deep:
        for name in ("independent_core.py", "check_chain_first.py"):
            subprocess.run([sys.executable, str(HERE / name)], cwd=ROOT, check=True)
    print("Complete stored five-loop census validation passed in %.3f seconds." %
          (perf_counter() - started), flush=True)


if __name__ == "__main__":
    main()
