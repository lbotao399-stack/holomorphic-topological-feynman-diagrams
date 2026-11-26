#!/usr/bin/env python3
"""Replay all four exact all-order holography certificates."""
from hashlib import sha256
from pathlib import Path
import json
import platform
import subprocess
import sys

import sympy


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CASES = (
    ("boundary/all_order_join.py", "boundary/all_order_join.json"),
    ("matching/ribbon_dictionary.py", "matching/ribbon_dictionary.json"),
    ("bulk/primary_tree_test.py", "bulk/primary_tree_test.json"),
    ("resummation/all_order_response.py", "resummation/all_order_response.json"),
)


def main():
    records = []
    for script, certificate in CASES:
        target = HERE / certificate
        before = target.read_bytes() if target.exists() else None
        result = subprocess.run(
            [sys.executable, str(HERE / script)], cwd=ROOT,
            capture_output=True, text=True, check=False,
        )
        if result.returncode:
            sys.stderr.write(result.stdout + result.stderr)
            raise SystemExit(f"Failed: {script}")
        after = target.read_bytes()
        if before is not None and before != after:
            raise SystemExit(f"Certificate changed on replay: {certificate}")
        data = json.loads(after)
        if data.get("status") != "PASS":
            raise SystemExit(f"No PASS status: {certificate}")
        records.append({
            "script": script,
            "script_sha256": sha256((HERE / script).read_bytes()).hexdigest(),
            "certificate": certificate,
            "certificate_sha256": sha256(after).hexdigest(),
            "saved_output_reproduced": before is not None,
            "status": "PASS",
        })
        print(f"PASS {script}", flush=True)
    output = {
        "status": "PASS",
        "python": platform.python_version(),
        "sympy": sympy.__version__,
        "independent_certificates": len(records),
        "checks": records,
        "scope": "Boundary all-order detector, ribbon topology, linear bulk states, and Ward targets; full nonlinear bulk amplitude remains open.",
    }
    (HERE / "replay.json").write_text(
        json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(f"PASS all {len(records)} exact certificates", flush=True)


if __name__ == "__main__":
    main()
