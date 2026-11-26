#!/usr/bin/env python3
"""Replay the five independent exact certificates for this research note.

Existing certificate JSON must reproduce byte for byte. The runner does not
reintegrate the inherited worldsheet periods or the original two-loop anomaly.
"""
from hashlib import sha256
from pathlib import Path
import json
import platform
import subprocess
import sys

import sympy


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CASES = (
    ("explore/sugra-11d-qme/developments_2026_11/further_twist_symbols.py",
     "explore/sugra-11d-qme/developments_2026_11/further_twist_symbols.json"),
    ("explore/twisted-holography-N4/developments_2026_11/delayed_bulk_response.py",
     "explore/twisted-holography-N4/developments_2026_11/delayed_bulk_response.json"),
    ("explore/iib-bcov/developments_2026_11/planar_ward_rigidity.py",
     "explore/iib-bcov/developments_2026_11/planar_ward_rigidity.json"),
    ("notes/twisted-sugra-bulk-tests/checks/compact_momentum_rank.py",
     "notes/twisted-sugra-bulk-tests/verification/compact_momentum_rank.json"),
    ("notes/twisted-sugra-bulk-tests/checks/cross_sector_checks.py",
     "notes/twisted-sugra-bulk-tests/verification/cross_sector_checks.json"),
)


def main():
    records = []
    for script, certificate in CASES:
        destination = ROOT / certificate
        before = destination.read_bytes() if destination.exists() else None
        completed = subprocess.run(
            [sys.executable, str(ROOT / script)], cwd=ROOT,
            text=True, capture_output=True, check=False,
        )
        if completed.returncode:
            sys.stderr.write(completed.stdout + completed.stderr)
            raise SystemExit(f"Certificate failed: {script}")
        after = destination.read_bytes()
        if before is not None and before != after:
            raise SystemExit(f"Saved certificate changed on replay: {certificate}")
        data = json.loads(after)
        if data.get("status") != "PASS":
            raise SystemExit(f"Certificate did not report PASS: {certificate}")
        records.append({
            "script": script,
            "script_sha256": sha256((ROOT / script).read_bytes()).hexdigest(),
            "certificate": certificate,
            "certificate_sha256": sha256(after).hexdigest(),
            "status": "PASS",
            "saved_output_reproduced": before is not None,
        })
        print(f"PASS {script}", flush=True)
    result = {
        "status": "PASS",
        "python": platform.python_version(),
        "sympy": sympy.__version__,
        "independent_certificates": len(records),
        "checks": records,
        "scope": "Exact algebraic and combinatorial replay; inherited periods are inputs.",
    }
    destination = HERE / "verification/replay.json"
    destination.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(f"PASS all {len(records)} exact certificates", flush=True)


if __name__ == "__main__":
    main()
