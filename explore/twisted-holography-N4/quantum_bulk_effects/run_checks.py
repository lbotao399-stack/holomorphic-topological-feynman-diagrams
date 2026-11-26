#!/usr/bin/env python3
"""Replay independent exact calculations and compare committed certificates.

Three constructions test boundary-derived predictions, one tests a classical
bulk kernel, and one tests spin selection and finite Hilbert complexes.
Passing these checks does not constitute a quantum bulk amplitude match.
"""

from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory


ROOT = Path(__file__).resolve().parent
JOBS = (
    ("branes/determinant_color.py", "branes/determinant_color.json"),
    ("branes/determinant_grassmann_check.py", "branes/determinant_grassmann_check.json"),
    ("branes/auxiliary_gauge_check.py", "branes/auxiliary_gauge_check.json"),
    ("bulk/flux_derived_kernel.py", "bulk/flux_derived_kernel.json"),
    ("spectrum/verify.py", "spectrum/results.json"),
)


def replay(job):
    script_name, certificate_name = job
    script = ROOT / script_name
    expected = (ROOT / certificate_name).read_bytes()
    source = script.read_bytes()
    with TemporaryDirectory(prefix="quantum_bulk_effects_") as scratch:
        isolated_script = Path(scratch) / script.name
        isolated_script.write_bytes(source)
        process = subprocess.run(
            [sys.executable, str(isolated_script)],
            cwd=scratch,
            text=True,
            capture_output=True,
            check=False,
        )
        if process.returncode:
            raise RuntimeError(
                f"{script_name} failed ({process.returncode})\n"
                f"{process.stdout}\n{process.stderr}"
            )
        generated = (Path(scratch) / Path(certificate_name).name).read_bytes()
        if generated != expected:
            raise AssertionError(f"Certificate changed: {certificate_name}")
        json.loads(generated)
    return {
        "script": script_name,
        "script_sha256": sha256(source).hexdigest(),
        "certificate": certificate_name,
        "certificate_sha256": sha256(expected).hexdigest(),
        "byte_for_byte_match": True,
    }


def main():
    with ThreadPoolExecutor(max_workers=4) as pool:
        results = list(pool.map(replay, JOBS))
    report = {
        "status": "PASS",
        "scope": "Exact boundary source predictions, classical bulk kernel, and spectral controls; no independent quantum bulk Q_n match.",
        "verification": results,
    }
    (ROOT / "replay.json").write_text(json.dumps(report, indent=2) + "\n")
    for result in results:
        print("PASS", result["script"], "committed certificate reproduced")
    print("PASS", len(results), "independent exact calculations")


if __name__ == "__main__":
    main()
