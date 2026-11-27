#!/usr/bin/env python3
"""Replay bounded exact manuscript controls and verify the existing PDF build.

No network, source reintegration, period-table regeneration, or score changes.
--math-only explicitly omits the PDF gate and cannot report a full replay pass.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import re
import signal
import subprocess
import sys
import time
import traceback


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PRIOR = ROOT / "explore/circuit-general"
JOBS = (
    ("generated_tables", "notes/methods-paper/generate_tables.py", None),
    ("printed_tables", "notes/methods-paper/checks/printed_tables.py", None),
    ("fork_control", "notes/methods-paper/checks/fork_control.py", None),
    ("wheel_residue", "notes/methods-paper/checks/wheel_residue.py", None),
    ("fixed_core_stabilization", "explore/circuit-general/core-period-stabilization/run_checks.py", None),
    ("wheel_closure", "explore/circuit-general/core-period-closure/verify.py", None),
    ("relative_engine_tests", "code/tests/test_core_sector_periods.py", 10),
    ("relative_certificates", "explore/circuit-general/core-sector-reduction/certify.py", None),
    ("supported_source_replay", "explore/circuit-general/forest-source-currents/run_checks.py", None),
    ("cube_expansion_tests", "code/tests/test_forest_cube_expansion.py", 7),
    ("complete_cube_cli_tests", "code/tests/test_solve_product_forest_cubes.py", 6),
)
PACKAGES = ("sympy", "python-flint", "pypdf")


def load(path):
    return json.loads(path.read_text())


def sha256(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def relative(path):
    return str(path.relative_to(ROOT))


def write_record(path, record):
    temporary = path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def source_paths():
    """Pin only inputs, never this runner's manifests, logs, or PDF outputs."""
    paths = {HERE / name for name in ("README.md", "run_checks.py", "build.py",
                                     "generate_tables.py", "paper.tex", "supplement.tex",
                                     "preamble.tex", "references.bib")}
    for directory, pattern in (("sections", "*.tex"), ("generated", "*.tex"), ("checks", "*.py")):
        paths.update((HERE / directory).glob(pattern))
    records = [HERE / "verification" / (name + ".json")
               for name in ("printed_tables", "fork_control", "wheel_residue")]
    for name, data in (
        ("core-period-stabilization", ("checks.json",)),
        ("core-period-closure", ("certificate.json",)),
        ("core-sector-reduction", ("reductions.json",)),
        ("forest-source-currents", ("certificate.json", "fresh_nonunit_frame.json",
                                   "degenerate_frame.json")),
    ):
        directory = PRIOR / name
        paths.update(directory.glob("*.py"))
        paths.update(directory.glob("*.md"))
        records.extend(directory / item for item in data)
    paths.update(records)
    paths.update((ROOT / "code/htk").glob("*.py"))
    paths.add(ROOT / "code/htk/CONVENTIONS.md")
    paths.add(ROOT / "code/solve_product_forest_cubes.py")
    paths.update(ROOT / script for _, script, _ in JOBS)
    paths.update(PRIOR / "forest-cube-expansion" / name
                 for name in ("README.md", "DERIVATION.md"))
    hash_keys = {"dependencies_sha256", "dependency_sha256", "verification_sha256",
                 "source_sha256", "dependencies"}

    def dependencies(value, parent):
        if isinstance(value, dict):
            for key, item in value.items():
                if key in hash_keys and isinstance(item, dict):
                    for name, digest in item.items():
                        if isinstance(digest, str) and re.fullmatch(r"[0-9a-f]{64}", digest):
                            paths.add((ROOT if "/" in name else parent) / name)
                else:
                    dependencies(item, parent)
        elif isinstance(value, list):
            for item in value:
                dependencies(item, parent)

    for path in records:
        dependencies(load(path), path.parent)
    for path in paths:
        if not path.is_file() or not path.resolve().is_relative_to(ROOT):
            raise FileNotFoundError("Missing or external replay input: " + str(path))
    return sorted(paths)


def versions():
    result = {"python": platform.python_version(), "python_implementation": platform.python_implementation(),
              "python_executable": sys.executable, "platform": platform.platform(), "packages": {}}
    for name in (*PACKAGES, "pytest"):
        try:
            result["packages"][name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            result["packages"][name] = None
    for program in ("pdflatex", "bibtex"):
        try:
            reply = subprocess.run([program, "--version"], capture_output=True, text=True, timeout=10)
            result[program] = reply.stdout.splitlines()[0] if reply.returncode == 0 else None
        except (OSError, subprocess.TimeoutExpired):
            result[program] = None
    return result


def counts(name, stdout):
    if name == "generated_tables":
        return {"complete_generated_files": len(list((HERE / "generated").glob("*.tex")))}
    if name == "printed_tables":
        return {key: value for key, value in load(HERE / "verification/printed_tables.json").items()
                if isinstance(value, int)}
    if name == "fork_control":
        data = load(HERE / "verification/fork_control.json")
        return {"cube_variables": len(data["variables"]),
                "complete_first_current_spinors": len(data["all_first_current_spinors_after_bracket_cubed"]),
                "polynomial_kirchhoff_identities": len(data["all_polynomial_kirchhoff_identities"]),
                "integrated_kirchhoff_identities": len(data["all_integrated_kirchhoff_identities"]),
                "independent_original_network_controls": len(data["original_network_controls"])}
    if name == "wheel_residue":
        return load(HERE / "verification/wheel_residue.json")["summary"]
    if name == "fixed_core_stabilization":
        return load(PRIOR / "core-period-stabilization/checks.json")["summary"]
    if name == "wheel_closure":
        data = load(PRIOR / "core-period-closure/certificate.json")
        return {key: len(data[key]) for key in ("finite_k4_base_range", "rational_witness_controls",
                                               "k4_embedding_controls", "spanning_controls")}
    if name == "relative_certificates":
        return {key: value for key, value in json.loads(stdout).items() if isinstance(value, int)}
    if name == "supported_source_replay":
        data = load(PRIOR / "forest-source-currents/certificate.json")
        return {"complete_original_spinors": data["complete_original_spinors"],
                "complete_original_scalar_entries": data["complete_original_scalar_entries"],
                "frame_coefficients": data["all_interpolation_entries_verified"],
                "interpolation_nodes": len(data["samples"]), "direct_nodes_stored": len(data["direct_points"]),
                "independent_complete_controls_stored": len(data["independent_controls"]),
                "residual_periods": len(data["residual_periods"])}
    if name == "pdf_build":
        return {key + "_pages": value["pages"]
                for key, value in load(HERE / "verification/build.json")["documents"].items()}
    return {}


def run_job(name, arguments, tests, logdir, timeout, blocked=None):
    output, error = logdir / (name + ".stdout.log"), logdir / (name + ".stderr.log")
    command = [sys.executable, *arguments]
    row = {"name": name, "command": command, "cwd": str(ROOT),
           "stdout": relative(output), "stderr": relative(error), "timeout_seconds": timeout}
    start = time.monotonic()
    env = dict(os.environ)
    env.pop("PYTHONOPTIMIZE", None)
    env.update(PYTHONDONTWRITEBYTECODE="1", PYTHONHASHSEED="0", PYTHONUNBUFFERED="1")
    with output.open("w") as out, error.open("w") as err:
        if blocked:
            err.write(blocked + "\n")
            row.update(status="blocked", returncode=None, error=blocked)
        else:
            process = subprocess.Popen(command, cwd=ROOT, env=env, stdout=out, stderr=err,
                                       start_new_session=True)
            try:
                returncode = process.wait(timeout=timeout)
                row.update(status="passed" if returncode == 0 else "failed", returncode=returncode)
            except (subprocess.TimeoutExpired, KeyboardInterrupt) as exception:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
                row.update(status="timeout" if isinstance(exception, subprocess.TimeoutExpired) else "interrupted",
                           returncode=process.returncode, error=type(exception).__name__)
    row["elapsed_seconds"] = round(time.monotonic() - start, 3)
    if row["status"] == "passed":
        try:
            if tests is not None:
                report = error.read_text() + output.read_text()
                match = re.search(r"Ran (\d+) tests?", report)
                if not match or int(match[1]) != tests or "skipped=" in report:
                    raise ArithmeticError("The complete expected unittest count did not pass without skips.")
                row["counts"] = {"tests_passed": int(match[1])}
            else:
                row["counts"] = counts(name, output.read_text())
        except Exception as exception:
            row.update(status="failed", error=str(exception))
    row["stdout_sha256"], row["stderr_sha256"] = sha256(output), sha256(error)
    print(name + ": " + row["status"] + " (" + str(row["elapsed_seconds"]) + " s)", flush=True)
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--math-only", action="store_true",
                        help="replay mathematics without claiming the complete PDF/build gate")
    parser.add_argument("--timeout", type=int, default=300, help="positive per-command timeout in seconds")
    args = parser.parse_args()
    if args.timeout <= 0:
        parser.error("--timeout must be positive")
    mode = "math_only" if args.math_only else "full"
    record_path = HERE / "verification" / ("math_checks.json" if args.math_only else "checks.json")
    logdir = HERE / "verification" / ("math_check_logs" if args.math_only else "check_logs")
    logdir.mkdir(parents=True, exist_ok=True)
    record = {"schema": "methods-paper-replay-v1", "mode": mode, "status": "running",
              "started_utc": datetime.now(timezone.utc).isoformat(),
              "full_replay_passed": False, "commands": [], "versions": versions(),
              "expensive_source_reintegration_performed": False,
              "prior_period_tables_regenerated": False,
              "scope": "Bounded exact controls, complete stored-certificate bindings, and existing PDF build verification; independent mathematical and visual review remain separate."}
    write_record(record_path, record)
    try:
        missing = [name for name in PACKAGES if record["versions"]["packages"][name] is None]
        if missing:
            raise RuntimeError("Missing required installed packages: " + ", ".join(missing))
        paths = source_paths()
        record["sources_sha256"] = {relative(path): sha256(path) for path in paths}
        for name, script, tests in JOBS:
            row = run_job(name, [script], tests, logdir, args.timeout)
            record["commands"].append(row)
            write_record(record_path, record)
            if row["status"] != "passed":
                raise RuntimeError("Replay failed: " + name + "; inspect the recorded logs.")
        if not args.math_only:
            required = [HERE / "verification/build.json", HERE / "output/pdf/paper.pdf",
                        HERE / "output/pdf/supplement.pdf"]
            absent = [relative(path) for path in required if not path.is_file()]
            row = run_job("pdf_build", ["notes/methods-paper/build.py", "--verify"], None, logdir,
                          args.timeout, "Missing required build artifacts: " + ", ".join(absent) if absent else None)
            record["commands"].append(row)
            if row["status"] != "passed":
                raise RuntimeError("Full replay requires both PDFs and a successful build.py --verify.")
            record["build_artifacts_sha256"] = {relative(path): sha256(path) for path in required}
        after = {relative(path): sha256(path) for path in source_paths()}
        if after != record["sources_sha256"]:
            record["sources_sha256_after"] = after
            raise RuntimeError("Replay input files changed during the run; rerun after they are frozen.")
        record["status"] = "math-only-passed" if args.math_only else "passed"
        record["full_replay_passed"] = not args.math_only
    except Exception as exception:
        record.update(status="failed", error=str(exception))
        traceback.print_exc()
    finally:
        record["finished_utc"] = datetime.now(timezone.utc).isoformat()
        record["commands_passed"] = sum(row["status"] == "passed" for row in record["commands"])
        record["commands_run"] = len(record["commands"])
        write_record(record_path, record)
    print("Manifest: " + relative(record_path), flush=True)
    print("Result: " + record["status"], flush=True)
    return 0 if record["status"] in ("passed", "math-only-passed") else 1


if __name__ == "__main__":
    sys.exit(main())
