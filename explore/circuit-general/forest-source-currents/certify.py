"""Reconstruct every supported-source current coefficient from exact samples."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
from time import perf_counter

from source_currents import (HERE, ROOT, EDGES, SOURCE_SUPPORT, ZERO_SOURCE_VERTICES,
    POINTS, DIRECT_POINTS, BASIS, OPTIONS, calculate_sources, canonical, digest,
    evaluate, interpolate, sources_at, swap_frame_record)

DEPENDENCIES = (
    "code/htk/forest_periods.py", "code/htk/forest_chain_periods.py",
    "code/htk/_forest_period_worker.py", "code/htk/treeflow.py",
    "code/htk/brute.py", "code/htk/series.py",
    "explore/circuit-general/symmetry/tree_minor.py",
    "explore/nonSP/verify/k4_periods/current_minor.py",
    "explore/circuit-general/fiveloop/per5.py",
    "explore/circuit-general/fiveloop/msys5.py",
    "explore/circuit-general/cores/gutil.py",
    "explore/circuit-general/cores/tables/prism_q7.json",
    "explore/circuit-general/assembly/tables/prism_q8.json",
    "explore/circuit-general/forest-source-currents/source_currents.py",
)


def dependency_hashes():
    return {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in DEPENDENCIES}


def verification_hashes():
    names = ("certify.py", "independent_controls.py", "verify.py", "render.py", "run_checks.py")
    return {str((HERE / name).relative_to(ROOT)): hashlib.sha256((HERE / name).read_bytes()).hexdigest()
            for name in names}


def atomic_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, sort_keys=True, indent=2) + "\n")
    temporary.replace(path)


def point_job(point, checkpoint):
    path = Path(checkpoint) / ("point_%d_%d.json" % point)
    hashes = dependency_hashes()
    if path.exists():
        value = json.loads(path.read_text())
        if value["dependency_sha256"] != hashes or value["point"] != list(point):
            raise ArithmeticError("an interpolation checkpoint has stale dependencies")
        return point, value["result"], 0
    start = perf_counter()
    result = calculate_sources(sources_at(*point))
    atomic_json(path, {"point": point, "dependency_sha256": hashes, "result": result})
    return point, result, perf_counter() - start


def assemble(samples):
    complete = {}
    for point, value in samples.items():
        complete[point] = value
        if point[0] != point[1]:
            complete[point[::-1]] = swap_frame_record(value)
    records = [complete[point] for point in POINTS]
    reconstruction = canonical(interpolate(records))
    for point, actual in zip(POINTS, records):
        expected = evaluate(reconstruction, sources_at(*point))
        for key in ("bare", "all_original_first_current_spinors"):
            if actual[key] != expected[key]:
                raise ArithmeticError("the complete polynomial missed an interpolation value")
    return {"schema_version": 1, "complete": True,
        "scope": "Complete bare kernel and every original first-current spinor for arbitrary balanced rational sources at vertices1,2,5,6; sources0,3,4 are zero. H=2,T=0; geometric kernel normalization.",
        "vertices": 7, "edges": EDGES, "marked": [5, 6],
        "source_support": SOURCE_SUPPORT, "structurally_zero_source_vertices": ZERO_SOURCE_VERTICES,
        "source_degree_bare": 10, "source_degree_first_current": 11,
        "frame_total_degree_bound": 6, "interpolation_monomials": BASIS,
        "interpolation_points": POINTS, "direct_points": DIRECT_POINTS,
        "complete_original_spinors": 11, "complete_original_scalar_entries": 22,
        "samples": records, "reconstruction": reconstruction,
        "independent_controls": {}, "all_values_finite": True, "residual_periods": [],
        "dependency_sha256": dependency_hashes(), "options": OPTIONS,
        "verification_sha256": verification_hashes(),
        "source_basis_meaning": ["lambda_1", "lambda_2", "lambda_5"],
        "bracket_monomial_meaning": "<lambda_1,lambda_2>^(5-i-j) <lambda_5,lambda_2>^i <lambda_1,lambda_5>^j",
        "all_interpolation_entries_verified": len(POINTS) * 23,
        "complete_reconstruction_sha256": digest(reconstruction)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--checkpoint", type=Path, default=HERE / "checkpoints")
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    samples = {}
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(point_job, point, args.checkpoint) for point in DIRECT_POINTS]
        for future in as_completed(futures):
            point, result, elapsed = future.result()
            samples[point] = result
            print("Completed full source sample %s: %d/%d, %.2f seconds." %
                  (point, len(samples), len(DIRECT_POINTS), elapsed), flush=True)
    certificate = canonical(assemble(samples))
    path = HERE / "certificate.json"
    if args.write:
        atomic_json(path, certificate)
    else:
        stored = json.loads(path.read_text())
        certificate["independent_controls"] = stored["independent_controls"]
        if certificate != stored:
            raise ArithmeticError("the complete supported-source certificate changed")
    print("All 28 nodes, 11 complete source-dependent spinors, and all source identities agree.", flush=True)


if __name__ == "__main__":
    main()
