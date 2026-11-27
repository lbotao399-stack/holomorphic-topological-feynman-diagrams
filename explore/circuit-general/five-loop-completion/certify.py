"""Compute complete supported-source bulk3 five-loop kernels with bounded jobs."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import signal
import time

import networkx as nx

from bulk3_periods import Bulk3Periods, HERE, ROOT, encoded, polynomial_record

DATA = ROOT / "explore/circuit-general/fiveloop/enum_d3.json"
ROWS = json.loads(DATA.read_text())["bulk3"]
DEPENDENCIES = (
    "explore/circuit-general/five-loop-completion/bulk3_periods.py",
    "explore/circuit-general/five-loop-completion/certify.py",
    "explore/circuit-general/fiveloop/enum_d3.json",
    "explore/circuit-general/fiveloop/skel5.json",
    "explore/circuit-general/fiveloop/per5.py",
    "explore/circuit-general/fiveloop/msys5.py",
    "explore/circuit-general/fiveloop/rel5.py",
    "explore/circuit-general/cores/gutil.py",
    "explore/circuit-general/cores/wheel.py",
    "explore/circuit-general/assembly/tables/prism_q8.json",
    "explore/circuit-general/assembly/tables/K33_q8.json",
    "explore/circuit-general/assembly/tables/K33_q9.json",
    "explore/nonSP/verify/k4_periods/current_minor.py",
    "code/htk/forest_periods.py", "code/htk/_forest_period_worker.py",
    "code/htk/series.py", "code/htk/treeflow.py", "code/htk/graph.py",
    "code/htk/brute.py", "code/htk/linalg.py", "code/htk/CONVENTIONS.md",
)


def hashes():
    return {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in DEPENDENCIES}


def sign(permutation):
    return (-1) ** sum(left > right for i, left in enumerate(permutation)
                      for right in permutation[i + 1:])


def zero_witness(row):
    graph = nx.Graph(row["edges"])
    vertices, sources = list(range(row["nv"])), row["S0"]
    for mapping in nx.algorithms.isomorphism.GraphMatcher(graph, graph).isomorphisms_iter():
        vertex_sign = sign([mapping[v] for v in vertices])
        reversals = [e for e, (left, right) in enumerate(row["edges"]) if mapping[left] > mapping[right]]
        source_permutation = [sources.index(mapping[v]) for v in sources]
        source_sign = sign(source_permutation)
        if vertex_sign * (-1) ** len(reversals) * source_sign == -1:
            return {"vertex_permutation": [mapping[v] for v in vertices],
                    "vertex_sign": vertex_sign, "reversed_original_edges": reversals,
                    "source_permutation": source_permutation, "source_determinant_sign": source_sign,
                    "complete_bare_multiplier": -1}
    return None


def atomic_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def calculate(index, currents, seconds):
    row = ROWS[index]
    start = time.monotonic()
    engine = None
    def timeout(*unused):
        raise TimeoutError("Complete graph exceeded %d seconds." % seconds)
    signal.signal(signal.SIGALRM, timeout)
    signal.alarm(seconds)
    result = {"schema": "five-loop-bulk3-supported-source-v1", "index": index,
              "ordered_oriented_edges": row["edges"], "vertices": row["nv"],
              "source_vertices": row["S0"], "core_classification": row["core"],
              "source_scope": "lambda_s0,lambda_s1,lambda_s2 arbitrary; lambda_s3=-sum; all other original sources zero",
              "bare_source_factor": "det(lambda_s0,lambda_s1,lambda_s2)^5",
              "first_current_meaning": "bare_source_factor times the complete linear combination of lambda_s0,lambda_s1,lambda_s2",
              "bare_status": "uncomputed", "first_current_status": "not_requested",
              "dependencies_sha256": hashes()}
    try:
        engine = Bulk3Periods(row)
        result.update(source_cycle_tree=list(engine.tree),
                      complete_bare_numerator=polynomial_record(engine._bare),
                      original_grounded_laplacian_density_checks=engine.density_checks,
                      complete_series_chains=[list(chain.edges) for chain in engine._chains],
                      scalar_core_edges=[list(edge) for edge in engine._core])
        bare = engine.integrate()
        result.update(bare_status="complete", bare_coefficient=encoded(bare))
        witness = zero_witness(row)
        if witness:
            if bare:
                raise ArithmeticError("The complete exact bare value violates its signed automorphism.")
            result["independent_zero_automorphism"] = witness
        if currents:
            result["first_current_status"] = "computing"
            def progress(edge):
                print("bulk3 %d: full original chord %d done; %.1f s" %
                      (index, edge, time.monotonic() - start), flush=True)
            values = engine.complete_first_currents(bare, progress=progress)
            result.update(first_current_status="complete",
                          all_original_first_current_coefficients=[[encoded(value) for value in item] for item in values],
                          integrated_original_source_identities=3 * row["nv"],
                          complete_affine_current_polynomial_identities=3 * len(row["edges"]))
    except Exception as error:
        result["error"] = type(error).__name__ + ": " + str(error)
        if result["first_current_status"] == "computing":
            result["first_current_status"] = "failed"
    finally:
        signal.alarm(0)
        if engine is not None:
            result["complete_finite_targets"] = engine.integrations
            engine.close()
        result["seconds"] = round(time.monotonic() - start, 3)
    return result


def selected(text):
    out = set()
    for item in text.split(","):
        if "-" in item:
            first, last = map(int, item.split("-"))
            out.update(range(first, last + 1))
        else:
            out.add(int(item))
    if not out or min(out) < 0 or max(out) >= len(ROWS):
        raise ValueError("Graph index outside the complete bulk3 census.")
    return sorted(out)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--indices", default="0-68")
    parser.add_argument("--currents", choices=("none", "lower", "all"), default="lower",
                        help="lower requests all first spinors for SP,K4,K33; all also tries the other selected cores")
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--seconds", type=int, default=600)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    if not 1 <= args.jobs <= 3 or args.seconds < 1:
        parser.error("Use one to three jobs and a positive graph budget.")
    tasks = []
    for index in selected(args.indices):
        currents = args.currents == "all" or (args.currents == "lower" and ROWS[index]["core"] in ("SP", "K4", "K33"))
        path = HERE / "cases" / ("bulk3_%03d.json" % index)
        if args.resume and path.exists():
            data = json.loads(path.read_text())
            if data["dependencies_sha256"] != hashes():
                raise ArithmeticError("A saved case has changed exact dependencies.")
            if data["bare_status"] == "complete" and (not currents or data["first_current_status"] == "complete"):
                print("Reused complete bulk3", index, flush=True)
                continue
        tasks.append((index, currents))
    failed = False
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        pending = {pool.submit(calculate, index, currents, args.seconds): (index, currents)
                   for index, currents in tasks}
        for future in as_completed(pending):
            index, currents = pending[future]
            result = future.result()
            atomic_json(HERE / "cases" / ("bulk3_%03d.json" % index), result)
            success = result["bare_status"] == "complete" and (not currents or result["first_current_status"] == "complete")
            failed |= not success
            print("bulk3 %d: bare %s, all first currents %s; %.3f s%s" %
                  (index, result["bare_status"], result["first_current_status"], result["seconds"],
                   "; " + result["error"] if "error" in result else ""), flush=True)
    raise SystemExit(1 if failed else 0)


if __name__ == "__main__":
    main()
