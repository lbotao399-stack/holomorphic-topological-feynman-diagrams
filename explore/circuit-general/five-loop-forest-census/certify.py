"""Evaluate all five-loop bare kernels and every original first-current spinor.

Only independent cotree currents need integration: all tree currents are
reconstructed from the exact original Kirchhoff equations, including the
complete affine source flow.  No original current or source is discarded.
The optional direct mode integrates every edge separately as a control.
Per-graph subprocess bounds never turn an unfinished graph into a value.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from time import perf_counter

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PREVIOUS = HERE.parent / "all-product-forests"
sys.path.insert(0, str(PREVIOUS))
from catalog import forest_graph, forests
from certify_integrals import add_records, finite_record, source_assignment
from htk.forest_reduction import ForestReduction
from htk.treeflow import tree_flow


def canonical(value):
    return json.loads(json.dumps(value, sort_keys=True))


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def check_conservation(record):
    edges = record["edges"]
    currents = record["all_original_first_current_spinors"]
    sources = record["full_balanced_sources"]
    if len(currents) != len(edges) or any(len(row) != 2 for row in currents):
        raise ArithmeticError("an original current spinor is incomplete")
    for vertex, source in enumerate(sources):
        for entry in range(2):
            actual = add_records((int(first == vertex) - int(second == vertex), currents[edge][entry])
                                 for edge, (first, second) in enumerate(edges))
            expected = add_records(((Fraction(source[entry]), record["bare"]),))
            if actual != expected:
                raise ArithmeticError("complete integrated source conservation failed")


def reconstruct_currents(graph, sources, bare, integrated):
    """Reconstruct every original current from a complete cotree and sources."""
    tree = tuple(graph.spanning_tree())
    cotree = tuple(graph.cotree(tree))
    if set(integrated) != set(cotree):
        raise ValueError("one complete spinor is required for every original cotree edge")
    cycles = graph.fundamental_cycles(tree)
    flow = tree_flow(graph, sources, tree)
    return [[add_records([(flow[edge][entry], bare)] +
             [(cycles[cycle][edge], integrated[original][entry])
              for cycle, original in enumerate(cotree)])
             for entry in range(2)] for edge in range(graph.ne)]


def solve(index, *, direct=False):
    graph, marked = forest_graph(forests(5)[index])
    sources, superpotential = source_assignment(graph, marked)
    with ForestReduction(graph, sources, marked, check=1, max_terms=4000000,
                         period_timeout=180, allow_series=False) as engine:
        if index in (42, 47, 50, 117, 120, 157):
            from htk.forest_chain_periods import ChainFirstPeriods
            engine._prepare()
            if engine._contractions.steps or engine._ears is not None:
                raise ArithmeticError("an irreducible census core unexpectedly acquired a local reduction")
            engine._core.close()
            engine._core = ChainFirstPeriods(graph, sources, marked, **engine.options)
        bare = finite_record(engine.moment())
        tree = tuple(graph.spanning_tree())
        cotree = tuple(graph.cotree(tree))
        selected = range(graph.ne) if direct else cotree
        integrated = {edge: [finite_record(engine.moment({(edge, entry): 1})) for entry in range(2)]
                      for edge in selected}
        if direct:
            currents = [integrated[edge] for edge in range(graph.ne)]
        else:
            currents = reconstruct_currents(graph, sources, bare, integrated)
        record = {"loops": 5, "index": index, "edges": graph.edges, "marked": marked,
                  "full_balanced_sources": sources, "superpotential_vertices": superpotential,
                  "bare": bare, "all_original_first_current_spinors": currents}
    check_conservation(record)
    return canonical(record)


def graph_job(index, checkpoint, timeout, direct):
    path = checkpoint / ("forest_%03d.json" % index)
    if path.exists():
        record = json.loads(path.read_text())
        check_conservation(record)
        return index, "reused", 0
    command = [sys.executable, str(Path(__file__).resolve()), "--worker", str(index),
               "--output", str(path)]
    if direct:
        command.append("--direct")
    started = perf_counter()
    run = subprocess.run(command, text=True, capture_output=True, timeout=timeout)
    if run.returncode:
        raise RuntimeError("forest %d did not complete: %s" % (index, run.stderr.strip()))
    record = json.loads(path.read_text())
    check_conservation(record)
    return index, "complete", perf_counter() - started


def assemble(checkpoint):
    records = []
    for index in range(len(forests(5))):
        record = json.loads((checkpoint / ("forest_%03d.json" % index)).read_text())
        graph, marked = forest_graph(forests(5)[index])
        sources, superpotential = source_assignment(graph, marked)
        expected = canonical({"loops": 5, "index": index, "edges": graph.edges, "marked": marked,
                              "full_balanced_sources": sources, "superpotential_vertices": superpotential})
        if {key: record[key] for key in expected} != expected:
            raise ArithmeticError("a census graph or complete source assignment changed")
        check_conservation(record)
        records.append(record)
    if len(records) != 200:
        raise ArithmeticError("the complete five-loop topology census changed")
    return {"schema_version": 1,
            "scope": "All 200 marked product-forest topologies at five loops; complete bare and every original-edge first-current spinor at the explicitly recorded admissible full source assignment. No all-source interpolation or Wick/color assembly is asserted.",
            "normalization": "htk H=2,T=0 original ordered and oriented geometric kernel",
            "topologies": len(records), "bare_kernels": len(records),
            "complete_first_current_entries": sum(2 * len(record["edges"]) for record in records),
            "exact_integrated_kirchhoff_identities": sum(2 * len(record["full_balanced_sources"]) for record in records),
            "reconstruction": "Every current equals its complete tree-supported source flow plus the full fundamental-cycle combination of independently integrated original cotree currents. The affine source term is retained.",
            "all_values_finite": True, "residual_periods": [],
            "records_sha256": digest(records), "records": records}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", type=int)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--checkpoint", type=Path)
    parser.add_argument("--indices")
    parser.add_argument("--jobs", type=int, default=3)
    parser.add_argument("--timeout", type=int, default=900)
    parser.add_argument("--direct", action="store_true")
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--assemble-only", action="store_true")
    options = parser.parse_args()
    if options.worker is not None:
        if options.output is None or not 0 <= options.worker < len(forests(5)):
            parser.error("a valid --worker requires --output")
        result = solve(options.worker, direct=options.direct)
        options.output.write_text(json.dumps(result, sort_keys=True, separators=(",", ":")) + "\n")
        return
    if options.checkpoint is None:
        parser.error("--checkpoint must name a work directory")
    if options.jobs < 1 or options.timeout < 1:
        parser.error("jobs and timeout must be positive")
    options.checkpoint.mkdir(parents=True, exist_ok=True)
    indices = list(range(len(forests(5)))) if options.indices is None else [int(x) for x in options.indices.split(",")]
    if any(not 0 <= index < len(forests(5)) for index in indices):
        parser.error("a selected forest index is outside the census")
    if not options.assemble_only:
        with ThreadPoolExecutor(max_workers=options.jobs) as pool:
            pending = {pool.submit(graph_job, index, options.checkpoint, options.timeout, options.direct): index
                       for index in indices}
            errors = []
            for future in as_completed(pending):
                try:
                    index, status, elapsed = future.result()
                    print("forest %d: %s; %.3f seconds" % (index, status, elapsed), flush=True)
                except Exception as error:
                    errors.append((pending[future], str(error)))
                    print("forest %d: incomplete; %s" % (pending[future], error), flush=True)
            if errors:
                raise RuntimeError("unfinished graphs retain no claimed finite values: %s" % errors)
    if options.indices is not None and not options.assemble_only:
        return
    result = assemble(options.checkpoint)
    path = HERE / "certificate.json"
    if options.write:
        path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    elif result != json.loads(path.read_text()):
        raise ArithmeticError("the complete finite census differs from its certificate")
    print("All %d bare kernels, %d first-current entries and %d integrated source identities verified." %
          (result["bare_kernels"], result["complete_first_current_entries"],
           result["exact_integrated_kirchhoff_identities"]), flush=True)


if __name__ == "__main__":
    main()
