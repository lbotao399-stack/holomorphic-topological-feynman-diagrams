"""Complete bare and first-current kernels on every forest through four loops.

The certificate fixes one admissible full source assignment on each graph.
Every original-edge first-current spinor is evaluated.  The new compositional
solver is checked against direct original-graph integration, with no passive
or sourced-ear precontraction in the reference.
"""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from time import perf_counter

from catalog import ROOT, forest_graph, forests
from htk.forest_model import recognize_product_forest
from htk.forest_periods import CoreCurrentPeriods, ExactPeriod

HERE = Path(__file__).resolve().parent


def source_assignment(graph, marked):
    structure = recognize_product_forest(graph, marked)
    support = {vertex for roots in structure.degree_two_vertices for vertex in roots}
    superpotential = []
    for component, higher in zip(structure.components, structure.higher_valence_vertices):
        # A higher superpotential at this vertex may leave an independent
        # output field even when its contracted graph degree is three.
        vertex = (higher or tuple(v for v in component if graph.degree(v) == 3) or component)[0]
        support.add(vertex)
        superpotential.append({"vertex": vertex, "degree": graph.degree(vertex) + 1})
    rows = [((vertex + 1), (-1) ** vertex * (vertex + 2)) if vertex in support else (0, 0)
            for vertex in range(graph.nv - 2)]
    rows.append((2 * graph.loops + 3, graph.loops + 3))
    rows.append(tuple(-sum(row[a] for row in rows) for a in range(2)))
    return tuple(rows), superpotential


def finite_record(value):
    if isinstance(value, (int, Fraction)):
        return {"1": str(value)} if value else {}
    if not isinstance(value, ExactPeriod):
        raise ArithmeticError("an enumerated four-loop kernel did not close as a finite exact expression")
    return {label: str(coefficient) for label, coefficient in sorted(value.items())}


def add_records(records):
    result = {}
    for weight, value in records:
        for label, coefficient in value.items():
            result[label] = result.get(label, Fraction()) + weight * Fraction(coefficient)
    return {label: str(coefficient) for label, coefficient in sorted(result.items()) if coefficient}


def one(job):
    loops, index, backend = job
    graph, marked = forest_graph(forests(loops)[index])
    sources, superpotential = source_assignment(graph, marked)
    start = perf_counter()
    if backend == "original":
        engine = CoreCurrentPeriods(graph, sources, marked=marked, check=0, max_terms=1000000)
        def moment(edge=None, entry=None):
            exponent = [0] * (2 * graph.ne)
            if edge is not None:
                exponent[2 * edge + entry] = 1
            return engine.integrate({tuple(exponent): Fraction(1)}, allow_series=False)
    else:
        from htk.forest_reduction import ForestReduction
        engine = ForestReduction(graph, sources, marked, check=1, max_terms=1000000,
                                 allow_series=False)
        def moment(edge=None, entry=None):
            return engine.moment(None if edge is None else {(edge, entry): 1})
    try:
        bare = finite_record(moment())
        currents = [[finite_record(moment(edge, entry)) for entry in range(2)] for edge in range(graph.ne)]
        # Full integrated Kirchhoff identities, separately in each stored
        # finite period channel.  No independence of zeta values is assumed.
        for vertex in range(graph.nv):
            for entry in range(2):
                actual = add_records((int(first == vertex) - int(second == vertex), currents[edge][entry])
                                     for edge, (first, second) in enumerate(graph.edges))
                expected = add_records(((sources[vertex][entry], bare),))
                if actual != expected:
                    raise ArithmeticError("a complete integrated source-conservation identity failed")
        return {"loops": loops, "index": index, "edges": graph.edges, "marked": marked,
                "full_balanced_sources": sources, "superpotential_vertices": superpotential,
                "bare": bare, "all_original_first_current_spinors": currents}, perf_counter() - start
    finally:
        engine.close()


def calculate(backend, jobs=3, last=4):
    tasks = [(loops, index, backend) for loops in range(1, last + 1) for index in range(len(forests(loops)))]
    records, start = [], perf_counter()
    with ProcessPoolExecutor(max_workers=jobs) as pool:
        for record, elapsed in pool.map(one, tasks):
            records.append(record)
            print("%s: loop %d forest %d, all %d first-current entries in %.3f s (%d/%d)." %
                  (backend, record["loops"], record["index"], 2 * len(record["edges"]), elapsed,
                   len(records), len(tasks)), flush=True)
    print("Complete %s integration finished in %.3f s." % (backend, perf_counter() - start), flush=True)
    return json.loads(json.dumps(records))


def build_certificate(records):
    counts = Counter(record["loops"] for record in records)
    periods = Counter(label for record in records for value in [record["bare"]] +
                      [entry for spinor in record["all_original_first_current_spinors"] for entry in spinor]
                      for label in value if label != "1")
    return {"schema_version": 1,
            "scope": "Complete bare kernel and every original-edge first-current spinor at the explicitly recorded admissible full source assignment on every marked forest through four loops; not an all-source interpolation or a full Wick/color assembly.",
            "convention": "htk H=2,T=0 original ordered/oriented geometric kernels before physical normalization or displacement factorials",
            "topologies": len(records), "counts": dict(sorted(counts.items())),
            "first_current_entries": sum(2 * len(record["edges"]) for record in records),
            "period_channels": dict(sorted(periods.items())), "records": records,
            "all_first_current_entries_retained": True,
            "original_graph_reference": "CoreCurrentPeriods on the complete original graph, every bare and first-current entry, no forest precontractions",
            "exact_integrated_kirchhoff_identities": sum(2 * (record["loops"] + 2) for record in records),
            "records_sha256": hashlib.sha256(json.dumps(records, sort_keys=True, separators=(",", ":")).encode()).hexdigest()}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--backend", choices=("original", "reduction"), default="reduction")
    parser.add_argument("--jobs", type=int, default=3)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--last", type=int, default=4)
    arguments = parser.parse_args()
    result = build_certificate(calculate(arguments.backend, arguments.jobs, arguments.last))
    path = HERE / "integral_certificate.json"
    if arguments.write:
        path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    elif json.loads(json.dumps(result)) != json.loads(path.read_text()):
        raise ArithmeticError("a complete original-edge integral differs from its certificate")
    print("All %d bare kernels and %d complete first-current entries verified." %
          (result["topologies"], result["first_current_entries"]), flush=True)
