"""Check the universal sourceful rule on every marked forest through six loops.

Every coefficient of the source-linear response is compared with the
original grounded-Laplacian inverse.  Independent cycle-form determinants
check the full current-area density.  No source or tensor entry is selected
to stand in for the complete response.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from time import perf_counter

from catalog import ROOT, count_table, forest_graph, forests
from check_catalog import check as check_catalog
from htk import brute
from htk.forest_kron import ForestKronPlan

HERE = Path(__file__).resolve().parent


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def full_sources(graph, salt):
    values = [(Fraction(vertex + 1 + salt), Fraction((-1) ** vertex * (vertex * vertex + 2 + salt)))
              for vertex in range(graph.nv - 1)]
    values.append(tuple(-sum(row[a] for row in values) for a in range(2)))
    return tuple(values)


def calculate(last=6):
    start = perf_counter()
    records, checked_coefficients, non_laman, maximum_boundary = [], 0, 0, 0
    for loops in range(1, last + 1):
        for index, forest in enumerate(forests(loops)):
            graph, marked = forest_graph(forest)
            plan = ForestKronPlan(graph, marked)
            non_laman += sum(step.remaining_edges != 2 * step.remaining_vertices - 3 for step in plan.steps)
            maximum_boundary = max(maximum_boundary, *(len(step.neighbors) for step in plan.steps))
            points = []
            for point in range(2):
                times = tuple(Fraction(2 + edge + (index % 7), 1 + (edge + point) % 3)
                              for edge in range(graph.ne))
                response = plan.response(times)
                reference = tuple(tuple(row) for row in brute.response(graph, times))
                if response.current_response != reference:
                    raise ArithmeticError("the complete source-linear current response differs")
                checked_coefficients += graph.ne * (graph.nv - 1)
                sources = full_sources(graph, point)
                currents = response.currents(sources)
                if not brute.check_kirchhoff(graph, times, sources, currents):
                    raise ArithmeticError("an original Kirchhoff identity failed")
                cycles = graph.cycle_matrix()
                if response.U != brute.symanzik_U(graph, times, cycles):
                    raise ArithmeticError("the local pivot product changed the Kirchhoff polynomial")
                density = response.density(sources)
                if density != brute.density_T0_det(graph, times, currents, cycles):
                    raise ArithmeticError("the complete current-area pullback differs from the original cycle form")
                if density != response.density(sources, gauge_edge=index % graph.ne):
                    raise ArithmeticError("the current-area density depends on its gauge edge")
                points.append({"times": list(map(str, times)),
                               "full_balanced_sources": [[str(v) for v in row] for row in sources],
                               "full_current_response_sha256": digest([[str(v) for v in row] for row in response.current_response]),
                               "kirchhoff_polynomial": str(response.U), "density": str(density)})
            records.append({"loops": loops, "index": index, "edges": graph.edges, "marked": marked,
                            "elimination_order": [step.vertex for step in plan.steps], "points": points})
        print("All %d topologies at %d loops passed full source-response and density checks (%.3f s)." %
              (len(forests(loops)), loops, perf_counter() - start), flush=True)
    paths = [Path(__file__), HERE / "catalog.py", HERE / "check_catalog.py",
             ROOT / "code/htk/forest_kron.py", ROOT / "code/htk/forest_model.py",
             ROOT / "code/htk/brute.py", ROOT / "code/htk/graph.py", ROOT / "code/htk/linalg.py"]
    return {"schema_version": 1,
            "scope": "Every marked product-forest topology through the listed loop number, all source-response coefficients, two positive rational resistance assignments, full density at two gauges.",
            "last_loop": last, "topologies": len(records), "resistance_assignments": 2 * len(records),
            "source_response_coefficients_compared": checked_coefficients,
            "non_laman_intermediate_stages": non_laman,
            "maximum_retained_boundary_vertices": maximum_boundary,
            "counts": count_table(last), "independent_catalog": check_catalog(),
            "records": records,
            "implementation_sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                                      for path in paths}}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--last", type=int, default=6)
    arguments = parser.parse_args()
    result = json.loads(json.dumps(calculate(arguments.last)))
    path = HERE / "kron_certificate.json"
    if arguments.write:
        path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    elif result != json.loads(path.read_text()):
        raise ArithmeticError("a complete census certificate field changed")
    print("Verified %d full source-response coefficients on %d forests." %
          (result["source_response_coefficients_compared"], result["topologies"]), flush=True)
