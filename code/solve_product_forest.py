"""Solve complete marked-forest current requests from an exact JSON input.

Example from the repository root:
  python3 code/solve_product_forest.py explore/circuit-general/all-product-forests/five_loop.json

Rational data are integers or fraction strings.  A tensor request retains
all of its free spinor indices and every stored tensor entry.  A missing
finite period is serialized as its complete convergent rational series,
including every numerator coefficient, never as an unevaluated integral.
"""
import argparse
from fractions import Fraction
import json
from pathlib import Path

from htk import Graph
from htk.forest_periods import ExactForestSeries, ExactPeriod
from htk.forest_reduction import ForestReduction


def powers(rows):
    result = {}
    for row in rows:
        if not isinstance(row, list) or len(row) != 3:
            raise ValueError("each current-power row must be [original edge, stored spinor entry, power]")
        edge, entry, power = row
        if any(isinstance(value, bool) or not isinstance(value, int) for value in row):
            raise ValueError("edge, stored spinor entry and power must be integers")
        if power < 0:
            raise ValueError("current powers must be nonnegative")
        result[edge, entry] = result.get((edge, entry), 0) + power
    return result


def polynomial(graph, rows):
    result = {}
    for term in rows:
        coefficient = term["coefficient"]
        if isinstance(coefficient, (float, bool)):
            raise ValueError("coefficients must be exact integers or fraction strings")
        coefficient = Fraction(coefficient)
        exponent = [0] * (2 * graph.ne)
        for (edge, entry), power in powers(term.get("powers", [])).items():
            if not 0 <= edge < graph.ne or entry not in (0, 1):
                raise ValueError("a current-power position lies outside the complete graph")
            exponent[2 * edge + entry] += power
        key = tuple(exponent)
        result[key] = result.get(key, Fraction()) + coefficient
    return {exponent: coefficient for exponent, coefficient in result.items() if coefficient}


def exact_value(value, order=None):
    if isinstance(value, (int, Fraction)):
        return {"type": "rational", "value": str(value)}
    if isinstance(value, ExactPeriod):
        return {"type": "finite_period", "coefficients":
                {label: str(coefficient) for label, coefficient in sorted(value.items())}}
    if not isinstance(value, ExactForestSeries):
        raise TypeError("unexpected exact result type")
    result = {"type": "complete_convergent_rational_series", "core_edges": value.edges,
              "rational_constant": str(value.constant),
              "components": [{"denominator_power": q,
                              "numerator": [{"exponents": exponent, "coefficient": str(coefficient)}
                                            for exponent, coefficient in sorted(numerator.items())]}
                             for q, numerator in value.components],
              "term_rule": "htk.forest_periods.ExactForestSeries.terms; default method gegenbauer",
              "finite_period_unavailable_reason": value.reason,
              "max_edges": value.max_edges, "max_terms": value.max_terms}
    if order is not None:
        rows = list(value.terms(order))
        result["certified_finite_terms"] = [{"order": row.order, "term": str(row.term),
                                            "partial_sum": str(row.partial_sum),
                                            "absolute_tail_bound": str(row.error_bound)} for row in rows]
    return result


def solve(document, *, finite_only=False, series_order=None):
    if isinstance(document["vertices"], bool) or not isinstance(document["vertices"], int):
        raise TypeError("the vertex count must be an integer")
    graph = Graph(document["vertices"], document["edges"], document.get("name"))
    if any(isinstance(vertex, bool) or not isinstance(vertex, int)
           for edge in document["edges"] for vertex in edge):
        raise TypeError("original vertex labels must be integers")
    requests = document.get("requests", [{"name": "bare", "powers": []}])
    if not isinstance(requests, list):
        raise TypeError("requests must be a list")
    options = dict(document.get("budgets", {}))
    if "allow_series" in options:
        raise ValueError("use --finite-only to set finite-period behavior")
    results = []
    with ForestReduction(graph, document["sources"], tuple(document["marked"]),
                         allow_series=not finite_only, **options) as reduction:
        for index, request in enumerate(requests):
            kinds = [key for key in ("powers", "polynomial", "tensor_edges") if key in request]
            if len(kinds) != 1:
                raise ValueError("each request must contain exactly one of powers, polynomial, tensor_edges")
            if kinds[0] == "tensor_edges":
                tensor = reduction.tensor(request["tensor_edges"], indices=request.get("indices"))
                value = {"type": "complete_spinor_tensor", "original_edges": tensor.edge_slots,
                         "free_lower_spinor_indices": tensor.indices,
                         "entries": [{"storage_indices": entry, "value": exact_value(coefficient, series_order)}
                                     for entry, coefficient in sorted(tensor.entries.items())]}
            elif kinds[0] == "polynomial":
                value = exact_value(reduction.integrate(polynomial(graph, request["polynomial"])), series_order)
            else:
                value = exact_value(reduction.moment(powers(request["powers"])), series_order)
            results.append({"name": request.get("name", str(index)), "result": value})
        return {"schema_version": 1, "graph": {"vertices": graph.nv, "edges": graph.edges,
                                                "marked": reduction.marked},
                "full_balanced_sources": [[str(value) for value in row] for row in reduction.sources],
                "convention": "htk H=2,T=0 geometric kernels with original edge order, before physical normalization and displacement factorials",
                "all_requests_complete": True, "results": results}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--finite-only", action="store_true")
    parser.add_argument("--series-order", type=int)
    arguments = parser.parse_args()
    if arguments.series_order is not None and arguments.series_order < 0:
        parser.error("--series-order must be nonnegative")
    result = solve(json.loads(arguments.input.read_text()), finite_only=arguments.finite_only,
                   series_order=arguments.series_order)
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.output:
        # Only complete successful results are written; budget failures leave
        # an existing output untouched.
        arguments.output.write_text(encoded)
    else:
        print(encoded, end="")


if __name__ == "__main__":
    main()
