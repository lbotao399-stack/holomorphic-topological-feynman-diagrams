"""Solve complete product-forest requests by direct exact cube expansions.

The JSON graph, source and request conventions agree with
``solve_product_forest.py``.  This additive entry point does not construct
the global Schwinger numerator or consult finite period tables.  Its result
is the complete specified infinite rational series; optional finite terms
and rational enclosures are certificates, never substituted for that result.

Examples, from the repository root:
  python3 code/solve_product_forest_cubes.py input.json
  python3 code/solve_product_forest_cubes.py input.json --series-order 2
  python3 code/solve_product_forest_cubes.py input.json --backend midpoint --width 1/2

``budgets`` contains the selected backend's explicit resource bounds.
``max_tensor_entries`` additionally bounds each complete tensor request.
The geometric backend uses ForestCubeExpansion; midpoint uses ForestCubes.
"""
import argparse
from fractions import Fraction
from itertools import product
import json
import os
from pathlib import Path
import tempfile

from htk import Graph
from htk.forest_cubes import ForestCubes
from solve_product_forest import polynomial, powers


def _integer(value, name, minimum=0):
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError("%s must be an integer at least %d" % (name, minimum))
    return value


def _rational(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, str, Fraction)):
        raise TypeError(name + " must be an exact integer or fraction string")
    return Fraction(value)


def _json_exact(value):
    """Serialize exact quantities without converting any rational to float."""
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _json_exact(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_exact(item) for item in value]
    if value is None or isinstance(value, (str, int, bool)):
        return value
    raise TypeError("an exact series specification contains an unsupported value")


def _row(row):
    result = {"order": row.order, "term": str(row.term),
              "partial_sum": str(row.partial_sum),
              "absolute_tail_bound": str(row.error_bound),
              "lower": str(row.lower), "upper": str(row.upper)}
    for name in ("evaluations", "polynomial_terms", "charts"):
        if hasattr(row, name):
            result[name] = getattr(row, name)
    return _json_exact(result)


def _current_polynomial(graph, request):
    if "polynomial" in request:
        return polynomial(graph, request["polynomial"])
    word = powers(request["powers"])
    exponent = [0] * (2 * graph.ne)
    for (edge, entry), power in word.items():
        if not 0 <= edge < graph.ne or entry not in (0, 1):
            raise ValueError("a current-power position lies outside the original graph")
        exponent[2 * edge + entry] = power
    return {tuple(exponent): Fraction(1)}


def _encode_polynomial(value):
    return [{"exponents": exponent, "coefficient": str(coefficient)}
            for exponent, coefficient in sorted(value.items())]


def solve(document, *, backend="geometric", series_order=None, width=None):
    """Return only complete exact series specifications for all requests.

    Every original source and current entry remains in the output.  Both
    optional evaluations may be requested together.  A failed tensor entry
    or exhausted resource budget aborts the whole document.
    """
    if backend not in ("geometric", "midpoint"):
        raise ValueError("backend must be geometric or midpoint")
    if series_order is not None:
        series_order = _integer(series_order, "series order")
    if width is not None:
        width = _rational(width, "requested interval width")
        if width <= 0:
            raise ValueError("requested interval width must be positive")
    vertices = _integer(document["vertices"], "vertex count", 2)
    edges = document["edges"]
    if not isinstance(edges, list) or any(not isinstance(edge, list) or len(edge) != 2 for edge in edges):
        raise ValueError("edges must contain ordered pairs of original vertex labels")
    for edge in edges:
        for vertex in edge:
            _integer(vertex, "original vertex label")
            if vertex >= vertices:
                raise ValueError("an original edge vertex lies outside the graph")
    graph = Graph(vertices, edges, document.get("name"))
    sources = document["sources"]
    if not isinstance(sources, list) or any(not isinstance(row, list) or len(row) != 2 for row in sources):
        raise ValueError("sources must contain complete two-entry spinor rows")
    sources = tuple(tuple(_rational(value, "source") for value in row) for row in sources)
    marked = document["marked"]
    if not isinstance(marked, list) or len(marked) != 2:
        raise ValueError("marked must contain the ordered vertices A,B")
    marked = tuple(_integer(vertex, "marked vertex") for vertex in marked)
    requests = document.get("requests", [{"name": "bare", "powers": []}])
    if not isinstance(requests, list):
        raise TypeError("requests must be a list")
    budgets = dict(document.get("budgets", {}))
    tensor_budget = _integer(budgets.pop("max_tensor_entries", 65536), "max_tensor_entries", 1)
    if backend == "geometric":
        from htk.forest_cube_expansion import ForestCubeExpansion
        engine = ForestCubeExpansion(graph, sources, marked, **budgets)
        atlas = engine.atlas
    else:
        engine = atlas = ForestCubes(graph, sources, marked, **budgets)

    def exact_series(current_polynomial):
        series = engine.series(polynomial=current_polynomial)
        result = {"type": "complete_convergent_rational_series", "backend": backend,
                  "original_current_polynomial": _encode_polynomial(current_polynomial),
                  "full_original_graph_and_sources": "enclosing_document",
                  "all_original_current_terms_retained": True,
                  "finite_certificate_is_not_the_exact_value": True}
        if backend == "geometric":
            result["series_specification"] = _json_exact(series.specification())
            result["term_rule"] = "htk.forest_cube_expansion.ExactCubeExpansion.terms"
        else:
            result["series_specification"] = {
                "dimension": atlas.dimension, "complete_cubes": atlas.chart_count,
                "spanning_trees": atlas.tree_count,
                "level_zero": "sum of the full oriented current-form midpoint values on every cube",
                "later_terms": "difference of complete dyadic midpoint sums at consecutive levels",
                "subdivisions_per_direction_at_level_k": "2**k",
                "remainder_rule": "rational interval arithmetic on every closed dyadic cell"}
            result["term_rule"] = "htk.forest_cubes.ExactCubeSeries.terms"
        if series_order is not None:
            result["certified_finite_terms"] = [_row(row) for row in series.terms(series_order)]
        if width is not None:
            result["requested_interval_width"] = str(width)
            result["certified_enclosure"] = _row(series.enclose(width))
        return result

    results = []
    for number, request in enumerate(requests):
        if not isinstance(request, dict):
            raise TypeError("each request must be an object")
        kinds = [key for key in ("powers", "polynomial", "tensor_edges") if key in request]
        if len(kinds) != 1:
            raise ValueError("each request must contain exactly one of powers, polynomial, tensor_edges")
        if kinds[0] == "tensor_edges":
            slots = tuple(_integer(edge, "tensor edge") for edge in request["tensor_edges"])
            if any(edge >= graph.ne for edge in slots):
                raise ValueError("a tensor edge lies outside the original graph")
            if len(slots) >= tensor_budget.bit_length():
                raise ValueError("the complete tensor exceeds max_tensor_entries")
            indices = tuple(request.get("indices", ["a%d" % (j + 1) for j in range(len(slots))]))
            if (len(indices) != len(slots) or len(set(indices)) != len(indices)
                    or any(not isinstance(index, str) or not index for index in indices)):
                raise ValueError("provide one distinct free-index name per tensor slot")
            entries = []
            for components in product((0, 1), repeat=len(slots)):
                exponent = [0] * (2 * graph.ne)
                for edge, component in zip(slots, components):
                    exponent[2 * edge + component] += 1
                value = exact_series({tuple(exponent): Fraction(1)})
                entries.append({"storage_indices": components, "value": value})
            result = {"type": "complete_spinor_tensor", "original_edges": slots,
                      "free_lower_spinor_indices": indices, "entries": entries}
        else:
            result = exact_series(_current_polynomial(graph, request))
        results.append({"name": request.get("name", str(number)), "result": result})
    return _json_exact({"schema_version": 1,
            "graph": {"vertices": graph.nv, "edges": graph.edges, "marked": marked},
            "full_balanced_sources": atlas.sources,
            "convention": "htk H=2,T=0 geometric kernels in the original edge order, before physical normalization and displacement factorials",
            "backend": backend, "budgets": dict(budgets, max_tensor_entries=tensor_budget),
            "all_requests_complete": True, "results": results})


def write_complete(path, result):
    """Atomically replace an output only after every request has succeeded."""
    path = Path(path)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=path.parent,
                                         prefix="." + path.name + ".", suffix=".tmp", delete=False) as output:
            temporary = Path(output.name)
            json.dump(result, output, indent=2, sort_keys=True)
            output.write("\n")
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary, path)
        temporary = None
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--backend", choices=("geometric", "midpoint"), default="geometric")
    parser.add_argument("--series-order", type=int)
    parser.add_argument("--width", help="positive exact rational maximum enclosure width, such as 1/100")
    arguments = parser.parse_args()
    result = solve(json.loads(arguments.input.read_text()), backend=arguments.backend,
                   series_order=arguments.series_order, width=arguments.width)
    if arguments.output:
        write_complete(arguments.output, result)
    else:
        print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
