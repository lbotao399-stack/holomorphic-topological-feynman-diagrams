"""Independent exact dispatcher, orientation and nonconstant-series controls."""
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "code"))
from htk.graph import Graph, theta
from htk.bulk_kernels import evaluate_bulk
from htk.bulk_series import ExactBulkSeries


def encoded(value):
    return {key: str(coefficient) for key, coefficient in value.coefficients.items()}


def pi_interval():
    def atan(inverse):
        terms = 30
        partial = sum((F((-1) ** k, (2 * k + 1) * inverse ** (2 * k + 1))
                       for k in range(terms)), F(0))
        adjacent = partial + F((-1) ** terms, (2 * terms + 1) * inverse ** (2 * terms + 1))
        return min(partial, adjacent), max(partial, adjacent)
    first, second = atan(5), atan(239)
    return 16 * first[0] - 4 * second[1], 16 * first[1] - 4 * second[0]


def finite_interval(answer):
    lo, hi = pi_interval()
    lower = upper = F(0)
    for label, coefficient in answer.coefficients.items():
        if label == "1":
            bounds = (F(1), F(1))
        elif label.startswith("pi^"):
            power = int(label[3:])
            assert power >= 0
            bounds = (lo ** power, hi ** power)
        else:
            raise AssertionError("this control expects a complete polynomial in pi")
        products = [coefficient * bound for bound in bounds]
        lower += min(products)
        upper += max(products)
    return lower, upper


def main():
    start = time.monotonic()
    original = ((0, 1), (0, 2), (1, 2), (3, 4), (3, 5), (4, 5), (0, 3), (1, 4), (2, 5))
    sources = ((2,), (-3,), (5,), (7,), (-13,), (2,))
    graph = Graph(6, original)
    permutation = (5, 2, 0, 4, 1, 3)
    order = (8, 2, 6, 4, 0, 7, 3, 1, 5)
    reversed_edges = {2, 6}
    transformed = []
    for edge in order:
        a, b = original[edge]
        if edge in reversed_edges:
            a, b = b, a
        transformed.append((permutation[a], permutation[b]))
    changed_sources = [None] * 6
    for old, new in enumerate(permutation):
        changed_sources[new] = sources[old]
    changed_graph = Graph(6, transformed)
    sign = (-1) ** sum(a > b for i, a in enumerate(order) for b in order[i + 1:])
    records = []
    for edge in (None,) + tuple(range(9)):
        word = () if edge is None else ((edge, 0),)
        changed_word = () if edge is None else ((order.index(edge), 0),)
        before = evaluate_bulk(graph, sources, H=1, T=1, word=word)
        after = evaluate_bulk(changed_graph, changed_sources, H=1, T=1, word=changed_word)
        factor = sign * (-1 if edge in reversed_edges else 1)
        assert dict(after.coefficients) == {key: factor * value for key, value in before.coefficients.items()}
        records.append({"original_word": word, "transformed_word": changed_word,
                        "complete_original_value": encoded(before), "complete_transformed_value": encoded(after),
                        "exact_orientation_factor": factor})
    forced = evaluate_bulk(changed_graph, changed_sources, H=1, T=1,
                           word=((order.index(2), 0),), route="periods", period_timeout=120)
    assert encoded(forced) == records[3]["complete_transformed_value"]
    print("all ten prism outputs retain source grounding, vertex labels, edge order and orientation", flush=True)

    graph3, _ = theta(3, 3, 3)
    sources3 = [[(v + 1) ** (a + 1) - (a + 2) for a in range(3)] for v in range(graph3.nv - 1)]
    word3 = ((0, 2), (8, 0))
    finite3 = evaluate_bulk(graph3, sources3, H=3, T=1, word=word3)
    periods3 = evaluate_bulk(graph3, sources3, H=3, T=1, word=word3,
                             route="periods", period_timeout=120)
    assert finite3.as_sympy() == periods3.as_sympy()
    assert finite3.coefficients
    print("three-component two-loop mixed current agrees by simplex and scalar periods", flush=True)

    graph2, _ = theta(2, 2, 1)
    sources2 = [[1], [2], [-1]]
    word2 = ((0, 0), (3, 0))
    known = evaluate_bulk(graph2, sources2, H=1, T=1, word=word2)
    series = evaluate_bulk(graph2, sources2, H=1, T=1, word=word2, route="series")
    assert isinstance(series, ExactBulkSeries)
    order2 = series.required_order(F(1, 1000))
    lower, upper = series.interval(order2)
    known_lower, known_upper = finite_interval(known)
    assert lower <= known_lower and known_upper <= upper
    assert upper - lower <= F(1, 1000)
    print("complete nonconstant mixed-current series encloses the exact finite answer", flush=True)
    paths = [ROOT / "code/htk" / name for name in
             ("bulk_kernels.py", "bulk_series.py", "topological_polynomial.py", "_bulk_period_worker.py")]
    paths += [Path(__file__)]
    result = {"schema": "independent-bulk-orientation-and-series-v1", "result": "passed",
              "all_ten_prism_relabeling_controls": records,
              "original_prism_graph": original, "complete_original_prism_sources": sources,
              "vertex_permutation": permutation, "edge_order": order, "reversed_original_edges": sorted(reversed_edges),
              "transformed_graph": transformed, "complete_transformed_sources": changed_sources,
              "forced_scalar_period_transformed_value": encoded(forced),
              "H3_T1_graph": graph3.edges, "H3_T1_sources": sources3, "H3_T1_word": word3,
              "H3_T1_simplex_value": encoded(finite3), "H3_T1_scalar_period_value": encoded(periods3),
              "mixed_current_series": {"original_graph": graph2.edges, "sources": sources2, "word": word2,
                                       "complete_finite_value": encoded(known), "order": order2,
                                       "interval": [str(lower), str(upper)],
                                       "known_interval_from_exact_Machin_identity": [str(known_lower), str(known_upper)],
                                       "requested_width": "1/1000", "complete_atlas_statistics": series.statistics},
              "implementation_sha256": {str(path.relative_to(ROOT)): sha256(path.read_bytes()).hexdigest() for path in paths},
              "seconds": round(time.monotonic() - start, 2)}
    (HERE / "bulk_independent_controls.json").write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")


if __name__ == "__main__":
    main()
