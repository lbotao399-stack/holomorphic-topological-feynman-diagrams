"""Replay exact geometric coefficients, tails, full source maps and spinors."""

import argparse
from dataclasses import asdict
from fractions import Fraction as F
import gzip
import hashlib
import json
from pathlib import Path
import sys
from unittest.mock import patch

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "code"))
from htk import brute
from htk.forest_cube_expansion import ForestCubeExpansion, _fraction, _q
from htk.graph import Graph
from solve_product_forest_cubes import solve


THREE = ((0, 1), (0, 3), (0, 4), (1, 2), (1, 4), (2, 3), (3, 4))
OUTPUT = HERE / "geometric_certificate.json.gz"
SUMMARY = HERE / "geometric_certificate.json"


def exact_json(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {str(key): exact_json(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [exact_json(item) for item in value]
    return value


def triangle():
    document = json.loads((HERE / "triangle_input.json").read_text())
    result = solve(document, series_order=6, width="1/1000")
    expected = {"bare": -F(1, 2), "mixed": -F(7, 12)}
    spinors = {"current_01": (-F(1, 6), F(1, 6)),
               "current_02": (-F(1, 3), -F(1, 6)),
               "current_12": (-F(1, 6), -F(1, 3))}
    for record in result["results"]:
        name, value = record["name"], record["result"]
        values = [(value, expected[name])] if name in expected else [
            (entry["value"], exact) for entry, exact in zip(value["entries"], spinors[name])]
        for expression, exact in values:
            for row in expression["certified_finite_terms"]+[expression["certified_enclosure"]]:
                if not F(row["lower"]) <= exact <= F(row["upper"]):
                    raise ArithmeticError("a complete triangle output fails its independent simplex value")
                if abs(F(row["partial_sum"])-exact) > F(row["absolute_tail_bound"]):
                    raise ArithmeticError("a triangle tail fails its independent exact value")
            row = expression["certified_enclosure"]
            if F(row["upper"])-F(row["lower"]) > F(1, 1000):
                raise ArithmeticError("a complete triangle tensor entry fails its requested interval width")
    graph = Graph(3, ((0, 1), (0, 2), (1, 2)))
    series = ForestCubeExpansion(graph, ((1, 0), (0, 1), (-1, -1)), marked=(1, 2)).series()
    all_terms = list(series.terms(20))
    for row in all_terms:
        term = F(0) if row.order % 2 else -F(3, 2**(row.order+3))
        actual_tail = F(1, 2**(2*(row.order//2)+3))
        if row.term != term or abs(row.partial_sum+F(1, 2)) != actual_tail:
            raise ArithmeticError("the universal expansion differs from the all-order triangle solution")
    return {"complete_cli_result": result,
            "all_order_triangle": {"even_term": "-3/2^(2*j+3)", "odd_term": "0",
                                   "exact_infinite_sum": "-1/2",
                                   "actual_tail_after_order": "1/2^(2*floor(order/2)+3)",
                                   "verified_terms_0_through_20": [asdict(row) for row in all_terms]},
            "independent_complete_first_current_spinors": spinors,
            "independent_mixed_simplex_value": -F(7, 12)}


def sourced_three_loop():
    sources = ((2, -1), (0, 0), (-1, 3), (4, 2), (-5, -4))
    engine = ForestCubeExpansion(Graph(5, THREE), sources, marked=(3, 4), check=0)
    # Certificate-only reuse of the exact already checked chart polynomials
    # across the bare value and both entries of the requested free spinor.
    original = engine._chart
    checked, cache = [], {}
    def checked_chart(cube):
        if cube.position not in cache:
            answer = original(cube)
            denominator, currents, numerator, power = answer
            point = tuple(F(j+2, 2*j+5) for j in range(cube.dimension))
            values = tuple(map(_q, point))
            divisor = _fraction(denominator(*values))
            times = cube.time_map(point)
            original_currents = brute.currents(engine.graph, times, sources)
            computed = tuple(_fraction(current(*values))/divisor for current in currents)
            if computed != tuple(value for row in original_currents for value in row):
                raise ArithmeticError("a geometric chart lost an original current or source")
            measure = F(1)
            for value in times:
                measure *= value
            for value in point:
                measure /= value
            expected = brute.density(engine.graph, times, sources).rational()*measure
            actual = _fraction(numerator(*values))/divisor**power
            if actual != expected:
                raise ArithmeticError("the augmented polynomial density differs from the original circuit density")
            checked.append({"cube": cube.position, "original_tree": cube.tree,
                            "original_cotree_order": cube.cotree_order,
                            "orientation": cube.orientation,
                            "point": point, "original_resistances": times,
                            "all_original_current_spinors": original_currents,
                            "normalized_denominator": divisor,
                            "verified_complete_density": actual,
                            "denominator_power_after_exact_cancellation": power})
            cache[cube.position] = answer
        return cache[cube.position]
    engine._chart = checked_chart
    references = (({}, {"1": F(2023, 96)}),
                  ({(1, 0): 1}, {"1": F(5477, 864), "z3": -F(304, 15)}),
                  ({(1, 1): 1}, {"1": -F(8205703, 302400), "z3": -F(6289, 2100)}))
    reference_path = HERE.parent / "forest-jet-identities" / "identity_certificate.json"
    archived = json.loads(reference_path.read_text())
    for powers, expected in references:
        target = next(record for record in archived["direct_period_and_rational_boundary_controls"]
                      if record["leg"] == 1 and record["degree"] == sum(powers.values())
                      and tuple(record["all_symmetric_spinor_powers"]) ==
                      (powers.get((1, 0), 0), powers.get((1, 1), 0)))
        if expected != {key: F(value) for key, value in target["complete_individual_graph_periods"][0].items()}:
            raise ArithmeticError("a complete independent finite-period reference differs")
    partial_zeta = sum((F(1, j**3) for j in range(1, 17)), F(0))
    zeta = (partial_zeta+F(1, 2*17**2), partial_zeta+F(1, 2*16**2))
    records = []
    with patch("htk.forest_periods.CoreCurrentPeriods._source_numerator",
               side_effect=AssertionError("global Schwinger numerator used")), \
         patch("htk.forest_periods._research_module",
               side_effect=AssertionError("finite period or earlier Hepp backend used")):
        for powers, expected in references:
            series = engine.series(powers)
            specification = series.specification()
            rows = list(series.terms(2))
            period_interval = sorted(expected["1"]+expected.get("z3", F(0))*endpoint for endpoint in zeta)
            for row in rows:
                if not row.lower <= period_interval[0] <= period_interval[1] <= row.upper:
                    raise ArithmeticError("a complete geometric remainder fails the independently known period")
            records.append({"original_current_word": [[edge, entry, power]
                                                      for (edge, entry), power in sorted(powers.items())],
                            "complete_finite_reference": expected,
                            "reference_rational_bracket": period_interval,
                            "complete_exact_series_specification": specification,
                            "verified_complete_orders_0_through_2": [asdict(row) for row in rows],
                            "sufficient_order_for_width_1_over_1000": series.required_order(F(1, 1000)),
                            "sufficient_order_coefficients_computed": False,
                            "statistics": series.last_stats})
    if len(checked) != engine.atlas.chart_count:
        raise ArithmeticError("a complete three-loop geometric atlas was not checked")
    return {"ordered_edges": THREE, "marked": (3, 4), "all_balanced_sources": sources,
            "complete_cubes": engine.atlas.chart_count,
            "every_original_chart_control": checked,
            "bare_and_complete_first_current_spinor": records,
            "zeta3_rational_bracket": zeta,
            "reference_certificate_sha256": hashlib.sha256(reference_path.read_bytes()).hexdigest(),
            "scope": "complete infinite rational specifications and first three exact coefficients; no high-precision three-loop enclosure is claimed"}


def certificate():
    triangle_data = triangle()
    print("Verified 21 all-order triangle coefficients and 8 complete bare/current/mixed enclosures.", flush=True)
    three_data = sourced_three_loop()
    print("Verified every original current and density on 144 sourced three-loop cubes; all 3 complete series specifications and 9 coefficients.", flush=True)
    dependencies = [ROOT / "code/htk/forest_cube_expansion.py", ROOT / "code/htk/forest_cubes.py",
                    ROOT / "code/htk/treeflow.py", ROOT / "code/htk/graph.py",
                    ROOT / "code/solve_product_forest_cubes.py"]
    return exact_json({"schema": 1,
                       "scope": "universal geometric rational expansion with complete exact polynomial data and certified finite tails",
                       "triangle": triangle_data, "sourced_three_loop": three_data,
                       "dependencies": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                                        for path in dependencies}})


def summary(result, encoded):
    three = result["sourced_three_loop"]
    records = three["bare_and_complete_first_current_spinor"]
    return {"schema": 1, "scope": result["scope"],
            "complete_certificate": {"file": OUTPUT.name,
                                     "encoding": "gzip-compressed complete UTF-8 JSON",
                                     "uncompressed_canonical_sha256": hashlib.sha256(encoded).hexdigest(),
                                     "uncompressed_bytes": len(encoded)},
            "triangle": {"all_order_even_term": "-3/2^(2*j+3)", "all_order_odd_term": "0",
                         "exact_sum": "-1/2", "complete_verified_terms": 21,
                         "complete_bare_first_current_and_mixed_enclosures": 8,
                         "requested_enclosure_width": "1/1000",
                         "complete_cli_result": "triangle_result.json"},
            "sourced_three_loop": {"ordered_edges": three["ordered_edges"],
                                   "marked": three["marked"], "all_balanced_sources": three["all_balanced_sources"],
                                   "complete_cubes_checked": three["complete_cubes"],
                                   "complete_original_current_entries_checked": 2*len(three["ordered_edges"])*three["complete_cubes"],
                                   "complete_original_density_values_checked": three["complete_cubes"],
                                   "series": [{"original_current_word": record["original_current_word"],
                                               "complete_rational_components": len(record["complete_exact_series_specification"]["components"]),
                                               "complete_finite_reference": record["complete_finite_reference"],
                                               "verified_complete_orders_0_through_2": record["verified_complete_orders_0_through_2"],
                                               "sufficient_order_for_width_1_over_1000": record["sufficient_order_for_width_1_over_1000"],
                                               "sufficient_order_coefficients_computed": False}
                                              for record in records]},
            "dependencies": result["dependencies"]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    arguments = parser.parse_args()
    result = certificate()
    encoded = json.dumps(result, separators=(",", ":"), sort_keys=True).encode("utf-8")
    compact = summary(result, encoded)
    if arguments.write:
        OUTPUT.write_bytes(gzip.compress(encoded, mtime=0))
        SUMMARY.write_text(json.dumps(compact, indent=2, sort_keys=True)+"\n")
        (HERE / "triangle_result.json").write_text(json.dumps(result["triangle"]["complete_cli_result"],
                                                             indent=2, sort_keys=True)+"\n")
        print("Wrote the complete compressed certificate, readable summary, and triangle_result.json.")
    else:
        if result != json.loads(gzip.decompress(OUTPUT.read_bytes())):
            raise ArithmeticError("the regenerated complete geometric certificate differs")
        if compact != json.loads(SUMMARY.read_text()):
            raise ArithmeticError("the regenerated certificate summary differs")
        if result["triangle"]["complete_cli_result"] != json.loads((HERE / "triangle_result.json").read_text()):
            raise ArithmeticError("the regenerated complete triangle CLI output differs")
        print("Complete geometric certificate and all original source/current results replay exactly.")
