"""Complete exact controls for the direct current-cube series backend.

Records every three-loop cube, every original current and Jacobian entry,
complete non-series-parallel series levels, achieved one-loop enclosures,
and original-graph five-loop controls. No expanded current numerator or
finite period worker is used to construct the new series.
"""

import argparse
from dataclasses import asdict
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys
from time import perf_counter
from unittest.mock import patch

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "code"))
from htk import brute
from htk.forest_cubes import ForestCubes
from htk.graph import Graph


THREE = ((0, 1), (0, 3), (0, 4), (1, 2), (1, 4), (2, 3), (3, 4))
FIVE = ((0, 1), (1, 2), (0, 3), (0, 4), (2, 5), (2, 6),
        (3, 5), (3, 6), (4, 5), (4, 6), (5, 6))
OUTPUT = HERE / "cube_certificate.json"


def jsonable(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {key: jsonable(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [jsonable(item) for item in value]
    return value


def sources(vertices):
    rows = tuple((F(vertex+1, vertex+2), F((vertex+1)**2+1, vertex+3))
                 for vertex in range(vertices-1))
    return rows + (tuple(-sum(row[a] for row in rows) for a in range(2)),)


def independent_jacobian(cube, point):
    graph = cube.atlas.G
    times = cube.time_map(point)
    currents = brute.currents(graph, times, cube.atlas.sources)
    incidence = graph.incidence()
    result = [[[F(0)]*cube.dimension for _ in range(2)] for _ in graph.edges]
    for direction in range(cube.dimension):
        logarithmic = [F(exponent[direction], point[direction]) for exponent in cube.time_exponents]
        changed_sources = [[sum(F(incidence[edge][vertex])*logarithmic[edge]*currents[edge][a]
                               for edge in range(graph.ne)) for a in range(2)]
                           for vertex in range(graph.nv-1)]
        changed_currents = brute.currents(graph, times, changed_sources)
        for edge in range(graph.ne):
            for a in range(2):
                result[edge][a][direction] = changed_currents[edge][a]-logarithmic[edge]*currents[edge][a]
    return tuple(tuple(tuple(row) for row in edge) for edge in result)


def chart_control(atlas, position):
    cube = atlas.chart(position)
    check = cube.check_original()
    point, times = check["point"], check["times"]
    found, inverse = atlas.locate(times)
    if (found.position, inverse) != (position, point):
        raise ArithmeticError("the complete original resistance point did not recover its cube")
    jacobian = cube.current_jacobian(point)
    if jacobian != independent_jacobian(cube, point):
        raise ArithmeticError("an original current Jacobian entry differs from the differentiated grounded equations")
    return {"chart_position": position, "original_tree_edges": cube.tree,
            "ordered_cotree_edges": cube.cotree_order,
            "tree_reference_cotrees": cube.tree_reference_cotrees,
            "original_time_monomial_exponents": cube.time_exponents,
            "all_normalized_tree_weight_exponents": cube.weight_exponents,
            "orientation": cube.orientation,
            "exact_cube_point": point, "all_original_resistances": times,
            "all_original_current_spinors": cube.currents(point),
            "all_original_current_jacobian_entries": jacobian,
            "complete_signed_density": check["density"]}


def non_sp_series():
    full = ((2, -1), (0, 0), (-1, 3), (4, 2), (-5, -4))
    atlas = ForestCubes(Graph(5, THREE), full, marked=(3, 4), check=0, max_evaluations=144)
    references = (({}, {"1": F(2023, 96)}),
                  ({(1, 0): 1}, {"1": F(5477, 864), "z3": -F(304, 15)}),
                  ({(1, 1): 1}, {"1": -F(8205703, 302400), "z3": -F(6289, 2100)}))
    reference_path = HERE.parent / "forest-jet-identities" / "identity_certificate.json"
    archived = json.loads(reference_path.read_text(encoding="utf-8"))
    for powers, expected in references:
        target = next(record for record in archived["direct_period_and_rational_boundary_controls"]
                      if record["leg"] == 1 and record["degree"] == sum(powers.values())
                      and tuple(record["all_symmetric_spinor_powers"]) ==
                      (powers.get((1, 0), 0), powers.get((1, 1), 0)))
        if expected != {key: F(value) for key, value in target["complete_individual_graph_periods"][0].items()}:
            raise ArithmeticError("the complete independent finite-period reference differs")
    # Integral comparison for a decreasing positive function gives a rigorous
    # rational zeta(3) bracket, without treating a decimal as an exact value.
    count = 16
    partial_zeta = sum((F(1, k**3) for k in range(1, count+1)), F(0))
    zeta = (partial_zeta+F(1, 2*(count+1)**2), partial_zeta+F(1, 2*count**2))
    independent_midpoints = [F(0)]*3
    midpoint = (F(1, 2),)*atlas.dimension
    for cube in atlas.charts():
        times = cube.time_map(midpoint)
        original_currents = brute.currents(atlas.G, times, full)
        measure = F(1)
        for value in times:
            measure *= value
        for value in midpoint:
            measure /= value
        original_density = brute.density(atlas.G, times, full).rational()*measure
        for index, (powers, expected) in enumerate(references):
            independent_midpoints[index] += original_density*brute.monomial(original_currents, powers)
    records = []
    with patch("htk.forest_periods.CoreCurrentPeriods._source_numerator", side_effect=AssertionError("global current numerator used")), \
         patch("htk.forest_periods._research_module", side_effect=AssertionError("period or old Hepp backend used")):
        for index, (powers, expected) in enumerate(references):
            row = next(atlas.series(powers).terms(0))
            if row.partial_sum != independent_midpoints[index] or row.evaluations != 144:
                raise ArithmeticError("a complete non-SP midpoint level omits a chart or has a sign error")
            coefficient = expected.get("z3", F(0))
            endpoint_values = tuple(expected["1"]+coefficient*value for value in zeta)
            known_interval = (min(endpoint_values), max(endpoint_values))
            if not row.lower <= known_interval[0] <= known_interval[1] <= row.upper:
                raise ArithmeticError("a complete cube interval fails the independent finite-period reference")
            records.append({"original_current_word": [{"edge": edge, "spinor_entry": a, "power": power}
                                                       for (edge, a), power in powers.items()],
                            "complete_known_period": expected,
                            "known_period_rational_bracket": known_interval,
                            "complete_level_zero": asdict(row),
                            "independent_original_graph_midpoint_sum": independent_midpoints[index]})
    return {"ordered_edges": THREE, "all_balanced_sources": full,
            "marked_vertices": (3, 4), "complete_atlas": atlas.stats,
            "zeta3_rational_bracket": zeta,
            "zeta3_bracket_proof": "sum(k^(-3),k=1..16)+1/(2*17^2) <= zeta(3) <= sum(k^(-3),k=1..16)+1/(2*16^2)",
            "complete_bare_and_both_first_current_entries": records,
            "scope": "complete exact series level and valid interval; its midpoint sum is not the exact period; initial determinant bounds are conservative",
            "reference_certificate_sha256": hashlib.sha256(reference_path.read_bytes()).hexdigest()}


def triangle_series():
    graph = Graph(3, ((0, 1), (0, 2), (1, 2)))
    full = ((1, 0), (0, 1), (-1, -1))
    atlas = ForestCubes(graph, full, marked=(1, 2), check=0)
    polynomial = {(0,)*6: F(1), (1, 0, 0, 0, 0, 0): F(1), (1, 0, 0, 1, 0, 0): F(-2)}
    expected = -F(7, 12)
    rows = list(atlas.series(polynomial=polynomial).terms(3))
    previous = None
    total = F(0)
    for row in rows:
        total += row.term
        if total != row.partial_sum or not row.lower <= expected <= row.upper:
            raise ArithmeticError("a complete triangle series or remainder differs from its independent Dirichlet value")
        if abs(expected-row.partial_sum) > row.error_bound:
            raise ArithmeticError("a reported exact tail bound fails")
        if previous is not None and not previous.lower <= row.lower <= row.upper <= previous.upper:
            raise ArithmeticError("the complete certified intervals do not nest")
        previous = row
    reached = atlas.series({(0, 0): 2, (2, 1): 1}).enclose(F(1, 2))
    if not reached.lower <= -F(1, 15) <= reached.upper or reached.upper-reached.lower > F(1, 2):
        raise ArithmeticError("the requested rational interval width was not achieved")
    return {"ordered_edges": graph.edges, "all_balanced_sources": full,
            "complete_original_current_polynomial": [{"exponents": powers, "coefficient": coefficient}
                                                      for powers, coefficient in polynomial.items()],
            "independent_simplex_value": expected,
            "independent_simplex_identity": "-1/2-1/6+2/24=-7/12",
            "all_levels_zero_through_three": [asdict(row) for row in rows],
            "requested_interval_word": [{"edge": 0, "spinor_entry": 0, "power": 2},
                                        {"edge": 2, "spinor_entry": 1, "power": 1}],
            "requested_rational_width": F(1, 2), "independent_inserted_value": -F(1, 15),
            "achieved_complete_enclosure": asdict(reached)}


def certify():
    atlas = ForestCubes(Graph(5, THREE), sources(5), marked=(3, 4), check=0)
    all_charts = [chart_control(atlas, position) for position in range(atlas.chart_count)]
    print("Checked every original current/Jacobian/density and inverse in all144 three-loop cubes.", flush=True)
    higher = ForestCubes(Graph(7, FIVE), sources(7), marked=(5, 6), check=0)
    higher_checks = [chart_control(higher, position) for position in (0, higher.chart_count//2, higher.chart_count-1)]
    print("Checked complete original five-loop maps in three representative cubes; full atlas has23400 cubes.", flush=True)
    result = {"schema": 1,
              "scope": "complete source-preserving current maps and certified rational series; no finite partial sum is identified with the exact period",
              "three_loop_full_atlas": {"ordered_edges": THREE, "all_balanced_sources": sources(5),
                                        "marked_vertices": (3, 4), "counts": atlas.stats,
                                        "original_spanning_trees_in_weight_order": atlas.trees,
                                        "all_exact_tree_flow_current_entries": atlas._flows,
                                        "every_cube_control": all_charts},
              "five_loop_atlas_controls": {"ordered_edges": FIVE, "all_balanced_sources": sources(7),
                                           "marked_vertices": (5, 6), "counts": higher.stats,
                                           "original_spanning_trees_in_weight_order": higher.trees,
                                           "all_exact_tree_flow_current_entries": higher._flows,
                                           "selected_complete_cube_controls": higher_checks},
              "complete_non_series_parallel_series_controls": non_sp_series(),
              "complete_triangle_series_and_achieved_width": triangle_series()}
    paths = (Path(__file__), ROOT / "code/htk/forest_cubes.py", ROOT / "code/tests/test_forest_cubes.py",
             ROOT / "code/htk/brute.py", ROOT / "code/htk/treeflow.py", ROOT / "code/htk/graph.py")
    result["implementation_sha256"] = {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                                       for path in paths}
    return jsonable(result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    start = perf_counter()
    result = certify()
    if args.write:
        OUTPUT.write_text(json.dumps(result, separators=(",", ":"))+"\n", encoding="utf-8")
        print("Wrote the complete cube/current/series certificate.")
    elif not OUTPUT.is_file() or json.loads(OUTPUT.read_text(encoding="utf-8")) != result:
        raise ArithmeticError("the complete regenerated cube certificate differs")
    else:
        print("Every complete cube map, interval, series term, and reference agrees.")
    print("Elapsed %.3f seconds." % (perf_counter()-start), flush=True)


if __name__ == "__main__":
    main()
