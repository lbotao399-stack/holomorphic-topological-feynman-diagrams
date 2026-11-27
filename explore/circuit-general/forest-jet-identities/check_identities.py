"""Independent exact controls for the complete source-simplex reduction.

All current entries, all 62 proper time subsets, and all balanced sources
are retained.  Finite sample zeros are recorded as probes, never used as
the proof of a universally vanishing face.  Original graph integrals use
CoreCurrentPeriods and do not use the boundary identity.
"""

import argparse
from fractions import Fraction
from functools import lru_cache
import hashlib
from itertools import combinations, permutations
import json
from pathlib import Path
import sys
from time import perf_counter

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "code"))
sys.path.insert(0, str(HERE.parent / "universal-product-jets"))

from htk import brute
from htk.forest_periods import CoreCurrentPeriods
from htk.graph import Graph
from htk.linalg import det
from htk.star import Dual
from htk.treeflow import tree_flow
from outer_faces import (CORE_EDGES, CHECK_POINT, SourceTriangleOuterFaces,
                         _fraction, _q, _time_map)
from boundary_faces import SourceTriangleTimeFaces


GRAPHS = (
    ((0, 1), (0, 3), (0, 4), (1, 2), (1, 4), (2, 3), (3, 4)),
    ((0, 1), (0, 3), (0, 4), (1, 2), (2, 3), (2, 4), (3, 4)),
    ((0, 1), (0, 3), (0, 4), (1, 2), (1, 3), (2, 4), (3, 4)),
)
VECTORS = ((Fraction(2), Fraction(-1)), (Fraction(-1), Fraction(3)),
           (Fraction(4), Fraction(2)))
OUTPUT = HERE / "identity_certificate.json"


def _jsonable(value):
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, dict):
        return {key: _jsonable(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_jsonable(item) for item in value]
    return value


def original_sources(vectors=VECTORS):
    P, Q, R = vectors
    B = tuple(-P[a]-Q[a]-R[a] for a in range(2))
    return ((P, (0, 0), Q, R, B), (P, Q, (0, 0), R, B), (P, (0, 0), Q, R, B))


@lru_cache(maxsize=None)
def _tree_data(vertices, edges):
    graph = Graph(vertices, edges)
    source_basis = [tuple(Fraction(i == j) for j in range(vertices-1)) for i in range(vertices-1)]
    rows = []
    for tree in combinations(range(len(edges)), vertices-1):
        try:
            graph._check_tree(tree)
        except ValueError:
            continue
        response = tree_flow(graph, source_basis, tree)
        rows.append((tree, tuple(edge for edge in range(len(edges)) if edge not in tree), response))
    return tuple(rows)


def tree_currents(vertices, edges, times, sources, shrinking=None):
    """Independent full-source tree average, optionally at a proper scale face."""
    trees = _tree_data(vertices, edges)
    if shrinking is not None:
        chosen = set(shrinking)
        degree = min(len(chosen.intersection(cotree)) for tree, cotree, response in trees)
        trees = tuple(row for row in trees if len(chosen.intersection(row[1])) == degree)
    zero = times[0] * 0
    denominator = zero
    numerator = [[zero for _ in range(2)] for edge in edges]
    for tree, cotree, response in trees:
        weight = zero + 1
        for edge in cotree:
            weight *= times[edge]
        denominator += weight
        for edge in range(len(edges)):
            for a in range(2):
                flow = sum((response[edge][vertex] * sources[vertex][a]
                            for vertex in range(vertices-1)), zero)
                numerator[edge][a] += weight * flow
    return tuple(tuple(value/denominator for value in row) for row in numerator)


def _simplex_sources(x, y, triangle, vectors=VECTORS):
    P, Q, R = vectors
    B = tuple(-P[a]-Q[a]-R[a] for a in range(2))
    if triangle == "TWA":
        return (tuple(P[a]+x*Q[a] for a in range(2)),
                tuple((1-x-y)*Q[a] for a in range(2)),
                tuple(R[a]+y*Q[a] for a in range(2)), B)
    return (P, tuple((1-x-y)*Q[a] for a in range(2)),
            tuple(R[a]+y*Q[a] for a in range(2)),
            tuple(B[a]+x*Q[a] for a in range(2)))


def radial_maps():
    records = []
    P, Q, R = VECTORS
    B = tuple(-P[a]-Q[a]-R[a] for a in range(2))
    point = (Fraction(2), Fraction(3), Fraction(5), Fraction(7), Fraction(11), Fraction(13, 17))
    parameters = tuple(Dual.var(value, j, 6) for j, value in enumerate(point))
    a, b, c, d, e, r = parameters
    maps = ((a, d, e, r*b, c, (1-r)*b, 1),
            ((1-r)*a, d, e, r*a, b, c, 1),
            (a, d, e, r*c, b, (1-r)*c, 1))
    core_sources = (
        (P, tuple((1-r)*Q[i] for i in range(2)), tuple(R[i]+r*Q[i] for i in range(2)), B),
        (tuple(P[i]+r*Q[i] for i in range(2)), tuple((1-r)*Q[i] for i in range(2)), R, B),
        (P, tuple((1-r)*Q[i] for i in range(2)), R, tuple(B[i]+r*Q[i] for i in range(2))),
    )
    for position, (edges, time_map, sources, moved_sources) in enumerate(
            zip(GRAPHS, maps, original_sources(), core_sources)):
        graph = Graph(5, edges)
        time_jacobian = det([value.g for value in time_map[:-1]])
        expected_jacobian = (-point[1], point[0], point[2])[position]
        if time_jacobian != expected_jacobian:
            raise ArithmeticError("an original series time-map orientation differs")
        times = tuple(value.v if isinstance(value, Dual) else Fraction(value) for value in time_map)
        original = brute.currents(graph, times, sources)
        core = tree_currents(4, CORE_EDGES, parameters[:5]+(Dual(1, n=6),), moved_sources)
        for old_edge, core_edge in zip((1, 2, 6), (3, 4, 5)):
            if original[old_edge] != [value.v for value in core[core_edge]]:
                raise ArithmeticError("a source-bearing radial map changed a common chord")
        volume = det([core[edge][component].g for edge in (3, 4, 5) for component in range(2)])
        expected_density = brute.density(graph, times, sources).rational() * abs(time_jacobian)
        signed_volume = (1, -1, -1)[position] * volume
        if expected_density != signed_volume:
            raise ArithmeticError("the full original density does not match its signed current volume")
        tree = (0, 3, 4, 5)
        cycles = graph.fundamental_cycles(tree)
        source_locations = (0, 2, 3) if position != 1 else (0, 1, 3)
        source_basis = [tuple(int(vertex == location) for location in source_locations) for vertex in range(4)]
        flow = tree_flow(graph, source_basis, tree)
        affine = tuple(tuple(flow[edge])+tuple(Fraction(cycles[j][edge]) for j in range(3))
                       for edge in range(7))
        basis = (P, Q, R, original[1], original[2], original[6])
        for edge in range(7):
            for component in range(2):
                if sum(affine[edge][j] * basis[j][component] for j in range(6)) != original[edge][component]:
                    raise ArithmeticError("an original current or source was lost in affine reconstruction")
        records.append({"graph_position": position+1, "ordered_edges": edges,
                        "all_balanced_sources": sources, "original_times": times,
                        "series_coordinate_jacobian": time_jacobian,
                        "complete_pulled_back_kernel_density": expected_density,
                        "full_original_current_affine_basis": ("lambda_W", "lambda_c", "lambda_A", "J_WA", "J_WB", "J_AB"),
                        "all_seven_original_current_coefficients": affine})
    return records


def outer_form_controls():
    records = []
    P, Q, R = VECTORS
    B = tuple(-P[a]-Q[a]-R[a] for a in range(2))
    for triangle in ("TWA", "TBA"):
        chart = SourceTriangleOuterFaces(VECTORS, triangle)
        check = chart.pointwise_check()
        parameters = tuple(Dual.var(value, j, 6) for j, value in enumerate(CHECK_POINT))
        times, ratio = _time_map(parameters, chart.path)
        if triangle == "TWA":
            sources = (tuple(P[a]+(1-ratio)*Q[a] for a in range(2)), (0, 0),
                       tuple(R[a]+ratio*Q[a] for a in range(2)), B)
        else:
            sources = (P, (0, 0), tuple(R[a]+ratio*Q[a] for a in range(2)),
                       tuple(B[a]+(1-ratio)*Q[a] for a in range(2)))
        exact = tree_currents(4, CORE_EDGES, times, sources)
        full_jacobian = det([exact[edge][component].g for edge in (3, 4, 5) for component in range(2)])
        polynomial_jacobian = _fraction(chart.jacobian(*map(_q, CHECK_POINT)))
        if full_jacobian != polynomial_jacobian:
            raise ArithmeticError("the complete corrected outer Jacobian differs from independent tree averaging")
        uncorrected = list(chart.currents)
        if triangle == "TWA":
            uncorrected[:2] = chart.source_half_current
        else:
            uncorrected[-2:] = tuple(-value for value in chart.source_half_current)
        half_jacobian = det([[_fraction(value.derivative(j)(*map(_q, CHECK_POINT))) for j in range(6)]
                             for value in uncorrected])
        if full_jacobian == half_jacobian:
            raise ArithmeticError("the chosen outer control did not detect the omitted affine-current correction")
        records.append({"triangle": triangle, "chart_variables": ("s1", "s2", "w1", "w2", "rX", "rQ"),
                        "exact_parameter_point": CHECK_POINT,
                        "coordinate_orientation": chart.orientation,
                        "coordinate_jacobian": chart.coordinate_jacobian_at_check,
                        "full_corrected_jacobian": full_jacobian,
                        "incorrect_half_edge_jacobian": half_jacobian,
                        "nonzero_required_correction": full_jacobian-half_jacobian,
                        "grounded_core_current_check": check})
    return records


def greedy_tree_controls():
    trees = _tree_data(4, CORE_EDGES)
    count = 0
    for order in permutations(range(6)):
        parent = list(range(4))
        def find(vertex):
            while parent[vertex] != vertex:
                vertex = parent[vertex]
            return vertex
        selected = []
        for edge in order:
            a, b = CORE_EDGES[edge]
            if find(a) != find(b):
                parent[find(b)] = find(a)
                selected.append(edge)
        prefixes = [set(order[:length]) for length in range(1, 6)]
        greedy = tuple(len(prefix.difference(selected)) for prefix in prefixes)
        constant_terms = 0
        for tree, cotree, response in trees:
            exponents = tuple(len(prefix.intersection(cotree)) for prefix in prefixes)
            if any(power < base for power, base in zip(exponents, greedy)):
                raise ArithmeticError("the greedy monomial does not divide a complete tree-current summand")
            constant_terms += exponents == greedy
        if constant_terms != 1:
            raise ArithmeticError("the residual Kirchhoff polynomial lacks its unique constant-one tree term")
        count += 1
    return {"all_ordered_Hepp_sectors": count, "all_core_spanning_trees": len(trees),
            "every_current_numerator_has_same_common_monomial": True,
            "residual_denominator_constant_coefficient": 1}


def proper_face_probes():
    records = []
    for triangle in ("TWA", "TBA"):
        faces = []
        for mask in range(1, 63):
            selected = tuple(edge for edge in range(6) if mask >> edge & 1)
            complement = tuple(edge for edge in range(6) if not mask >> edge & 1)
            independent = selected[:-1]+complement[:-1]
            times = [Dual(1, n=6) for _ in range(6)]
            for j, edge in enumerate(independent):
                times[edge] = Dual.var(Fraction(edge+2, edge+1), j, 6)
            x, y = Dual.var(Fraction(1, 5), 4, 6), Dual.var(Fraction(2, 7), 5, 6)
            sources = _simplex_sources(x, y, triangle)
            currents = tree_currents(4, CORE_EDGES, times, sources, selected)
            jacobian = det([currents[edge][component].g for edge in (3, 4, 5) for component in range(2)])
            faces.append({"shrinking_edge_subset": selected, "time_face_coordinate_edges": independent,
                          "complete_six_current_jacobian": jacobian})
        nonzero = [face["shrinking_edge_subset"] for face in faces if face["complete_six_current_jacobian"]]
        if nonzero != ([(2,)] if triangle == "TWA" else [(0,)]):
            raise ArithmeticError("the complete generic proper-face probe changed")
        records.append({"triangle": triangle, "source_simplex_point": (Fraction(1, 5), Fraction(2, 7)),
                        "scope": "all 62 exact pointwise leading-tree probes; sample zeros are not a universal vanishing proof",
                        "every_proper_subset": faces})
    return records


def _multiply(left, right):
    out = {}
    for a, ca in left.items():
        for b, cb in right.items():
            powers = tuple(i+j for i, j in zip(a, b))
            out[powers] = out.get(powers, Fraction(0)) + ca*cb
    return {powers: coefficient for powers, coefficient in out.items() if coefficient}


def _word(leg, first, second):
    answer = {(0,)*6: Fraction(1)}
    for component, power in enumerate((first, second)):
        direct = [0]*6
        direct[component] = 1
        factor = {tuple(direct): Fraction(1)}
        if leg == 2:
            ab = [0]*6
            ab[4+component] = 1
            factor = {tuple(ab): Fraction(1), tuple(direct): Fraction(-1),
                      (0,)*6: -VECTORS[2][component]}
        for _ in range(power):
            answer = _multiply(answer, factor)
    return answer


def _period(value):
    if isinstance(value, (Fraction, int)):
        return {"1": Fraction(value)} if value else {}
    return dict(value)


def direct_period_controls():
    engines = [CoreCurrentPeriods(Graph(5, edges), sources, marked=(3, 4), check=1)
               for edges, sources in zip(GRAPHS, original_sources())]
    outer = [SourceTriangleOuterFaces(VECTORS, triangle) for triangle in ("TWA", "TBA")]
    time = [SourceTriangleTimeFaces(VECTORS, triangle) for triangle in ("TWA", "TBA")]
    records = []
    try:
        for degree in range(4):
            for first in range(degree+1):
                for leg, edges in enumerate(((1, 1, 1), (5, 4, 4)), 1):
                    word = _word(leg, first, degree-first)
                    values = []
                    for engine, edge in zip(engines, edges):
                        exponents = [0]*14
                        exponents[2*edge] = first
                        exponents[2*edge+1] = degree-first
                        values.append(_period(engine.integrate({tuple(exponents): Fraction(1)}, allow_series=False)))
                    pairs = []
                    for partner in (1, 2):
                        summed = {label: values[0].get(label, Fraction(0))+values[partner].get(label, Fraction(0))
                                  for label in set(values[0]) | set(values[partner])}
                        summed = {label: value for label, value in summed.items() if value}
                        if set(summed)-{"1"}:
                            raise ArithmeticError("a nonrational coefficient survived the complete direct graph pair")
                        exterior = outer[partner-1].integrate(word)
                        boundary = time[partner-1].integrate(word)
                        if summed.get("1", Fraction(0)) != exterior-boundary:
                            raise ArithmeticError("the complete rational-boundary algorithm differs from direct original graph integration")
                        pairs.append({"original_graph_positions": (1, partner+1), "outer_complete_current_form": exterior,
                                      "all_time_faces_sum": boundary, "complete_rational_pair": exterior-boundary})
                    records.append({"degree": degree, "leg": leg, "all_symmetric_spinor_powers": (first, degree-first),
                                    "common_current_polynomial": [{"powers": powers, "coefficient": coefficient}
                                                                  for powers, coefficient in sorted(word.items())],
                                    "complete_individual_graph_periods": values, "both_exact_boundary_identities": pairs})
            print("Checked all complete direct tensor entries through degree", degree, flush=True)
    finally:
        for engine in engines:
            engine.close()
    return records


def certify():
    result = {"schema": 1, "scope": "complete original-current/source maps, smooth Hepp charts, every proper time face, and direct original-graph controls",
              "independent_full_source_spinors": VECTORS,
              "radial_original_graph_maps": radial_maps(),
              "complete_outer_affine_current_correction": outer_form_controls(),
              "all_greedy_tree_sector_checks": greedy_tree_controls(),
              "all_proper_time_face_probes": proper_face_probes(),
              "direct_period_and_rational_boundary_controls": direct_period_controls(),
              "complete_rational_pair_checks": 40,
              "all_degree_proof": "DERIVATION.md; finite samples are independent controls, not extrapolation"}
    paths = (HERE / "check_identities.py", HERE / "outer_faces.py",
             HERE.parent / "universal-product-jets" / "boundary_faces.py",
             ROOT / "code" / "htk" / "brute.py", ROOT / "code" / "htk" / "treeflow.py",
             ROOT / "code" / "htk" / "forest_periods.py")
    result["implementation_sha256"] = {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                                       for path in paths}
    return _jsonable(result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    start = perf_counter()
    result = certify()
    if args.write:
        OUTPUT.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
        print("Wrote the complete exact source-simplex identity certificate.")
    elif not OUTPUT.is_file() or json.loads(OUTPUT.read_text(encoding="utf-8")) != result:
        raise ArithmeticError("the regenerated complete identity certificate differs")
    else:
        print("Verified all original maps, both corrected outer forms, 720 sectors, 124 proper-face probes, and 40 complete period-pair identities.")
    print("Elapsed %.3f seconds." % (perf_counter()-start), flush=True)


if __name__ == "__main__":
    main()
