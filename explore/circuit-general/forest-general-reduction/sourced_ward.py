"""Complete four/five-loop Ward identities with genuinely sourced ear vertices.

The three higher-loop quotients are integrated on their original graphs by
the new ear engine. Their boundary is computed separately: the proven
three-loop all-word Ward boundary is multiplied by the full ear-current
brackets and integrated with independent one-dimensional Lagrange rules.
No source-moment weights from forest_ears are used for that boundary.

Run: python3 explore/circuit-general/forest-general-reduction/sourced_ward.py
Use --write to regenerate the complete deterministic certificate.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
from itertools import combinations, product
import json
from math import factorial
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "code"))
sys.path.insert(0, str(HERE.parent / "forest-leibniz"))

from flint import fmpq_poly, fmpq_mpoly_ctx
from htk.graph import Graph
from htk.linalg import perm_sign
from htk.forest_ears import ForestEarMoments
from ward_three import (EDGES, MARKED, GENERIC_SOURCES, PHYSICAL_SOURCES,
                        boundary_value, coefficients, _fraction, _q)


def fixture(ears, sources):
    edges = list(EDGES)
    rows = [list(map(F, row)) for row in sources]
    added = []
    for i in range(ears):
        vertex = len(rows)
        source = (F(i + 2), F(3 * i - 3))
        added.append(source)
        edges.extend(((1, vertex), (vertex, 3)))
        rows.append(list(source))
        for a in range(2):
            rows[4][a] -= source[a]
    return Graph(len(rows), edges), tuple(map(tuple, rows)), tuple(added)


def quotient(graph, sources, position):
    left, right = graph.edges[position]
    retained = tuple(v for v in range(graph.nv) if v != right)
    mapping = {v: retained.index(left if v == right else v) for v in range(graph.nv)}
    edges = tuple((mapping[a], mapping[b]) for e, (a, b) in enumerate(graph.edges) if e != position)
    rows = [[F(0), F(0)] for _ in retained]
    for vertex, row in enumerate(sources):
        for a in range(2):
            rows[mapping[vertex]][a] += row[a]
    return Graph(len(retained), edges), tuple(map(tuple, rows)), mapping


def overfull(nv, edges):
    for size in range(2, nv):
        for vertices in combinations(range(nv), size):
            members = set(vertices)
            inside = tuple(e for e, edge in enumerate(edges) if set(edge) <= members)
            if len(inside) > 2 * size - 3:
                return {"vertices": vertices, "edge_indices": inside,
                        "edge_count": len(inside), "allowed_count": 2 * size - 3}
    return None


def face_census(graph, sources, ears):
    tree = []
    for position in range(graph.ne):
        reduced, _, _ = quotient(graph, sources, position)
        tree.append({"position": position, "edge": graph.edges[position], "sign": (-1) ** position,
                     "quotient_edges": tuple(reduced.edges), "zero_witness": overfull(reduced.nv, reduced.edges)})
    if tuple(row["position"] for row in tree if row["zero_witness"] is None) != (0, 1, 6):
        raise ArithmeticError("an additional nonzero tree face appeared")
    proper = []
    for size in range(3, graph.nv):
        for cut in combinations(range(graph.nv), size):
            members = set(cut)
            inside = tuple(e for e, edge in enumerate(graph.edges) if set(edge) <= members)
            if len(inside) != 2 * size - 3:
                continue
            rest = tuple(e for e in range(graph.ne) if e not in inside)
            retained = tuple(v for v in range(graph.nv) if v not in members) + (cut[0],)
            mapping = {v: retained.index(cut[0] if v in members else v) for v in range(graph.nv)}
            inner = tuple((cut.index(graph.edges[e][0]), cut.index(graph.edges[e][1])) for e in inside)
            outer = tuple((mapping[a], mapping[b]) for e, (a, b) in enumerate(graph.edges) if e in rest)
            proper.append({"cut": cut, "inner_original_edges": inside, "outer_original_edges": rest,
                           "inner_edges": inner, "outer_edges": outer,
                           "inner_zero_witness": overfull(size, inner),
                           "outer_zero_witness": overfull(len(retained), outer),
                           "sign": -(-1) ** ((size - 2) * (graph.loops - size + 2)) * perm_sign(inside + rest)})
    nonzero = [row for row in proper if row["inner_zero_witness"] is None and row["outer_zero_witness"] is None]
    if len(nonzero) != 1 or nonzero[0]["cut"] != tuple(v for v in range(graph.nv) if v not in (0, 5)):
        raise ArithmeticError("the higher-loop proper-face structure changed")
    if nonzero[0]["sign"] != (-1) ** ears:
        raise ArithmeticError("the sourced-ear proper-face orientation changed")
    return {"tree_faces": tuple(tree), "proper_faces": tuple(proper)}


@lru_cache(maxsize=None)
def interval_rule(degree):
    """Independent exact integral of each Lagrange polynomial at integer nodes."""
    x = fmpq_poly([0, 1])
    values = []
    for node in range(degree + 1):
        polynomial = fmpq_poly([1])
        for other in range(degree + 1):
            if other != node:
                polynomial *= (x - other) / (node - other)
        weight = sum(_fraction(polynomial[power]) / (power + 1) for power in range(degree + 1))
        values.append((F(node), weight))
    for power in range(degree + 1):
        if sum(weight * node ** power for node, weight in values) != F(1, power + 1):
            raise ArithmeticError("an exact interval moment was lost")
    return tuple(values)


def independent_boundary(ears, sources):
    """Bare sourced-ear lift of the full rational diamond--triangle boundary."""
    _, full, added = fixture(ears, sources)
    context = fmpq_mpoly_ctx.get(("j0", "j1"))
    j0, j1 = context.gens()
    polynomial = context.constant(1)
    for source in added:
        polynomial *= _q(source[0]) * j1 - _q(source[1]) * j0
    integral = F(0)
    for point in product(interval_rule(6 + ears), repeat=ears):
        ratios, weights = zip(*point)
        rows = [list(row) for row in full[:6]]
        weight = F(1)
        for source, ratio, local_weight in zip(added, ratios, weights):
            weight *= local_weight
            for a in range(2):
                rows[1][a] += (1 - ratio) * source[a]
                rows[3][a] += ratio * source[a]
        value = F(0)
        for exponent, coefficient in polynomial.terms():
            word = {(2, component): int(power) for component, power in enumerate(exponent) if power}
            value += _fraction(coefficient) * boundary_value(rows, word)
        integral += weight * value
    sign = (-1) ** (ears * (3 + ears) - ears * (ears - 1) // 2)
    return sign * integral / factorial(ears + 1)


def verify(ears, sources):
    graph, full, added = fixture(ears, sources)
    census = face_census(graph, full, ears)
    terms, total = [], {}
    for position in (0, 1, 6):
        reduced, rows, mapping = quotient(graph, full, position)
        marked = tuple(mapping[v] for v in MARKED)
        with ForestEarMoments(reduced, rows, marked, check=1) as engine:
            value = coefficients(engine.moment(allow_series=False))
            statistics, steps = dict(engine.last_stats), engine.steps
        sign = (-1) ** position
        for label, scalar in value.items():
            total[label] = total.get(label, F(0)) + sign * scalar
        terms.append({"position": position, "contracted_edge": graph.edges[position], "sign": sign,
                      "edges": tuple(reduced.edges), "sources": rows, "marked": marked,
                      "vertex_map": tuple(mapping[v] for v in range(graph.nv)),
                      "coefficients": value, "signed_coefficients": {label: sign * scalar for label, scalar in value.items()},
                      "ear_steps": steps, "stats": statistics})
    total = {label: scalar for label, scalar in total.items() if scalar}
    boundary = independent_boundary(ears, sources)
    if total != coefficients(boundary):
        raise ArithmeticError("the complete sourced higher-loop Ward expression disagrees with its rational boundary")
    if not any(term["coefficients"].get("z3", 0) for term in terms):
        raise ArithmeticError("the nontrivial period-cancellation control became zero")
    return {"loops": graph.loops, "ear_count": ears, "edges": tuple(graph.edges),
            "sources": full, "marked": MARKED, "added_nonzero_sources": added,
            "original_current_word": (), "face_census": census, "terms": tuple(terms),
            "complete_signed_sum": total, "independent_rational_boundary": boundary,
            "boundary_source_degree_bound": 6 + ears,
            "boundary_method": "independent tensor Lagrange integration of the exact base all-word spider boundary"}


def serial(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {str(key): serial(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [serial(item) for item in value]
    return value


def certificate():
    cases = []
    for name, sources in (("generic", GENERIC_SOURCES), ("physical_base", PHYSICAL_SOURCES)):
        for ears in (1, 2):
            record = verify(ears, sources)
            record["base_sources_name"] = name
            cases.append(record)
            print("%s, %d sourced loops: %s" % (name, record["loops"], record["independent_rational_boundary"]), flush=True)
    paths = (Path(__file__), ROOT / "code/htk/forest_ears.py", ROOT / "code/htk/forest_periods.py",
             ROOT / "code/htk/_forest_period_worker.py", HERE.parent / "forest-leibniz/ward_three.py")
    return serial({"scope": "Complete graph-level sourced-ear Ward identities; no full Wick/color/operator assembly is asserted.",
                   "records": cases, "max_terms": 1000000,
                   "source_hashes": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    actual = certificate()
    path = HERE / "sourced_ward_certificate.json"
    if args.write:
        path.write_text(json.dumps(actual, sort_keys=True, indent=2) + "\n")
        print("Wrote the complete sourced-ear Ward certificate.", flush=True)
    else:
        if actual != json.loads(path.read_text()):
            raise ArithmeticError("the sourced-ear Ward certificate changed")
        print("Every sourced-ear Ward certificate field replays.", flush=True)


if __name__ == "__main__":
    main()
