"""Exact full first/second current moments of the (5,0) HT tadpole.

The bare density is derived by a polynomial Pfaffian identity, not fitted.
Every original edge is retained.  The root-edge submatrix is the complete
order-two input required by the 3d supercharge.  Periods use the existing
exact half-integer scalar evaluator, after summing full graph orbits.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations
import json
from math import prod
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path[:0] = [str(ROOT / "code"), str(ROOT / "explore/circuit-general/cores")]
from flint import fmpq_mpoly_ctx
from htk.graph import Graph
from htk import brute
from htk.linalg import det
import numer
import cores
import gutil
from scalar_cache import fixed_parameter_cache

EDGES = tuple(sorted(((0, 1), (1, 2), (2, 3), (3, 4), (0, 4),
                      (1, 5), (2, 5), (3, 5), (4, 5))))
SOURCES = ((1,), (0,), (0,), (0,), (0,), (-1,))


def rational(x):
    if hasattr(x, "to_dict"):
        coefficients = x.to_dict()
        assert all(not any(e) for e in coefficients)
        return sum((rational(c) for c in coefficients.values()), Fraction(0))
    if hasattr(x, "p"):
        return Fraction(int(x.p), int(x.q))
    return Fraction(x)


def sign(p):
    return (-1) ** sum(p[i] > p[j] for i in range(len(p))
                       for j in range(i + 1, len(p)))


def determinant(a):
    n = len(a)
    return sum(sign(p) * prod(a[i][p[i]] for i in range(n))
               for p in permutations(range(n)))


def serialize_polynomial(p):
    return [[list(map(int, e)), str(rational(c))] for e, c in sorted(p.to_dict().items())]


def atomic_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def derive_density():
    """Return P,U,(U*y_e), with D=P/U^(9/2), by exact polynomial division."""
    context = fmpq_mpoly_ctx.get(tuple("t%d" % i for i in range(9)))
    times = context.gens()
    graph = Graph(6, list(EDGES))
    cycles = graph.cycle_matrix()
    matrix = [[sum(cycles[i][e] * cycles[j][e] * times[e]
                   for e in range(9)) for j in range(4)] for i in range(4)]
    kirchhoff = determinant(matrix)
    adjugate = [[(-1) ** (i + j) * determinant(
        [[matrix[a][b] for b in range(4) if b != i]
         for a in range(4) if a != j]) for j in range(4)] for i in range(4)]
    assert all(sum(matrix[i][k] * adjugate[k][j] for k in range(4))
               == (kirchhoff if i == j else 0)
               for i in range(4) for j in range(4))
    response = [[sum(cycles[a][i] * adjugate[a][b] * cycles[b][j]
                     for a in range(4) for b in range(4))
                 for j in range(9)] for i in range(9)]
    tree_currents = numer.tree_current_polys(graph, SOURCES)
    currents = [context.from_dict(tree_currents[e][0]) for e in range(9)]
    incidence = graph.incidence()
    assert all(sum(incidence[e][v] * currents[e] for e in range(9))
               == SOURCES[v][0] * kirchhoff for v in range(5))
    assert all(sum(cycles[a][e] * times[e] * currents[e] for e in range(9)) == 0
               for a in range(4))
    # Extract dt_0...dt_7 from topological four-form wedge holomorphic
    # four-form.  In the four topological edge slots Wick's three pairings
    # have the same sign; each off-diagonal covariance is -response/(2U).
    numerator = context.constant(0)
    partitions = []
    for topological in combinations(range(8), 4):
        holomorphic = tuple(i for i in range(8) if i not in topological)
        coefficient = (sign(topological + holomorphic)
                       * int(det([[r[e] for e in topological] for r in cycles]))
                       * int(det([[r[e] for e in holomorphic] for r in cycles])))
        if not coefficient:
            continue
        i, j, k, l = topological
        wick = (response[i][j] * response[k][l]
                + response[i][k] * response[j][l]
                + response[i][l] * response[j][k])
        numerator += coefficient * wick * prod(currents[e] for e in holomorphic)
        partitions.append([list(topological), list(holomorphic), coefficient])
    # Before cancellation D=numerator/(4*t_8*U^(15/2)).  The full exact
    # division proves the lower pole order, including every coefficient.
    polynomial, remainder = divmod(numerator, 4 * times[8] * kirchhoff ** 3)
    assert not remainder, "Pfaffian numerator did not divide exactly"
    assert {sum(e) for e in polynomial.to_dict()} == {9}
    assert len(partitions) == 16
    controls = []
    for values in (tuple(range(1, 10)), (2, 5, 3, 7, 11, 13, 17, 19, 23)):
        substitutions = dict(zip(tuple("t%d" % i for i in range(9)), values))
        p = rational(polynomial.subs(substitutions))
        u = rational(kirchhoff.subs(substitutions))
        actual = brute.density(graph, values, SOURCES, T=1, check_basic=True)
        assert actual.radicand == u
        assert actual.coeff * u ** 4 == p
        y = brute.currents(graph, values, SOURCES)
        assert all(rational(currents[e].subs(substitutions)) == y[e][0] * u
                   for e in range(9))
        controls.append({"times": list(values), "U": str(u), "P": str(p),
                         "density_rational_coefficient": str(actual.coeff),
                         "all_currents": [str(row[0]) for row in y]})
    return context, polynomial, kirchhoff, currents, partitions, controls


def dependencies():
    paths = [Path(__file__), HERE / "scalar_cache.py", ROOT / "code/htk/brute.py", ROOT / "code/htk/graph.py",
             ROOT / "code/htk/linalg.py"]
    paths += sorted((ROOT / "explore/circuit-general/cores").glob("*.py"))
    return {str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in paths}


def _compute(output, root_only=False):
    start = time.monotonic()
    ctx, polynomial, kirchhoff, currents, partitions, controls = derive_density()
    print("exact density", len(polynomial.to_dict()), "terms; U", len(kirchhoff.to_dict()),
          "terms", round(time.monotonic()-start, 2), "s", flush=True)
    base = {"schema": "tadpole50-complete-moments-v1", "edges": EDGES,
            "sources": SOURCES, "H": 1, "T": 1, "pole": "9/2",
            "bare_numerator": serialize_polynomial(polynomial),
            "kirchhoff": serialize_polynomial(kirchhoff),
            "current_numerators": [serialize_polynomial(p) for p in currents],
            "pfaffian_partitions": partitions, "original_graph_controls": controls,
            "implementation_sha256": dependencies(), "complete": False, "moments": []}
    atomic_json(output, base)
    selected = [i for i, e in enumerate(EDGES) if e[1] == 5] if root_only else list(range(9))
    words = [()] + [(e,) for e in selected] + list(combinations(selected, 2)) + [(e, e) for e in selected]
    words.sort(key=lambda w: (len(w), w))
    automorphisms = gutil.edge_automorphisms(EDGES)
    period_cache = {}
    orbit_cache = {}
    for word in words:
        p = polynomial * prod(currents[e] for e in word)
        orbit = {}
        for exponent, coefficient in p.to_dict().items():
            alpha = tuple(int(k) + 1 for k in exponent)
            if alpha not in orbit_cache:
                orbit_cache[alpha] = gutil.orbit_rep(alpha, automorphisms)
            representative = orbit_cache[alpha]
            orbit[representative] = orbit.get(representative, 0) + rational(coefficient)
        value = {}
        pole = Fraction(9, 2) + len(word)
        for alpha, coefficient in orbit.items():
            if not coefficient:
                continue
            key = (pole, alpha)
            if key not in period_cache:
                period_cache[key] = cores.period(EDGES, alpha, pole)
                if len(period_cache) % 50 == 0:
                    print("scalar periods", len(period_cache), "word", word,
                          "seconds", round(time.monotonic()-start, 2), flush=True)
            for label, scalar in period_cache[key].items():
                value[label] = value.get(label, 0) + coefficient * scalar
        value = {cores.enc(label): str(v) for label, v in value.items() if v}
        record = {"edges": list(word), "pole": str(pole),
                  "numerator_terms": len(p.to_dict()),
                  "nonzero_scalar_orbits": sum(bool(c) for c in orbit.values()),
                  "value": value}
        base["moments"].append(record)
        atomic_json(output, base)
        print("moment", word, value, "seconds", round(time.monotonic()-start, 2), flush=True)
    base["complete"] = True
    base["scope"] = "all root-edge words of degree at most two" if root_only else "all original-edge words of degree at most two"
    base["scalar_periods_evaluated"] = len(period_cache)
    atomic_json(output, base)


def compute(output, root_only=False):
    with fixed_parameter_cache() as stats:
        _compute(output, root_only)
    data = json.loads(Path(output).read_text())
    data["exact_memoization"] = stats
    atomic_json(output, data)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, default=HERE / "tadpole_moments.json")
    ap.add_argument("--root-only", action="store_true")
    args = ap.parse_args()
    compute(args.output, args.root_only)
