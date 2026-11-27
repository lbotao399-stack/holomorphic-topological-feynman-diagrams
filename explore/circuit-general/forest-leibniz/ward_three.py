"""Complete current-word Ward relation with a rational spider boundary.

The ordered sliding graph is
    01, 05, 13, 14, 23, 24, 25, 34,
with marked insertions A=3, B=4. Its three nonzero edge contractions
occur at positions 0,1,6 with signs +,-,+. Its only proper tight face
with a nonparallel quotient is {1,2,3,4}: a two-loop diamond composed
with an outer triangle. Both are rational spider graphs.

Every finite original-edge current word is supported. When the word hits
the contracted edge a->b, its limiting current is recovered from the
Kirchhoff equation at a, including the complete source lambda_a.
The proper face is integrated independently using its polynomial spider
chart and exact Dirichlet moments. No source row or period term is removed.

Run: python3 explore/circuit-general/forest-leibniz/ward_three.py
"""
import argparse
from collections.abc import Mapping
from fractions import Fraction as F
from functools import lru_cache
import hashlib
from itertools import combinations
import json
from math import factorial
from numbers import Integral
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "code"))

from flint import fmpq, fmpq_mpoly_ctx
from htk.graph import Graph
from htk.forest_periods import CoreCurrentPeriods


EDGES = ((0, 1), (0, 5), (1, 3), (1, 4), (2, 3), (2, 4), (2, 5), (3, 4))
MARKED = (3, 4)
INNER_EDGES = ((0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
INNER_ORIGINAL = (2, 3, 4, 5, 7)
OUTER_EDGES = ((0, 2), (0, 1), (2, 1))
OUTER_ORIGINAL = (0, 1, 6)
GENERIC_SOURCES = ((1, 2), (3, -1), (-2, 4), (5, 3), (-9, -4), (2, -4))
PHYSICAL_SOURCES = ((1, 2), (0, 0), (0, 0), (5, 3), (-8, -1), (2, -4))


def _integer(value, name):
    if isinstance(value, bool) or not isinstance(value, Integral) or value < 0:
        raise ValueError(name + " must be a nonnegative integer")
    return int(value)


def _sources(sources):
    rows = tuple(tuple(row) for row in sources)
    if len(rows) != 6 or any(len(row) != 2 for row in rows):
        raise ValueError("six complete two-component source rows are required")
    if any(isinstance(value, float) for row in rows for value in row):
        raise TypeError("source values must be exact rational numbers")
    rows = tuple(tuple(F(value) for value in row) for row in rows)
    if any(sum(row[i] for row in rows) for i in range(2)):
        raise ValueError("the complete sources must sum to zero")
    return rows


def _word(powers):
    if powers is not None and not isinstance(powers, Mapping):
        raise TypeError("current powers must map (original edge, component) to integers")
    out = []
    for location, value in ({} if powers is None else powers).items():
        if not isinstance(location, tuple) or len(location) != 2:
            raise ValueError("current locations are (edge index, component)")
        edge, component = (_integer(location[0], "edge"), _integer(location[1], "component"))
        power = _integer(value, "current power")
        if edge >= len(EDGES) or component >= 2:
            raise ValueError("the current location lies outside the sliding graph")
        if power:
            out.append((edge, component, power))
    return tuple(sorted(out))


def _q(value):
    value = F(value)
    return fmpq(value.numerator, value.denominator)


def _fraction(value):
    return F(int(value.p), int(value.q))


def _add(*vectors):
    return tuple(sum(vector[i] for vector in vectors) for i in range(2))


def _scale(scalar, vector):
    return tuple(scalar * value for value in vector)


def _bracket(left, right):
    return left[0] * right[1] - left[1] * right[0]


def coefficients(value):
    if isinstance(value, (Integral, F)):
        return {"1": F(value)} if value else {}
    return {key: F(scalar) for key, scalar in value.items() if scalar}


def _sum_periods(groups):
    out = {}
    for group in groups:
        for key, value in group.items():
            out[key] = out.get(key, F(0)) + value
    return {key: value for key, value in out.items() if value}


def proper_faces():
    """All proper tight faces, including the independently vanishing quotients."""
    result = []
    for size in range(3, 6):
        for cut in combinations(range(6), size):
            inner = tuple(i for i, edge in enumerate(EDGES) if set(edge) <= set(cut))
            if len(inner) != 2 * size - 3:
                continue
            rest = tuple(i for i in range(len(EDGES)) if i not in inner)
            vertices = tuple(v for v in range(6) if v not in cut) + (cut[0],)
            mapping = {v: vertices.index(cut[0] if v in cut else v) for v in range(6)}
            outer = tuple((mapping[a], mapping[b]) for i, (a, b) in enumerate(EDGES) if i in rest)
            parallel = len(set(map(frozenset, outer))) != len(outer)
            order = inner + rest
            parity = (-1) ** sum(a > b for i, a in enumerate(order) for b in order[i + 1:])
            inner_loops, outer_loops = size - 2, 5 - size
            result.append(dict(cut=cut, inner_edges=inner, outer_edges=outer, parallel_quotient=parallel,
                               boundary_sign=-(-1) ** (inner_loops * outer_loops) * parity))
    surviving = [row for row in result if not row["parallel_quotient"]]
    if len(surviving) != 1 or surviving[0]["cut"] != (1, 2, 3, 4) or surviving[0]["boundary_sign"] != 1:
        raise ArithmeticError("the specified sliding graph has an unexpected proper-face decomposition")
    return tuple(result)


def diamond_first_spinor(left_source, right_source, A_source):
    """Universal two-loop diamond first A-spoke moment in spinor brackets."""
    u = _bracket(left_source, right_source)
    v = _bracket(left_source, A_source)
    w = _bracket(right_source, A_source)
    return _scale(F(1, 288), _add(
        _scale(3 * (2 * v + u) * (2 * w - u), A_source),
        _scale((3 * v + 2 * u) * (2 * w - u), right_source),
        _scale(-3 * (2 * v + u) * (3 * w - 2 * u), left_source)))


def boundary_first_spinor(sources):
    """Universal rational first-spoke boundary, using spinor brackets only.

    The three vectors in ``vertices`` are the affine shifted first source
    at the outer triangle's corners. The formulas below are the exact first,
    second, and third Dirichlet moments of those affine vectors; consequently
    this evaluates the integrated formula without expanding any coordinate.
    """
    sources = _sources(sources)
    first, second, A = sources[0], sources[5], sources[3]
    total = _add(sources[1], sources[2], first, second)
    vertices = (_add(sources[1], first), sources[1], _add(sources[1], first, second))
    q = tuple(_bracket(vector, _add(total, _scale(2, A))) for vector in vertices)
    r = tuple(_bracket(vector, _add(_scale(2, total), _scale(3, A))) for vector in vertices)
    c = _bracket(total, A)
    vertex_sum, q_sum, r_sum = _add(*vertices), sum(q), sum(r)
    qr_sum = sum(a * b for a, b in zip(q, r))
    Pq_sum = _add(*(_scale(a, vector) for a, vector in zip(q, vertices)))
    Pr_sum = _add(*(_scale(b, vector) for b, vector in zip(r, vertices)))
    Pqr_sum = _add(*(_scale(a * b, vector) for a, b, vector in zip(q, r, vertices)))
    int_q, int_r = q_sum / 6, r_sum / 6
    int_qq = (q_sum * q_sum + sum(a * a for a in q)) / 24
    int_qr = (q_sum * r_sum + qr_sum) / 24
    int_Pq = _scale(F(1, 24), _add(_scale(q_sum, vertex_sum), Pq_sum))
    int_Pr = _scale(F(1, 24), _add(_scale(r_sum, vertex_sum), Pr_sum))
    int_Pqr = _scale(F(1, 120), _add(_scale(q_sum * r_sum + qr_sum, vertex_sum),
                                      _scale(q_sum, Pr_sum), _scale(r_sum, Pq_sum),
                                      _scale(2, Pqr_sum)))
    return _scale(_bracket(first, second) / 288, _add(
        _scale(6 * c * int_q - 3 * int_qq, A),
        _scale(2 * c * int_r - int_qr, total),
        _scale(4, int_Pqr), _scale(-2 * c, int_Pr), _scale(-9 * c, int_Pq)))


@lru_cache(maxsize=256)
def _diamond_polynomial(word, max_terms):
    """Integrate the rational spider chart, retaining all six source coordinates."""
    names = ("p0", "p1", "q0", "q1", "a0", "a1")
    context = fmpq_mpoly_ctx.get(names + ("s", "t", "r", "u"))
    p0, p1, q0, q1, a0, a1, s, t, r, u = context.gens()
    P, Q, A = (p0, p1), (q0, q1), (a0, a1)
    total = _add(A, _scale(s, P), _scale(t, Q))
    currents = (_add(_scale(s, P), _scale(-r, total)),
                _add(_scale(1 - s, P), _scale(r, total)),
                _add(_scale(t, Q), _scale(-u, total)),
                _add(_scale(1 - t, Q), _scale(u, total)),
                _scale(1 - r - u, total))
    # The remaining simplex factor (1-r-u) is integrated in the weight.
    polynomial = -_bracket(P, total) * _bracket(Q, total)
    if len(polynomial) > max_terms:
        raise RuntimeError("the exact diamond polynomial exceeded max_terms")
    for edge, component, power in word:
        polynomial *= currents[edge][component] ** power
        if len(polynomial) > max_terms:
            raise RuntimeError("the exact diamond polynomial exceeded max_terms")
    output = {}
    for exponent, coefficient in polynomial.terms():
        first, second, third, fourth = map(int, exponent[6:])
        weight = fmpq(factorial(third) * factorial(fourth),
                      (first + 1) * (second + 1) * factorial(third + fourth + 3))
        key = exponent[:6]
        output[key] = output.get(key, fmpq(0)) + coefficient * weight
    result = fmpq_mpoly_ctx.get(names).from_dict(output)
    if len(result) > max_terms:
        raise RuntimeError("the exact diamond source polynomial exceeded max_terms")
    return result


def boundary_value(sources, powers=None, *, max_terms=1000000):
    """Complete rational proper-face composition for any finite original word."""
    sources, word = _sources(sources), _word(powers)
    max_terms = _integer(max_terms, "max_terms")
    if not max_terms:
        raise ValueError("max_terms must be positive")
    proper_faces()
    inner_word = tuple((INNER_ORIGINAL.index(edge), component, power)
                       for edge, component, power in word if edge in INNER_ORIGINAL)
    polynomial = _diamond_polynomial(inner_word, max_terms)
    context = fmpq_mpoly_ctx.get(("x", "z"))
    x, z = context.gens()
    first, second = tuple(map(_q, sources[0])), tuple(map(_q, sources[5]))
    outer = (_add(_scale(1 - x, first), _scale(z, second)),
             _add(_scale(x, first), _scale(-z, second)),
             _add(_scale(-x, first), _scale(z - 1, second)))
    inner_sources = (_add(tuple(map(_q, sources[1])), outer[0]),
                     _add(tuple(map(_q, sources[2])), _scale(-1, outer[2])),
                     tuple(context.constant(_q(v)) for v in sources[3]))
    mapped = polynomial.compose(*(value for row in inner_sources for value in row))
    if len(mapped) > max_terms:
        raise RuntimeError("the exact triangle composition exceeded max_terms")
    for edge, component, power in word:
        if edge in OUTER_ORIGINAL:
            mapped *= outer[OUTER_ORIGINAL.index(edge)][component] ** power
            if len(mapped) > max_terms:
                raise RuntimeError("the exact triangle composition exceeded max_terms")
    value = F(0)
    for exponent, coefficient in mapped.terms():
        a, b = map(int, exponent)
        value += _fraction(coefficient) * F(factorial(a) * factorial(b), factorial(a + b + 2))
    return _bracket(sources[0], sources[5]) * value


class WardThree:
    """Reusable exact three-loop Ward calculation, retaining all source rows."""
    def __init__(self, sources, *, check=1, max_terms=1000000):
        self.sources = _sources(sources)
        self.max_terms = _integer(max_terms, "max_terms")
        self.check_count = _integer(check, "check")
        if not self.max_terms:
            raise ValueError("max_terms must be positive")
        self.quotients = []
        self._results = {}
        context = fmpq_mpoly_ctx.get(tuple("y%d" % i for i in range(14)))
        variables = context.gens()
        for position, (left, right) in enumerate(EDGES):
            vertices = tuple(v for v in range(6) if v != right)
            mapping = {v: vertices.index(left if v == right else v) for v in range(6)}
            rest = tuple(i for i in range(8) if i != position)
            edges = tuple((mapping[a], mapping[b]) for i, (a, b) in enumerate(EDGES) if i != position)
            rows = [[F(0), F(0)] for _ in vertices]
            for vertex, row in enumerate(self.sources):
                for component in range(2):
                    rows[mapping[vertex]][component] += row[component]
            lifted = []
            for edge, (a, b) in enumerate(EDGES):
                if edge != position:
                    index = rest.index(edge)
                    lifted.append(variables[2 * index:2 * index + 2])
                else:
                    lifted.append(tuple(context.constant(_q(self.sources[left][component]))
                                        - sum((int(a0 == left) - int(b0 == left)) * variables[2 * index + component]
                                              for index, (a0, b0) in enumerate(EDGES[j] for j in rest))
                                        for component in range(2)))
            engine = CoreCurrentPeriods(Graph(5, edges), rows, check=self.check_count, max_terms=self.max_terms)
            self.quotients.append(dict(position=position, contracted_edge=EDGES[position], sign=(-1) ** position,
                                       edges=edges, sources=tuple(map(tuple, rows)), lifted=tuple(lifted),
                                       vertex_map=tuple(mapping[v] for v in range(6)),
                                       marked=tuple(mapping[v] for v in MARKED),
                                       context=context, engine=engine))

    def close(self):
        for quotient in self.quotients:
            quotient["engine"].close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, traceback):
        self.close()

    def contracted_terms(self, powers=None):
        word = _word(powers)
        if word not in self._results:
            terms = []
            for quotient in self.quotients:
                polynomial = quotient["context"].constant(1)
                for edge, component, power in word:
                    polynomial *= quotient["lifted"][edge][component] ** power
                    if len(polynomial) > self.max_terms:
                        raise RuntimeError("a contracted current polynomial exceeded max_terms")
                values = {tuple(map(int, exponent)): _fraction(coefficient)
                          for exponent, coefficient in polynomial.terms() if coefficient}
                value = coefficients(quotient["engine"].integrate(values, allow_series=False))
                terms.append({key: quotient[key] for key in ("position", "contracted_edge", "sign", "edges", "sources",
                                                             "vertex_map", "marked")}
                             | {"coefficients": value,
                                "signed_coefficients": {key: quotient["sign"] * v for key, v in value.items()}})
            self._results[word] = tuple(terms)
        return self._results[word]

    def boundary_value(self, powers=None):
        return boundary_value(self.sources, powers, max_terms=self.max_terms)

    def verify(self, powers=None):
        terms = self.contracted_terms(powers)
        left = _sum_periods(row["signed_coefficients"] for row in terms)
        right = self.boundary_value(powers)
        if left != ({"1": right} if right else {}):
            raise ArithmeticError("the complete edge-contraction sum differs from its rational spider boundary")
        return dict(powers=_word(powers), terms=terms, signed_sum=left, rational_boundary=right)


def contracted_terms(sources, powers=None, **options):
    with WardThree(sources, **options) as engine:
        return engine.contracted_terms(powers)


def certificate():
    records = []
    for label, sources in (("generic", GENERIC_SOURCES), ("physical", PHYSICAL_SOURCES)):
        with WardThree(sources) as engine:
            bare = engine.verify()
            first = []
            for component in range(2):
                result = engine.verify({(2, component): 1})
                first.append(result)
                print(label, "physical A-spoke coordinate", component, flush=True)
                for term in result["terms"]:
                    if term["coefficients"]:
                        print(" ", term["position"], term["sign"], term["coefficients"], flush=True)
                print("  complete rational boundary", result["rational_boundary"], flush=True)
            spinor = boundary_first_spinor(sources)
            if spinor != tuple(row["rational_boundary"] for row in first):
                raise ArithmeticError("the universal bracket formula disagrees with direct rational face integration")
            records.append(dict(name=label, full_balanced_sources=sources,
                                bare_identity=bare, first_A_spoke_spinor=first,
                                universal_bracket_boundary=spinor))
    files = (Path(__file__).resolve(), ROOT / "code/htk/forest_periods.py",
             ROOT / "code/htk/_forest_period_worker.py")
    return dict(ordered_sliding_edges=EDGES, marked_vertices=MARKED, original_A_spoke=2,
                proper_faces=proper_faces(),
                identity="sum_e (-1)^e K_(G/e)[complete lifted current word] = rational diamond-triangle face",
                mathematical_scope="graph-level Ward identity; no Wick/color assembly is claimed",
                cancellation="complete exact coefficient sums; no independence of zeta values is assumed",
                implementation_sha256={str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                                       for path in files},
                source_cases=records)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="write the complete exact Ward certificate")
    arguments = parser.parse_args()
    result = json.loads(json.dumps(certificate(), default=str))
    path = HERE / "ward_three_certificate.json"
    if arguments.write:
        path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        print("Wrote the complete exact graph-level Ward certificate.", flush=True)
    else:
        if result != json.loads(path.read_text(encoding="utf-8")):
            raise ArithmeticError("the stored exact Ward certificate differs")
        print("Every exact Ward certificate field replays.", flush=True)


if __name__ == "__main__":
    main()
