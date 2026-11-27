"""Complete rational outer source faces, including the affine chord correction.

The six current coordinates are (J_WA[0], J_WA[1], J_WB[0], J_WB[1],
J_AB[0], J_AB[1]).  A polynomial is a dictionary from full six-entry
nonnegative exponent tuples to exact rational coefficients.  Source input is
(lambda_W, lambda_c, lambda_A), each a full two-entry exact spinor.

The integration orientation is the positive coordinate order
(R_TW, R_TA, R_TB, R_WA, R_WB, source_path_fraction), with R_AB=1.
No ordinary subdivided-graph kernel is substituted for the complete form.
"""

from fractions import Fraction
from math import factorial
from pathlib import Path
import sys

from flint import fmpq, fmpq_mpoly_ctx

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
if str(ROOT / "code") not in sys.path:
    sys.path.insert(0, str(ROOT / "code"))

from htk.brute import currents as brute_currents
from htk.graph import Graph
from htk.linalg import det
from htk.star import Dual


CORE_EDGES = ((1, 0), (1, 2), (1, 3), (0, 2), (0, 3), (2, 3))
CORE = Graph(4, CORE_EDGES)
VARIABLES = ("s1", "s2", "w1", "w2", "rX", "rQ")
CHECK_POINT = (Fraction(2, 5), Fraction(3, 7), Fraction(1, 5),
               Fraction(1, 4), Fraction(2, 7), Fraction(3, 8))


def _fraction(value):
    if isinstance(value, (float, complex, bool)):
        raise TypeError("exact rational inputs are required")
    return Fraction(str(value)) if isinstance(value, fmpq) else Fraction(value)


def _q(value):
    value = _fraction(value)
    return fmpq(value.numerator, value.denominator)


def _validated_word(polynomial):
    if polynomial is None:
        return {(0,) * 6: Fraction(1)}
    out = {}
    for powers, coefficient in polynomial.items():
        if len(powers) != 6 or any(type(power) is not int or power < 0 for power in powers):
            raise ValueError("each current exponent must be a six-tuple of nonnegative integers")
        coefficient = _fraction(coefficient)
        if coefficient:
            out[tuple(powers)] = out.get(tuple(powers), Fraction(0)) + coefficient
    return {powers: coefficient for powers, coefficient in out.items() if coefficient}


def _poly_det(matrix, context, max_terms):
    """Fraction-free elimination, checking each exact polynomial division."""
    rows = [list(row) for row in matrix]
    size = len(rows)
    previous = context.constant(1)
    sign = 1
    for column in range(size - 1):
        pivot_row = next((row for row in range(column, size) if rows[row][column]), None)
        if pivot_row is None:
            return context.constant(0)
        if pivot_row != column:
            rows[column], rows[pivot_row] = rows[pivot_row], rows[column]
            sign = -sign
        pivot = rows[column][column]
        for row in range(column + 1, size):
            for other in range(column + 1, size):
                numerator = pivot * rows[row][other] - rows[row][column] * rows[column][other]
                quotient, remainder = divmod(numerator, previous)
                if remainder:
                    raise ArithmeticError("a fraction-free determinant division was not exact")
                if len(quotient) > max_terms:
                    raise RuntimeError("the exact outer-face determinant exceeded max_terms")
                rows[row][other] = quotient
            rows[row][column] = context.constant(0)
        previous = pivot
    return sign * rows[-1][-1]


def _bounded_power(value, exponent, context, max_terms):
    result = context.constant(1)
    while exponent:
        if exponent & 1:
            result *= value
            if len(result) > max_terms:
                raise RuntimeError("the exact outer-face insertion exceeded max_terms")
        exponent >>= 1
        if exponent:
            value *= value
            if len(value) > max_terms:
                raise RuntimeError("the exact outer-face insertion exceeded max_terms")
    return result


def _time_map(parameters, path):
    s1, s2, w1, w2, rX, rQ = parameters
    w0 = 1 - w1 - w2
    reduced_uw = 1 / w0
    reduced_ux = rX / w2
    reduced_wx = (1 - rX) / w2
    delta_uw = reduced_uw
    delta_ux = reduced_ux / (1 - s1)
    delta_wx = reduced_wx / (1 - s2)
    total = delta_uw + delta_ux + delta_wx
    star_u = delta_uw * delta_ux / total
    star_w = delta_uw * delta_wx / total
    star_x = delta_ux * delta_wx / total
    outer_uw = 1 / w1
    outer_ux = reduced_ux / s1
    outer_wx = reduced_wx / s2
    if path == "WA":
        return (star_u, star_w, star_x, outer_uw, outer_ux, outer_wx), rQ
    if path == "BA":
        return (star_x, star_w, star_u, outer_wx, outer_ux, outer_uw), rQ
    raise ValueError("outer path must be WA or BA")


class SourceTriangleOuterFaces:
    """Exact corrected outer face for either source-simplex identity.

    ``integrate(polynomial)`` returns a Fraction.  No maximum jet degree is
    imposed; a transparent polynomial term budget bounds the computation.
    """

    def __init__(self, vectors, triangle="TWA", max_terms=1_000_000):
        if triangle not in ("TWA", "TBA"):
            raise ValueError("triangle must be TWA or TBA")
        if type(max_terms) is not int or max_terms <= 0:
            raise ValueError("max_terms must be a positive integer")
        vectors = tuple(tuple(_fraction(value) for value in row) for row in vectors)
        if len(vectors) != 3 or any(len(row) != 2 for row in vectors):
            raise ValueError("three complete two-entry source spinors are required")
        self.vectors = vectors
        self.triangle = triangle
        self.path = "WA" if triangle == "TWA" else "BA"
        self.max_terms = max_terms
        self.context = fmpq_mpoly_ctx.get(VARIABLES, ordering="lex")
        s1, s2, w1, w2, rX, rQ = self.context.gens()
        P, Q, R = vectors
        B = tuple(-P[a] - Q[a] - R[a] for a in range(2))
        source_u, source_x = (P, B) if self.path == "WA" else (B, P)
        edge_uw, edge_ux, edge_wx = [], [], []
        half_edge = []
        for a in range(2):
            total = _q(source_u[a]) + (1-rQ) * _q(Q[a]) + (1-rX) * _q(source_x[a])
            edge_uw.append(w1 * total)
            edge_ux.append(s1 * (w2 * total - (1-rX) * _q(source_x[a])))
            edge_wx.append(-s2 * (w2 * total + rX * _q(source_x[a])))
            half_edge.append(w1 * total - (1-rQ) * _q(Q[a]))
        rows = (edge_uw, edge_ux, edge_wx) if self.path == "WA" else (
            tuple(-value for value in edge_wx), tuple(-value for value in edge_ux),
            tuple(-value for value in edge_uw))
        self.currents = tuple(value for row in rows for value in row)
        self.source_half_current = tuple(half_edge)
        self.effective_source_chord = tuple(edge_uw)
        self.jacobian = _poly_det([[current.derivative(j) for j in range(6)]
                                   for current in self.currents], self.context, max_terms)
        parameters = tuple(Dual.var(value, j, 6) for j, value in enumerate(CHECK_POINT))
        times, ratio = _time_map(parameters, self.path)
        gauge_coordinates = tuple(value / times[-1] for value in times[:-1]) + (ratio,)
        self.coordinate_jacobian_at_check = det([value.g for value in gauge_coordinates])
        if not self.coordinate_jacobian_at_check:
            raise ArithmeticError("the positive outer-face coordinate map is singular")
        self.orientation = 1 if self.coordinate_jacobian_at_check > 0 else -1
        self.stats = {"variables": 6, "jacobian_terms": len(self.jacobian),
                      "orientation": self.orientation, "max_terms": max_terms}

    def pointwise_check(self, parameters=CHECK_POINT):
        """Independent grounded-Laplacian check of all common chord currents."""
        parameters = tuple(_fraction(value) for value in parameters)
        times, ratio = _time_map(parameters, self.path)
        P, Q, R = self.vectors
        B = tuple(-P[a] - Q[a] - R[a] for a in range(2))
        if self.path == "WA":
            sources = (tuple(P[a] + (1-ratio)*Q[a] for a in range(2)), (0, 0),
                       tuple(R[a] + ratio*Q[a] for a in range(2)), B)
        else:
            sources = (P, (0, 0), tuple(R[a] + ratio*Q[a] for a in range(2)),
                       tuple(B[a] + (1-ratio)*Q[a] for a in range(2)))
        exact = brute_currents(CORE, times, sources)
        expected = tuple(exact[e][a] for e in (3, 4, 5) for a in range(2))
        got = tuple(_fraction(current(*map(_q, parameters))) for current in self.currents)
        if got != expected:
            raise ArithmeticError("a corrected common chord differs from the original core solve")
        for a in range(2):
            difference = self.effective_source_chord[a] - self.source_half_current[a]
            wanted = (1-self.context.gen(5)) * _q(Q[a])
            if difference != wanted or difference.derivative(5) != -_q(Q[a]):
                raise ArithmeticError("the complete affine chord correction was not retained")
        return {"core_times": tuple(str(value) for value in times),
                "core_sources": tuple(tuple(str(value) for value in row) for row in sources),
                "common_currents": tuple(str(value) for value in got),
                "affine_chord_and_differential_checked": True}

    def integrate(self, polynomial=None):
        polynomial = _validated_word(polynomial)
        substitution = self.context.constant(0)
        powers_cache = {}
        for exponents, coefficient in polynomial.items():
            term = self.context.constant(_q(coefficient))
            for edge_component, exponent in enumerate(exponents):
                if exponent:
                    key = edge_component, exponent
                    if key not in powers_cache:
                        powers_cache[key] = _bounded_power(self.currents[edge_component], exponent,
                                                          self.context, self.max_terms)
                    term *= powers_cache[key]
                    if len(term) > self.max_terms:
                        raise RuntimeError("the exact outer-face insertion exceeded max_terms")
            substitution += term
            if len(substitution) > self.max_terms:
                raise RuntimeError("the exact outer-face insertion exceeded max_terms")
        integrand = self.jacobian * substitution
        if len(integrand) > self.max_terms:
            raise RuntimeError("the exact outer-face integrand exceeded max_terms")
        answer = Fraction(0)
        for exponents, coefficient in integrand.to_dict().items():
            a, b, c, d, e, f = map(int, exponents)
            moment = Fraction(factorial(c) * factorial(d), factorial(c+d+2))
            moment /= (a+1) * (b+1) * (e+1) * (f+1)
            answer += _fraction(coefficient) * moment
        self.stats.update(insertion_terms=len(substitution), integrand_terms=len(integrand))
        return self.orientation * answer


__all__ = ["SourceTriangleOuterFaces", "CORE_EDGES", "VARIABLES", "CHECK_POINT"]
