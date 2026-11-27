"""Exact local obstruction to polynomial interpolation for a sourced star.

The signed SingleMove chart contains a logarithmic interval integral even
though the complete original two-loop kernel is the rational number -26.
All sources and both current coordinates are retained in the checks.
"""
from fractions import Fraction
from pathlib import Path
import sys

import sympy as sp

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "code"))

from htk import brute
from htk.forest_ears import ForestEarMoments
from htk.forest_periods import CoreCurrentPeriods
from htk.graph import Graph
from htk.star import SingleMove


EDGES = ((0, 1), (0, 2), (0, 3), (1, 3), (2, 3))
SOURCES = ((1, 2), (3, -1), (-2, 4), (-2, -5))
MARKED = (2, 3)


def check():
    graph = Graph(4, EDGES)
    assert all(sum(row[coordinate] for row in SOURCES) == 0 for coordinate in range(2))
    move = SingleMove(graph, 0, 1, 2, 3)
    retained_times = tuple(map(Fraction, (2, 3, 5)))
    first, second = sp.symbols("r1 r2")
    denominator = 2 * first * second - 7 * first - 5 * second + 10
    expected = (first * second
                * (2 * first * second - 11 * first - 9 * second + 18)
                * (20 * first * second - 119 * first - 71 * second + 170)
                / (500 * denominator**2))

    # Substitute the complete original numerator at once, using the exact
    # U = 150/(r1*r2*denominator) and signed star Jacobian.  This small
    # polynomial identity does not use the claimed closed density formula.
    original = CoreCurrentPeriods(graph, SOURCES, marked=MARKED, check=2)
    try:
        original._source_numerator()
        terms = []
        for powers, coefficient in original._bare.terms():
            a, b, c, d, e = powers
            assert a + b + c + d + e == 3
            assert max(powers) <= 1
            term = sp.Rational(int(coefficient.p), int(coefficient.q))
            term *= 6**a * 10**b * 15**c * 3**d * 5**e
            term *= ((1 - second)**a * (1 - first)**b
                     * first**(2 - d) * second**(2 - e)
                     * denominator**(3 - a - b - c))
            terms.append(term)
        transformed = -sp.Rational(450, 150**4) * sp.Add(*terms) / denominator**2
        assert sp.cancel(transformed - expected) == 0
        complete_original = original.integrate({(0,) * (2 * graph.ne): Fraction(1)})
        assert complete_original == -26
    finally:
        original.close()

    # Independent direct current/density comparisons at three positive
    # rational chart points, preserving all original edge orientations.
    points = ((Fraction(2, 5), Fraction(3, 7)),
              (Fraction(4, 9), Fraction(5, 8)),
              (Fraction(1, 3), Fraction(1, 3)))
    for ratios in points:
        times = move.tmap(ratios, retained_times)
        assert move.currents(ratios, retained_times, SOURCES) == brute.currents(graph, times, SOURCES)
        jacobian = move.jacobian(ratios, retained_times)
        assert jacobian < 0
        direct = brute.density(graph, times, SOURCES).rational() * jacobian
        exact = expected.subs({first: sp.Rational(ratios[0].numerator, ratios[0].denominator),
                               second: sp.Rational(ratios[1].numerator, ratios[1].denominator)})
        assert direct == Fraction(int(sp.numer(exact)), int(sp.denom(exact)))

    with ForestEarMoments(graph, SOURCES, marked=MARKED, check=2) as ears:
        assert ears.moment() == complete_original

    # A real primitive on 0 <= r1 <= 1, where 25 - 19*r1 is positive.
    line = sp.factor(expected.subs(second, sp.Rational(1, 3)))
    primitive = (sp.Rational(10447, 1083000) * first**2
                 - sp.Rational(6089, 2572125) * first
                 - sp.Rational(30781, 9774075) * sp.log(25 - 19 * first)
                 + sp.Rational(112, 130321) / (19 * first - 25))
    assert sp.cancel(sp.diff(primitive, first) - line) == 0
    integral = primitive.subs(first, 1) - primitive.subs(first, 0)
    exact_integral = (sp.Rational(147541, 20577000)
                      + sp.Rational(30781, 9774075) * (sp.log(25) - sp.log(6)))
    assert sp.simplify(integral - exact_integral) == 0
    return {
        "signed_density": str(sp.factor(expected)),
        "fixed_r2_interval_integral": str(exact_integral),
        "complete_original_kernel": str(complete_original),
        "all_original_current_coordinates_checked": 2 * graph.ne * len(points),
    }


if __name__ == "__main__":
    for label, value in check().items():
        print(label + ": " + str(value))
    print("The exact sourced-star local obstruction and complete rational control pass.")
