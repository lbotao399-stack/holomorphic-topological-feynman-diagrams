"""Independent complete 9-, 13-, and 19-loop sourced-comb controls.

Run the analytic transfer followed by the existing exact SP integrator:
    python3 explore/circuit-general/sourced-forests/comb_reference.py

Also compare against the new public product-forest implementation:
    python3 explore/circuit-general/sourced-forests/comb_reference.py --compare-public

The reference calculation never calls ForestContractions, CoreCurrentPeriods,
or ProductForestMoments. The optional public comparison is a separate step.

For k nested passive paired moves, the deepest A-spoke insertion becomes

    2*(-1)^k L^a <L,Y>^k <H,Y>^k
    / (k! * ((k+1)!)^2 * (k+2)! * (2*k+1)!).

Indeed the next paired map substitutes
    (r1 L^a-r3 Y^a) * r1^k*r2^k*s^(2*k),  s=1-r3-r4,
and supplies -(r1-r2)*s*<L,Y>*<H,Y>. The Y^a term vanishes
by interchange of r1,r2. The surviving square and triangle moments are

    1/((k+1)*(k+2)^2*(k+3)),  1/((2*k+2)*(2*k+3)).

Their product gives the stated factorial coefficient by induction. The
remaining graph, all five nonzero source rows, and all edge conventions
are explicitly specified below, independently of the contraction planner.
"""
import argparse
from fractions import Fraction
from math import factorial
from numbers import Integral
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "code"))

from htk.graph import Graph
from htk.spred import plan_generic, sp_reduce


CORE_EDGES = ((0, 3), (0, 4), (3, 4), (0, 1), (0, 2), (1, 3), (2, 4))
CORE_SOURCES = ((1, 2), (3, -1), (-2, 4), (5, 3), (-7, -8))
EXPECTED = {
    3: Fraction(-11549895079003417, 586642839552000000),
    5: Fraction(-653729304337869467963209, 3257778542909926806650880000000),
    8: Fraction(241018376797614391720336915629690895727,
                23196061207306780193313894655375109127143424000000000000),
}


def reference(pairs):
    """Exact analytic comb transfer, then the independently sourced SP density."""
    if isinstance(pairs, bool) or not isinstance(pairs, Integral) or pairs < 0:
        raise ValueError("the number of nested paired moves must be a nonnegative integer")
    pairs = int(pairs)
    graph = Graph(5, CORE_EDGES)
    reduction = sp_reduce(graph, D=2, plan=plan_generic(graph))
    density = reduction.density(CORE_SOURCES, T=0, method="generic", check=2)
    left, right, marked = reduction.currents(CORE_SOURCES)[:3]
    local = (left[0], left[0] * marked[1] - left[1] * marked[0],
             right[0] * marked[1] - right[1] * marked[0])
    scalar = Fraction(2 * (-1) ** pairs,
                      factorial(pairs) * factorial(pairs + 1) ** 2
                      * factorial(pairs + 2) * factorial(2 * pairs + 1))
    factors = []
    # Each factor is simplified before taking its repeated power. The old
    # SP block integrator then expands only factors incident to that block.
    for polynomial, degree in ([(factor, 1) for factor in density.factors]
                               + list(zip(local, (1, pairs, pairs)))):
        coefficient, irreducible = polynomial.factor_list()
        scalar *= Fraction(str(coefficient)) ** degree
        factors.extend(factor for factor, power in irreducible
                       for _ in range(int(power) * degree))
    integral = density._integrate_polynomial(factors)
    if any(any(exponent) for exponent in integral):
        raise ArithmeticError("the independent SP reference retained a ratio variable")
    return scalar * Fraction(str(integral.LC))


def compare_public(pairs, expected):
    """A separate comparison; the reference above does not call this route."""
    from htk.product_forests import ProductForestMoments
    from fixtures import higher_paired_comb
    graph, sources, marked, powers = higher_paired_comb(pairs)
    with ProductForestMoments(graph, sources, marked, check=1,
                              max_terms=1000000, allow_series=False) as engine:
        actual = engine.moment(powers)
    if actual != expected:
        raise ArithmeticError("public and independent analytic comb moments disagree")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compare-public", action="store_true")
    arguments = parser.parse_args()
    for pairs, expected in EXPECTED.items():
        actual = reference(pairs)
        if actual != expected:
            raise ArithmeticError("the complete independent comb reference changed")
        if arguments.compare_public:
            compare_public(pairs, actual)
        print("%d loops: %s" % (3 + 2 * pairs, actual), flush=True)
    print("All three complete exact comb moments agree.", flush=True)


if __name__ == "__main__":
    main()
