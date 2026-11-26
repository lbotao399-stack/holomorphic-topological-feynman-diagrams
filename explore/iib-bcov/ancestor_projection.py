"""Cotangent-line pullback and rational-tail fusions on elementary fivefold probes."""
from collections import Counter
from fractions import Fraction
from itertools import combinations
from math import comb, prod
from pathlib import Path
import sys

from elliptic_descendants import vertex
from longitudinal_theta import primary_coefficient

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'twisted-holography-N4' / 'code'))
from soft_quartic_two_betti import descendant_vertex, pixton_quartic


def fusions(retained, forgotten, level, count):
    """Fuse a subset into the first retained mark, keeping the other waves labeled."""
    for subset in combinations(range(len(forgotten)), count):
        merged = tuple(sum(row) for row in zip(retained[0], *(forgotten[i] for i in subset)))
        yield ([merged, *retained[1:], *(p for i, p in enumerate(forgotten) if i not in subset)],
               level - count)


def check_one_loop():
    # Each forgotten coordinate is saturated by one elliptic mark.
    for n in (2, 3, 4, 5):
        hard = [(i+1, i*i, 2*i-1, -3*i, 2) for i in range(n-1)]
        for strengths in ((1, 1, 1, 1), (2, -1, 3, -2)):
            forgotten = [tuple(a * int(i == j) for j in range(5))
                         for i, a in enumerate(strengths, 1)]
            balance = tuple(-sum(row) for row in zip(*hard, *forgotten))
            retained = [balance, *hard]
            level = n-1
            expected = prod(a*a for a in strengths) * vertex(
                [(p[0],) for p in retained], [level] + [0]*(n-1))
            result = sum((-1)**count * vertex(points, [k]+[0]*(len(points)-1))
                         for count in range(level+1)
                         for points, k in fusions(retained, forgotten, level, count))
            assert result == expected and expected > 0
    print('Eight elliptic ancestor projections agree with the full descendant partition sum.')


def check_two_loops():
    root = (1, 0, 0, 0, 0)
    vector = (-2, 1, 0, 0, 0)
    polarization = (1, 2, 0, 0, 0)
    retained = [(1, -1, 0, 0, 0), root, vector]
    forgotten = [tuple(sign * int(i == j) for j in range(5))
                 for i in range(1, 5) for sign in (1, -1)]
    charges = [p[0] for p in retained]
    assert sum(x*y for x,y in zip(vector, polarization)) == 0
    assert sum(x*y for x,y in zip(root, polarization)) == 1
    base = descendant_vertex(2, [(p,) for p in charges], (4, 0, 0))
    assert base == pixton_quartic(charges, (4, 0, 0)) == Fraction(91, 5760)
    ancestor = 2**4 * base
    counts = Counter()
    for points, level in fusions(retained, forgotten, 4, 4):
        assert level == 0 and not any(map(sum, zip(*points)))
        scalar = [points[0], *points[3:]]
        counts[-primary_coefficient(scalar, vector, polarization)] += 1
    assert counts == {0:42, Fraction(1,31104):12, Fraction(29,155520):8, Fraction(53,155520):8}
    primary = sum(weight * n for weight, n in counts.items())
    assert ancestor+primary == Fraction(10007, 38880)
    print('Labeled two-loop fusions by cardinality:', [comb(8, j) for j in range(5)])
    print('Fivefold ancestor:', ancestor, 'primary tails:', primary, 'descendant lower bound:', ancestor+primary)


if __name__ == '__main__':
    check_one_loop()
    check_two_loops()
