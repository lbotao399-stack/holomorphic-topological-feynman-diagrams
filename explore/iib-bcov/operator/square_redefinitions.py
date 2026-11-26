"""Marked scalar-square redefinitions at two and three loops on the repeated waves."""
from collections import Counter
from fractions import Fraction
from itertools import combinations
from flint import fmpz_mat


def covering_minors(momenta, degree):
    """Ordered products of squared five-row minors involving every scalar slot."""
    slots = range(len(momenta))
    minors = []
    for indices in combinations(slots, 5):
        value = int(fmpz_mat([momenta[i] for i in indices]).det())
        if value:
            minors.append((sum(1 << i for i in indices), value*value))
    weights = Counter({0: 1})
    for _ in range(degree):
        expanded = Counter()
        for covered, weight in weights.items():
            for mask, square in minors:
                expanded[covered | mask] += weight*square
        weights = expanded
    return weights[(1 << len(momenta)) - 1]


def tree_images(n):
    axes = (0,)*n + (1,)*n + (2,3)
    basis = tuple(tuple(int(i == j) for j in range(5)) for i in range(5))
    changes = [0]*4
    for wave, axis in enumerate(axes):
        momenta = [basis[i] for i in axes]
        momenta[wave] = tuple(x+y for x,y in zip(momenta[wave],basis[4]))
        if axis in (2,3):
            assert fmpz_mat(momenta).rank() == 4
        changes[axis] += covering_minors(momenta, n)
    return tuple(changes)


def hexagon_mark():
    # Of the five active waves, four must be scalars and the fifth the vector.
    scalars = (0,0,0,1,1,1,2,3)
    contributing = 0
    for active in combinations(range(8), 4):
        if len({scalars[i] for i in active}) < 4:
            continue
        spectator = {scalars[i] for i in range(8) if i not in active}
        assert spectator == {0,1} and len(spectator) + 2 < 5
        contributing += 1
    assert contributing == 9


def quadratic_source_check():
    """Compare the polarized and Taylor-divided planar Hessians to the full two-loop cell sum."""
    from itertools import product
    from math import factorial
    import sys
    from pathlib import Path
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
    from longitudinal_theta import primary_coefficient, scalar_square
    from planar_bulk import mass
    axes = [tuple(int(i==j) for i in range(5)) for j in range(5)]
    rank_one = primary_coefficient([axes[i] for i in (0,0,1,2,3)],axes[4],(0,1,1,1,0))
    assert rank_one/2 == Fraction(1,20736)
    polarization = (0,0,1,1,0)
    def evaluate(plane):
        waves = [(*p,0,0,0) for p in plane] + axes[2:4]
        return scalar_square(waves,axes[4],polarization)
    for plane in [((1,0),(2,1),(3,-1),(1,-2)), ((1,1),)*4, ((1,0),(0,1))*2]:
        assert evaluate(plane)==Fraction(mass(plane),4860)
    # Independent Taylor expansion with a non-diagonal, indefinite Hessian.
    directions = [(1,0),(1,1),(2,-1)]
    strengths = (2,3,-1)
    h00,h01,h11 = [sum(c*p[i]*p[j] for c,p in zip(strengths,directions))
                   for i,j in ((0,0),(0,1),(1,1))]
    from math import prod
    from functools import cache
    integral = cache(evaluate)
    result = sum(prod(strengths[i] for i in ids)*integral(tuple(directions[i] for i in sorted(ids)))
                 for ids in product(range(3),repeat=4))/factorial(4)
    assert result == Fraction((h00*h11-h01*h01)**2,4860)


if __name__ == '__main__':
    assert tree_images(1) == (0,0,0,0)
    assert tree_images(2) == (4,4,0,0)
    assert tree_images(3) == (108,108,0,0)
    hexagon_mark()
    quadratic_source_check()
    for n in (1,2,3):
        print('Order', n, 'point counterterm components', tree_images(n))
    print('Two-loop protected component:', Fraction(1,405)/3)
    print('Three-loop protected rational and zeta(3) coefficients:',
          Fraction(411065203,326918592000)/4, -Fraction(3617,3503500)/4)
