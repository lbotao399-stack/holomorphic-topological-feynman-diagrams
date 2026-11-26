"""Planar genus-three descendant zeta(3) coefficients and one flat fivefold scalar probe."""
from collections import Counter
from fractions import Fraction
from itertools import combinations
from math import factorial, prod

from mixed_sextic_vertex import partitions, sextic


def tetrahedron_numerator(charges):
    """485100 times the planar tetrahedron coefficient, over the supplied scalar ring."""
    if len(charges) != 4 or any(len(p) != 2 for p in charges) or any(map(sum, zip(*charges))):
        raise ValueError("Supply four conserved planar momenta.")
    pairs = {}
    for i, j in combinations(range(4), 2):
        x, y = charges[i]
        u, v = charges[j]
        pairs[i, j] = (x*v-y*u)**2
    stars = sum(prod(value for edge, value in pairs.items() if i in edge) for i in range(4))
    matchings = (pairs[0, 1]*pairs[2, 3] + pairs[0, 2]*pairs[1, 3] + pairs[0, 3]*pairs[1, 2])
    return stars + sum(pairs.values())*matchings


def descendant_zeta(charges, levels):
    """Sum the positive four-vertex contributions with the genus-zero cotangent weights."""
    if (len(charges) != len(levels) or any(len(p) != 2 for p in charges)
            or any(k < 0 for k in levels) or sum(levels) != len(levels)
            or any(map(sum, zip(*charges)))):
        raise ValueError("Supply conserved planar charges and complementary genus-three levels.")
    total = Fraction()
    for groups in partitions(len(levels)):
        if len(groups) != 4 or any(sum(levels[i] for i in g) != len(g) for g in groups):
            continue
        momenta = [tuple(sum(charges[i][j] for i in g) for j in range(2)) for g in groups]
        weight = Fraction(prod(factorial(len(g)) for g in groups), prod(map(factorial, levels)))
        total += weight*tetrahedron_numerator(momenta)/485100
    return total


def verify_polynomial():
    from flint import fmpz_mpoly_ctx
    from mixed_k4_interior import EDGES, polynomial
    x = fmpz_mpoly_ctx.get(tuple('x'+str(i) for i in range(6))).gens()
    momenta = (x[:2], x[2:4], x[4:6], tuple(-sum(x[k::2]) for k in (0, 1)))
    weights = [-sum(a*b for a, b in zip(momenta[i], momenta[j])) for i, j in EDGES]
    reference = 0
    for powers, coefficient in polynomial()[1].items():
        integer = 1940400*coefficient
        assert integer.denominator == 1
        reference += integer.numerator*prod(w**k for w, k in zip(weights, powers) if k)
    assert reference == 4*tetrahedron_numerator(momenta)
    print('All six independent planar components agree as an exact polynomial.')


def fivefold_probe(levels, labeled=False):
    """Forget the two zero planar primary marks by the string equation."""
    charges = ((1, 0), (-1, 0), (0, 1), (0, -1))
    representatives = ((1, 1, 1, 1), (0, 2, 1, 1), (1, 2, 0, 1), (2, 2, 0, 0), (0, 2, 2, 0))
    terms = Counter({levels: 1})
    for _ in range(2):
        reduced = Counter()
        for current, weight in terms.items():
            for i, k in enumerate(current):
                if k:
                    reduced[current[:i] + (k-1,) + current[i+1:]] += weight
        terms = reduced
    result = [Fraction(), Fraction()]
    for current, weight in terms.items():
        rational, zeta = sextic(3, charges, current)
        assert zeta == descendant_zeta(charges, current)
        if labeled and current in representatives:
            from mixed_sextic_vertex import labeled_sum
            assert labeled_sum(3, charges, current) == (rational, zeta)
        result[0] += weight*rational
        result[1] += weight*zeta
    assert result[1] == descendant_zeta(charges+((0, 0), (0, 0)), levels+(0, 0))
    return tuple(factorial(3)**3*x for x in result)


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--labeled', action='store_true',
                        help='Independently sum every marked assignment on the distinct symmetry representatives.')
    args = parser.parse_args()
    verify_polynomial()
    for levels, expected in (((2, 2, 1, 1), (Fraction(4093949, 11642400), Fraction(96, 13475))),
                             ((2, 2, 2, 0), (Fraction(11147, 63000), 0))):
        value = fivefold_probe(levels, args.labeled)
        assert value == expected
        print('Signed fivefold Q3 scalar measurement', levels, '(rational, zeta(3)):', *value)
