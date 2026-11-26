"""Integrate two core lengths of an automorphism-averaged four-loop placement."""
from argparse import ArgumentParser
from fractions import Fraction as F
from functools import cache
from itertools import combinations, permutations
from math import comb, factorial, lcm, prod
from pathlib import Path
import gzip
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'nonSP'))
from graphs import PRISM


def edge_pair_orbit(edges, first, second):
    index = {tuple(sorted(edge)): i for i, edge in enumerate(edges)}
    orbit = set()
    for perm in permutations(range(6)):
        image = [index.get(tuple(sorted((perm[i], perm[j])))) for i, j in edges]
        if None not in image:
            orbit.add((image[first], image[second]))
    return sorted(orbit)


def kirchhoff(edges, first, second):
    terms = [0] * 4
    for tree in combinations(range(9), 5):
        vertices = {0}
        for _ in range(5):
            for e in tree:
                i, j = edges[e]
                if i in vertices or j in vertices:
                    vertices.update((i, j))
        if len(vertices) == 6:
            terms[int(first not in tree) + 2 * int(second not in tree)] += 1
    return tuple(terms)


def numerator(source, orbit):
    bins = [0] * 144
    scale = 1
    facts = [factorial(i) for i in range(13)]
    with gzip.open(source, 'rt') as stream:
        for line in stream:
            key, coefficient = line.split()
            key = int(key)
            powers = [((key >> shift) & 15) - 1 for shift in range(32, -1, -4)]
            assert sum(powers) == 43 and 0 <= min(powers) <= max(powers) < 12
            value = F(coefficient) * F(facts[12], prod(facts[i] for i in powers))
            new = lcm(scale, value.denominator)
            if new != scale:
                bins = [v * (new // scale) for v in bins]
                scale = new
            weight = value.numerator * (scale // value.denominator)
            for i, j in orbit:
                bins[12 * powers[i] + powers[j]] += weight
    return {(i // 12, i % 12): F(v, scale * len(orbit)) for i, v in enumerate(bins) if v}


def pole_coefficients(j, u, v, p, q):
    # Taylor coefficients of y^j/(y+v)^q at y=-u, before division by (y+u)^p.
    return [sum(comb(j, r) * (-u)**(j-r) * (-1)**(m-r)
                * comb(q+m-r-1, m-r) / (v-u)**(q+m-r)
                for r in range(min(j, m)+1)) for m in range(p)]


@cache
def moment(i, j, kirchhoff_coefficients):
    a, b, c, d = map(F, kirchhoff_coefficients)
    u, v = a/c, b/d
    p, q = 12-i, i+1
    scale = F(factorial(i) * factorial(11-i), factorial(12)) / c**p / d**q
    if u == v:
        return scale * u**(j-12) * F(factorial(j)*factorial(11-j), factorial(12)), F(0)
    left = pole_coefficients(j, u, v, p, q)
    right = pole_coefficients(j, v, u, q, p)
    assert left[-1] + right[-1] == 0
    rational = sum(row[-k] * root**(1-k) / (k-1)
                   for row, root in ((left, u), (right, v)) for k in range(2, len(row)+1))
    return scale * rational, scale * left[-1]


def sign_bound(rational, log_weight, argument, terms):
    # 2 atanh((x-1)/(x+1)); bound the omitted tail by a geometric series.
    z = (argument-1)/(argument+1)
    center = 2*sum(z**(2*k+1)/(2*k+1) for k in range(terms))
    error = 2*abs(z)**(2*terms+1)/((2*terms+1)*(1-z*z))
    center = rational + log_weight*center
    error *= abs(log_weight)
    return center-error, center+error


def evaluate(source, edges, pair):
    orbit = edge_pair_orbit(edges, *pair)
    coeffs = kirchhoff(edges, *pair)
    assert all(kirchhoff(edges, i, j) == coeffs for i, j in orbit)
    a, b, c, d = coeffs
    poly = numerator(source, orbit)
    unit = sum(poly.values()) / sum(coeffs)**13
    rational = log_weight = F(0)
    for (i, j), value in poly.items():
        r, l = moment(i, j, coeffs)
        r1, l1 = moment(j, i, (a, c, b, d))
        assert (r, l) == (r1, l1)
        rational += value * r
        log_weight += value * l
    argument = F(b*c, a*d)
    return unit, (rational, log_weight, argument), sign_bound(rational, log_weight, argument, 120)


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('bank', type=Path)
    parser.add_argument('core', choices=('prism', 'bipartite'))
    parser.add_argument('edges', type=int, nargs=2, choices=range(9))
    args = parser.parse_args()
    assert args.edges[0] != args.edges[1]
    edges = PRISM if args.core == 'prism' else [(i, 3+j) for i in range(3) for j in range(3)]
    unit, value, bounds = evaluate(args.bank, edges, args.edges)
    print('Unit density', unit)
    print('Rational coefficient, logarithmic coefficient, log argument', *value)
    print('Bounds in units of 10^-20', (bounds[0]*10**20)//1, -((-bounds[1]*10**20)//1))
