"""Exact current moments for the theta166 flow with c at N and gamma at B1."""
from fractions import Fraction
from functools import cache
from itertools import product
from math import comb, factorial, prod
from pathlib import Path
import sys

from sympy import ZZ
from sympy.polys.rings import ring

CURRENT_RING, *CURRENT_VARIABLES = ring('b1,b2,b3,b4,b5,c1,c2,c3,c4,c5', ZZ)
B = CURRENT_VARIABLES[:5]
C = CURRENT_VARIABLES[5:]
IDENTITY = tuple(range(5))
SWAPPED = (0, 1, 3, 2, 4)


def numerator():
    area = (B[2] - B[4]) * (C[3] - C[4]) - (B[3] - B[4]) * (C[2] - C[4])
    return (1 + B[0])**3 * (1 - C[0])**4 * (2 * B[0] + 1) * (B[0] + 2 * C[0] - 1) * area


@cache
def rising(start, degree):
    return prod((start + i for i in range(degree)), start=Fraction(1))


@cache
def weight_moment(a, u, v, d):
    """Expectation of (wB-1)^a (-wB)^u wC^v (1-wC)^d."""
    alpha = Fraction(1, 2)
    total = Fraction(13, 2)
    return sum(
        (-1)**(a - i + u + j) * comb(a, i) * comb(d, j)
        * rising(alpha, i + u) * rising(alpha, j + v) / rising(total, i + u + j + v)
        for i in range(a + 1) for j in range(d + 1)
    )


@cache
def path_moment(exponents):
    """Expectation of a monomial in 1 > t1 > ... > t5 > 0, density 5!."""
    tail = 0
    denominator = 1
    for j, exponent in enumerate(reversed(exponents), 1):
        tail += exponent
        denominator *= j + tail
    return Fraction(factorial(len(exponents)), denominator)


@cache
def choices(powers):
    return [(choice, sum(choice), prod(comb(n, k) for n, k in zip(powers, choice)))
            for choice in product(*(range(n + 1) for n in powers))]


def mean(polynomial, order):
    answer = Fraction(0)
    for exponents, coefficient in polynomial.items():
        k, l = exponents[:5], exponents[5:]
        for uvec, u, cu in choices(k):
            for vvec, v, cv in choices(l):
                b_path = tuple(k[i] - uvec[i] + vvec[i] for i in IDENTITY)
                c_path = tuple(l[i] - vvec[i] + uvec[i] for i in order)
                answer += (int(coefficient) * cu * cv * path_moment(b_path) * path_moment(c_path)
                           * weight_moment(sum(k) - u, u, v, sum(l) - v))
    return answer


def reference_setup(order):
    verify = Path(__file__).resolve().parents[1] / 'verify'
    sys.path.insert(0, str(verify))
    from vcommon import theta_graph
    from v_integrand import flows_m1
    nv, edges, paths = theta_graph(1, 6, 6)
    basis = [[int(i == j) for i in IDENTITY] for j in IDENTITY]
    lam = {0: [0] * 5, 11: [0] * 5}
    lam.update({i + 1: basis[i] for i in IDENTITY})
    lam.update({i + 6: [-x for x in basis[j]] for i, j in enumerate(order)})
    external = {v: dict(lam=lam[v], u=basis[0], vg=basis[1]) for v in range(nv)}
    return nv, edges, paths, lam, external, flows_m1(nv, edges)


def main():
    polynomial = numerator()
    for order, determinant, expected in (
        (SWAPPED, 1, Fraction(-3419207, 5558054638500)),
        (IDENTITY, -1, Fraction(0)),
    ):
        nv, edges, paths, lam, external, flows = reference_setup(order)
        from vcommon import D_theta, kirchhoff
        from v_integrand import grouped_integrand
        flow, types = next((f, t) for f, t in flows if t[0] == 'c' and all(h != 1 for _, h in f))
        assert D_theta(lam, paths) == determinant
        for shift in (4, 6, 11):
            times = [1 + (3 * e + shift) % 17 for e in range(13)]
            y = kirchhoff(nv, edges, times, lam)
            exact = grouped_integrand(nv, edges, paths, flow, types, y, external)
            currents = dict(zip(polynomial.ring.symbols, y[1] + y[7]))
            assert polynomial.as_expr().subs(currents) == exact != 0
        averaged = mean(polynomial, order)
        assert averaged == expected
        print(order, 'signed mean', averaged, 'kernel/pi', determinant * averaged / 158400)


if __name__ == '__main__':
    main()
