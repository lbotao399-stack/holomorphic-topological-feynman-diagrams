"""Two-loop contact and counterterm with a single scalar Poisson cubic."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations, permutations
from math import prod

import sympy as sp
from sympy.polys.rings import ring

from theta_moyal import Theta, permuted_monomial, ribbon


def poisson_graph(shape):
    theta = Theta(shape)
    u, v, w = theta.brackets
    # The scalar trivalent Poisson vertex has coefficient [q1 q2].
    integrand = u * (v+w) * prod(theta.vertices)
    moments = defaultdict(F)
    for exponents, integer in integrand.items():
        moments[exponents[:3]] += int(integer) * theta.moment(exponents[3:])
    factor = ribbon.orientation_sign(6, theta.edges, theta.ext) * theta.prefactor / 4
    answer = defaultdict(F)
    for order in permutations(range(4)):
        sign = ribbon.perm_sign(order)
        for powers, coefficient in moments.items():
            for permuted, integer in permuted_monomial(powers, order).items():
                answer[permuted] += factor * sign * coefficient * int(integer)
    return {k: c for k, c in answer.items() if c}


def shuffle(first_degree, second_degree):
    algebra, u, v, w = ring("u v w", sp.ZZ)
    brackets = {(0,1):u, (0,2):v, (0,3):-u-v, (1,2):w, (1,3):u-w, (2,3):v+w}
    result = algebra.zero
    for i, j in combinations(range(4), 2):
        k, l = [s for s in range(4) if s not in (i, j)]
        result += (-1)**(i+j-1) * brackets[i,j]**first_degree * brackets[k,l]**second_degree
    return {powers: int(c) for powers, c in result.items()}


def probe(cochain):
    # (p0,p1,p2,p3) = ((1,0),(0,1),(-3,2),(2,-3)), hence (u,v,w) = (1,2,3).
    return sum(c * 2**v * 3**w for (_,v,w), c in cochain.items())


def check():
    first = poisson_graph((1,3,3))
    second = poisson_graph((2,2,3))
    tree = shuffle(1,7)
    assert not shuffle(1,1)
    assert probe(first) == F(4823,56100)
    assert probe(second) == -F(1483,336600)
    assert probe(tree) == 65520
    coefficient = F(323,259459200)
    for powers in first.keys() | second.keys() | tree.keys():
        assert first.get(powers,0) + second.get(powers,0) == coefficient * tree.get(powers,0)
    print("Theta probes:", probe(first), probe(second), "sum", probe(first) + probe(second))
    print("Full polynomial verified: (323/259459200) * { [pq], [pq]^7 }")


if __name__ == "__main__":
    check()
