"""Symbolic path-ratio integrals in both multi-trace Moyal channels."""
from collections import defaultdict
from itertools import combinations, permutations
from math import factorial, prod

import sympy as sp
from sympy.polys.fields import field

from probe_laurent import DOUBLE, PAIRS, TRIPLE, RibbonProbe


# Polynomial-exponential coefficients have keys (frequency, power of e);
# the common e^-4 from the four path-ratio integrations is kept outside.
K, a, b = field("a,b", sp.QQ)
weights = (a, b, 1-a-b)


def divided_difference(nodes):
    """Coefficients of the divided difference of exp(e*x), including repeated nodes."""
    unique = set(nodes)
    if len(unique) == 1:
        degree = len(nodes)-1
        return {(nodes[0],degree): K.one/factorial(degree)}
    if len(unique) == len(nodes):
        return {(node,0): 1/prod(node-v for v in nodes if v != node) for node in nodes}
    nodes = list(nodes)
    j = next(j for j,v in enumerate(nodes) if v != nodes[0])
    nodes[j],nodes[-1] = nodes[-1],nodes[j]
    result = defaultdict(lambda: K.zero)
    for key,c in divided_difference(nodes[1:]).items():
        result[key] += c/(nodes[-1]-nodes[0])
    for key,c in divided_difference(nodes[:-1]).items():
        result[key] -= c/(nodes[-1]-nodes[0])
    return result


def integrate_ratios(phase, shape, shares):
    coefficients = defaultdict(lambda: K.zero)
    for (i,j,*powers), integer in phase:
        assert sum(powers) <= 1 and i+j <= 1
        coefficients[tuple(powers)] += shares[0]**i * shares[1]**j * integer / 2
    constant = coefficients[(0,0,0,0)]
    slopes = [coefficients[tuple(int(k == l) for k in range(4))] for l in range(4)]
    result = {(constant,0): K.one}
    pos = 0
    for length in shape:
        simplex = divided_difference([K.zero]+slopes[pos:pos+length-1])
        pos += length-1
        terms = defaultdict(lambda: K.zero)
        for (x,m),c in result.items():
            for (y,n),d in simplex.items():
                terms[x+y,m+n] += c*d*factorial(length-1)
        result = terms
    return result


def phase_integral(momenta, word, expected_rank):
    result = defaultdict(lambda: K.zero)
    for shape in ((1,3,3),(2,2,3)):
        probe = RibbonProbe(shape, momenta, word)
        prefactor = probe.theta.prefactor / (4*6)
        for order in permutations(weights):
            density = 5*order[0]**2 if shape == (1,3,3) else 15*order[0]*order[1]
            for (rank,phase),orientation in probe.phases.items():
                assert rank == expected_rank
                for key,c in integrate_ratios(phase,shape,order).items():
                    result[key] += orientation * prefactor.numerator * density * c / prefactor.denominator
    return result


def check():
    result = phase_integral(PAIRS, TRIPLE, 0)
    expected = defaultdict(lambda: K.zero)
    for sign in (-1,1):
        for size in range(4):
            for subset in combinations(weights,size):
                expected[sign*sum(subset,K.zero),0] += (-1)**(3-size) / prod(weights)
    assert all(result[key] == expected[key] for key in result.keys() | expected.keys())
    print("The three-trace exponential identity holds over Q(w_A,w_B).")

    momenta = ((1,0),(0,1),(0,1),(-1,-2))
    result = phase_integral(momenta, DOUBLE, 1)
    expected = defaultdict(lambda: K.zero)
    for i, share in enumerate(weights):
        others = [w for j,w in enumerate(weights) if j != i]
        density = sum(w**2 for w in others) / (3*share**2*(1-share)**2)
        first = [(share,K.one/2),(-share,K.one/2),(K.zero,-K.one)]
        second = [(1-share,K.one/2),(share-1,K.one/2),(K.zero,-K.one)]
        for x,c in first:
            for y,d in second:
                expected[x+y,1] += c*d*density
    assert all(result[key] == expected[key] for key in result.keys() | expected.keys())
    print("The two-trace exponential identity holds over Q(w_A,w_B).")


if __name__ == "__main__":
    check()
