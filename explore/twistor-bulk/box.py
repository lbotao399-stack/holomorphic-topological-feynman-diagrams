"""Directed twistor box and the physical one-loop anomaly coefficients."""
from itertools import combinations, permutations
from math import factorial, prod

import sympy as sp

T = sp.symbols('t:4')


def wheel(order, gravity, brackets):
    """Dirichlet moment after factoring out the common holomorphic determinant."""
    value = sp.Poly(1, *T)
    for vertex in gravity:
        row = [sum(brackets[vertex, other] for other in order[:j+1])
               for j in range(4)]
        local = sum(t*(row[order.index(vertex)]-v) for t,v in zip(T,row))
        value *= sp.Poly(local,*T)
    return sp.expand(sum(c*prod(map(factorial,k)) for k,c in value.terms()) /
                     factorial(3+len(gravity)))


def trace(order, brackets):
    return prod(brackets[i,j] for i,j in zip(order, (*order[1:],order[0])))


def source_brackets():
    brackets = sp.zeros(4)
    for (i,j), symbol in zip(combinations(range(3),2),sp.symbols('u v w')):
        brackets[i,j], brackets[j,i] = symbol, -symbol
    for i in range(3):
        brackets[i,3] = -sum(brackets[i,j] for j in range(3))
        brackets[3,i] = -brackets[i,3]
    return brackets


def verify():
    b = source_brackets()
    orders = [(0,*rest) for rest in permutations((1,2,3))]
    for order in orders:
        assert wheel(order,[],b) == sp.Rational(1,6)
    # Distinct color words remain separate when there are three gauge insertions.
    for next_label in (2,3):
        same = [o for o in orders if o[1:][(o[1:].index(1)+1)%3] == next_label]
        assert len(same) == 3
        assert sp.expand(sum(wheel(o,[0],b) for o in same)) == 0
    for i,j in combinations(range(4),2):
        mixed = sum(wheel(o,(i,j),b) for o in orders)
        assert sp.expand(mixed+2*trace((i,j),b)/24) == 0
    gravity = sum(wheel(o,range(4),b) for o in orders)
    assert sp.expand(gravity-sum(trace(o,b) for o in orders)/120) == 0
    print('all six cyclic orders; arbitrary conserved momenta: pure gauge 1/6, mixed -1/24, gravity 1/120')


if __name__ == '__main__':
    verify()
