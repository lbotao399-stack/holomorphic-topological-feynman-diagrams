"""Evaluate the matrix-unit BV shuffle directly on the commuting detector."""
from itertools import combinations, product
from math import prod
from fractions import Fraction

from invariant_ranks import balance, bracket
from poisson_center import columns, invariants


def trace(x):
    return sum(x[i][i] for i in range(len(x)))


def pair(x,y):
    return sum(x[i][j]*y[j][i] for i,j in product(range(len(x)), repeat=2))


def matrix(n, i, j):
    return [[int((a,b) == (i,j)) for b in range(n)] for a in range(n)]


def mixed(xs, ps):
    n = len(xs[0])
    t = [trace(x) for x in xs]
    return bracket(ps[0],ps[1]) * (sum(t[i]*pair(xs[i-1], xs[i-2])
                                      for i in range(3))-Fraction(3,n)*prod(t))


def shuffle(xs, ps, central_momentum):
    n = len(xs[0])
    total = 0
    for i,j in combinations(range(len(xs)), 2):
        rest = [s for s in range(len(xs)) if s not in (i,j)]
        P = tuple(x+y for x,y in zip(ps[i],ps[j]))
        polys = central_momentum([P]+[ps[s] for s in rest])
        sign = (-1)**(i+j-1)
        for a,b in product(range(n), repeat=2):
            left = mixed([xs[i],xs[j],matrix(n,a,b)], [ps[i],ps[j],tuple(-x for x in P)])
            right = trace(matrix(n,b,a))*prod(trace(xs[s]) for s in rest)*polys
            total += sign*left*right
    return total


for n in (2,3,5):
    identity = [[int(i==j) for j in range(n)] for i in range(n)]
    traceless = [[int(i==j==0)-int(i==j==1) for j in range(n)] for i in range(n)]
    for w in (4,6,8):
        ps = balance([(1,0),(0,2),(2,1)])
        value = shuffle([traceless,traceless,identity,identity], ps,
                        lambda q: bracket(q[1],q[2])**(w-1))
        assert value == 2*n**3*bracket(ps[0],ps[1])*bracket(ps[2],ps[3])**(w-1) != 0
    for w in range(3,10):
        ps = balance([(1,0),(0,2),(2,1),(-3,5)])
        P = tuple(a+b for a,b in zip(*ps[:2]))
        for powers in columns(w):
            def central(q):
                el,v = invariants(q)
                return prod(x**k for x,k in zip(el,powers))*(v if w%2==0 else 1)
            value = shuffle([traceless,traceless,identity,identity,identity], ps, central)
            assert value == 2*n**4*bracket(ps[0],ps[1])*central([P]+ps[2:]) != 0
    print('rank', n, 'mixed shuffle and five-input commuting detector pass')
