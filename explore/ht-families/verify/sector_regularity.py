"""Exact sector-frame, curvature and Kirchhoff boundary checks for T=1."""
from fractions import Fraction
from itertools import combinations
from math import prod
from random import Random

import sympy as sp

from vcore import Net, pfaffian_forms, rho_forms


def pfaffian(matrix_values, n):
    # A pair of tangent directions selects one coefficient of each curvature two-form.
    forms = [[{} for _ in range(n)] for _ in range(n)]
    for (a,b), matrix in matrix_values.items():
        mask = (1 << a) | (1 << b)
        for i,j in combinations(range(n), 2):
            if matrix[i,j]:
                forms[i][j][mask] = Fraction(matrix[i,j])
    return pfaffian_forms(forms).get((1 << n)-1, 0)


def check(net, order, rng):
    cycle = sp.Matrix(net.cycle_basis()).T[order,:]
    edges, n = cycle.shape
    basis = list(cycle.T.rref()[1])
    coordinates = cycle * cycle[basis,:].inv()
    for edge in range(edges):
        for j,i in enumerate(basis):
            assert not coordinates[edge,j] or i <= edge

    bases = [(J, cycle[J,:].det()) for J in combinations(range(edges), n)]
    bases = [(J,minor) for J,minor in bases if minor]
    assert all(all(i<=j for i,j in zip(basis,J)) for J,_ in bases)
    ratios = [sp.Rational(rng.randint(1,4),5) for _ in range(edges-1)]
    # Use all tangent directions through n random linear combinations so a
    # local zero minor is not mistaken for an orientation check.
    slopes = [[rng.randrange(-4,5) for _ in range(n)] for _ in ratios]
    roots = [prod(ratios[:e]) for e in range(edges)]
    times = [root**2 for root in roots]
    frame = sp.Matrix(edges,n,lambda e,j: coordinates[e,j]*prod(ratios[basis[j]:e]))
    gram = frame.T * frame
    complement = sp.eye(edges) - frame * gram.inv() * frame.T
    derivatives = []
    for a in range(n):
        derivatives.append(sp.Matrix(edges,n,lambda e,j:
            frame[e,j]*sum(slopes[k][a]/ratios[k] for k in range(basis[j],e))))
    curvature = {}
    for a,b in combinations(range(n),2):
        raw = derivatives[a].T * complement * derivatives[b]
        curvature[a,b] = raw-raw.T
    pf = pfaffian(curvature,n)
    dts = [{1<<a:2*times[e]*sum(slopes[k][a]/ratios[k] for k in range(e))
            for a in range(n)} for e in range(edges)]
    original,U,_ = rho_forms(cycle.T.tolist(),times,dts)
    minor = cycle[basis,:].det() * prod(roots[i] for i in basis)
    assert U == gram.det()*minor**2
    assert original.get((1<<n)-1,0) == 2**(n//2)*minor*pf

    source_flow = sp.Matrix([rng.randrange(-7,8) for _ in range(edges)])
    tree_currents = [(J,d,-cycle[J,:].inv()*source_flow[list(J),:]) for J,d in bases]

    def average(values):
        weights = [d**2*prod(prod(values[i:j])**2 for i,j in zip(basis,J))
                   for J,d,_ in tree_currents]
        assert sum(weights) >= 1
        return sum((w*ell for w,(_,_,ell) in zip(weights,tree_currents)),sp.zeros(n,1))/sum(weights)

    loop = average(ratios)
    assert cycle.T*sp.diag(*times)*(source_flow+cycle*loop) == sp.zeros(n,1)
    boundary = list(ratios)
    for i in range(0,len(boundary),2):
        boundary[i] = 0
    average(boundary)
    corner = [0]*len(ratios)
    assert average(corner) == -cycle[basis,:].inv()*source_flow[basis,:]
    return pf != 0


if __name__=="__main__":
    rng = Random(291)
    fixtures = [
        Net(4,[(0,1),(0,2),(0,3),(1,2),(1,3)]),
        Net(6,[(0,1),(0,2),(0,4),(0,5),(1,2),(1,3),(1,4),(3,4),(4,5)]),
        Net(7,[(i,j) for i in range(2) for j in range(2,7)]+[(2,3),(4,5)]),
    ]
    for graph in fixtures:
        count = 0
        for trial in range(3):
            order = list(range(graph.N))
            rng.shuffle(order)
            count += check(graph,order,rng)
        assert count, "all curvature comparisons had zero top coefficient"
        print(f"{graph.n} loops, 3 sectors: regular coordinates and tree means, {count} nonzero Pfaffian checks")
