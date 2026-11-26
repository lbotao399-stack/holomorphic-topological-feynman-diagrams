"""Exact checks of chamber grouping, both expansions and both numerator signs."""
from fractions import Fraction as F
from itertools import permutations
from math import prod

from period import HeppSeries


def direct_chambers(engine, a, beta, cotree=None):
    total=F(0)
    for order in permutations(range(engine.n)):
        # No determinant or subset recurrence is used for this sum. In the
        # chamber, reversed lexicographic order picks the largest cotree product.
        largest = max(engine.cotrees, key=lambda I:tuple((I >> i) & 1 for i in reversed(order)))
        if cotree is not None and cotree!=largest:
            continue
        margin=0
        factors=[]
        for e in order[:-1]:
            margin += a[e]+1-beta*((largest>>e)&1)
            factors.append(margin)
        total += F(1)/prod(factors)
    return total


def arctan_interval(n, depth):
    partial=sum(F((-1)**k,(2*k+1)*n**(2*k+1)) for k in range(depth))
    other=partial+F((-1)**depth,(2*depth+1)*n**(2*depth+1))
    return min(partial,other),max(partial,other)


def check():
    # The holomorphic one-loop triangle in two complex dimensions. The simplex
    # chart has U=1 and constant density 1 for unit source determinant.
    triangle=HeppSeries([(0,1),(1,2),(2,0)],3,{(0,0,0):1})
    # An asymmetric signed numerator tests the dominant-cone decomposition.
    signed=HeppSeries([(0,1),(1,2),(2,0)],5,{(2,0,0):3,(0,1,1):-2})
    # Simplex Dirichlet moments 2!/4! and 1/4!, respectively.
    for engine, target in ((triangle,F(1,2)),(signed,F(1,6))):
        for method in ('binomial','gegenbauer'):
            seq=list(engine.terms(15,method=method))
            symmetric=list(engine.terms(15,method=method,permutations='graph'))
            for (_,t,p,b),(_,t2,p2,b2) in zip(seq,symmetric):
                assert t==t2 and p==p2 and b2<=b and abs(F(str(p))-target)<=F(str(b))
            cone_sequences=[list(engine.terms(3,I,method)) for I in engine.cotrees]
            assert all(sum(s[k][2] for s in cone_sequences)==seq[k][2] for k in range(4))
    banana=HeppSeries([(0,1)]*3,F(3,2),{(0,0,0):1})
    # Integral over x>0,y>0 of (xy+x+y)^(-3/2): integrate y first to get
    # 2/((x+1)*sqrt(x)). Machin's formula encloses its integral 2*pi.
    a,b=arctan_interval(5,14);c,d=arctan_interval(239,4)
    lower,upper=32*a-8*d,32*b-8*c
    _,_,p,bound=list(banana.terms(30,permutations='graph'))[-1]
    p,bound=F(str(p)),F(str(bound))
    assert p-bound<lower<upper<p+bound and bound<F(1,10**13)
    from moment_bounds import enclose
    mass=[m for k,m,b in banana.moments(3,permutations='graph')]
    hermite_lo,hermite_hi=enclose(mass,banana.alpha,banana.tree_count)
    assert hermite_lo<lower<upper<hermite_hi
    k4=HeppSeries([(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)],2,{(0,)*6:1})
    a=(0,)*6
    H=k4.moment_function(F(2),1)
    assert F(str(H(k4.full,0)))==direct_chambers(k4,a,F(2))
    for I in k4.cotrees:
        H=k4.moment_function(F(2),1,I)
        assert F(str(H(k4.full,0)))==direct_chambers(k4,a,F(2),I)
    print('Triangle 1/2 and signed Dirichlet moment 1/6 verified at all 16 orders.')
    print('Three-edge banana 2*pi enclosed at order 30; error bound',bound)
    print('K4 720 direct chambers agree with global DP and all 16 cotree cones.')


if __name__=='__main__':
    check()
