"""Independent degree-four cubature of the original 13-edge form at fixed core times."""
from fractions import Fraction as F
from functools import cache
from itertools import product
from math import comb, factorial, prod
from random import Random
from flint import fmpq

from chain_moments import Reduction, graphs, vcore


@cache
def binomial_coefficients(n):
    # Coefficients of binomial(10 x - 1, n) in the ordinary power basis.
    a = [1]
    for j in range(1,n+1):
        a = [-j*x+10*y for x,y in zip(a+[0],[0]+a)]
    return [F(c,factorial(n)) for c in a]


def cubature(red):
    weights = {}
    for a in product(range(5),repeat=4):
        if sum(a)>4:
            continue
        avg = F(0)
        for m in product(*(range(v+1) for v in a)):
            moment = F(prod(factorial(v) for v in m),
                       prod(factorial(len(g)+sum(m[i] for i in g)) for g in red.groups))
            avg += moment*prod(binomial_coefficients(v)[j] for v,j in zip(a,m))
        for k in product(*(range(v+1) for v in a)):
            pt = tuple(F(i+1,10) for i in k)
            weights[pt] = weights.get(pt,0)+avg*(-1)**(sum(a)-sum(k))*prod(
                comb(v,j) for v,j in zip(a,k))
    return {point:weight for point,weight in weights.items() if weight}


def raw(red, R, p, x):
    t = []
    for e,g in enumerate(red.groups):
        shares = [x[i] for i in g]
        assert sum(shares)<1
        t.extend(R[e]*v for v in [1-sum(shares)]+shares)
    val,U,_ = vcore.density_slice(red.net,red.expanded_cycles,t,[(0,0)]*6+p[:3])
    return val*prod(R[e]**len(g) for e,g in enumerate(red.groups)),U


if __name__=='__main__':
    rng = Random(584)
    for i,g in enumerate(graphs):
        red = Reduction(g)
        R = [rng.randrange(2,10) for e in range(9)]
        p = [(rng.randrange(-5,6),rng.randrange(-5,6)) for _ in range(3)]
        p += [tuple(-sum(col) for col in zip(*p))]
        x = [F(k+1,12) for k in range(4)]
        assert red.density(R,p,at=[fmpq(t.numerator,t.denominator) for t in x]) == raw(red,R,p,x)
        if i in (0,3,11,39):
            tested = {pt:raw(red,R,p,pt)[0] for pt in cubature(red)}
            weight = cubature(red)
            assert sum(weight[pt]*v for pt,v in tested.items()) == red.density(R,p)[0]
            print('cubature',i,len(tested),'points passed',flush=True)
    print('all 50 original-graph signed comparisons passed')
