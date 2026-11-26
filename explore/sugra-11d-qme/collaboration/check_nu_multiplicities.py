"""Nonconstant nu modes on the theta166 endpoint component and a two-species check."""
from fractions import Fraction
from itertools import product
from math import comb, factorial, prod

from check_endpoint_moment import direct_component, mass
from nu_multiplicities import continue_multiplicities, falling_moment, rising
from vcommon import D_theta, theta_graph


def compositions(n, length):
    if length==1:
        yield (n,)
        return
    for first in range(n+1):
        for rest in compositions(n-first,length-1):
            yield (first,)+rest


SHAPE = (1,6,6)
ZERO = [[0]*n for n in (3,*SHAPE)]
MASS = mass(SHAPE,ZERO)
MOMENTUM = (1,2,0,0,0)


def endpoint(occupations):
    """Cubic gamma at N. All eleven b vertices may receive external nu.
    Compensate the third external vector's momentum at vertex 3."""
    loads = list(occupations)
    loads[2] -= sum(loads)
    # Matrix determinant lemma in the two five-column path blocks B,C.
    det = -(1+sum(k*a for k,a in zip(MOMENTUM,loads[:5]))) * (
        1-sum(k*a for k,a in zip(reversed(MOMENTUM),loads[5:10])))
    # Uniform path fractions: the j-th interior source enters (6-j)/6 of the mean.
    lb,lc = (sum((5-i)*a for i,a in enumerate(path)) for path in (loads[:5],loads[5:10]))
    wedge = (5+lb)*(-2+2*lc) - (4+2*lb)*(-1+lc)
    return MASS*det*Fraction(11,13)*Fraction(wedge,36)


def direct(occupations):
    loads = list(occupations)
    loads[2] -= sum(loads)
    sources = {v+1: [int(i==v) if v<5 else -int(i==9-v) for i in range(5)]
               for v in range(10)}
    for v in range(10):
        sources[v+1] = [x+loads[v]*k for x,k in zip(sources[v+1],MOMENTUM)]
    _,_,paths = theta_graph(*SHAPE)
    prefix = []
    for verts,_ in paths:
        seq = [[0]*5]
        for v in verts:
            seq.append([x+y for x,y in zip(seq[-1],sources[v])])
        prefix.append(seq)
    return D_theta(sources,paths)*direct_component(prefix,SHAPE,ZERO,0,1)


def tests():
    for loads in ([0]*11,[1,0,2,0,0,0,1,0,0,1,1],[0,2,0,1,2,0,0,1,0,0,2]):
        assert direct(loads) == endpoint(loads)
    observed = {k:factorial(k)*sum(map(endpoint,compositions(k,11))) for k in range(7)}
    samples = {(k,): observed[k] for k in range(4)}
    for k,value in observed.items():
        assert continue_multiplicities(11,(k,),samples) == value
    normalized = [samples[(k,)]/rising(11,k) for k in range(4)]
    diffs = []
    while normalized:
        diffs.append(normalized[0])
        normalized = [b-a for a,b in zip(normalized,normalized[1:])]
    assert diffs[0]==Fraction(1,1123200)
    print("theta166 endpoint coefficient / pi divided by (11)_K, in binomial(K,j) order:",diffs)
    # Track four distinguishable labels, of two species, at two different sites.
    orders = {(0,0):1,(0,1):1,(1,0):1}
    for sizes in ((1,2),(2,2),(2,3),(3,2)):
        total = 0
        for assignment in product(range(3),repeat=sum(sizes)):
            left,right = assignment[:sizes[0]],assignment[sizes[0]:]
            total += (left.count(0)*left.count(1)*right.count(0)
                      *prod(factorial(assignment.count(i)) for i in range(3)))
        assert falling_moment(3,sizes,orders) == Fraction(total,rising(3,sum(sizes)))
    grid = {(a,b): rising(3,a+b)*falling_moment(3,(a,b),orders)
            for a in range(4) for b in range(4-a)}
    assert continue_multiplicities(3,(7,8),grid)==rising(3,15)*falling_moment(3,(7,8),orders)
    print("All 12,376 weak compositions at K=0,...,6 and four distinguishable-label checks passed.")


if __name__=="__main__":
    tests()
