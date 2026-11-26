"""Integrated scalar cubic at two loops, including every marked vertex."""
from collections import Counter
from fractions import Fraction as F
from math import factorial,prod
from functools import cache
from marked_theta import MOMENTA,partitions,square,green_moments,conductance

# Scalar momenta are e1,e1,e2,e2,e3,e4,e5. The vector wave balances them.
# Project onto e_1-2e_3, transverse to its momentum.
PROJECTED=(1,1,0,0,-2,0,0,0)


@cache
def linear(k,j):
    return F(k[j]+1,(k[j]+2)*prod(v+1 for v in k))


def run(shape,momenta=MOMENTA,projected=PROJECTED):
    ds=[6-k for k in shape]
    base=conductance(ds)
    answer=[F(),F()]
    for groups in partitions(shape):
        path={i:h for h,g in enumerate(groups) for i in g}
        for k,c in square(groups,momenta).items():
            top=F()
            phi=[F() for _ in path]
            for j,q in enumerate(projected):
                if not q:
                    continue
                top+=q*linear(k,j)
                for i in range(8):
                    xy,bridge=green_moments(k,i,j)
                    phi[i]+=q*xy
                    if path[i]==path[j]:
                        phi[i]+=q*F(9,ds[path[i]])*bridge
            # N' grounded; the vector is at label 7. A mark on that vertex
            # would put a vector in the scalar cubic and thus gives zero.
            answer[0]+=c*base*(2*phi[7]-top)
            answer[1]+=c*base*(8*phi[7]-sum(phi))
    divisor=2*prod(factorial(n) for n in Counter(shape).values())
    return tuple(a/divisor for a in answer)


if __name__=='__main__':
    totals=[F(),F()]
    for shape in ((1,5,5),(2,4,5),(3,3,5),(3,4,4)):
        values=run(shape)
        print(shape,*values,'total',sum(values))
        totals=[x+y for x,y in zip(totals,values)]
    print('trivalent, bivalent, combined',*totals,sum(totals))
