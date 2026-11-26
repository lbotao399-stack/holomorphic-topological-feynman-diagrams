"""Four-simplex order cubatures for every marking of the integrated cubic."""
from collections import Counter
from fractions import Fraction as F
from functools import cache
from math import factorial,prod
import numpy as np
from flint import fmpz_mat
from theta_cubature import cube,UNIT
from marked_theta import partitions,conductance
from integrated_cubic import PROJECTED,run


@cache
def green(paths):
    linear=[F() for i in range(8)]
    product=[[F() for j in range(8)] for i in range(8)]
    bridge=[[F() for j in range(8)] for i in range(8)]
    zero=[0]*3
    cols=[(v if h!=2 else zero)+(v if h!=0 else zero) for v,h in zip(UNIT,paths)]
    # The first four independent tail times drop out of the determinant.
    nodes=cube()
    for d in (5,7,9,11):
        xlist=[t for t,den in nodes if den==d]
        j=(11-d)//2
        w=F((-1)**j*d**7,64*factorial(j)*factorial(11-j))
        multipliers=[]
        for t in xlist:
            v=[a-t[3] for a in t[:3]]
            mult=nodes[t,d]/w
            assert mult.denominator==1
            multipliers.append(mult.numerator*int(fmpz_mat([v+zero,zero+v,*cols]).det())**2)
        weights=np.array(multipliers,dtype=np.int64)
        vals=[None]*4+list(np.array(xlist,dtype=np.int64).T)
        def mean(row):
            return int(weights.sum()*row if np.isscalar(row) else weights @ row)*w/(12*d**6)
        for i,x in enumerate(vals):
            linear[i]+=mean(6*d*d if x is None else 12*d*x)
            for k,y in enumerate(vals):
                if x is None and y is None:
                    xy=(4 if i==k else 3)*d*d
                    brown=(2 if i==k else 1)*d*d
                elif x is None or y is None:
                    active=y if x is None else x
                    xy=6*d*active
                    brown=6*active*(d-active)
                else:
                    xy=12*x*y
                    brown=12*d*np.minimum(x,y)-xy
                product[i][k]+=mean(xy)
                bridge[i][k]+=mean(brown)
    return linear,product,bridge


def integrate(shape):
    answer=[F(),F()]
    deg=[6-i for i in shape]
    for groups in partitions(shape):
        path={i:h for h,g in enumerate(groups) for i in g}
        if path[0]==path[1] or path[2]==path[3]:
            continue
        linear,product,bridge=green(tuple(path[i] for i in range(4,8)))
        base=conductance(deg)
        phi=[F() for i in path]
        north=F()
        for j,q in enumerate(PROJECTED):
            if q:
                north+=q*linear[j]*base
                for i in range(8):
                    phi[i]+=q*product[i][j]*base
                    if path[i]==path[j]:
                        powers=list(deg)
                        powers[path[i]]-=1
                        phi[i]+=q*bridge[i][j]*conductance(powers)
        answer[0]+=2*phi[7]-north
        answer[1]+=8*phi[7]-sum(phi)
    divisor=2*prod(factorial(n) for n in Counter(shape).values())
    return tuple(x/divisor for x in answer)


if __name__=='__main__':
    total=F()
    for shape in ((1,5,5),(2,4,5),(3,3,5),(3,4,4)):
        exact=integrate(shape)
        assert exact==run(shape)
        total+=sum(exact)
    assert total==F(1,1215)
    print('Eight integrated shape/marking coefficients independently checked, sum',total)
