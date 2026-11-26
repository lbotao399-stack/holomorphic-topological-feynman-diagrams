"""Two-loop marked scalar square in minimal BCOV on the fivefold."""
from collections import Counter
from fractions import Fraction as F
from functools import cache
from itertools import combinations
from math import factorial,prod
from flint import fmpz_mat

# The last two labels are the transverse vector and the point insertion.
MOMENTA=tuple(tuple(int(i==j) for j in range(5)) for i in (0,0,1,1,2,3,4))+((-2,-2,-1,-1,-1),)
CURRENT=(1,1,1,1,1,1,0,-6)


def partitions(shape):
    labels=set(range(sum(shape)-3))
    for a in combinations(sorted(labels),shape[0]-1):
        rest=labels-set(a)
        for b in combinations(sorted(rest),shape[1]-1):
            yield a,b,tuple(sorted(rest-set(b)))


def square(groups,momenta=MOMENTA):
    cols=[]
    zero=(0,)*len(momenta[0])
    for path,group in enumerate(groups):
        for i in group:
            p=momenta[i]
            cols.append((p if path!=2 else zero)+(p if path!=0 else zero))
    # Determinant multilinearity in the first two columns, at L=sum(s_i p_i).
    quad=Counter()
    for i,p in enumerate(momenta):
        for j,q in enumerate(momenta):
            value=int(fmpz_mat([p+zero,zero+q,*cols]).det())
            if value:
                quad[tuple(sorted((i,j)))]+=value
    out=Counter()
    for a,x in quad.items():
        for b,y in quad.items():
            if x*y:
                monomial=a+b
                out[tuple(monomial.count(i) for i in range(len(momenta)))]+=x*y
    return {k:v for k,v in out.items() if v}


@cache
def green_moments(k,i,j):
    # Integrate s_i s_j and min(s_i,s_j)-s_i s_j on the uniform cube.
    cube=F(1,prod(a+1 for a in k))
    x,y=k[i],k[j]
    if i==j:
        product=cube*F(x+1,x+3)
        minimum=cube*F(x+1,x+2)
    else:
        product=cube*F((x+1)*(y+1),(x+2)*(y+2))
        minimum=cube*F((x+1)*(y+1),x+y+3)*(F(1,x+2)+F(1,y+2))
    return product,minimum-product


def conductance(degrees):
    return F(prod(factorial(k) for k in degrees),factorial(2+sum(degrees)))


def integrate(shape):
    degrees=tuple(6-j for j in shape)
    base=conductance(degrees)
    local=[conductance([v-int(h==i) for h,v in enumerate(degrees)]) for i in range(3)]
    value=F()
    active=0
    for groups in partitions(shape):
        polynomial=square(groups)
        if not polynomial:
            continue
        active+=1
        path={i:h for h,group in enumerate(groups) for i in group}
        reduced=F()
        for k,c in polynomial.items():
            for at,sign in ((6,1),(7,-1)):
                for j,source in enumerate(CURRENT):
                    if source:
                        xy,brownian=green_moments(k,at,j)
                        weight=base*xy
                        if path[at]==path[j]:
                            weight+=local[path[j]]*brownian
                        reduced+=sign*c*source*weight
        assert reduced>0
        value+=reduced
    aut=2*prod(factorial(n) for n in Counter(shape).values())
    return active,value/aut


if __name__=='__main__':
    result=F()
    for shape in ((1,5,5),(2,4,5),(3,3,5),(3,4,4)):
        active,value=integrate(shape)
        print(shape,active,value,flush=True)
        result+=value
    print('sum',result)
