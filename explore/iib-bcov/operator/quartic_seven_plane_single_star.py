"""Rational independent integration after a source-free tetrahedron star elimination."""
from collections import defaultdict
from fractions import Fraction as F
from functools import cache
from math import factorial, prod
from itertools import combinations
import sys
from flint import fmpq, fmpq_mpoly_ctx
from quartic_seven_plane_single_projection import Q, integrated_weights
from quartic_seven_plane_single import cuts, orbit, tetrahedron_bank


@cache
def parallel_pair(a,b):
    assert 0<a<Q and 0<b<Q and a+b>Q
    return F(factorial(Q-a-1)*factorial(Q-b-1)*factorial(a+b-Q-1),factorial(2*Q-a-b-1))


def star_integral(raw):
    # Old 01,02,12 survive. The three new edges opposite old 03,13,23 are 12,02,01.
    ctx=fmpq_mpoly_ctx.get(tuple(f't{i}' for i in range(6)))
    t=ctx.gens()
    s=t[3]+t[4]+t[5]
    terms=defaultdict(dict)
    for p,w in raw.items():
        m=[x-1 for x in p]
        a,b,c=(m[i] for i in (2,4,5))
        mon=(m[0],m[1],m[3],1+b+c,1+a+c,1+a+b)
        v=w*F(factorial(Q-1),prod(map(factorial,m)))
        terms[a+b+c][mon]=fmpq(v.numerator,v.denominator)
    total=ctx.constant(0)
    for degree,row in terms.items():
        poly=ctx.from_dict(row)
        n=Q-3-degree
        if n<0:
            divisor=s**(-n)
            poly,rem=divmod(poly,divisor)
            assert not rem
        else:
            poly*=s**n
        total+=poly
    assert not divmod(total,s)[1]
    value=F()
    for mon,coefficient in total.terms():
        p=tuple(1+int(k) for k in mon)
        assert sum(p)==4*Q
        moment=prod(parallel_pair(p[i],p[j]) for i,j in ((0,5),(1,4),(2,3)))/factorial(Q-1)
        value+=F(int(coefficient.numerator),int(coefficient.denominator))*moment
    return value


def star_moment(powers):
    """Direct rational normalized tetrahedron cut at any vertex of degree sum <= Q."""
    edges=list(combinations(range(4),2))
    weight=dict(zip(edges,powers))
    for center in range(4):
        others=[j for j in range(4) if j!=center]
        arms=[weight[tuple(sorted((center,j)))] for j in others]
        excess=Q-sum(arms)
        if excess<0:
            continue
        opposite=[weight[tuple(x for x in others if x!=j)] for j in others]
        ans=F()
        for i in range(excess+1):
            for j in range(excess-i+1):
                h=(i,j,excess-i-j)
                term=F(factorial(excess),prod(map(factorial,h)))
                for k,b in enumerate(opposite):
                    term*=parallel_pair(sum(arms)-arms[k]+h[k],b)
                ans+=term
        return ans/prod(factorial(i-1) for i in powers)
    return None


def verify_rational_cells():
    count=0
    for powers,(a,b) in tetrahedron_bank().items():
        value=star_moment(powers)
        if value is None:
            assert b
        else:
            assert (a,b)==(value,0)
            count+=1
    assert count==1088
    return count


if __name__=='__main__':
    if '--moments' in sys.argv:
        print('rational tetrahedron cells',verify_rational_cells(),flush=True)
    marks=(0,0,0,4,6)
    raw=integrated_weights(marks,fold=False)
    v=star_integral(raw)
    reduced=defaultdict(F)
    for p,w in raw.items():
        reduced[orbit(p)]+=w
    assert cuts(reduced)==(v,0)
    assert v==-F(238826266733,22279187700000000)
    print('star-triangle integral',v)
