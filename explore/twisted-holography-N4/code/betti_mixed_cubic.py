"""Polarized cubic and quartic Betti intersections for the mixed D3 bound."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement, product
from math import comb, factorial, prod
from pathlib import Path
import sys

from betti_flow_polytope import weak
from betti_power_graphs import LINKS, graph_terms, scalar_power, vector_coefficient, weighted_vertex


def triangle_subdivision(charges):
    """Exact constant-in-m limit of the up/down toric localization sum on x+y+z=1."""
    import sympy as sp
    i,j,m=sp.symbols('i j m',integer=True,positive=True)
    x,y=sp.symbols('x y')
    a,b,c=charges
    h=(a*x*y+b*y*(1-x-y)+c*(1-x-y)*x)/2
    X,Y,Z=1,2,-3
    slope=sp.diff(h,x)*X+sp.diff(h,y)*Y
    up=slope+(sp.diff(h,x,2)*X+sp.diff(h,y,2)*Y)/(2*m)
    down=up+(X+Y)*sp.diff(h,x,y)/m
    def triangle_sum(poly,offset):
        lattice=sp.Poly(poly.subs({x:i/m,y:j/m}),i,j)
        return sum(v*sp.summation(i**k*sp.summation(j**l,(j,0,m-offset-i)),(i,0,m-offset))
                   for (k,l),v in lattice.terms())
    moment=sp.expand(-(triangle_sum(up**3,1)-triangle_sum(down**3,2))/(6*m))
    assert sp.degree(sp.cancel(moment*m**4),m)<=4
    faces=(a**3*X*Y*(X+Y)+b**3*Y*Z*(Y+Z)+c**3*X*Z*(X+Z))/sp.Integer(240)
    return sp.factor((moment.coeff(m,0)-faces)/(X*Y*Z))


def bound(retained):
    degree=len(retained)
    if degree not in (3,4,5):
        raise ValueError('Retain between three and five actual Betti forms.')
    remaining=Counter({(i,j):m for i,j,m in LINKS})
    remaining.subtract(retained)
    if min(remaining.values())<0: raise ValueError('Invalid retained multiplicity')
    degrees=Counter({(0,)*5:1})
    for (i,j),m in remaining.items():
        expanded=Counter()
        for powers,weight in degrees.items():
            for k in range(m+1):
                v=list(powers)
                v[i]+=k; v[j]+=m-k
                expanded[tuple(v)]+=weight*comb(m,k)
        degrees=expanded
    # Polarize the sum of actual forms, putting each retained charge in its own target coordinate.
    stencil=Counter()
    for count in range(1,degree+1):
        for pairs,multiplicity in Counter(combinations(retained,count)).items():
            waves=tuple(tuple(int(i==a)-int(i==b) for a,b in pairs) for i in range(5))
            for g,f,k,w in graph_terms(degree,waves):
                stencil[g,f,k]+=(-1)**(degree-count)*multiplicity*w
    total=F()
    for (groups,flags,extra),weight in stencil.items():
        if not weight: continue
        value=0
        for powers,mu in degrees.items():
            shifted=tuple(x+y for x,y in zip(extra,powers))
            for group,branch in zip(groups,flags):
                mu*=weighted_vertex(group,branch,shifted)
                if not mu: break
            value+=mu
        total+=value*weight
    return 2*total/(2**(14-degree)*prod(factorial(m) for _,_,m in LINKS))


def check_vertices():
    import sympy as sp
    sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'iib-bcov'))
    from elliptic_descendants import vertex
    a,b,c=sp.symbols('a b c')
    assert sp.expand(triangle_subdivision((a,b,c))+(a+b+c)**3/960
                     -(a+b+c)*(a*b+b*c+c*a)/120+7*a*b*c/240)==0
    for target in (3,4):
        marks=[tuple(int(i==j) for j in range(target)) for i in range(target)]+[(-1,)*target]
        for n in range(target+1,target+4):
            if n>target+1:
                v=(n-target,2,-1,1)[:target]
                marks[-1]=tuple(x-y for x,y in zip(marks[-1],v))
                marks.append(v)
            for levels in weak(n-target,n):
                assert scalar_power(1,target,tuple(marks),levels)==vertex(marks,levels)
            print('Genus-one target dimension',target,'on',n,'marks verified',flush=True)
            graph_terms.cache_clear()


def check_bounds():
    pairs=tuple((i,j) for i,j,m in LINKS)
    for degree in (3,4):
        values={}
        for keep in combinations_with_replacement(pairs,degree):
            if keep.count((1,4))>2: continue
            reflected=tuple(sorted((5-i,j) if (i,j) in ((2,4),(3,4)) else (i,j) for i,j in keep))
            if reflected in values:
                values[keep]=values[reflected]
            else:
                values[keep]=bound(keep)
        print('Best bound at degree',degree,'from',len(values),'choices:',
              min(values,key=values.get),min(values.values()),flush=True)
        expected = F(66530581,10146860236800) if degree==3 else F(1584203,684913065984)
        assert {keep for keep,v in values.items() if v==min(values.values())}=={
            ((2,4),)*degree, ((3,4),)*degree}
        assert min(values.values())==expected
        mixed={keep:value for keep,value in values.items() if len(set(keep))>1}
        print('Best with distinct retained pairs:',min(mixed,key=mixed.get),min(mixed.values()),flush=True)
        graph_terms.cache_clear()


if __name__=='__main__':
    check_vertices()
    check_bounds()
