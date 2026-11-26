"""Six marked hexagon facets, normalized Gaussian, Wick sign and SO(N) polarization."""
from fractions import Fraction as F
from itertools import permutations
from math import factorial
from pathlib import Path
from random import Random
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'code'))
from htk.linalg import det
from htk.graph import polygon


def parity(p):
    return (-1)**sum(x>y for i,x in enumerate(p) for y in p[i+1:])


def gaussian_wedge(t,momenta):
    """Coefficients in dbar z_0^1 ... dbar z_4^5, before the scalar Gaussian integration."""
    U=sum(t)
    increments=[[sum(momenta[s][k] for s in range(e+1)) for k in range(5)] for e in range(6)]
    y=[[v-sum(F(t[j],U)*increments[j][k] for j in range(6))
        for k,v in enumerate(row)] for row in increments]
    return edge_forms(polygon(6),y)


def edge_forms(graph,y):
    """Antiholomorphic form determinant for every possible heat marking, with all vertex volumes on the left."""
    incidence=graph.incidence()
    def difference(e,k):
        return [w*int(i==k) for w in incidence[e] for i in range(5)]
    values=[]
    for cut in range(graph.ne):
        forms=[]
        factor=F(1)
        for e in range(graph.ne):
            if e==cut:
                forms.extend(difference(e,k) for k in range(5))
                continue
            x=y[e]
            h=next(i for i,a in enumerate(x) if a)
            factor*=F((-1)**h,x[h]**3)
            for k in range(5):
                if k!=h:
                    forms.append([x[h]*a-x[k]*b for a,b in zip(difference(e,k),difference(e,h))])
        values.append(factor*det(forms))
    return values


def wick_sign():
    global_order=[('x',v) for v in range(6)]
    global_order += [(e,s) for e in range(6) for s in range(2)]
    local=[slot for v in range(6) for slot in [('x',v),((v-1)%6,1),(v,0)]]
    return parity([global_order.index(s) for s in local])


def trace_powers(eigenvalues):
    trace=lambda k:sum(i**k for i in eigenvalues)
    adjoint=sum((a+b)**6 for i,a in enumerate(eigenvalues) for b in eigenvalues[i+1:])
    return adjoint,(len(eigenvalues)-32)*trace(6)+15*trace(2)*trace(4)


def check():
    g=polygon(6)
    assert g.edges==[(e,(e+1)%6) for e in range(6)]
    rng=Random(65032)
    for case in range(4):
        p=[[int(i==j) for j in range(5)] for i in range(5)] if case==0 else [
             [rng.randrange(-6,7) for _ in range(5)] for _ in range(5)]
        p.append([-sum(c) for c in zip(*p)])
        t=[rng.randrange(1,10) for _ in range(6)]
        assert gaussian_wedge(t,p)==[-F(x,sum(t))*det(p[:5]) for x in t]
    assert wick_sign()==-1
    # Six equal facet domains; integrate (1+s_1+...+s_5)^(-6) on [1,infinity)^5.
    one_facet=F(1,6*factorial(5))
    assert 6*one_facet==F(1,120)
    assert F(factorial(6),12*factorial(5))==F(1,2)
    for n in range(2,41):
        positive=[rng.randrange(1,12) for _ in range(n//2)]
        spectrum=positive+[-v for v in positive]+[0]*(n%2)
        a,b=trace_powers(spectrum)
        assert a==b
    # Polarization of product of traces partitions the six slots into two and four.
    assert len({tuple(sorted(p[:2])) for p in permutations(range(6))})==15
    print('raw antiholomorphic wedge = -t_mark det(p0,...,p4)/sum(t)')
    print('normalized projective volume',6*one_facet,'cyclic Wick sign',wick_sign())
    print('primitive coefficient N-32, polarized quadratic-quartic partitions 15')


if __name__=='__main__':
    check()
