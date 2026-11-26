"""Compactified closed-string ring contractions and marked threefold kernels."""
from fractions import Fraction as F
from itertools import permutations,combinations
from math import factorial,prod
from collections import Counter
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'kodaira-spencer'))
from compactification import torus_product,k3_product,cubic


def glued_ring(multiply,top,parity):
    labels=range(top+1)
    dual=[]
    handle=[0]*(top+1)
    for i in labels:
        j,= [j for j in labels if cubic(multiply,top,i,j,0)]
        dual.append((j,cubic(multiply,top,i,j,0)))
        sign,k=multiply(i,j)
        handle[k]+=(-1)**parity(i)*dual[-1][1]*sign
    # Expanding each internal tree edge in dual basis reproduces multiplication.
    for a in labels:
        for b in labels:
            ans=multiply(a,b)
            for c in labels:
                for d in (0,top):
                    result=sum(cubic(multiply,top,a,b,i)*sg*cubic(multiply,top,j,c,d)
                               for i,(j,sg) in enumerate(dual))
                    # This ordering places the dual of the first contraction first
                    # in the second trace. The even external tests used here
                    # need no additional odd interchange.
                    if any(parity(x) for x in (a,b,c,d)):
                        continue
                    assert result==ans[0]*cubic(multiply,top,ans[1],c,d)
    assert not any(handle[:-1])
    assert multiply(top,top)[0]==0
    return handle[-1]


def square():
    from theta_periods import blocks

    sector=[]
    for shape in ((1,3,3),(2,2,3)):
        mat=blocks(shape)
        raw=sum((mat[p[2],p[j]]-mat[p[3],p[j]])*charge for p in permutations(range(4))
                for j,charge in enumerate((1,1,0,-2)))
        sector.append(F(raw)/4)
    assert sector==[F(1,10080),F(29,272160)]
    return sum(sector)


def third_derivative():
    total=F()
    for shape in ((1,3,3),(2,2,3)):
        ds=[4-length for length in shape]
        mass=F(prod(factorial(k) for k in ds),factorial(2+sum(ds)))
        for a in combinations(range(4),shape[0]-1):
            for b in combinations(sorted(set(range(4))-set(a)),shape[1]-1):
                groups=[a,b,sorted(set(range(4))-set(a)-set(b))]
                h={i:j for j,g in enumerate(groups) for i in g}
                if h[0]==h[1] or h[2]==h[3]:
                    continue
                order=prod(factorial(m) for m in Counter(shape).values())
                along=h[3]
                total+=mass*(F(7,48)+F((2+sum(ds))*(shape[along]-2),120*ds[along]))/order
    return total


def extremal(n):
    return F(1,2**n*factorial(n)*factorial(2*n+1))


if __name__=='__main__':
    hk=glued_ring(k3_product,23,lambda i:0)
    ht=glued_ring(torus_product,15,lambda i:i.bit_count()%2)
    assert (hk,ht)==(24,0)
    two=square()
    assert (hk*two,ht*two)==(F(2,405),0)
    print('Scalar-square two-loop coefficients on unit labels:',hk*two,ht*two)
    assert third_derivative()==extremal(2)==F(1,960)
    print('Tree after removing the mark, loops one to four:',[extremal(n) for n in range(1,5)])
