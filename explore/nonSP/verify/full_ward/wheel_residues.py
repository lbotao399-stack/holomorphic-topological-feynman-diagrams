"""All channels of an integer-power wheel period from one univariate partial fraction."""
from collections import Counter
from fractions import Fraction
from functools import cache
from math import comb,factorial,prod
from operator import index

from flint import fmpq,fmpq_poly
from wheel_angular import multiplicity,radial_master,X


def indices(spoke,rim,q):
    spoke=tuple(map(index,spoke))
    rim=tuple(map(index,rim))
    q=index(q)
    n=len(rim)
    if len(spoke)!=n or n<3 or sum(spoke)+sum(rim)!=n*q:
        raise ValueError('Expected n spoke and n rim powers of total n*q.')
    if any(not 1<=v<q for v in (*spoke,*rim)) or sum(rim)<=q:
        raise ValueError('Divergent projective wheel powers.')
    c=0
    nodes=[]
    for i in range(n):
        c+=2*spoke[i]+rim[i-1]+rim[i]-2*q
        d=rim[i]-1
        nodes.extend(range(c-d,c+d+1,2))
    assert c==0
    minimum=min(nodes)
    assert all((v-minimum)%2==0 for v in nodes)
    counts=Counter((v-minimum)//2 for v in nodes)
    length=max(counts)+1
    if length>=q:
        raise ValueError('Divergent wheel radial gap.')
    histogram=tuple(counts[j] for j in range(length))
    factor=prod(factorial(q-p-1)*factorial(r-1) for p,r in zip(spoke,rim))
    return min(histogram,histogram[::-1]),factor


@cache
def summation(q,h,u):
    """Summing M(s)/(s-h)^u from s=q-1; polynomial parts cancel in the complete sum."""
    coefficients=multiplicity(q)(X+h)
    zeta=[fmpq(0) for _ in range(u+1)]
    for j in range(u):
        weight=u-j
        zeta[weight]=coefficients[j]
        zeta[0]-=coefficients[j]*sum((fmpq(1,k**weight) for k in range(1,q-1-h)),fmpq(0))
    return tuple(zeta)


@cache
def period(q,histogram):
    total=sum(histogram)
    parts=[]
    for a,n in enumerate(histogram):
        if not n:
            continue
        # Taylor expand the regular part of P at a, where P(w)=prod_j (w-j)^(-N_j).
        denominator=prod(((X+a-j)**m for j,m in enumerate(histogram) if j!=a and m),
                         start=fmpq_poly([1]))
        coeff=[]
        for i in range(n):
            value=(int(i==0)-sum(denominator[j]*coeff[i-j] for j in range(1,i+1)))/denominator[0]
            coeff.append(value)
            parts.append((a,n-i,value))
    pairs=Counter()
    for a,k,x in parts:
        for b,l,y in parts:
            pairs[a-b,k+l-1]+=(-1)**(total+l)*comb(k+l-2,k-1)*x*y
    out=Counter()
    for (h,u),value in pairs.items():
        if value:
            for w,c in enumerate(summation(q,h,u)):
                if c:
                    out[w]+=c*value
    assert all(not v or w==0 or w>=3 and w%2 for w,v in out.items())
    return tuple((w,Fraction(str(v))) for w,v in sorted(out.items()) if v)


def integrate(spoke,rim,q):
    histogram,factor=indices(spoke,rim,q)
    values={w:factor*c for w,c in period(q,histogram)}
    assert all(w<=2*len(spoke)-3 for w in values)
    return values


def integrate_terms(terms,q):
    """Sum (powers, coefficient) rows, spoke powers followed by cyclic rim powers."""
    by_histogram=Counter()
    for powers,coefficient in terms:
        if coefficient:
            n=len(powers)//2
            h,scale=indices(powers[:n],powers[n:],q)
            by_histogram[h]+=Fraction(coefficient)*scale
    result=Counter()
    for h,coefficient in by_histogram.items():
        if coefficient:
            for weight,value in period(q,h):
                result[weight]+=value*coefficient
    return dict(sorted((w,v) for w,v in result.items() if v))


if __name__=='__main__':
    from wheel_angular import cut_rows,integrate as by_sectors
    for n in range(3,11):
        for q in range(2,n):
            assert integrate([q-1]*n,[1]*n,q)==radial_master(n,q)
    for n in range(3,13):
        spoke=[1]+[n+1]*(n-1)
        rim=[2]*n
        h,_=indices(spoke,rim,n+2)
        assert h==(1,)+(2,)*(n-1)+(1,)
        sharp=integrate(spoke,rim,n+2)
        assert list(sharp)==[0] and sharp[0]>0
    for spoke,rim,vals in cut_rows():
        actual=integrate(spoke,rim,6)
        assert [actual.get(w,0) for w in (0,3,5)]==vals
    for spoke,rim in [([5,5,4,5,5],[3,2,2,2,2]),([6,4,5,5,4],[2,2,1,3,3]),
                      ([3,5,5,5,3],[3,3,4,2,2]),([4,4,4,4,4],[3,3,3,3,3])]:
        assert integrate(spoke,rim,7)==by_sectors(spoke,rim,7)
    raw=[((5,5,4,5,5,3,2,2,2,2),Fraction(2,3)),((4,4,4,4,4,3,3,3,3,3),-2)]
    weighted=Counter()
    for powers,c in raw:
        for w,v in integrate(powers[:5],powers[5:],7).items():
            weighted[w]+=c*v
    assert dict(weighted)==integrate_terms(raw,7)
    print('All channels agree for all 1421 four-loop cuts, four five-loop sector integrals and radial masters through ten loops.')
