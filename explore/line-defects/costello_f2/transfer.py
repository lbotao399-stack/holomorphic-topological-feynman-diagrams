"""Polynomial replay of the rotation-weight polydisc contraction used in the PBW projection."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations, product


def apply(forms, operation):
    result = defaultdict(F)
    for monomial,value in forms.items():
        for key,factor in operation(*monomial):
            result[key] += value*factor
    return {k:v for k,v in result.items() if v}


def dbar(z,bar,slots):
    for i,n in enumerate(bar):
        if n and i not in slots:
            b = list(bar)
            b[i] -= 1
            sign = (-1)**sum(j<i for j in slots)
            yield (z,tuple(b),tuple(sorted((*slots,i)))),sign*n


def homotopy(z,bar,slots):
    # The radial integral on a monomial is exactly 1/(|beta|+q).
    for j,i in enumerate(slots):
        b = list(bar)
        b[i] += 1
        yield (z,tuple(b),slots[:j]+slots[j+1:]),F((-1)**j,sum(bar)+len(slots))


def dz(z,bar,slots):
    if z[0]:
        yield ((z[0]-1,*z[1:]),bar,slots),z[0]


def verify():
    for z in product(range(3),repeat=2):
        for bar in product(range(3),repeat=2):
            for q in range(3):
                for slots in combinations(range(2),q):
                    form = {(z,bar,slots):F(1)}
                    a = apply(apply(form,homotopy),dbar)
                    b = apply(apply(form,dbar),homotopy)
                    actual = {k:a.get(k,0)+b.get(k,0) for k in a.keys()|b.keys()}
                    assert {k:v for k,v in actual.items() if v}==(form if q+sum(bar) else {})
                    assert apply(apply(form,homotopy),dz)==apply(apply(form,dz),homotopy)
    print('324 monomial forms verify the homotopy identity and holomorphic derivative commutation.')


if __name__=='__main__':
    verify()
