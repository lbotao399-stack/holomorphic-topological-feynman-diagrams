"""Exact one-interior cyclic disk, including all descendant powers and the divergence term."""
from collections import Counter
from fractions import Fraction as F
from functools import cache
from itertools import combinations, product
from math import factorial


@cache
def expand(arguments):
    result=Counter()
    for terms in product(*arguments):
        s=tuple(i for i,_ in terms)
        if len(set(s))!=len(s):
            continue
        coefficient=(-1)**sum(i>j for index,i in enumerate(s) for j in s[index+1:])
        for _,sign in terms:
            coefficient*=sign
        result[tuple(sorted(s))]+=coefficient
    return {key:value for key,value in result.items() if value}


def matching(m,slots):
    """Coefficient in (sum dtheta_i wedge dtheta_{i+1})^ell divided by ell!."""
    remaining=sorted(set(range(m))-set(slots))
    return all(b==a+1 for a,b in zip(remaining[::2],remaining[1::2]))


@cache
def cyclic_average(k,ell):
    m=k+2*ell
    # Ordered slots 0,...,m-1 after the spectator. Label -1 is the incoming closed momentum.
    args=[((i,1),) for i in range(m)]
    spectator=tuple((i,-1) for i in range(-1,m))
    all_inputs=[spectator]+args
    result=Counter()
    coefficient=F(factorial(ell),factorial(m)*(m+1))
    for selected in combinations(range(m),k):
        if not matching(m,selected):
            continue
        for rotation in range(m+1):
            ordered=tuple(all_inputs[(1+j+rotation)%(m+1)] for j in selected)
            factor=coefficient*(-1)**(m*rotation)
            for word,sign in expand(ordered).items():result[word]+=factor*sign
    return {key:value for key,value in result.items() if value}


@cache
def operation(k,ell):
    m=k+2*ell
    unit=F(factorial(ell),factorial(m))
    symbol={}
    for slots in combinations(range(m),k):
        if matching(m,slots):
            symbol[slots]=unit
    if k:
        for slots in combinations(range(m),k-1):
            weight=sum((-1)**sum(s<j for s in slots) for j in range(m)
                       if j not in slots and matching(m,(*slots,j)))
            if weight:
                symbol[(-1,*slots)]=unit*F(weight,m+1)
    return symbol


def substitute(symbol,arguments):
    out=Counter()
    for word,value in symbol.items():
        for key,c in expand(tuple(((-1,1),) if i==-1 else arguments[i] for i in word)).items():
            out[key]+=value*c
    return {key:value for key,value in out.items() if value}


def hochschild(k,ell):
    m=k+2*ell
    args=[((i,1),) for i in range(m+1)]
    symbol=operation(k,ell)
    out=Counter(substitute(symbol,args[1:]))
    for i in range(m):
        values=substitute(symbol,args[:i]+[args[i]+args[i+1]]+args[i+2:])
        for key,value in values.items():out[key]+=(-1)**(i+1)*value
    for key,value in substitute(symbol,args[:-1]).items():out[key]+=(-1)**(m+1)*value
    return {key:value for key,value in out.items() if value}


def verify(k,ell):
    assert operation(k,ell)==cyclic_average(k,ell)
    left=hochschild(k,ell)
    if not k:
        assert not left
        return
    right={(-1,*key):value for key,value in operation(k-1,ell+1).items() if -1 not in key}
    ratios={left.get(key,0)/value for key,value in right.items()}
    assert set(left)==set(right)
    assert ratios=={1}


if __name__=='__main__':
    for ell in range(4):
        for k in range(6):
            verify(k,ell)
    print("24 primary/descendant symbols agree with the cyclic projection and the linear disk Ward identity")
