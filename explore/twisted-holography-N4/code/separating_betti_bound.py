"""Retain a separating trace after four actual single-target Betti forms in the D3 bound."""
from collections import Counter
from fractions import Fraction as F
from functools import cache
from itertools import combinations, product
from math import comb, factorial, prod

from betti_flow_polytope import weak
from betti_power_graphs import retained_bound, scalar_power
from soft_quartic_psi_bound import cut_intersection, intersection, upper_bound
from soft_quartic_two_betti import LINKS, add_at, four_times_betti_pair


@cache
def component(genus, power, charges, levels, independent):
    if not power:
        return intersection(genus,levels)
    if not any(charges):
        return F()
    if not independent:
        return scalar_power(genus,power,charges,levels)
    first,second = [i for i,v in enumerate(charges) if v]
    if power==1:
        return (sum(intersection(genus,add_at(levels,i)) for i in (first,second))
                -cut_intersection(genus,levels,first,second,False))/2
    if power==2:
        return four_times_betti_pair(genus,levels,(first,second),(first,second),False)/8
    raise ValueError('The divisor and two-Betti control supports degree at most two.')


def boundary(genus, degree, charges, levels, pair, independent=False):
    total = F()
    spare = tuple(i for i in range(len(levels)) if i not in pair)
    for count in range(len(spare)+1):
        for subset in combinations(spare,count):
            left = (pair[0],*subset)
            right = tuple(i for i in range(len(levels)) if i not in left)
            local = []
            for side in (left,right):
                q=tuple(charges[i] for i in side)
                local.append(((*q,-sum(q)),tuple(levels[i] for i in side)+(0,)))
            for h in range(genus+1):
                gside=(h,genus-h)
                powers=tuple(3*g-3+len(k)-sum(k) for g,(_,k) in zip(gside,local))
                if sum(powers)!=degree or any(d<0 or d>g for d,g in zip(powers,gside)):
                    continue
                total+=prod(component(g,d,q,k,independent) for g,d,(q,k) in zip(gside,powers,local))
    return total


def bounds(degree=4, retained=(2,4)):
    genus=4
    charges=tuple(int(i==retained[0])-int(i==retained[1]) for i in range(5))
    multiplicities=Counter({(i,j):m for i,j,m in LINKS})
    multiplicities[retained]-=degree
    normalization=F(2*factorial(degree),2**(14-degree)*prod(factorial(m) for _,_,m in LINKS))
    results={}
    for pair,m in multiplicities.items():
        if not m: continue
        counts=multiplicities.copy()
        counts[pair]-=1
        naive=correction=F()
        for choice in product(*(range(k+1) for k in counts.values())):
            levels=[0]*5
            weight=1
            for ((i,j),k),x in zip(counts.items(),choice):
                levels[i]+=x
                levels[j]+=k-x
                weight*=comb(k,x)
            naive+=weight*sum(component(genus,degree,charges,tuple(add_at(levels,i)),False) for i in pair)
            correction+=weight*boundary(genus,degree,charges,levels,pair)
        assert correction>=0
        results[pair]=normalization*naive,normalization*correction,normalization*(naive-correction)
    assert len({row[0] for row in results.values()})==1
    return results


def check():
    # Degree zero reproduces the direct one-divisor computation with the
    # separating genus sum oriented by the two distinguished marks.
    for pair,(naive,correction,bound) in bounds(0).items():
        links=tuple(sorted(LINKS,key=lambda l:l[:2]!=pair))
        assert (naive,bound)==upper_bound(links,False)
    # In degree two the independent compact square contains the extra
    # two-parallel-edge term, including when both charged marks lie on one side.
    tested=0
    charges=(0,1,-1,0)
    for g in (2,3,4):
        for levels in weak(3*g-3+4-3,4):
            for cut in ((0,1),(0,3),(1,2)):
                value=boundary(g,2,charges,levels,cut)
                assert value==boundary(g,2,charges,levels,cut,True)
                tested+=1
    print('Independent degree-zero and',tested,'degree-two separating tests passed.',flush=True)
    all_values={}
    for retained in ((0,1),(2,4),(3,4)):
        result=bounds(4,retained)
        previous=retained_bound(retained,4)
        for pair,row in result.items():
            assert row[0]==previous
            all_values[retained,pair]=row[-1]
            print('Retained pair, fifth pair, previous bound, correction, improved bound:',
                  retained,pair,*row,flush=True)
    best=F(1172347,570760888320)
    assert min(all_values.values())==best
    assert {k for k,v in all_values.items() if v==best}=={
        ((2,4),(1,4)),((3,4),(1,4))}


if __name__=='__main__':
    check()
