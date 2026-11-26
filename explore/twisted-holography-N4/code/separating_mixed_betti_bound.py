"""One separating trace after four copies of a Betti form and one independent copy."""
from collections import Counter
from fractions import Fraction as F
from functools import cache
from itertools import combinations, product
from math import comb, factorial, prod

from betti_first_variation import charges, stencil
from betti_power_graphs import LINKS, scalar_power, weighted_vertex
from separating_betti_bound import component


@cache
def mixed_rows(h,p,q):
    return tuple((groups,flags,extra,w) for (groups,flags,extra),w in
                 stencil(h+1,p,q).items() if w)


@cache
def mixed_component(h,p,q,levels):
    if not h or not any(q) or not any(p) or q==p or q==tuple(-x for x in p):
        return F()
    total=F()
    for groups,flags,extra,w in mixed_rows(h,p,q):
        powers=tuple(x+y for x,y in zip(levels,extra))
        for group,branch in zip(groups,flags):
            w*=weighted_vertex(group,branch,powers)
            if not w: break
        total+=w
    return total


def boundary(levels,p,q,cut):
    n=len(levels)
    free=tuple(i for i in range(n) if i not in cut)
    total=F()
    for size in range(len(free)+1):
        for group in combinations(free,size):
            left=(cut[0],*group)
            sides=(left,tuple(i for i in range(n) if i not in left))
            local=[]
            for side in sides:
                triples=[(p[i],q[i],levels[i]) for i in side]
                triples.append((-sum(p[i] for i in side),-sum(q[i] for i in side),0))
                local.append(tuple(zip(*sorted(triples))))
            for h in range(5):
                genera=(h,4-h)
                for active in (0,1):
                    if not all(sum(k)+g+int(i==active)==3*g-3+len(k)
                               for i,(g,(_,_,k)) in enumerate(zip(genera,local))):
                        continue
                    j=1-active
                    pa,qa,ka=local[active]
                    pj,_,kj=local[j]
                    passive=component(genera[j],genera[j],pj,kj,False)
                    if passive:
                        total+=passive*mixed_component(genera[active],pa,qa,ka)
    return total


def corrected(base,fifth,previous):
    remaining=Counter({(i,j):m for i,j,m in LINKS})
    remaining[base]-=4
    remaining[fifth]-=1
    denominator=2**9*prod(factorial(m) for _,_,m in LINKS)
    result={}
    for cut,multiple in remaining.items():
        if not multiple: continue
        counts=remaining.copy()
        counts[cut]-=1
        correction=F()
        for choices in product(*(range(m+1) for m in counts.values())):
            k=[0]*5
            weight=1
            for ((i,j),m),r in zip(counts.items(),choices):
                k[i]+=r
                k[j]+=m-r
                weight*=comb(m,r)
            correction+=weight*boundary(k,charges(base),charges(fifth),cut)
        result[cut]=previous-F(2*factorial(4),denominator)*correction
        assert 0<result[cut]<=previous
    return result


def check():
    from betti_flow_polytope import weak
    from betti_first_variation import bound
    # Finite-difference polarization of the unexpanded vector height provides an
    # independent control of each component function needed below through genus three.
    for h in (1,2,3):
        n=4
        p=(1,0,0,-1)
        q=(0,1,0,-1)
        d=h+1
        weights=[sum(F(-1,j*j) for j in range(1,d+1))]+[
            F(1,j*j)*prod(F(-s*s,j*j-s*s) for s in range(1,d+1) if s!=j) for j in range(1,d+1)]
        for k in list(weak(3*h-3+n-d,n))[::3]:
            independent=sum(w*scalar_power(h,d,tuple(zip(p,(j*x for x in q))),k)
                            for j,w in enumerate(weights))
            assert independent==mixed_component(h,p,q,k)
    previous={((0,1),(1,4)):F(171793,152882380800),
              ((0,1),(2,4)):F(3629,1651507200),
              ((2,4),(0,1)):F(3226883,1712282664960),
              ((2,4),(1,4)):F(12053471,8561413324800),
              ((2,4),(3,4)):F(2102297,1426902220800)}
    assert bound((0,1),(1,4))==previous[(0,1),(1,4)]
    results={}
    for (p,q),before in previous.items():
        for cut,value in corrected(p,q,before).items():
            results[p,q,cut]=value
            print('Fourfold pair, fifth, separating pair, bound:',p,q,cut,value,flush=True)
    print('Overall minimum:',min(results,key=results.get),min(results.values()),flush=True)
    best=F(113203,107017666560)
    assert len(results)==15 and min(results.values())==best
    assert {key for key,w in results.items() if w==best}=={
        ((0,1),(1,4),(2,4)),((0,1),(1,4),(3,4))}


if __name__=='__main__':
    check()
