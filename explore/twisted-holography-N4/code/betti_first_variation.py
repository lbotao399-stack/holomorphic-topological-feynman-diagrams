"""First variation of the scalar Abel height in a second independent target coordinate."""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from math import comb, factorial, prod

from flint import fmpq_mat
from betti_flow_polytope import charge_cuts, normal_coefficient, simplices, simplex_moment
from betti_power_graphs import LINKS, graph_terms, vector_coefficient, weighted_vertex


@lru_cache(maxsize=8192)
def facet_moment(powers,bad,cut):
    """Derivative in the specified capacity bound; all original cut normals have entries 0 or 1."""
    dim=len(powers)
    if len(bad)==1:
        (only,bound),=bad
        assert only==cut and len(cut)==dim
        degree=dim+sum(powers)
        return F(bound**(degree-1)*prod(map(factorial,powers)),factorial(degree-1))
    bound=dict(bad)[cut]
    forgotten=cut[0]
    total=F()
    for _,simplex in simplices(dim,bad):
        boundary=tuple(p for p in simplex if sum(p[i] for i in cut)==bound)
        if len(boundary)!=dim:
            continue
        first,*others=boundary
        jac=abs(F(str(fmpq_mat([[str(p[i]-first[i]) for i in range(dim) if i!=forgotten]
                                for p in others]).det())))
        total+=simplex_moment(powers,jac,boundary)
    return total


@lru_cache(maxsize=32768)
def normal_derivative(edges,base,tangent,powers):
    if not edges or not any(tangent): return F()
    primary=charge_cuts(edges,base)
    secondary=dict(charge_cuts(edges,tangent))
    # A bridge of zero primary charge first moves in energy a, where its
    # single-edge coefficient has derivative -p_B²/2 and no branch power.
    missing=[c[0] for c in secondary if len(c)==1 and c not in dict(primary)]
    if missing:
        if len(missing)>1: return F()
        i,=missing
        if powers[i]: return F()
        reduced=tuple((tuple(j-int(j>i) for j in c),s) for c,s in primary)
        return -F(secondary[(i,)]**2,2)*normal_coefficient(reduced,powers[:i]+powers[i+1:])
    all_edges=set(range(len(edges)))
    if not primary or set().union(*(set(c) for c,_ in primary))!=all_edges: return F()
    velocities={c:F(secondary[c]**2,2*s) for c,s in primary if c in secondary}
    if not velocities: return F()
    if len(velocities)==len(primary):
        ratios={v/dict(primary)[c] for c,v in velocities.items()}
        if len(ratios)==1:
            return 2*(len(edges)+sum(powers))*ratios.pop()*normal_coefficient(primary,powers)
    exponents=tuple(2*k+1 for k in powers)
    total=F()
    for size in range(1,len(primary)+1):
        for bad in combinations(primary,size):
            if set().union(*(set(c) for c,_ in bad))!=all_edges: continue
            total+=(-1)**size*sum(velocities[c]*facet_moment(exponents,bad,c)
                                  for c,_ in bad if c in velocities)
    return total/prod(2**k*factorial(k) for k in powers)


def charges(pair):
    return tuple(int(i==pair[0])-int(i==pair[1]) for i in range(5))


def stencil(d,p,q):
    terms=Counter()
    for g,f,k,w in graph_terms(d,p,4,q):
        terms[g,f,k]+=w
    for g,f,k,w in graph_terms(d-1,p,4):
        for i,v in enumerate(q):
            if v:
                shifted=k[:i]+(k[i]+1,)+k[i+1:]
                terms[g,f,shifted]+=w*F(v*v,2)
    return terms


def bound(base,fifth,d=5):
    remaining=Counter({(i,j):m for i,j,m in LINKS})
    remaining[base]-=d-1
    remaining[fifth]-=1
    assert min(remaining.values())>=0
    degrees=Counter({(0,)*5:1})
    for (i,j),m in remaining.items():
        expanded=Counter()
        for powers,weight in degrees.items():
            for k in range(m+1):
                v=list(powers)
                v[i]+=k
                v[j]+=m-k
                expanded[tuple(v)]+=weight*comb(m,k)
        degrees=expanded
    total=F()
    for (groups,flags,extra),weight in stencil(d,charges(base),charges(fifth)).items():
        if not weight: continue
        for powers,mu in degrees.items():
            shifted=tuple(x+y for x,y in zip(extra,powers))
            for group,branch in zip(groups,flags):
                mu*=weighted_vertex(group,branch,shifted)
                if not mu: break
            total+=weight*mu
    return 2*factorial(d-1)*total/(2**(14-d)*prod(factorial(m) for _,_,m in LINKS))


def check():
    from betti_power_graphs import shapes
    from betti_flow_polytope import weak
    tested=0
    for e in range(1,5):
        for v,edges,*_ in shapes(e):
            p=tuple(range(1,v)) + (-v*(v-1)//2,)
            q=(1,)+(0,)*(v-2)+(-1,)
            for normal in weak(4-e,e):
                # Degree d in t. Values at t=j² give actual integer vector currents.
                d=e+sum(normal)
                values=[vector_coefficient(edges,tuple(zip(p,(j*x for x in q))),normal) for j in range(d+1)]
                weights=[sum(F(-1,j*j) for j in range(1,d+1))]
                weights.extend(F(1,j*j)*prod(F(-s*s,j*j-s*s) for s in range(1,d+1) if s!=j)
                               for j in range(1,d+1))
                assert normal_derivative(edges,p,q,normal)==sum(a*b for a,b in zip(weights,values)),(edges,normal)
                tested+=1
    print('Independent vector polynomial derivative tests:',tested,flush=True)
    from betti_mixed_cubic import bound as polarized
    for p,q in (((0,1),(1,4)),((2,4),(3,4))):
        assert bound(p,q,4)==polarized((p,p,p,q))
    for p in ((2,4),(0,1)):
        for q in ((0,1),(1,4),(2,4),(3,4)):
            if p==q or (p==(0,1) and q==(3,4)): continue
            print('Five retained',p,q,bound(p,q),flush=True)


if __name__=='__main__':
    check()
