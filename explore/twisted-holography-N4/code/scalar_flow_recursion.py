"""Recover the scalar graph polynomial from second derivatives along a spanning tree."""
from fractions import Fraction
from functools import lru_cache
from time import perf_counter

from flint import fmpq_mpoly_ctx


@lru_cache(maxsize=4096)
def tree_coordinates(n,edges):
    """For each tree edge return its index and the vertex charge sum on its tail side."""
    parts=list(range(n))
    tree=[]
    for i,(u,v,_) in enumerate(edges):
        if parts[u]==parts[v]:
            continue
        old,new=parts[v],parts[u]
        parts=[new if part==old else part for part in parts]
        tree.append(i)
    if len(tree)!=n-1: return ()
    cuts=[]
    for i in tree:
        reach={edges[i][0]}
        for _ in range(n-1):
            for j in tree:
                if j==i: continue
                u,v,_=edges[j]
                if u in reach or v in reach: reach.update((u,v))
        cuts.append((i,tuple(sorted(reach))))
    return tuple(cuts)


def ring(n):
    return fmpq_mpoly_ctx.get(tuple('s'+str(i) for i in range(max(1,n-1))))


def substitute(edges,charges,ctx):
    n=len(charges)
    if not edges: return ctx.constant(int(n==1))
    coordinates=tree_coordinates(n,edges)
    if len(coordinates)!=n-1 or any(u==v for u,v,_ in edges): return ctx.constant(0)
    p=graph_polynomial(n,edges)
    return p.compose(*(sum((charges[j] for j in side),ctx.constant(0)) for _,side in coordinates),ctx=ctx)


@lru_cache(maxsize=4096)
def graph_polynomial(n,edges):
    """Use oriented tree currents as independent variables for conserved charge space."""
    ctx=ring(n)
    if not edges: return ctx.constant(int(n==1))
    tree=tree_coordinates(n,edges)
    if len(tree)!=n-1 or any(u==v for u,v,_ in edges): return ctx.constant(0)
    charges=[ctx.constant(0) for _ in range(n)]
    for x,(edge,_) in zip(ctx.gens(),tree):
        u,v,_=edges[edge]
        charges[u]+=x
        charges[v]-=x
    terms={}
    for j,(i,_) in enumerate(tree):
        u,v,m=edges[i]
        deleted=edges[:i]+edges[i+1:]
        if m==1:
            rhs=substitute(deleted,charges,ctx)
            contracted=tuple(k-int(k>v) if k!=v else u for k in range(n))
            new_edges=tuple(sorted((*sorted((contracted[a],contracted[b])),k) for a,b,k in deleted))
            new_charges=[ctx.constant(0) for _ in range(n-1)]
            for a,p in zip(contracted,charges):
                new_charges[a]+=p
            rhs-=substitute(new_edges,new_charges,ctx)
        else:
            lowered=tuple(sorted((*deleted,(u,v,m-1))))
            rhs=(2*m-1)*substitute(lowered,charges,ctx)
        for powers,c in rhs.terms():
            # Every degree-2d monomial has a first variable occurring twice.
            if any(k>=2 for k in powers[:j]): continue
            k=powers[j]
            key=(*powers[:j],k+2,*powers[j+1:])
            terms[key]=c/((k+1)*(k+2))
    return ctx.from_dict(terms)


@lru_cache(maxsize=32768)
def coefficient(edges,charges,powers):
    if sum(charges) or len(edges)!=len(powers) or any(k<0 for k in powers):
        raise ValueError('Supply conserved charges and one nonnegative branch power per edge.')
    decorated=tuple(sorted((*sorted((u,v)),k+1) for (u,v),k in zip(edges,powers)))
    n=len(charges)
    if n==1: return Fraction(int(not edges))
    tree=tree_coordinates(n,decorated)
    if len(tree)!=n-1: raise ValueError('Use a connected graph.')
    x=tuple(sum(charges[j] for j in side) for _,side in tree)
    return Fraction(str(graph_polynomial(n,decorated)(*x)))


def check():
    from betti_flow_polytope import charge_cuts,normal_coefficient,weak
    from betti_power_graphs import shapes
    count=0
    for d in range(1,5):
        for e in range(1,d+1):
            for n,edges,*_ in shapes(e):
                p=(*range(1,n),-n*(n-1)//2)
                for powers in weak(d-e,e):
                    assert coefficient(edges,p,powers)==normal_coefficient(charge_cuts(edges,p),powers)
                    count+=1
    print('Exact cross-checks on positive and negative cut contributions:',count,flush=True)
    from itertools import combinations
    for n in (4,5):
        edges=tuple(combinations(range(n),2))
        p=(*range(1,n),-n*(n-1)//2)
        start=perf_counter()
        actual=coefficient(edges,p,(0,)*len(edges))
        if n==4: assert actual==normal_coefficient(charge_cuts(edges,p),(0,)*len(edges))
        print('Complete graph',n,'value',actual,'elapsed',round(perf_counter()-start,3),flush=True)
    from betti_power_graphs import scalar_power
    from one_target_dr import homogeneous_dr
    count=0
    for genus in (5,6):
        for p in ((1,-1),(2,-3,1)):
            expected=homogeneous_dr(genus,p).to_dict()
            for levels in weak(2*genus-3+len(p),len(p)):
                assert scalar_power(genus,genus,p,levels)==Fraction(str(expected.get(levels,0)))
                count+=1
    assert count==179
    print('Genus five and six descendant comparisons:',count,flush=True)
    # Tadpoles and changes of vertex labels, including a decorated bridge.
    assert coefficient(((0,1),(1,1)),(2,-2),(1,0))==0
    edges=((0,1),(1,2),(1,2),(0,2),(2,3))
    p=(1,-3,4,-2)
    powers=(1,0,2,0,1)
    actual=coefficient(edges,p,powers)
    perm=(2,0,3,1)
    renamed=tuple((perm[u],perm[v]) for u,v in edges)
    assert actual==coefficient(renamed,tuple(p[perm.index(i)] for i in range(4)),powers)


if __name__=='__main__':
    check()
