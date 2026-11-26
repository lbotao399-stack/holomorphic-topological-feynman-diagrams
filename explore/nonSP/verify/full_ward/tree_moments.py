"""Exact polarization of the outer graph integral in its spanning-tree weights."""
from functools import cache
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement, product
from math import comb, factorial, lcm, prod
import numpy as np
from flint import fmpq_mat, fmpz_mat

from lower_faces import (G, Graph, NotSeriesParallel, EDGES, checked_cut, orbit,
                         plan_generic, skeleton_chains, sp_monomial, tree_flows)


def compositions(total,parts):
    if parts==1:
        yield (total,)
    else:
        for i in range(total+1):
            for rest in compositions(total-i,parts-1):
                yield (i,)+rest


def subsets_with_multiplicity(indices):
    from collections import Counter
    counts=sorted(Counter(indices).items())
    for ex in product(*(range(k+1) for _,k in counts)):
        key=tuple(i for (i,_),n in zip(counts,ex) for _ in range(n))
        w=prod(comb(k,n) for (_,k),n in zip(counts,ex))
        yield key,(-1)**(len(indices)-len(key))*w


def monomial_integrator(graph,q):
    try:
        plan=plan_generic(Graph(graph.nv,graph.edges))
    except NotSeriesParallel:
        plan=None
    if plan is not None:
        return cache(lambda ex: (sp_monomial(ex,q,plan),F()))
    cubic,chains=skeleton_chains(graph.edges,graph.nv)
    paths={tuple(sorted((cubic.index(v[0]),cubic.index(v[-1])))):es for v,es in chains}
    assert set(paths)==set(EDGES)
    @cache
    def integrate(ex):
        powers=tuple(sum(1+ex[e] for e in paths[p]) for p in EDGES)
        weight=F(prod(map(factorial,ex)),factorial(q-1))
        return tuple(weight*c for c in checked_cut(orbit(powers)))
    return integrate


def bare_table(graph, trees):
    """Bare Schwinger numerator on a unisolvent lattice of source components."""
    m=graph.n
    axis=np.array(list(compositions(m,m+1)),dtype=np.int64)
    flows=np.array(list(trees.values()),dtype=np.int64)
    cycles=np.array(graph.C,dtype=np.int64)
    keys=list(trees)
    minors=[]
    for subset in combinations(range(2*m),m):
        rest=tuple(j for j in range(2*m) if j not in subset)
        sign=(-1)**sum(i-j for j,i in enumerate(subset))
        coef=sign*int(fmpz_mat(cycles[:,np.array(subset)+1].tolist()).det())
        coef*=int(fmpz_mat(cycles[:,np.array(rest)+1].tolist()).det())
        if coef:
            minors.append((np.array(subset)+1,np.array(rest)+1,coef))
    shape=(m,)*(2*m)
    values=np.empty(shape+(len(axis),len(axis)),dtype=np.int64)
    for idx in np.ndindex(shape):
        t=np.array((1,)+tuple(v+1 for v in idx))
        weights=np.prod(t[np.array(keys)],axis=1)
        u=int(weights.sum())
        currents=(np.tensordot(weights,flows,axes=1) @ axis.T).T
        det=0
        # Each m-current product fits int64; the two blocks can require more.
        bound=int(np.max(np.abs(currents)))**m
        assert bound<2**63
        dtype=object if bound**2*sum(abs(c) for _,_,c in minors)>=2**63 else np.int64
        for cols,other,c in minors:
            x=np.prod(currents[:,cols],axis=1).astype(dtype)
            y=np.prod(currents[:,other],axis=1).astype(dtype)
            det+=c*np.outer(x,y)
        assert not np.any(det % u**m)
        values[idx]=det//u**m
    vand=fmpz_mat([[v**k for k in range(m)] for v in range(1,m+1)]).inv()
    d=lcm(*(int(x.q) for x in vand.entries()))
    inv=np.array([[int(vand[i,j]*d) for j in range(m)] for i in range(m)],dtype=np.int64)
    for j in range(2*m):
        view=np.moveaxis(values,j,0)
        assert int(np.max(np.abs(view)))*int(np.abs(inv).sum(axis=1).max())<2**63
        converted=(inv @ view.reshape(m,-1)).reshape(view.shape)
        assert not np.any(converted%d)
        values=np.moveaxis(converted//d,0,j)
    exponents=[]
    coeffs=[]
    for idx in np.ndindex(shape):
        if np.any(values[idx]):
            first=m*m-1-sum(idx)
            assert 0<=first<m
            exponents.append((first,)+idx)
            coeffs.append(values[idx].ravel())
    return axis,exponents,fmpz_mat(np.array(coeffs).T.tolist())


def integer_matrix(rows):
    den=lcm(*(v.denominator for row in rows for v in row))
    return den,fmpz_mat([[int(v*den) for v in row] for row in rows])


@cache
def moment_rule(edges):
    graph=G(max(max(e) for e in edges)+1,edges)
    m=graph.n
    k=2*(4-m)
    source=np.vstack([np.eye(m+1,dtype=int),-np.ones(m+1,dtype=int)])
    trees=tree_flows(graph,source.tolist())
    axis,exponents,bare=bare_table(graph,trees)
    keys=list(trees)
    alpha=list(combinations_with_replacement(range(len(keys)),k))
    beta=sorted({u for a in alpha for u,_ in subsets_with_multiplicity(a) if u})
    positions={b:i for i,b in enumerate(beta)}
    difference=fmpq_mat(len(alpha),len(beta))
    from collections import Counter
    for j,a in enumerate(alpha):
        factorials=prod(factorial(n) for n in Counter(a).values())
        for b,v in subsets_with_multiplicity(a):
            if b:
                difference[j,positions[b]]=f'{v}/{factorials}'
    # Each bare numerator monomial is convergent. Each additional tree
    # complement supplies at least the cycle rank on every scaling subset.
    integrate=monomial_integrator(graph,m+2+k)
    shifts=[tuple(sum(e in keys[i] for i in a) for e in range(graph.N)) for a in alpha]
    periods=[[integrate(tuple(x+y for x,y in zip(ex,s))) for s in shifts] for ex in exponents]
    result=[]
    for channel in (0,1):
        den,table=integer_matrix([[p[channel] for p in row] for row in periods])
        rational=fmpq_mat(bare*table)*difference/den
        denom=lcm(*(int(x.q) for x in rational.entries()))
        array=np.array([[int(rational[i,j]*denom) for j in range(len(beta))]
                        for i in range(len(axis)**2)],dtype=object)
        result.append((denom,array))
    multiplicities=np.zeros((len(beta),len(keys)),dtype=np.int64)
    for i,b in enumerate(beta):
        for j in b:
            multiplicities[i,j]+=1
    return axis,np.array(list(trees.values())),multiplicities,result
