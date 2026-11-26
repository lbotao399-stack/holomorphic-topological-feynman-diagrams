"""Unprojected proper tight-face compositions for seven-vertex sliding graphs on C²."""
from collections import defaultdict
from fractions import Fraction as F
from functools import cache
from itertools import combinations, combinations_with_replacement
from math import factorial, lcm, prod
from pathlib import Path
import sys
from random import Random
import json

import networkx as nx
from flint import fmpz_mat, fmpz_mpoly_ctx, nmod_mat

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[3]/'code'))
sys.path.insert(0,str(HERE.parent))
from four_mark_relations import edge_order, sign, types
sys.path.insert(0,str(HERE.parent/'k4_periods'))
from kernel import G, Graph, NotSeriesParallel, numerator, plan_generic, skeleton_chains, sp_monomial, tree_flows
from core import EDGES, checked_cut, orbit
from htk.catalog import tight_graphs


def integrate(graph, polynomial, q):
    if not polynomial:
        return F(),F()
    try:
        plan = plan_generic(Graph(graph.nv,graph.edges))
    except NotSeriesParallel:
        plan = None
    if plan is not None:
        return sum((int(c)*sp_monomial(tuple(map(int,e)),q,plan) for e,c in polynomial.terms()),F()),F()
    cubic,chains = skeleton_chains(graph.edges,graph.nv)
    assert len(cubic)==4
    paths={tuple(sorted((cubic.index(path[0]),cubic.index(path[-1])))):es for path,es in chains}
    assert set(paths)==set(EDGES)
    grouped=defaultdict(int)
    for ex,c in polynomial.terms():
        ex=tuple(map(int,ex))
        powers=tuple(sum(1+ex[k] for k in paths[edge]) for edge in EDGES)
        grouped[orbit(powers)]+=int(c)*prod(map(factorial,ex))
    out=[F(),F()]
    for powers,c in grouped.items():
        if c:
            for i,x in enumerate(checked_cut(powers)):
                out[i]+=F(c,factorial(q-1))*x
    return tuple(out)


def tight_faces(edges):
    for k in (3,4,5):
        for cut in combinations(range(7),k):
            inner=[e for e in edges if e[0] in cut and e[1] in cut]
            if len(inner)!=2*k-3:
                continue
            rest=[e for e in edges if e not in inner]
            vertices=[v for v in range(7) if v not in cut]+[cut[0]]
            collapse={v:vertices.index(cut[0] if v in cut else v) for v in range(7)}
            mapped=[(collapse[i],collapse[j]) for i,j in rest]
            if len(set(map(frozenset,mapped)))!=len(mapped):
                continue
            inner_edges=[(cut.index(i),cut.index(j)) for i,j in inner]
            parity=sign([edges.index(e) for e in inner+rest])
            yield cut,inner_edges,rest,collapse,mapped,parity


@cache
def standard_basis(n):
    pairs=list(combinations(range(n+1),2))
    return tuple(entries for entries in combinations_with_replacement(pairs,n)
                 if all(a[1]<=b[1] for a,b in zip(entries,entries[1:])))


def brackets(momenta, basis):
    determinants={(a,b):momenta[a][0]*momenta[b][1]-momenta[a][1]*momenta[b][0]
                  for a,b in {p for row in basis for p in row}}
    return [prod(determinants[p] for p in row) for row in basis]


@cache
def inner_polynomial(edges):
    n=(len(edges)-1)//2
    graph=G(n+2,edges)
    basis=standard_basis(n)
    rng=Random(2719+n)
    probes=[]
    rows=[]
    for _ in range(2*len(basis)+12):
        p=[tuple(rng.randint(-2,2) for _ in range(2)) for _ in range(n+1)]
        p.append(tuple(-sum(x) for x in zip(*p)))
        rows.append(brackets(p,basis))
        probes.append(p)
    reduced,rank=nmod_mat(rows,65521).transpose().rref()
    assert rank==len(basis)
    indices=[]
    for i in range(rank):
        start=indices[-1]+1 if indices else 0
        indices.append(next(j for j in range(start,len(rows)) if reduced[i,j]))
    indices.extend([i for i in range(len(rows)) if i not in indices][-2:])
    ctx=fmpz_mpoly_ctx.get(tuple('t'+str(i) for i in range(graph.N)))
    results=[]
    for i in indices:
        momenta=probes[i]
        raw=numerator(graph,momenta,tree_flows(graph,momenta))
        results.append(integrate(graph,ctx.from_dict(raw),n+2))
    matrix=fmpz_mat([rows[i] for i in indices[:rank]])
    denom=lcm(*(v.denominator for row in results for v in row))
    data=[[int(v*denom) for v in row] for row in results]
    coeff=matrix.solve(fmpz_mat(data[:rank]))
    assert fmpz_mat([rows[i] for i in indices])*coeff==fmpz_mat(data)
    return tuple(tuple(F(str(coeff[i,j]))/denom for i in range(rank)) for j in range(2))


@cache
def canonical_inner(edges):
    n=(len(edges)-1)//2
    original=nx.Graph(edges)
    for ref in tight_graphs(2,n):
        target=tuple(ref.edges)
        matcher=nx.isomorphism.GraphMatcher(original,nx.Graph(target))
        if matcher.is_isomorphic():
            mapping=tuple(matcher.mapping[i] for i in range(n+2))
            reordered=[tuple(sorted((mapping[a],mapping[b]))) for a,b in edges]
            return target,mapping,sign([target.index(e) for e in reordered])
    raise ValueError('No catalog representative for the inner graph.')


def lower_faces(edges,source):
    if len(source)!=7 or any(len(p)!=2 for p in source) or any(map(sum,zip(*source))):
        raise ValueError('Expected seven conserved two-dimensional sources.')
    den=lcm(*(F(v).denominator for row in source for v in row))
    source=[[int(den*v) for v in row] for row in source]
    output=[]
    for cut,inner,rest,mapping,outer,parity in tight_faces(edges):
        n=len(cut)-2
        m=4-n
        graph=G(m+2,outer)
        momenta=[[0,0] for _ in range(m+2)]
        for v,vec in enumerate(source):
            for j,p in enumerate(vec):
                momenta[mapping[v]][j]+=p
        trees=tree_flows(graph,momenta)
        ctx=fmpz_mpoly_ctx.get(tuple('t'+str(i) for i in range(graph.N)))
        exponents={key:tuple(int(i in key) for i in range(graph.N)) for key in trees}
        U=ctx.from_dict({ex:1 for ex in exponents.values()})
        density=ctx.from_dict(numerator(graph,momenta,trees))
        target,reorder,inner_sign=canonical_inner(tuple(inner))
        tau=[None]*len(cut)
        for k,v in enumerate(cut):
            row=[]
            for j in range(2):
                shifts={exponents[key]:sum((int(b==v)-int(a==v))*currents[e][j]
                                          for e,(a,b) in enumerate(rest)) for key,currents in trees.items()}
                row.append(source[v][j]*U+ctx.from_dict(shifts))
            tau[reorder[k]]=row
        assert all(not sum(v[j] for v in tau) for j in range(2))
        values=brackets(tau,standard_basis(n))
        result=[F(),F()]
        prefactor=F(-(-1)**(m*n)*parity*inner_sign,den**8)
        for j,coeffs in enumerate(inner_polynomial(target)):
            if not any(coeffs):
                continue
            d=lcm(*(c.denominator for c in coeffs))
            poly=sum(int(c*d)*v for c,v in zip(coeffs,values) if c)
            parts=integrate(graph,density*poly,m+2+2*n)
            if j:
                assert not parts[1]
            for k,val in enumerate(parts):
                if val:
                    result[j+k]+=prefactor*val/d
        output.append((cut,tuple(result)))
    return output


SOURCES=((1,2),(-2,1),(3,-1),(-1,-2),(2,3),(-4,1),(1,-4))


def certificate():
    data=[]
    for i,g in enumerate(types(7)):
        edges=edge_order(g)
        faces=lower_faces(edges,SOURCES)
        total=tuple(sum((row[d] for _,row in faces),F()) for d in range(2))
        record={'graph':i,'edges':edges,'faces':[(c,list(map(str,v))) for c,v in faces],
                'rhs':list(map(str,total))}
        data.append(record)
        print(i, ' '.join(record['rhs']),flush=True)
    return {'sources':SOURCES,'basis':['1','zeta(3)'],'graphs':data}


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--save',action='store_true')
    args=parser.parse_args()
    artifact=HERE/'lower_faces.json'
    data=json.dumps(certificate(),separators=(',',':'))+'\n'
    if args.save:
        artifact.write_text(data)
    else:
        assert data==artifact.read_text()
