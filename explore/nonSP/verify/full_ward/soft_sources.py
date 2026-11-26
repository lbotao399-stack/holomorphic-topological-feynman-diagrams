"""Exact source restrictions on homogeneous four-loop Ward solutions."""
from itertools import combinations
from pathlib import Path

import networkx as nx
import numpy as np
from flint import fmpz_mat

from contractions import record, scaled_coefficients, substitution
from homogeneous import exact_kernel
from four_mark_relations import edge_order,sign
def zero_sum_cut(cut):
    # Eliminate one source on either side of the zero-total partition.
    other=[v for v in range(6) if v not in cut]
    images=np.zeros((6,4),dtype=np.int64)
    count=0
    for group in [list(cut),other]:
        for v in group[:-1]:
            images[v,count]=1
            images[group[-1],count]=-1
            count+=1
    return substitution(tuple(tuple(int(v) for v in row) for row in images))


def flips(graph,targets):
    edges=edge_order(graph)
    stubs=np.concatenate((np.eye(5,dtype=np.int64),-np.ones((1,5),dtype=np.int64)))
    for u,v in combinations(range(6),2):
        interior=graph.copy()
        interior.remove_nodes_from((u,v))
        comps=list(nx.connected_components(interior))
        for mask in range(1,(1<<len(comps))-1):
            chosen=set().union(*(part for i,part in enumerate(comps) if (mask>>i)&1))
            def switch(a,b):
                if b not in chosen: return a
                if a==u:return v
                if a==v:return u
                return a
            changed=[(switch(a,b),switch(b,a)) for a,b in edges]
            transformed=nx.Graph(changed)
            images=stubs.copy()
            shift=stubs[list(chosen)].sum(axis=0)
            images[list(chosen)]*=-1
            images[u]+=shift;images[v]+=shift
            for i,target in enumerate(targets):
                matcher=nx.isomorphism.GraphMatcher(transformed,target)
                if matcher.is_isomorphic():
                    p=matcher.mapping
                    reordered=[tuple(sorted((p[a],p[b]))) for a,b in changed]
                    sorted_edges=edge_order(target)
                    eps=sign([sorted_edges.index(e) for e in reordered])
                    where=np.argsort([p[j] for j in range(6)])
                    yield i,eps,substitution(tuple(tuple(int(k) for k in row) for row in images[where])),chosen
                    break
            else:
                raise ValueError('Whitney flip missing from catalog.')


def edge_cuts(g):
    for n in (1,2,3):
        for group in combinations(range(6),n):
            if n==3 and 0 not in group:continue
            if sum(1 for _ in nx.edge_boundary(g,group))==2:
                yield group


def check(cache_path=None):
    data,graphs=record()
    if cache_path and Path(cache_path).is_file():
        kernels=np.load(cache_path)
    else:
        invariant,_,_,lift=exact_kernel()
        offsets=np.cumsum([0]+[x.shape[1] for x in invariant])
        kernels=np.stack([x@lift[a:b] for x,a,b in zip(invariant,offsets,offsets[1:])])
        if cache_path:
            np.save(cache_path,kernels)
    _,coeffs=scaled_coefficients(data)
    rows=[]
    for i,g in enumerate(graphs):
        terms=[]
        for group in edge_cuts(g):
            matrix=zero_sum_cut(group)
            assert not np.any(matrix.toarray().astype(object)@coeffs[i])
            part=matrix@kernels[i]
            terms.append((group,fmpz_mat(part.tolist()).rank()))
            rows.append(part)
        print(i,terms,flush=True)
    combined=fmpz_mat(np.concatenate(rows).tolist())
    assert combined.rank()==0
    print('all cut source ranks',combined.rank(),'remaining',15-combined.rank(),flush=True)

    flip_rows,one_rows=[],[]
    for i,g in enumerate(graphs):
        terms=[]
        for j,eps,matrix,chosen in flips(g,graphs):
            assert not np.any(coeffs[i]-eps*(matrix.toarray().astype(object)@coeffs[j]))
            diff=kernels[i]-eps*(matrix@kernels[j])
            terms.append((j, len(chosen),fmpz_mat(diff.tolist()).rank()))
            flip_rows.append(diff)
            if len(chosen)==1:one_rows.append(diff)
        print('flips',i,terms,flush=True)
    for kind,arr in [('one',one_rows),('all',flip_rows)]:
        rank=fmpz_mat(np.concatenate(arr).tolist()).rank()
        assert rank==0
        print(kind,'flip rank',rank,'remaining',15-rank)

    series_rows=[]
    stubs=np.concatenate((np.eye(5,dtype=np.int64),-np.ones((1,5),dtype=np.int64)))
    for i,g in enumerate(graphs):
        ranks=[]
        for w,deg in g.degree:
            if deg!=2:
                continue
            half=[]
            for v in g.neighbors(w):
                mat=2*stubs
                mat[v]+=stubs[w]
                mat[w]=stubs[w]
                half.append(substitution(tuple(tuple(int(x) for x in row) for row in mat)))
            rule=half[0]+half[1]
            assert not np.any(256*coeffs[i]-rule.toarray().astype(object)@coeffs[i])
            terms=256*kernels[i]-rule@kernels[i]
            series_rows.append(terms)
            ranks.append((w,fmpz_mat(terms.tolist()).rank()))
        print('subdivisions',i,ranks,flush=True)
    r=fmpz_mat(np.concatenate(series_rows).tolist()).rank()
    assert r==0
    print('all subdivision rank',r,'remaining',15-r,flush=True)


if __name__=='__main__':
    import sys
    check(sys.argv[1] if len(sys.argv)>1 else None)
