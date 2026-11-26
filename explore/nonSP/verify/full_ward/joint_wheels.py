"""Common positive-resistance wheel fibers, including planar inversion on the wheel."""
from fractions import Fraction
from itertools import product
from math import gcd, prod
from pathlib import Path
import sys

import networkx as nx
import numpy as np
from flint import fmpz_mat, nmod_mat

from contractions import edge_order, record
from homogeneous import exact_kernel
from joint_soft import remaining
from lower_faces import G,tree_flows
from sampling import evaluation,sampling


WHEEL=nx.wheel_graph(5)
CORE=edge_order(WHEEL)
ORBIT=np.array([[CORE.index(tuple(sorted((p[a],p[b])))) for a,b in CORE]
               for p in nx.isomorphism.GraphMatcher(WHEEL,WHEEL).isomorphisms_iter()])
# 0 is the outer dual face. Bounded face i contains vertices 0,i,i+1 cyclically.
DUAL=[]
for a,b in CORE:
    if a==0:
        edge=tuple(sorted((((b-2)%4)+1,b)))
    else:
        edge=(0,4 if (a,b)==(1,4) else a)
    DUAL.append(CORE.index(edge))

WEIGHTS=[460,-954,2196,-2459,2196,-954,460]
assert all(sum(w*Fraction(j,8)**d for j,w in enumerate(WEIGHTS,1))==Fraction(945,d+1)
           for d in range(7))


def presentation(graph):
    edges=edge_order(graph)
    w,=[v for v,d in graph.degree if d==2]
    u,v=graph.neighbors(w)
    collapsed=graph.copy()
    collapsed.remove_node(w)
    collapsed.add_edge(u,v)
    matcher=nx.isomorphism.GraphMatcher(collapsed,WHEEL)
    assert matcher.is_isomorphic()
    mapping=matcher.mapping
    labels=[]
    for edge in edges:
        ends=(u,v) if w in edge else edge
        labels.append(CORE.index(tuple(sorted(mapping[k] for k in ends))))
    a,b=[j for j,e in enumerate(edges) if w in e]
    g=G(6,edges)
    sources=np.concatenate((np.eye(5,dtype=np.int64),-np.ones((1,5),dtype=np.int64)))
    trees=tree_flows(g,sources)
    return labels,a,b,np.array(g.C),np.array(list(trees)),np.array(list(trees.values())).reshape(len(trees),-1)


def fiber(wheel,t,probes):
    labels,a,b,cycles,ct,trees=wheel
    times=t[ORBIT]
    series=times[:,labels[a]]
    original=8*np.repeat(times[:,labels],7,axis=0)
    original[:,a]=np.repeat(series,7)*np.tile(np.arange(1,8),len(ORBIT))
    original[:,b]=8*np.repeat(series,7)-original[:,a]
    products=original[:,ct].prod(axis=2)
    U=products.sum(axis=1)
    # Four tree times <=48, at most 126 spanning trees, and bounded source probes.
    assert original.max()<=48 and np.max(np.abs(probes))<=2
    currents=(products@trees).reshape(-1,9,5)
    y=(currents[:,None,:,:]@probes[None,:,:5,:]).transpose(0,1,3,2)
    matrices=(y[:,:,:,None,1:]*cycles[None,None,:,1:]).reshape(len(original),len(probes),8,8)
    out=[0]*len(probes)
    for j,pack in enumerate(matrices):
        denom=int(original[j,0])*int(U[j])**4
        scale=int(series[j//7])*WEIGHTS[j%7]
        for k,matrix in enumerate(pack):
            c,rem=divmod(int(fmpz_mat(matrix.tolist()).det()),denom)
            assert rem==0
            out[k]+=scale*c
    return out


def pivot_columns(rows):
    prime=1009
    mat,rank=nmod_mat([[int(v)%prime for v in row] for row in rows],prime).rref()
    cols=[]
    for i in range(rank):
        start=cols[-1]+1 if cols else 0
        cols.append(next(j for j in range(start,mat.ncols()) if mat[i,j]))
    return cols


def check(kernels):
    probes=np.array(sampling()[0][:88],dtype=np.int64)
    kept=remaining(kernels)
    kvalues=fmpz_mat(evaluation(probes).tolist())
    rhs=[]
    for i in (8,10):
        rhs.append((kvalues*fmpz_mat(kernels[i].tolist())*kept).transpose().tolist())
    target=[a+b for a,b in zip(*rhs)]
    graphs=record()[1]
    wheels=[presentation(graphs[i]) for i in (8,10)]
    rng=np.random.default_rng(171235)
    sample=[]
    for _ in range(148):
        t=rng.choice((1,2,3,6),size=8)
        reflected=6//t[DUAL]
        jac=prod(map(int,t))**4
        row=[]
        for g in wheels:
            straight=fiber(g,t,probes)
            reverse=fiber(g,reflected,probes)
            row.extend(6**16*a+jac*b for a,b in zip(straight,reverse))
        scale=gcd(*row)
        sample.append([x//scale for x in row] if scale else row)
    indices=pivot_columns(zip(*sample))
    sample=[sample[i] for i in indices]
    cols=pivot_columns(sample)
    assert len(indices)==len(cols)==136
    minor=fmpz_mat([[r[i] for i in cols] for r in sample])
    right=fmpz_mat([[r[i] for i in cols] for r in target])
    combination=minor.transpose().solve(right.transpose()).transpose()
    assert combination*fmpz_mat(sample)==fmpz_mat(target)
    print('Joint wheel fibers: 136 independent rows; all ten directions lie in their exact rational span.')


if __name__=='__main__':
    if len(sys.argv)>1:
        k=np.load(sys.argv[1])
    else:
        invariant,_,_,lift=exact_kernel()
        offsets=np.cumsum([0]+[g.shape[1] for g in invariant])
        k=np.stack([g@lift[a:b] for g,a,b in zip(invariant,offsets,offsets[1:])])
    check(k)
