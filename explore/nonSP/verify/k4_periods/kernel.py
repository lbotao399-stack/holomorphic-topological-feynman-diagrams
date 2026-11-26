"""Exact K4-subdivision kernel with arbitrary rational momenta and edge-current insertions."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations
from math import comb, factorial, lcm, prod
from pathlib import Path
import sys

import networkx as nx
import numpy as np
from flint import fmpz_mat, fmpz_mpoly_ctx

from core import EDGES, checked_cut, orbit

sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from graphs import Graph, K4SUB, skeleton_chains
from htk.spred import NotSeriesParallel, plan_generic
from vk import G, density_exact


def lower_grid(size, remaining, degree):
    if size==0:
        yield ()
        return
    for i in range(min(degree,remaining)+1):
        for rest in lower_grid(size-1,remaining-i,degree):
            yield (i,)+rest


def tree_flows(graph,source):
    rows = {}
    for missing in combinations(range(graph.N),graph.n):
        tree = nx.Graph((u,v,{"index":j}) for j,(u,v) in enumerate(graph.edges) if j not in missing)
        if len(tree)!=graph.nv or not nx.is_tree(tree):
            continue
        current = [[0]*len(source[0]) for _ in graph.edges]
        for v,path in nx.single_source_shortest_path(tree,graph.root).items():
            for u,w in zip(path,path[1:]):
                e = tree[u][w]["index"]
                sign = 1 if graph.edges[e]==(w,u) else -1
                for j,m in enumerate(source[v]):
                    current[e][j] += sign*m
        rows[missing] = current
    return rows


def interpolate(samples,degree_limit,total):
    # Newton differences on a downward closed set, followed by expansion of
    # (t-1)(t-2). All integer divisions are checked.
    order = sorted(samples,key=sum)
    for axis in range(len(order[0])):
        for ex in reversed(order):
            degree = ex[axis]
            if not degree:
                continue
            def at(k):
                return samples[ex[:axis]+(k,)+ex[axis+1:]]
            v = sum((-1)**(degree-j)*comb(degree,j)*at(j) for j in range(degree+1))
            samples[ex],r = divmod(v,factorial(degree))
            assert r==0
    # Integer monomial coefficients of the falling factorials based at time one.
    falling = [[1]]
    for k in range(1,degree_limit+1):
        previous = falling[-1]
        falling.append([a-k*b for a,b in zip([0]+previous,previous+[0])])
    for axis in range(len(order[0])):
        for ex in order:
            degree=ex[axis]
            if degree==degree_limit:
                continue
            def at(k):
                return samples.get(ex[:axis]+(k,)+ex[axis+1:],0)
            samples[ex] += sum(at(j)*falling[j][degree] for j in range(degree+1,degree_limit+1))
    result = {(total-sum(ex),)+ex:c for ex,c in samples.items() if c}
    assert all(0<=ex[0]<=degree_limit for ex in result)
    return result


def numerator(graph,source,trees):
    h = len(source[0])
    keys = list(trees)
    bound = len(keys)*graph.n**graph.n*max(1,max(sum(abs(v[j]) for v in source) for j in range(h)))
    dtype = object if bound >= 2**62 else np.int64
    tree_array = np.array(list(trees.values()),dtype=dtype).reshape(len(keys),-1)
    cols = np.array(keys)
    cycles = np.asarray(graph.C)
    samples = {}
    def batch(points):
        times = np.array([(1,)+tuple(1+k for k in row) for row in points],dtype=np.int64)
        products = times[:,cols].prod(axis=2)
        numerators = (products@tree_array).reshape(-1,graph.N,h)
        for ex,u,y in zip(points,products.sum(axis=1),numerators):
            matrix = (y.T[:,None,1:]*cycles[None,:,1:]).reshape(graph.N-1,graph.N-1)
            det = int(fmpz_mat(matrix.tolist()).det())
            samples[ex],r = divmod(det,int(u)**(graph.n*(h-1)))
            assert r==0
    group = []
    for point in lower_grid(graph.N-1,graph.n**2-1,graph.n-1):
        group.append(point)
        if len(group)==256:
            batch(group)
            group=[]
    if group:
        batch(group)
    poly = interpolate(samples,graph.n-1,graph.n**2-1)
    # A separate grounded linear solve at positive times outside the interpolation grid.
    times=[4+(i%3) for i in range(graph.N)]
    u=sum(prod(times[i] for i in miss) for miss in trees)
    assert sum(c*prod(t**k for t,k in zip(times,ex)) for ex,c in poly.items())==(
        u**(h+graph.n)*density_exact(graph,times,source[:-1]))
    return poly


def sp_monomial(exponents,q,plan):
    degree = [k+1 for k in exponents]
    result = F(1)
    def beta(powers):
        if len(powers)==1:
            return F(1)
        assert min(powers)>0
        return F(prod(factorial(i-1) for i in powers),factorial(sum(powers)-1))
    for _,_,chains in plan:
        child = []
        for chain in chains:
            powers = [degree[e] for e in chain]
            child.append(sum(powers))
            result *= beta(powers)
        result *= beta([q-k for k in child])
        degree.append(sum(child)-(len(child)-1)*q)
    assert degree[-1]==0
    return result


def evaluate(edges,source,insertions=()):
    """An insertion (e,i) multiplies the integrand by component i on oriented edge e."""
    n = max(max(e) for e in edges)+1
    graph = G(n,edges)
    h,r = divmod(len(edges)-1,3)
    if graph.n!=3 or r or h<1 or len(source)!=n or any(len(v)!=h for v in source):
        raise ValueError('Expected a three-loop holomorphic graph and one momentum per vertex.')
    source = [list(map(F,v)) for v in source]
    if any(map(sum,zip(*source))):
        raise ValueError('Momenta must sum to zero.')
    denom = lcm(*(x.denominator for v in source for x in v))
    source = [[int(x*denom) for x in v] for v in source]
    trees = tree_flows(graph,source)
    raw = numerator(graph,source,trees)
    if insertions:
        ctx = fmpz_mpoly_ctx.get(tuple(f't{i}' for i in range(graph.N)))
        poly = ctx.from_dict(raw)
        for e,j in insertions:
            if not 0<=e<graph.N or not 0<=j<h:
                raise ValueError('Unknown oriented edge or source component.')
            poly *= ctx.from_dict({tuple(int(i in m) for i in range(graph.N)):row[e][j]
                                   for m,row in trees.items() if row[e][j]})
        raw = {tuple(map(int,ex)):int(c) for ex,c in poly.terms()}
    q = h+3+len(insertions)
    try:
        plan = plan_generic(Graph(n,edges))
    except NotSeriesParallel:
        plan = None
    if plan is not None:
        value = sum(c*sp_monomial(ex,q,plan) for ex,c in raw.items())
        return F(value,denom**(3*h+len(insertions))),F(0)
    cubic,chains = skeleton_chains(edges,n)
    assert len(cubic)==4
    paths = {tuple(sorted((cubic.index(path[0]),cubic.index(path[-1])))):es for path,es in chains}
    assert set(paths)==set(EDGES)
    grouped = defaultdict(int)
    for ex,c in raw.items():
        powers = tuple(sum(1+ex[k] for k in paths[edge]) for edge in EDGES)
        grouped[orbit(powers)] += c*prod(map(factorial,ex))
    normalization = factorial(q-1)*denom**(3*h+len(insertions))
    value = [F(),F()]
    for powers,c in grouped.items():
        if c:
            for j,v in enumerate(checked_cut(powers)):
                value[j] += F(c,normalization)*v
    return tuple(value)


if __name__=='__main__':
    import json
    reference = json.loads((Path(__file__).resolve().parents[2]/'results_nonSP_C3.json').read_text())['graphs']
    for name,edges in K4SUB:
        vertices=[v for v,k in nx.Graph(edges).degree() if k==2]
        source=[(0,0,0)]*8
        for v,p in zip(sorted(vertices),((1,0,0),(0,1,0),(0,0,1),(-1,-1,-1))):
            source[v]=p
        ans=evaluate(edges,source)
        print(name,ans,flush=True)
        for basis in reference[name]['exact']:
            assert ans==tuple(F(basis[k]) for k in ('c_rational_part','c_zeta3_part'))
