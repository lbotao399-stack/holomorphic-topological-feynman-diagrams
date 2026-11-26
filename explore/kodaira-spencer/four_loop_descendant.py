"""Integrate the four bivalent marked sites of the genus-four C3 descendant."""
from argparse import ArgumentParser
from collections import defaultdict
from fractions import Fraction
from itertools import combinations
import json
from math import factorial, prod
from pathlib import Path
import sys

import networkx as nx
import numpy as np
from flint import fmpq

IIB = Path(__file__).resolve().parents[1]/'iib-bcov'
sys.path.insert(0,str(IIB))
from bulk_three import forest_data, integrate_core
from current_minor import numerator
from quartic_master_weights import AUTOMORPHISMS, PERMS, orbit_groups, unpack
from quartic_planar_cores import cubic_cores, graph_shapes
from quartic_planar_periods import (G, checked_cut, dirichlet_values, kite_powers, ordered_integrals,
                                   reduced_powers, reduction, skeleton_chains, tree_flows)
from point_graph import solve_graph
from scalar_fan import cofactor


def diagram(core_name,index):
    distinct_edges = {'two-double':7,'one-double':8,'prism':9,'bipartite':9}[core_name]
    core, = (c for c in cubic_cores() if len(set(c))==distinct_edges
             and (distinct_edges<9 or nx.is_bipartite(nx.Graph(c))==(core_name=='bipartite')))
    edges = list(graph_shapes(core))[index]
    graph = nx.Graph(edges)
    branch,chains = skeleton_chains(edges,len(graph))
    slots = [v for v,d in graph.degree() if d==2]
    paths = []
    for verts,es in chains:
        i,j = branch.index(verts[0]),branch.index(verts[-1])
        if i>j:
            verts,es = verts[::-1],es[::-1]
        paths.append((tuple(sorted((i,j))),verts,es))
    paths.sort(key=lambda x:(x[0],tuple(x[2])))
    core = tuple(q for q,_,_ in paths)
    ctx,u,forest = forest_data(core,len(branch))
    sources = [(0,0,0)]*len(graph)
    for v,p in zip(slots,((1,0,0),(0,1,0),(0,0,1),(-1,-1,-1))):
        sources[v] = p
    g = G(len(graph),edges)
    current = numerator(g,sources,next(iter(tree_flows(g,sources).values())))
    big,full_u,full_forest = forest_data(edges,len(graph))
    resistance = sum(full_forest[v][w] for v,w in combinations(slots,2))
    orders = [(vs,es) for _,vs,es in paths]
    _,integrated = ordered_integrals(current,orders,core,ctx,u,forest,slots,4,retained_only=True)
    # Independent simplex cubature makes no use of interpolation between endpoints.
    times = [1+(k*k+2*k)%17 for k in range(9)]
    assert dirichlet_values(current,resistance,[es for _,es in orders],times)[1]==integrated(*times)
    full_times = [2+(k*k+3*k)%13 for k in range(len(edges))]
    y,_ = solve_graph(edges,sources,full_times,0)
    assert abs(Fraction(str(current(*full_times))))==int(full_u(*full_times))**4*cofactor(y,g.C,full_times)
    aut = sum(1 for _ in nx.isomorphism.GraphMatcher(graph,graph).isomorphisms_iter())
    return core,2*integrated/aut


def integrate(name,index,output=None,builder=diagram):
    core,poly = builder(name,index)
    print(name,index,"integrated numerator monomials",len(poly),flush=True)
    if name=='two-double':
        return integrate_core(core,poly,4,12)
    moves,kind,keep = reduction(core)
    if name=='one-double':
        assert kind=='kite'
        weights = defaultdict(Fraction)
        for i in range(len(poly)):
            ex = tuple(map(int,poly.monomial(i)))
            powers,weight = reduced_powers(ex,12,moves,keep)
            key,multiple = kite_powers(powers)
            weights[key] += weight*multiple*Fraction(str(poly.coefficient(i)))/factorial(11)
        return tuple(sum(c*checked_cut(key)[j] for key,c in weights.items() if c) for j in range(2))
    if output is None:
        raise ValueError("Specify an output file for primitive moment weights")
    assert not moves and kind==('prism' if name=='prism' else 'K33')
    perms = np.asarray(PERMS if name=='prism' else AUTOMORPHISMS,dtype=np.intp)
    weights = orbit_groups(poly,perms,keep)
    del poly
    import gzip
    with gzip.open(output,'wt') as f:
        for key,value in sorted(weights.items()):
            powers = unpack(key)
            value *= fmpq(prod(factorial(i-1) for i in powers),factorial(11))
            print(*powers,value,0,file=f)
    print(len(weights),"projective orbit weights",flush=True)


def verify_periods():
    from descendant_ladder import ladder_weights
    from quartic_prism_physical_ibp import zeta_interval
    data = json.loads(Path(__file__).with_name('four_loop_descendant.json').read_text())
    bounds = [zeta_interval(s,128) for s in (3,5)]
    total = [Fraction(0)]*3
    for name,count in [('two-double',3),('one-double',8),('prism',28),('bipartite',11)]:
        diagrams = data['graphs'][name]
        assert len(diagrams)==count
        for index,row in enumerate(diagrams):
            a,b,c = map(Fraction,row)
            assert a+sum(min(k*lo,k*hi) for k,(lo,hi) in zip((b,c),bounds))>0
            if name=='two-double' and index==0:
                assert (a,b,c)==(ladder_weights(4)[2],0,0)
            if name in ('two-double','one-double'):
                assert c==0
            total = [v+w for v,w in zip(total,(a,b,c))]
    assert total==list(map(Fraction,data['sum']))
    print('D4 in the 1, zeta(3), zeta(5) basis:',*(str(x) for x in total))


if __name__=='__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('core',nargs='?',choices=['two-double','one-double','prism','bipartite'])
    parser.add_argument('shape',type=int,nargs='?')
    parser.add_argument('--weights',type=Path)
    args = parser.parse_args()
    if args.core is None and args.weights is None:
        verify_periods()
    else:
        if args.core is None or args.shape is None or args.shape<0:
            parser.error('Specify a core and a nonnegative graph index')
        result = integrate(args.core,args.shape,args.weights)
        if result is not None:
            print("Integrated",*(str(x) for x in result))
