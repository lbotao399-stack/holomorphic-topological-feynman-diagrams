"""Unlabelled tight cubic contact graphs in D=3, loop number four."""
from collections import Counter
from itertools import combinations, combinations_with_replacement
from math import factorial, prod
import json

import networkx as nx


def automorphisms(graph):
    match = nx.algorithms.isomorphism.categorical_edge_match('multiplicity', 1)
    return list(nx.algorithms.isomorphism.GraphMatcher(graph, graph, edge_match=match).isomorphisms_iter())


def cores():
    for graph in nx.graph_atlas_g():
        if len(graph) != 6 or not nx.is_biconnected(graph):
            continue
        degree = dict(graph.degree())
        if set(degree.values()) - {2, 3}:
            continue
        missing = {i for i in degree if degree[i] == 2}
        edges = [e for e in graph.edges if set(e) <= missing]
        seen = []
        for double in combinations(edges, len(missing)//2):
            if {i for e in double for i in e} != missing:
                continue
            g = graph.copy()
            for e in graph.edges:
                g.edges[e]['multiplicity'] = 1 + int(e in double)
            match = nx.algorithms.isomorphism.categorical_edge_match('multiplicity', 1)
            if not any(nx.is_isomorphic(g, h, edge_match=match) for h in seen):
                seen.append(g)
                yield g


def is_tight(edges):
    # It suffices to check induced sets; the empty edge set imposes no bound.
    masks = [(1<<a) | (1<<b) for a,b in edges]
    for subset in range(1<<10):
        edge_count = sum(subset & m == m for m in masks)
        if edge_count and 2*edge_count > 3*subset.bit_count()-4:
            return False
    return True


def graph_table():
    for g in cores():
        base = sorted(tuple(sorted(e)) for e in g.edges)
        index = {e:i for i,e in enumerate(base)}
        mult = [g.edges[e]['multiplicity'] for e in base]
        core = [e for e,m in zip(base,mult) for _ in range(m)]
        actions = [tuple(index[tuple(sorted((f[i], f[j])))] for i,j in base)
                   for f in automorphisms(g)]
        seen = set()
        for split in combinations_with_replacement(range(9),4):
            lengths = [1+split.count(j) for j in range(9)]
            decoration = tuple(tuple(sorted(l for e,l in zip(core,lengths) if e==pair))
                               for pair in base)
            if decoration in seen or any(t.count(1)>1 for t in decoration):
                continue
            orbit = {tuple(decoration[i] for i in s) for s in actions}
            seen.update(orbit)
            # Group the parallel edges in a fixed order for reproducible replay.
            n = 6
            edges = []
            chains = []
            for (i,j), group in zip(base,decoration):
                for length in group:
                    nodes = [i]+list(range(n,n+length-1))+[j]
                    n += length-1
                    edges.extend(zip(nodes,nodes[1:]))
                    chains.append(nodes)
            if not is_tight(edges):
                continue
            stabilizer = len(actions)//len(orbit)
            parallel = prod(factorial(v) for ds in decoration for v in Counter(ds).values())
            graph = nx.Graph(edges)
            assert len(graph)==10 and len(graph.edges)==13 and nx.is_biconnected(graph)
            automorphism = stabilizer*parallel
            yield dict(chains=chains, aut=automorphism)


if __name__=='__main__':
    graphs = list(graph_table())
    print(json.dumps(graphs, separators=(',',':')))
