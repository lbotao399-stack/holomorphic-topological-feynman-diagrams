"""Detect every odd triangle in the stored H=T=1 census and replay its current law."""
from itertools import combinations
from pathlib import Path
import json

import networkx as nx
import sympy as sp

CENSUS = Path(__file__).resolve().parents[1]/"c_pi_survey.json"


def partitions(edges):
    graph = nx.Graph(edges)
    for vertices in combinations(sorted(graph),3):
        pairs = list(zip(vertices,vertices[1:]+vertices[:1]))
        blocks = {frozenset((u,v)): [(u,v)] if graph.has_edge(u,v) else [] for u,v in pairs}
        interior = graph.subgraph(set(graph)-set(vertices))
        for component in nx.connected_components(interior):
            adjacent = frozenset(w for v in component for w in graph[v] if w in vertices)
            if len(adjacent) != 2:
                break
            blocks[adjacent].extend(row for row in edges if row[0] in component or row[1] in component)
        else:
            groups = [blocks[frozenset(pair)] for pair in pairs]
            ranks = tuple(len(group)-len(set(v for edge in group for v in edge))+1 for group in groups)
            if all(group and nx.is_connected(nx.Graph(group)) for group in groups) and all(r%2 for r in ranks):
                yield vertices,groups,ranks


def network(vertices,edges,times):
    root = vertices[-1]
    columns = {v:i for i,v in enumerate(vertices[:-1])}
    B = sp.zeros(len(edges),len(columns))
    for i,(u,v) in enumerate(edges):
        if u != root:
            B[i,columns[u]] = 1
        if v != root:
            B[i,columns[v]] = -1
    conductance = sp.diag(*(1/t for t in times))
    response = (B.T*conductance*B).inv()
    return B, conductance*B*response, response


def check_currents(vertices,edges,groups):
    # A terminal unit-current vector in each piece has squared norm R
    # in the resistance metric, and is orthogonal to every internal cycle.
    times = {tuple(sorted(edge)): sp.Integer(1+i%5) for i,edge in enumerate(edges)}
    rescaled = {}
    flow = {}
    for (u,v),group in zip(zip(vertices,vertices[1:]+vertices[:1]),groups):
        local = sorted({w for edge in group for w in edge})
        g_edges = [tuple(sorted(e)) for e in group]
        local_times = [times[e] for e in g_edges]
        B, currents, response = network(local,g_edges,local_times)
        source = sp.Matrix([int(i==u)-int(i==v) for i in local[:-1]])
        h = currents*source
        R = (source.T*response*source)[0]
        assert R == (h.T*sp.diag(*local_times)*h)[0]
        assert all((c.T*sp.diag(*local_times)*h)[0]==0 for c in B.T.nullspace())
        for e,current in zip(g_edges,h):
            flow[e] = current
            rescaled[e] = times[e]/R
    all_vertices = sorted({v for edge in edges for v in edge})
    oriented = [tuple(sorted(e)) for e in edges]
    membership = {tuple(sorted(e)):i for i,g in enumerate(groups) for e in g}
    def currents_at(scales):
        t = [rescaled[e]*scales[membership[e]] for e in oriented]
        return network(all_vertices,oriented,t)[1]
    samples = [(1,1,1),(2,1,1),(1,2,1)]
    currents = [currents_at(point) for point in samples]
    loop = sp.Matrix([flow[e] for e in oriented])
    variations = sp.Matrix.hstack(loop,currents[1]-currents[0],currents[2]-currents[0])
    assert variations.rank()==1
    # Affine interpolation in resistance shares on all unit-source columns at once.
    triangle = sp.Matrix([[sp.Rational(q,sum(p)) for q in p] for p in samples])
    for target in ((1,2,3),(4,3,2),(2,7,3)):
        weights = triangle.T.inv()*sp.Matrix([sp.Rational(q,sum(target)) for q in target])
        predicted = sum((w*c for w,c in zip(weights,currents)),sp.zeros(*currents[0].shape))
        assert predicted == currents_at(target)


if __name__=="__main__":
    census = json.loads(CENSUS.read_text())
    for n,expected in [(2,[]),(4,[2]),(6,[6,13,16,25,29,30])]:
        hits = []
        for index,_,edges,powers in census[f"2,{n}"]:
            cases = list(partitions(edges))
            assert bool(cases)==(powers["small"]>n//2)==(powers["large"]>n//2)
            if not cases:
                continue
            hits.append(index)
            for vertices,groups,ranks in cases:
                check_currents(vertices,edges,groups)
                print(n,index,vertices,ranks,"all source columns checked")
        assert hits == expected
