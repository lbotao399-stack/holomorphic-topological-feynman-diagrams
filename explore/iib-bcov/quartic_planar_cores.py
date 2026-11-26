"""Cubic four-loop cores for the two-plane bulk reduction and the strict sign of J."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement, permutations

import networkx as nx
from flint import fmpq_mat

from fourfold_bulk import symmetries
from scalar_fan import cofactor
from point_graph import solve_graph
from square_three import G


def cubic_cores():
    # At six vertices a bridgeless cubic graph has no loops or triple edges.
    # The endpoints of its double edges form a matching of the support's bivalent vertices.
    cores = set()
    for graph in nx.graph_atlas_g():
        if len(graph)!=6 or not nx.is_biconnected(graph):
            continue
        vertices = set(graph)
        degrees = dict(graph.degree())
        if min(degrees.values())<2 or max(degrees.values())>3:
            continue
        bivalent = {i for i in vertices if degrees[i]==2}
        for doubled in combinations(graph.edges,len(bivalent)//2):
            if {i for edge in doubled for i in edge}!=bivalent:
                continue
            edges = [*graph.edges,*doubled]
            images = (tuple(sorted(tuple(sorted((g[a],g[b]))) for a,b in edges))
                      for g in permutations(range(6)))
            cores.add(min(images))
    return sorted(cores)


def constraints(core):
    tests = []
    for size in range(2,len(core)+1):
        for indices in combinations(range(len(core)),size):
            subs = nx.MultiGraph([core[e] for e in indices])
            if not nx.is_connected(subs):
                continue
            # A tree appendage can only increase length at fixed cycle rank.
            if min(dict(subs.degree()).values())<2:
                continue
            genus = size-len(subs)+1
            deficit = 3*genus+1-size
            if deficit>0:
                tests.append((indices,deficit))
    return tests


def graph_shapes(core):
    tests = constraints(core)
    symmetry = symmetries(core)
    seen = set()
    for marks in combinations_with_replacement(range(9),4):
        if marks in seen:
            continue
        orbit = {tuple(sorted(g[e] for e in marks)) for g in symmetry}
        seen.update(orbit)
        counts = Counter(marks)
        if any(sum(counts[e] for e in indices)<limit for indices,limit in tests):
            continue
        edges = []
        next_vertex = 6
        for i,(a,b) in enumerate(core):
            path = [a,*range(next_vertex,next_vertex+counts[i]),b]
            next_vertex += counts[i]
            edges.extend(zip(path,path[1:]))
        assert nx.number_of_selfloops(nx.Graph(edges))==0 and len(set(tuple(sorted(e)) for e in edges))==13
        yield edges


def metric_graph_factor(edges,times):
    # Compute the same resistance transfer from a grounded conductance matrix
    # and from the Gram matrix of the cycle projection and midpoint evaluations.
    n = len(edges)-max(max(e) for e in edges)
    columns = [[int(a==v)-int(b==v) for v in range(1,len(edges)+1-n)] for a,b in edges]
    incidence = fmpq_mat(columns)
    conductance = fmpq_mat([[str(F(int(i==j),times[i])) for j in range(len(edges))]
                            for i in range(len(edges))])
    green = (incidence.transpose()*conductance*incidence).inv()
    transfer = conductance-conductance*incidence*green*incidence.transpose()*conductance
    mid = fmpq_mat([[int(a==v)+int(b==v) for v in range(1,len(columns)+1-n)] for a,b in edges])/2
    midpoint = mid*green*mid.transpose()
    def resistance(a,b):
        def entry(i,j):
            return green[i-1,j-1] if i and j else 0
        return entry(a,a)+entry(b,b)-2*entry(a,b)
    probabilities = [times[i]*transfer[i,i]/n for i in range(len(edges))]
    assert sum(probabilities)==1
    difference = gram = diagonal = F(0)
    for i,(a,b) in enumerate(edges):
        for j,(c,d) in enumerate(edges):
            covariance = F(str(times[i]*times[j]*transfer[i,j]**2/n-probabilities[i]*probabilities[j]))
            pair = F(str(sum(resistance(v,w) for v in (a,b) for w in (c,d))))/4
            if i==j:
                assert covariance==(n-1)*F(str(probabilities[i]))**2
                pair -= F(times[i],6)
                diagonal += covariance*times[i]/6
            difference += covariance*pair
            gram += covariance*F(str(midpoint[i,j]))
    assert difference == -2*gram-diagonal and gram>=0 and diagonal>0
    return n*difference/2


def certificate():
    source = [(0,0,0)]*6+[(1,0,0),(0,1,0),(0,0,1),(-1,-1,-1)]
    expected = {
        16: (3,F(-53200,13089)), 48: (0,F(-97835,22896)),
        8: (8,F(-375221,135288)), 12: (28,F(-522493,250232)),
        72: (11,F(-149341,76634)),
    }
    checked = set()
    for core in cubic_cores():
        count = ok = 0
        for edges in graph_shapes(core):
            graph = G(10,edges)
            count += 1
            for offset in range(3):
                times = [1+(7*k*k+17*k+31*offset*k)%101 for k in range(13)]
                current,_ = solve_graph(edges,source,times,0)
                if cofactor(current,graph.C,times):
                    ok += 1
                    break
        net = nx.Graph(core)
        reduced = net.copy()
        while len(reduced)>2 and min(dict(reduced.degree()).values())<3:
            v=next(v for v,k in reduced.degree if k<3)
            nb=list(reduced[v])
            reduced.remove_node(v)
            if len(nb)==2:reduced.add_edge(*nb)
        assert ok==count,(core,ok,count)
        sign = metric_graph_factor(core,list(range(1,10)))
        order = len(symmetries(core))
        assert order not in checked and (count,sign)==expected[order]
        checked.add(order)
        print('edges',core,'automorphisms',order, 'simple reduced V',len(reduced),
              'tight shapes',count,'J fixed-time transfer',sign,flush=True)
    assert checked==set(expected)


if __name__=='__main__':
    certificate()
