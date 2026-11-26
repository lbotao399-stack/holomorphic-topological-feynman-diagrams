"""Three-loop point scalar square by two cotree sums and exact threefold periods."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations, permutations
from math import factorial, prod
from pathlib import Path
import sys

import networkx as nx
from flint import fmpz_mpoly_ctx

from point_graph import solve_graph
from scalar_fan import cofactor

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'nonSP/verify/k4_periods'))
from kernel import (EDGES, G, Graph, NotSeriesParallel, checked_cut, numerator,
                    orbit, plan_generic, skeleton_chains, sp_monomial, tree_flows)
from htk.catalog import tight_graphs


def summed_forest(edges, slots):
    """Numerator of the sum of all six pair resistances, with denominator U."""
    terms = {}
    for missing in combinations(range(len(edges)), 4):
        forest = nx.Graph()
        forest.add_nodes_from(range(8))
        forest.add_edges_from(edge for i, edge in enumerate(edges) if i not in missing)
        components = list(nx.connected_components(forest))
        if len(components) == 2:
            weight = prod(len(set(slots) & part) for part in components)
            if weight:
                terms[tuple(int(i in missing) for i in range(len(edges)))] = weight
    return terms


def check_subdivision(graph, times, slots, cotrees, resistance_sum):
    """Check the reduction before any edge integration on an original fivefold graph."""
    cuts = (cotrees[1], cotrees[-2])
    edges, split_times, parent, colors = [], [], [], []
    for i, (left, right) in enumerate(graph.edges):
        chain = [left]
        for color, cut in enumerate(cuts):
            if i in cut:
                chain.append(graph.nv + len(colors))
                colors.append(color)
        chain.append(right)
        for j, edge in enumerate(zip(chain, chain[1:]), 1):
            edges.append(edge)
            parent.append(i)
            split_times.append(F(2*j*times[i], len(chain)*(len(chain)-1)))
    cycles = [[row[i] for i in parent] for row in graph.C]
    voltage_sum = [F(0)]*4
    factors = set()
    for first, second, vector, mark in permutations(slots):
        source = [(0, 0, 0)] * graph.nv
        source[first], source[second], source[vector] = (1,0,0), (0,1,0), (0,0,1)
        source[mark] = (-1,-1,-1)
        flow, _ = solve_graph(graph.edges, source, times, mark)
        small = cofactor(flow, graph.C, times)
        expanded = [(0,0,*p) for p in source]
        expanded[mark] = (-3,-3,*source[mark])
        expanded.extend((int(color == 0),int(color == 1),0,0,0) for color in colors)
        full_flow, phi = solve_graph(edges, expanded, split_times, mark)
        assert cofactor(full_flow, cycles, split_times) == small
        factors.add(small)
        for axis in range(4):
            voltage_sum[axis] += phi[vector][axis]
    assert len(factors) == 1 and next(iter(factors)) > 0
    assert voltage_sum == [2*m*resistance_sum for m in (3,3,1,1)]


def integral(graph, polynomial):
    """Integrate a degree-twenty numerator over the ten-edge simplex with U^10."""
    try:
        plan = plan_generic(Graph(graph.nv, graph.edges))
    except NotSeriesParallel:
        plan = None
    if plan is not None:
        return sum(int(c)*sp_monomial(tuple(map(int, ex)), 10, plan)
                   for ex, c in polynomial.terms()), F(0)
    core, chains = skeleton_chains(graph.edges, graph.nv)
    paths = {tuple(sorted((core.index(path[0]), core.index(path[-1])))): indices
             for path, indices in chains}
    assert set(paths) == set(EDGES)
    grouped = defaultdict(int)
    for ex, coefficient in polynomial.terms():
        ex = tuple(map(int, ex))
        powers = tuple(sum(1+ex[i] for i in paths[e]) for e in EDGES)
        grouped[orbit(powers)] += int(coefficient)*prod(map(factorial, ex))
    value = [F(0), F(0)]
    for powers, coefficient in grouped.items():
        if coefficient:
            for i, v in enumerate(checked_cut(powers)):
                value[i] += F(coefficient, factorial(9))*v
    return tuple(value)


def evaluate(edges):
    network = nx.Graph(edges)
    slots = sorted(v for v, degree in network.degree() if degree == 2)
    assert len(slots) == 4 and max(dict(network.degree()).values()) == 3
    graph = G(8, edges)
    sources = [(0, 0, 0)]*8
    for v, p in zip(slots, ((1,0,0), (0,1,0), (0,0,1), (-1,-1,-1))):
        sources[v] = p
    trees = tree_flows(graph, sources)
    poly = numerator(graph, sources, trees)
    forest = summed_forest(edges, slots)
    times = [2+i%7 for i in range(len(edges))]
    u = sum(prod(times[i] for i in missing) for missing in trees)
    resistance_sum = F(sum(c*prod(t**k for t,k in zip(times,ex)) for ex,c in forest.items()), u)
    # Verify the forest identity against 24 grounded, labeled potential evaluations.
    voltage_sum = F(0)
    for scalar1, scalar2, vector, mark in permutations(slots):
        source = [(0,)]*8
        source[scalar1] = source[scalar2] = (1,)
        source[mark] = (-2,)
        _, phi = solve_graph(edges, source, times, mark)
        voltage_sum += phi[vector][0]
    assert voltage_sum == 4*resistance_sum
    check_subdivision(graph, times, slots, list(trees), resistance_sum)
    ctx = fmpz_mpoly_ctx.get(tuple(f't{i}' for i in range(len(edges))))
    terms = ctx.from_dict(poly)**2 * ctx.from_dict(forest)
    automorphisms = sum(1 for _ in nx.isomorphism.GraphMatcher(network, network).isomorphisms_iter())
    return tuple(F(4, automorphisms)*x for x in integral(graph, terms))


def check_bulk_insertion():
    # The other four vertices of the hexagon each carry a single external wave.
    axes = (0,0,0,1,1,1,2,3,4)
    subsets = 0
    for singles in combinations(range(9), 4):
        directions = {axes[i] for i in singles}
        if len(directions) < 4:
            continue
        grouped = {axes[i] for i in range(9) if i not in singles}
        assert len(grouped) + 1 <= 4  # only one independent incident edge current
        subsets += 1
    return subsets


if __name__ == '__main__':
    assert check_bulk_insertion() == 33
    total = [F(0), F(0)]
    count = 0
    for graph in tight_graphs(3, 3):
        if max(dict(nx.Graph(graph.edges).degree()).values()) != 3:
            continue
        value = evaluate(graph.edges)
        print(graph.name, *value, flush=True)
        for i in (0, 1):
            total[i] += value[i]
        count += 1
    assert count == 9
    assert total == [F(411065203,47076277248000), -F(3617,504504000)]
    print('Threefold (rational, zeta(3)):', *total)
    print('Fivefold (rational, zeta(3)):', *(factorial(3)**2 * 4 * x for x in total))
