"""Admissible epsilon and Zhang phi inserted into the threefold four-wave graph measure."""
from fractions import Fraction as F
from itertools import combinations
import networkx as nx
from flint import fmpq

from zero_wave_descendant import exact_numerator,tree_map
from bulk_three import canonical_resistance,forest_data,integrate_core,pushforward
from square_three import G,skeleton_chains,tight_graphs


def metric_numerators(core, ctx, u, forest):
    """Kirchhoff determinant times each admissible metric-graph invariant."""
    vertices = len(forest)
    g = len(core)-vertices+1
    rho,rem = divmod(canonical_resistance(core,ctx,u,forest,0),u)
    assert not rem
    epsilon = (2*(g-1)*rho+sum(forest[i][j] for i,j in combinations(range(vertices),2)))/g
    phi = -sum(ctx.gens())*u/4+3*rho/2+epsilon/2
    return epsilon,phi


def graph_periods(edges):
    net = nx.Graph(edges)
    graph = G(len(net),edges)
    g = graph.n
    branch,chains = skeleton_chains(edges,len(net))
    ordered = sorted((tuple(sorted((branch.index(verts[0]),branch.index(verts[-1])))),es)
                     for verts,es in chains)
    core,paths = map(tuple,zip(*ordered))
    ctx,u,forest = forest_data(core,len(branch))
    epsilon,phi = metric_numerators(core,ctx,u,forest)
    slots = [v for v,d in net.degree if d==2]
    sources = [(0,0,0)]*len(net)
    for i,p in zip(slots,((1,0,0),(0,1,0),(0,0,1),(-1,-1,-1))):
        sources[i] = p
    current = exact_numerator(graph,tree_map(graph),sources)
    aut = sum(1 for _ in nx.isomorphism.GraphMatcher(net,net).isomorphisms_iter())
    base = pushforward(current**2,paths,ctx)*fmpq(24,aut)
    return [integrate_core(core,base*weight,g,2*g+4) for weight in (epsilon,phi)]


if __name__=='__main__':
    for g in (2,3):
        totals = [[F(),F()],[F(),F()]]
        for diagram in tight_graphs(3,g):
            if max(dict(nx.Graph(diagram.edges).degree()).values())!=3:
                continue
            for row,weight in zip(totals,graph_periods(diagram.edges)):
                for j,c in enumerate(weight):
                    row[j] += c
        print(g,'epsilon and phi',totals)
