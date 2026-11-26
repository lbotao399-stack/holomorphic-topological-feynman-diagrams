"""A transverse three-loop bulk potential difference on the repeated-momentum source."""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from math import factorial, prod
from flint import fmpq, fmpq_mpoly_ctx, fmpz_mat
import networkx as nx

from square_three import (EDGES, G, Graph, NotSeriesParallel, checked_cut, numerator,
                          orbit, plan_generic, skeleton_chains, sp_monomial, tight_graphs,
                          tree_flows)


def forest_data(edges, nv):
    ctx = fmpq_mpoly_ctx.get(tuple(f'x{i}' for i in range(len(edges))))
    x = ctx.gens()
    u = ctx.constant(0)
    resistances = [[ctx.constant(0) for _ in range(nv)] for _ in range(nv)]
    for count in (nv-1, nv-2):
        for keep in combinations(range(len(edges)), count):
            graph = nx.MultiGraph()
            graph.add_nodes_from(range(nv))
            graph.add_edges_from(edges[j] for j in keep)
            parts = list(nx.connected_components(graph))
            if len(parts) != nv-count:
                continue
            weight = prod(x[j] for j in range(len(edges)) if j not in keep)
            if len(parts) == 1:
                u += weight
            else:
                for i in parts[0]:
                    for j in parts[1]:
                        resistances[i][j] += weight
                        resistances[j][i] += weight
    return ctx, u, resistances


def pushforward(poly, chains, core_ctx):
    # Integrate the path fractions in a monomial at fixed total path resistance.
    values = defaultdict(F)
    for ex, coefficient in poly.terms():
        ex = tuple(map(int, ex))
        powers = tuple(len(path)-1+sum(ex[i] for i in path) for path in chains)
        values[powers] += F(str(coefficient))*F(prod(map(factorial, ex)),prod(map(factorial, powers)))
    return core_ctx.from_dict({ex:fmpq(value.numerator,value.denominator)
                               for ex,value in values.items() if value})


@lru_cache(maxsize=8)
def core_plan(edges,n):
    try:
        return plan_generic(Graph(2*n-2,edges))
    except NotSeriesParallel:
        if edges!=EDGES:
            raise
        return None


@lru_cache(maxsize=16384)
def core_moment(edges,n,q,ex):
    steps = core_plan(edges,n)
    if steps is not None:
        return sp_monomial(ex,q,steps),F(0)
    powers = orbit(tuple(k+1 for k in ex))
    factor = F(prod(map(factorial,ex)),factorial(q-1))
    return tuple(factor*v for v in checked_cut(powers))


def integrate_core(edges, polynomial, n, q):
    edges = tuple(edges)
    values = [F(0),F(0)]
    for ex,c in polynomial.terms():
        weight = F(int(c.numerator),int(c.denominator))
        for i,v in enumerate(core_moment(edges,n,q,tuple(map(int,ex)))):
            if v:
                values[i] += weight*v
    return tuple(values)


def canonical_resistance(edges, ctx, u, forest, at):
    # U^2 times the mean resistance from a fixed point to the canonical measure.
    x = ctx.gens()
    branch = range(len(forest))
    result = -u*sum(forest[v][at] for v in branch)/2
    for i,(left,right) in enumerate(edges):
        excluded = x[i]*u.derivative(i)
        result += excluded*(3*(forest[left][at]+forest[right][at])
                            +x[i]*u-forest[left][right])/6
    return result


def direct_response(edges, slots, u, forest, times):
    # Independently average the two randomly drawn cotree sources. Their
    # Brownian bridges contribute only if they occupy the very same edge.
    tree = int(u(*times))
    resistance = [[F(int(x(*times)),tree) for x in row] for row in forest]
    n = len(edges)-len(forest)+1
    probs = [F(times[i]*int(u.derivative(i)(*times)),n*tree) for i in range(len(edges))]
    assert sum(probs)==1
    mu = [F(0)]*len(forest)
    for (v,w),weight in zip(edges,probs):
        mu[v] += weight/2
        mu[w] += weight/2
    mixed = sum(m*resistance[i][j] for i,m in enumerate(mu) for j in slots)/4
    same = sum(x*y*resistance[i][j] for i,x in enumerate(mu) for j,y in enumerate(mu))
    two = sum(resistance[i][j] for i,j in combinations(slots,2))/6
    return F(n,2)*(2*mixed-same-two+sum(t*p*p for t,p in zip(times,probs))/6)


def evaluate(edges):
    net = nx.Graph(edges)
    nv = len(net)
    n = len(edges)-nv+1
    slots = sorted(v for v,k in net.degree if k==2)
    assert len(slots)==4
    branch, paths = skeleton_chains(edges,nv)
    ordered = sorted((tuple(sorted((branch.index(p[0]),branch.index(p[-1])))),es) for p,es in paths)
    core_edges, chains = zip(*ordered)
    cctx, core_u, core_f = forest_data(core_edges,len(branch))
    rho = canonical_resistance(core_edges,cctx,core_u,core_f,0)
    assert all(canonical_resistance(core_edges,cctx,core_u,core_f,i)==rho for i in range(1,len(branch)))
    ctx,u,f = forest_data(edges,nv)
    rest = sum(f[b][t] for b in branch for t in slots)/8
    rest -= sum(f[b][c] for b,c in combinations(branch,2))/(4*n)
    rest -= n*sum(f[v][w] for v,w in combinations(slots,2))/12
    for times in ([2+i%7 for i in range(len(edges))], [5+(3*i*i)%11 for i in range(len(edges))]):
        totals = [sum(times[i] for i in path) for path in chains]
        tree = int(u(*times))
        assert tree == core_u(*totals)
        combined = F(str(rho(*totals)))/(2*n*tree**2) + F(str(rest(*times)))/tree
        assert combined == direct_response(edges,slots,u,f,times)
    g = G(nv,edges)
    source = [(0,0,0)]*nv
    for v,p in zip(slots,((1,0,0),(0,1,0),(0,0,1),(-1,-1,-1))):
        source[v] = p
    raw = ctx.from_dict(numerator(g,source,tree_flows(g,source)))**2
    poly = pushforward(raw,chains,cctx)*rho/(2*n) + pushforward(raw*rest,chains,cctx)*core_u
    aut = sum(1 for _ in nx.isomorphism.GraphMatcher(net,net).isomorphisms_iter())
    return tuple(F(24*factorial(n)**2,aut)*v for v in integrate_core(core_edges,poly,n,2*n+5))


def check_theta():
    from bulk_eight import MOMENTA, potentials
    total = sum(evaluate(graph.edges)[0] for graph in tight_graphs(3,2)
                if max(dict(nx.Graph(graph.edges).degree()).values())==3)
    phi = potentials(MOMENTA)
    assert total == phi[0][1]-phi[4][1] == -F(7,77760)
    for a,b in (((2,-3),(3,1)), ((0,2),(4,-1))):
        p = [tuple(v*s for v in row) for row,s in zip(MOMENTA,(*a,*b,1,1,1))]
        p.append(tuple(-sum(col) for col in zip(*p)))
        fields = [list(row) for row in potentials(p)]
        factor = total*prod((*a,*b))**2/2
        for i in range(2):
            fields[i][1] -= sum(b)*factor
            fields[i+2][0] -= sum(a)*factor
        common = [fields[5 if k==2 else 4][k] for k in range(5)]
        for q,row in zip(p,fields):
            delta = [x-y for x,y in zip(row,common)]
            assert all(delta[i]*q[j]==delta[j]*q[i] for i,j in combinations(range(5),2))
    print('Two-loop transverse bulk check:', total, flush=True)


def bulk_contact_test():
    p = [tuple(int(i==j) for i in range(5)) for j in (0,0,0,1,1,1,2,3,4)]
    p.append(tuple(-sum(column) for column in zip(*p)))
    full_rank = 0
    for singles in combinations(range(10),5):
        if fmpz_mat([p[i] for i in singles]).det():
            group = [p[i] for i in range(10) if i not in singles]
            assert fmpz_mat(group).rank() <= 3
            full_rank += 1
    assert full_rank == 42
    def trace(i,j):
        fused = tuple(x+y for x,y in zip(p[i],p[j]))
        q = [row for k,row in enumerate(p) if k not in (i,j)] + [fused]
        return sum(int(fmpz_mat(rows).det())**6 for rows in combinations(q,5))
    pairing = 3*(trace(0,3)-trace(0,9)-trace(6,3)+trace(6,9))
    assert pairing == -1200
    print('Sixth-power determinant trace boundary:',pairing,flush=True)


if __name__=='__main__':
    check_theta()
    bulk_contact_test()
    total = [F(0),F(0)]
    for graph in tight_graphs(3,3):
        if max(dict(nx.Graph(graph.edges).degree()).values()) != 3:
            continue
        values = evaluate(graph.edges)
        print(graph.name, *values, flush=True)
        for j,v in enumerate(values):
            total[j] += v
    assert total == [-F(4059837269,70614415872000), F(1285913,27243216000)]
    print('Bare transverse three-loop bulk potential:', *total, flush=True)
