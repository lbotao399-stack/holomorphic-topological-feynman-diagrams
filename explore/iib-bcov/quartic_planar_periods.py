"""Integrate the four-loop planar BCOV periods on series-parallel cores."""
from collections import defaultdict
from fractions import Fraction as F
from functools import cache
from itertools import combinations
from math import factorial, prod
from pathlib import Path
import sys

import networkx as nx
from flint import fmpq, fmpq_mpoly_ctx

from bulk_three import (canonical_resistance, direct_response, forest_data,
                        integrate_core, pushforward)
from planar_bulk import transfer
from path_moments import PathMoments
from fourfold_bulk import symmetries
from quartic_planar_cores import cubic_cores, graph_shapes, metric_graph_factor
from scalar_fan import cofactor
from point_graph import solve_graph
from square_three import G, checked_cut, orbit, skeleton_chains, tree_flows
from current_minor import numerator

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'nonSP/verify/prism_forcer'))
from four_loop import reduced_powers, reduction


@cache
def core_data(edges):
    nv = max(max(edge) for edge in edges)+1
    n = len(edges)-nv+1
    ctx,u,forest = forest_data(edges,nv)
    rho = canonical_resistance(edges,ctx,u,forest,0)
    assert all(rho==canonical_resistance(edges,ctx,u,forest,i) for i in range(1,nv))
    pair = transfer(edges,ctx,u,forest,n)
    times = list(range(1,len(edges)+1))
    check = n*F(str(pair(*times)))/(2*F(str(u(*times)))**3)
    assert check==metric_graph_factor(edges,times)
    short_rho,r = divmod(rho,u)
    short_pair,s = divmod(pair,u**2)
    assert not r and not s and short_rho*u==rho and short_pair*u**2==pair
    return ctx,u,short_rho,short_pair


def ordered_integrals(current,paths,core,ctx,u,forest,slots,n,*,retained_only=False):
    """Integrate the current square against the retained-slot resistances using four positions."""
    space = PathMoments(current,[es for _,es in paths],n-1,ctx=ctx)
    sites = {v:(i,j) for i,(verts,_) in enumerate(paths)
             for j,v in enumerate(verts[1:-1],1)}
    def inner(left,right):
        return space.volume*space.inner(left,right)
    base = inner(space.base,space.base)
    modes = {v:space.position(*sites[v]) for v in slots}
    branch = sorted({v for edge in core for v in edge})
    # The threefold descendant uses just the sum over the six retained pairs.
    pair_weight = fmpq(1) if retained_only else fmpq(-n,12)
    branch_weight = fmpq(0) if retained_only else fmpq(1,8)
    constant = ctx.constant(0) if retained_only else -sum(forest[a][b] for a,b in combinations(branch,2))/(4*n)
    linear = {}
    quadratic = {}
    for v in slots:
        path,position = sites[v]
        a,b = core[path]
        constant += sum(forest[h][a] for h in branch)*branch_weight
        linear[v] = sum(forest[h][b]-forest[h][a] for h in branch)*branch_weight
        # Three retained pairs meet each slot.
        t = ctx.gen(path)
        variance = (t*u-forest[a][b])*(len(branch)*branch_weight+3*pair_weight)
        linear[v] += variance
        quadratic[v,v] = -variance
    for v,w in combinations(slots,2):
        i,s = sites[v]
        j,t = sites[w]
        a,b = core[i]
        c,d = core[j]
        constant += pair_weight*forest[a][c]
        linear[v] += pair_weight*(forest[b][c]-forest[a][c])
        linear[w] += pair_weight*(forest[a][d]-forest[a][c])
        quadratic[v,w] = pair_weight*(forest[b][d]-forest[b][c]-forest[a][d]+forest[a][c])
        if i==j:
            bridge = -2*pair_weight*ctx.gen(i)*u
            linear[v if s<t else w] += bridge
            quadratic[v,w] -= bridge
    result = base*constant
    for v,weight in linear.items():
        result += inner(space.base,modes[v])*weight
    for (v,w),weight in quadratic.items():
        if weight:
            result += inner(modes[v],modes[w])*weight
    return base,result


def dirichlet_values(current,rest,chains,times):
    """Direct four-coordinate Dirichlet integration at specified total core lengths."""
    small = fmpq_mpoly_ctx.get(tuple(f's{i}' for i in range(4)))
    x = small.gens()
    segments = [None]*13
    offset = 0
    groups = []
    for path,t in zip(chains,times):
        ids = range(offset,offset+len(path)-1)
        groups.append(ids)
        for i,j in zip(path,ids):
            segments[i] = t*x[j]
        segments[path[-1]] = t*(1-sum((x[j] for j in ids),small.constant(0)))
        offset += len(path)-1
    assert offset==4
    square = rest.context().from_dict(current.to_dict()).compose(*segments)**2
    jac = prod(t**len(ids) for t,ids in zip(times,groups))
    values = []
    for p in (square,square*rest.compose(*segments)):
        answer = 0
        for ex,c in p.terms():
            ex = tuple(map(int,ex))
            den = prod(factorial(len(ids)+sum(ex[j] for j in ids)) for ids in groups)
            answer += c*jac*prod(factorial(k) for k in ex)/den
        values.append(answer)
    return values


def polynomials(edges,compare=False,components=False,unit_only=False,samples=()):
    net = nx.Graph(edges)
    n = len(edges)-len(net)+1
    branch,paths = skeleton_chains(edges,len(net))
    from_branch = {v:i for i,v in enumerate(branch)}
    ordered = []
    for verts,es in paths:
        if from_branch[verts[0]]>from_branch[verts[-1]]:
            verts,es = verts[::-1],es[::-1]
        pair = from_branch[verts[0]],from_branch[verts[-1]]
        ordered.append((pair,verts,es))
    ordered.sort(key=lambda row:(row[0],tuple(row[2])))
    core = tuple(pair for pair,verts,es in ordered)
    chains = tuple(es for pair,verts,es in ordered)
    ctx,u,rho,pair = core_data(core)
    big,full_u,forest = forest_data(edges,len(net))
    slots = sorted(v for v,degree in net.degree if degree==2)
    source = [(0,0,0)]*len(net)
    for v,p in zip(slots,((1,0,0),(0,1,0),(0,0,1),(-1,-1,-1))):
        source[v] = p
    graph = G(len(net),edges)
    current = numerator(graph,source,next(iter(tree_flows(graph,source).values())))
    rest = sum(forest[b][v] for b in branch for v in slots)/8
    rest -= sum(forest[b][c] for b,c in combinations(branch,2))/(4*n)
    rest -= n*sum(forest[v][w] for v,w in combinations(slots,2))/12
    for seed in (3,7):
        times = [1+(i*i*seed+17*i)%47 for i in range(len(edges))]
        totals = [sum(times[i] for i in path) for path in chains]
        denom = F(str(full_u(*times)))
        assert denom==F(str(u(*totals)))
        flow,_ = solve_graph(edges,source,times,0)
        assert abs(F(str(current(*times))))==denom**n*cofactor(flow,graph.C,times)
        mixed = (F(str(rho(*totals)))/(2*n)+F(str(rest(*times))))/denom
        assert mixed==direct_response(edges,slots,full_u,forest,times)
    aut = sum(1 for _ in nx.isomorphism.GraphMatcher(net,net).isomorphisms_iter())
    if unit_only or samples:
        values = []
        for times in ([[1]*9] if unit_only else samples):
            bz,rz = dirichlet_values(current,rest,chains,times)
            values.append((24*factorial(n)**2*(bz*rho(*times)/(2*n)+rz)/aut,
                           12*n*bz*pair(*times)/aut))
        return values[0] if unit_only else values
    _,_,core_forest = forest_data(core,len(branch))
    base,integrated = ordered_integrals(current,[(verts,es) for pair,verts,es in ordered],
                                       core,ctx,u,core_forest,slots,n)
    times = [1+((i+3)**2)%17 for i in range(9)]
    assert dirichlet_values(current,rest,chains,times)==[base(*times),integrated(*times)]
    if compare:
        square = big.from_dict(current.to_dict())**2
        assert base==pushforward(square,chains,ctx)
        assert integrated==pushforward(square*rest,chains,ctx)
    if components:
        return core,base/aut,integrated/aut
    k = (base*rho/(2*n)+integrated)*(24*factorial(n)**2)/aut
    j = (base*pair)*(12*n)/aut
    return core,((k,2*n+4),(j,2*n+4))


def series_parallel():
    core = next(core for core in cubic_cores() if len(set(core))==7)
    expected = (
        (F(-1453,136914624000000),F(-541,8762535936000000)),
        (F(-4877,410743872000000),F(-1889,26287607808000000)),
        (F(-73,22819104000000),F(-269,13143803904000000)),
    )
    total = [F(0),F(0)]
    whole_k = 0
    shapes = list(graph_shapes(core))
    assert len(shapes)==len(expected)
    for edges,row in zip(shapes,expected):
        reduced,terms = polynomials(edges,compare=True)
        whole_k += terms[0][0]
        if edges==shapes[0]:
            k = terms[0][0]
            assert F(str(k(*([1]*9))))==F(-834,49)
            assert k(1,7,7,1,7,7,1,1,7)==757810515778611603
        for i,(poly,q) in enumerate(terms):
            value,zeta = integrate_core(reduced,poly,4,q)
            assert zeta==0 and value==row[i]
            total[i] += value
            print('series-parallel',i,'terms',len(poly),'value',value,flush=True)
    assert total==[-F(211,8214877440000),-F(1,6490767360000)]
    times = (1,1,1,10,10,1,10,1,1)
    group = symmetries(core)
    average = sum(whole_k(*(times[i] for i in p)) for p in group)/len(group)
    assert F(str(average))==F(65625600000,49)
    assert F(str(whole_k(*([1]*9))))==F(-3408,49)
    print('series-parallel K4, J4',*total,flush=True)
    return total


@cache
def kite_powers(powers):
    return orbit(powers),prod(factorial(k-1) for k in powers)


def one_double():
    core = next(core for core in cubic_cores() if len(set(core))==8)
    columns = [defaultdict(F),defaultdict(F)]
    for i,edges in enumerate(graph_shapes(core)):
        reduced,terms = polynomials(edges)
        moves,name,keep = reduction(reduced)
        assert name=='kite'
        for column,(poly,q) in zip(columns,terms):
            for ex,coefficient in poly.terms():
                powers,weight = reduced_powers(tuple(map(int,ex)),q,moves,keep)
                assert min(powers)>0 and max(powers)<q and sum(powers)==3*q
                representative,multiple = kite_powers(powers)
                column[representative] += weight*multiple*F(int(coefficient.numerator),
                                                              int(coefficient.denominator))/factorial(q-1)
        print('one double',i+1,'reduced monomials',*(len(p) for p,q in terms),flush=True)
    total = [[F(0),F(0)] for _ in columns]
    keys = sorted(set.union(*(set(c) for c in columns)))
    assert len(keys)==1258
    print('tetrahedral moment orbits',len(keys),flush=True)
    for i,key in enumerate(keys):
        values = checked_cut(key)
        for column,row in zip(columns,total):
            for j,v in enumerate(values):
                row[j] += column[key]*v
        if (i+1)%500==0:
            print('completed moments',i+1,flush=True)
    print('one double K4, J4',total,flush=True)
    assert total==[
        [F(11144483303,16895050672500000),-F(186919,340540200000)],
        [F(197793221359,36248733480960000000),-F(507937,111878807040000)],
    ]
    return total


if __name__=='__main__':
    if sys.argv[1:]==['kite']:
        one_double()
    elif not sys.argv[1:]:
        series_parallel()
    else:
        raise SystemExit('Use no arguments for the series-parallel sum, or kite for the one-doubled core.')
