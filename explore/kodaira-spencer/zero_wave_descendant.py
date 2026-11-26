"""The C3 vertex with a zero wave carrying two cotangent-line insertions."""
from fractions import Fraction as F
from itertools import permutations
from math import comb
from pathlib import Path
import sys

import networkx as nx
from flint import fmpq

IIB = Path(__file__).resolve().parents[1]/'iib-bcov'
sys.path.insert(0,str(IIB))
from bulk_three import forest_data, integrate_core, pushforward
from point_graph import solve_graph
from scalar_fan import cofactor
from square_three import G, skeleton_chains, tight_graphs
from current_minor import numerator


def tree_map(graph):
    edges = {tuple(sorted(edge)):(i,edge) for i,edge in enumerate(graph.edges)}
    paths = nx.single_source_shortest_path(nx.Graph(graph.edges),0)
    return [[(edges[tuple(sorted((u,v)))][0],1 if edges[tuple(sorted((u,v)))][1]==(v,u) else -1)
             for u,v in zip(path,path[1:])] for _,path in sorted(paths.items())]


def particular_flow(graph,paths,sources):
    currents = [[0]*3 for _ in graph.edges]
    for path,source in zip(paths,sources):
        for edge,sign in path:
            for j,x in enumerate(source):
                currents[edge][j] += sign*x
    return currents


def exact_numerator(graph,paths,sources):
    return numerator(graph,sources,particular_flow(graph,paths,sources))


def change(sources,w,k,i,t):
    result = [list(p) for p in sources]
    result[w][i] += t
    result[k][i] -= t
    return result


def literal_derivative(edges,slots,sources,times):
    """Differentiate the unsimplified flux at fixed metric using an exact central stencil."""
    graph = G(len(sources),edges)
    g = graph.n
    order = g+1
    total = F()
    data = {}
    for i,h,k in permutations(range(4),3):
        # Express the four momenta in the frame with i=e^2, h=e^1, k balancing.
        frame = [(0,0,0)]*len(sources)
        frame[slots[i]],frame[slots[h]],frame[slots[k]] = (0,1,0),(1,0,0),(-1,-1,-1)
        rest, = set(range(4))-{i,h,k}
        frame[slots[rest]] = (0,0,1)
        for w in range(len(sources)):
            for j in [k for k in range(-order,order+1) if k]:
                # All 24 frames have the same squared determinant at fixed i,k,w,j
                # and give different voltages depending on the choice of h.
                key = (i,k,w,j)
                if key not in data:
                    moved = change(frame,w,slots[k],1,j)
                    currents,potential = solve_graph(edges,moved,times,0)
                    data[key] = cofactor(currents,graph.C,times)**2,[potential[v][1] for v in range(len(sources))]
                square,v = data[key]
                coefficient = F((-1)**(abs(j)+1)*comb(2*order,order-abs(j)),j*comb(2*order,order))
                total += coefficient*square*(v[w]-v[slots[h]])
    return total


def graph_polynomial(edges,check=False):
    net = nx.Graph(edges)
    slots = [v for v,d in net.degree if d==2]
    branch = [v for v,d in net.degree if d==3]
    assert len(slots)==4
    graph = G(len(net),edges)
    times = [1+(h*h+2*h)%13 for h in range(len(edges))]
    ctx,tree,dist = forest_data(edges,len(net))
    paths = tree_map(graph)
    sources = [(0,0,0)]*len(net)
    for i,p in zip(slots,((1,0,0),(0,1,0),(0,0,1),(-1,-1,-1))):
        sources[i] = p
    def polynomial(src):
        return ctx.from_dict(exact_numerator(graph,paths,src).to_dict())

    n = graph.n
    base = polynomial(sources)
    if not base:
        return ctx.constant(0)
    total = ctx.constant(0)
    pivot = slots[-1]
    for w in [*slots,*branch]:
        if w in slots:
            derivatives = [n*base if i==w else 0 for i in slots]
        else:
            derivatives = []
            for i in range(3):
                diff = polynomial(change(sources,w,pivot,i,1))-polynomial(change(sources,w,pivot,i,-1))
                if n>2:
                    diff = (8*diff-polynomial(change(sources,w,pivot,i,2))
                            +polynomial(change(sources,w,pivot,i,-2)))/6
                derivatives.append(diff/2)
            derivatives.append(n*base-sum(derivatives))
        for i,h,k in permutations(range(4),3):
            iw,hw,kw = slots[i],slots[h],slots[k]
            voltage = (dist[w][kw]-dist[w][iw]+dist[hw][iw]-dist[hw][kw])/2
            variation = (dist[w][kw]+dist[hw][w]-dist[hw][kw])/2
            total += 2*derivatives[i]*voltage+base*variation
        if check and w in branch:
            # A different origin for the simplex tests the fourth derivative.
            lifted = [list(p) for p in sources]
            for v,sign in ((w,1),(slots[0],-1)):
                for j in range(3):
                    lifted[v][j] += sign*sources[slots[-1]][j]
            y,phi = solve_graph(edges,lifted,times,0)
            assert abs(F(str(polynomial(lifted)(*times)))) == cofactor(y,graph.C,times)*F(str(tree(*times)))**n
            def evaluate(t):
                moved = [list(p) for p in sources]
                for v,sign in ((w,1),(slots[0],-1)):
                    for j in range(3):
                        moved[v][j] += sign*t*sources[slots[-1]][j]
                return F(str(polynomial(moved)(*times)))
            value = (8*(evaluate(1)-evaluate(-1))-(evaluate(2)-evaluate(-2)))/12
            assert value==F(str(derivatives[3](*times)))
    aut = sum(1 for _ in nx.isomorphism.GraphMatcher(net,net).isomorphisms_iter())
    result = base*total/fmpq(aut)
    if check:
        actual = literal_derivative(edges,slots,sources,times)
        assert F(str(result(*times)))*aut == actual*F(str(tree(*times)))**(2*n+1)
    return result


def reduced_diagram(edges,check=False):
    net = nx.Graph(edges)
    branch,chains = skeleton_chains(edges,len(net))
    ordered = sorted((tuple(sorted((branch.index(verts[0]),branch.index(verts[-1])))),es)
                     for verts,es in chains)
    core,paths = map(tuple,zip(*ordered))
    ctx,_,_ = forest_data(core,len(branch))
    return core,pushforward(graph_polynomial(edges,check),paths,ctx)


if __name__=='__main__':
    for genus in (2,3):
        result = [F(),F()]
        for diagram in tight_graphs(3,genus):
            edges = diagram.edges
            if max(dict(nx.Graph(edges).degree()).values())!=3:
                continue
            core,poly = reduced_diagram(edges,check=not any(result))
            value = integrate_core(core,poly,genus,2*genus+4)
            result = [x+y for x,y in zip(result,value)]
            print(genus,diagram.name,*value,flush=True)
        print('Two cotangent lines at one zero wave',genus,*result,flush=True)
        expected = {2:(F(59,77760),F()),3:(F(84137363,4707627724800),-F(3617,252252000))}
        assert tuple(result)==expected[genus]
