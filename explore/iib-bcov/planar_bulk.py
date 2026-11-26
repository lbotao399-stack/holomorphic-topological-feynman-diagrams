"""Two planar source directions with arbitrary waves on the fivefold."""
from fractions import Fraction as F
from itertools import combinations
from math import factorial
from functools import cache
import networkx as nx
from bulk_three import forest_data, integrate_core, pushforward
from square_three import G, numerator, skeleton_chains, tight_graphs, tree_flows


def bracket(p,q):
    return p[0]*q[1]-p[1]*q[0]


@cache
def mass(p):
    """Coefficient of all distinct variables in det(sum_i x_i p_i p_i^t)^n."""
    if not p:
        return 1
    return len(p)//2*sum(bracket(p[0],p[j])**2*mass(p[1:j]+p[j+1:]) for j in range(1,len(p)))


def response(p,k,j):
    n = len(p)//2
    z = mass(p)
    out = [[k*z*sum(row[d] for row in p)/(n*factorial(n)**2) for d in range(2)] for _ in p]
    for i,h in combinations(range(2*n),2):
        others = tuple(q for pos,q in enumerate(p) if pos not in (i,h))
        b = F(z-n*n*bracket(p[i],p[h])**2*mass(others),n*(n-1))
        for d in range(2):
            out[i][d] += j*b*p[h][d]
            out[h][d] += j*b*p[i][d]
    return [r+[0,0,0] for r in out]+[[0]*5 for _ in range(4)]


def transfer(edges, ctx, u, f, n):
    # U^3 times the mean-resistance difference. Both edge distributions have
    # the same marginals; their interpolation variances cancel.
    x = ctx.gens()
    deriv = [u.derivative(i) for i in range(len(edges))]
    out = ctx.constant(0)
    for i,(a,b) in enumerate(edges):
        for j,(c,d) in enumerate(edges[:i+1]):
            weight = x[i]*x[j]*((n-1)*deriv[i]*deriv[j]-n*u*deriv[i].derivative(j))/n**2
            distances = (f[a][c]+f[a][d]+f[b][c]+f[b][d])/4
            if i==j:
                distances -= x[i]*u/6
            out += (1 if i==j else 2)*weight*distances
    return out


def evaluate(edges):
    graph = nx.Graph(edges)
    n = len(edges)-len(graph)+1
    branch,paths = skeleton_chains(edges,len(graph))
    ordered = sorted((tuple(sorted((branch.index(p[0]),branch.index(p[-1])))),es) for p,es in paths)
    core_edges,chains = zip(*ordered)
    ctx,u,f = forest_data(core_edges,len(branch))
    factor = transfer(core_edges,ctx,u,f,n)
    g = G(len(graph),edges)
    source = [(0,0,0)]*len(graph)
    slots = sorted(i for i,val in graph.degree if val==2)
    for i,p in zip(slots,((1,0,0),(0,1,0),(0,0,1),(-1,-1,-1))):
        source[i] = p
    big,_,_ = forest_data(edges,len(graph))
    raw = big.from_dict(numerator(g,source,tree_flows(g,source)))**2
    poly = pushforward(raw,chains,ctx)*factor
    aut = sum(1 for _ in nx.isomorphism.GraphMatcher(graph,graph).isomorphisms_iter())
    return tuple(F(12*n,aut)*v for v in integrate_core(core_edges,poly,n,2*n+6))


def check_two(j):
    from bulk_eight import potentials
    for plane in (((1,0),(2,1),(0,1),(1,-2)), ((2,1),(-1,4),(3,-2),(2,3))):
        p = [(*q,0,0,0) for q in plane]+[(0,0,1,0,0),(0,0,0,1,0),(0,0,0,0,1)]
        p.append(tuple(-sum(col) for col in zip(*p)))
        actual = potentials(p)
        prediction = response(plane,-F(7,77760),j)
        diff = [[a-b for a,b in zip(x,y)] for x,y in zip(actual,prediction)]
        common = [diff[5 if k==2 else 4][k] for k in range(5)]
        for row,q in zip(diff,p):
            delta = [x-y for x,y in zip(row,common)]
            assert all(delta[k]*q[h]==delta[h]*q[k] for k,h in combinations(range(5),2))
    print('Independent eight-primary theta responses agree for four planar momenta.',flush=True)


if __name__=='__main__':
    for n in (2,3):
        total = [F(0),F(0)]
        for graph in tight_graphs(3,n):
            if max(dict(nx.Graph(graph.edges).degree()).values()) == 3:
                value = evaluate(graph.edges)
                print(graph.name,*value,flush=True)
                for i,m in enumerate(value):
                    total[i] += m
        print('Planar transfer coefficient at loop',n,*total,flush=True)
        if n==2:
            assert total == [-F(1,20736),F(0)]
            check_two(total[0])
        else:
            assert total == [-F(2078212919,635529742848000),F(131371,49037788800)]
