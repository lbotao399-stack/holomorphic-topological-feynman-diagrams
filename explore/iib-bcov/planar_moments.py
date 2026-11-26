"""Direct placement determinants for the planar first and second edge moments."""
from collections import Counter
from fractions import Fraction as F
from itertools import product
from math import prod
from flint import fmpz_mat
from planar_bulk import bracket, forest_data, mass
from square_three import EDGES, G, cofactor, solve_graph


def verify(edges,nv,p):
    n = len(p)//2
    graph = G(nv,edges)
    times = tuple(2+i for i in range(len(edges)))
    columns = [[[q[d]*row[e] for d in range(2) for row in graph.C] for e in range(len(edges))]
               for q in p]
    total = 0
    first, second = Counter(), Counter()
    for order in product(range(len(edges)),repeat=len(p)):
        if max(Counter(order).values()) > 2:
            continue
        det = int(fmpz_mat([col[e] for col,e in zip(columns,order)]).det())
        coeff = det**2*prod(times[e] for e in order)
        total += coeff
        first[order[0]] += coeff
        second[order[:2]] += coeff
    _,u,_ = forest_data(edges,nv)
    tree = int(u(*times))
    z = mass(p)
    pair = bracket(p[0],p[1])**2*mass(p[2:])
    b = F(z-n*n*pair,n*(n-1))
    a = b+pair
    assert total==tree**2*z
    for e,t in enumerate(times):
        ue = u.derivative(e)
        assert first[e]==F(tree*z*t*int(ue(*times)),n)
        for h,s in enumerate(times):
            d = int(ue(*times)*u.derivative(h)(*times))
            assert second[e,h] == t*s*(a*d-b*(d-tree*int(ue.derivative(h)(*times))))
    print(n,'loops: mass and all one-edge and ordered edge-pair weights verified.',flush=True)


def square_family(limit):
    edges = [(0,1),(1,2),(2,3),(3,0)]
    left,right = 1,2
    for n in range(1,limit+1):
        nv = 2*n+2
        slots = [i for i in range(nv) if sum(i in e for e in edges)==2]
        sources = [(0,0,0)]*nv
        for i,q in zip(slots,((1,0,0),(0,1,0),(0,0,1),(-1,-1,-1))):
            sources[i] = q
        times = [1+(7*i+3*i*i)%29 for i in range(len(edges))]
        flow,_ = solve_graph(edges,sources,times,0)
        assert cofactor(flow,G(nv,edges).C,times)>0
        edges.extend(((left,nv),(nv,nv+1),(nv+1,right)))
        left,right = nv,nv+1
    print('Four retained sources give nonzero threefold minors through loop',limit)


if __name__=='__main__':
    verify(((0,1),(0,1),(0,1)),2,((1,0),(2,1),(0,1),(1,-2)))
    verify(EDGES,4,((1,0),(0,1),(1,1),(1,-1),(2,1),(1,2)))
    square_family(8)
