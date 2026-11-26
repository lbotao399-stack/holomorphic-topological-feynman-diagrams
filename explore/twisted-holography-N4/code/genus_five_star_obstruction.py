"""First unreduced full electrical face for three planar insertions at genus five."""
from functools import cache
from fractions import Fraction
from itertools import combinations

import sympy as sp


def edge(u,v):
    return tuple(sorted((u,v)))


@cache
def reductions(edges,marked):
    """Parallel pairs are identified after each positive series or empty-star move."""
    vertices = set().union(*map(set,edges))
    if len(vertices)<=2:
        return True
    for v in vertices:
        neighbors = {j for e in edges if v in e for j in e if j!=v}
        if len(neighbors)>3 or len(neighbors)==3 and v in marked:
            continue
        reduced = frozenset(e for e in edges if v not in e) | frozenset(
            edge(i,j) for i,j in combinations(neighbors,2))
        loaded = marked-{v}
        if v in marked:
            loaded = loaded|neighbors
        if reductions(reduced,loaded):
            return True
    return False


def connected(vertices,edges):
    reached = {min(vertices)}
    for _ in vertices:
        for i,j in edges:
            if i in reached: reached.add(j)
            if j in reached: reached.add(i)
    return reached==set(vertices)


def check_smaller():
    count = 0
    for n in (4,5):
        vertices = set(range(n))
        for e in range(n+2,min(7,n*(n-1)//2)+1):
            for graph in combinations(combinations(range(n),2),e):
                if any(not connected(vertices-{v},[p for p in graph if v not in p]) for v in vertices):
                    continue
                for m in range(4):
                    for marked in combinations(range(n),m):
                        assert reductions(frozenset(graph),frozenset(marked))
                        count += 1
    return count


def check_obstruction():
    graph = tuple(e for e in combinations(range(5),2) if e not in ((0,1),(0,2)))
    assert not reductions(frozenset(graph),frozenset((0,1,2)))
    # A genus-one degree-two marked vertex, two rational trivalent marked
    # vertices, two rational tetravalent empty vertices. Both psi flags
    # sit on 34, so the branch binomial cancels the order-two automorphism.
    degrees = [sum(i in e for e in graph) for i in range(5)]
    marks = (1,1,1,0,0)
    genus = (1,0,0,0,0)
    levels = (3,1,1,0,0)
    gap = [3*h-3+d+m-l for h,d,m,l in zip(genus,degrees,marks,levels)]
    assert gap == [0,0,0,1,1] and len(graph)-len(degrees)+1+sum(genus)==5
    weight = Fraction(1,24)*2/2
    assert weight==Fraction(1,24)
    # Exact derivative of the energy vector in logarithmic conductances.
    incidence = sp.zeros(4,8)
    for e,(u,v) in enumerate(graph):
        if u<4: incidence[u,e]=1
        if v<4: incidence[v,e]=-1
    conductance = sp.diag(*range(1,9))
    inverse = (incidence*conductance*incidence.T).inv()
    charge = sp.Matrix([[1,0],[0,1],[-1,-1],[0,0]])
    current = conductance*incidence.T*inverse*charge
    response = sp.eye(8)-conductance*incidence.T*inverse*incidence
    gram = current*current.T
    jacobian = response.multiply_elementwise(gram)
    assert incidence*current==charge and all(gram[i,i]>0 for i in range(8))
    assert jacobian*sp.ones(8,1)==sp.zeros(8,1)
    scaled_minors = [conductance[i,i]*jacobian.minor_submatrix(i,i).det() for i in range(8)]
    assert set(scaled_minors)=={sp.Rational(139094453210825,42716824170539300917583554944)}
    print('Stable weight',weight,'; eight nonzero projected Jacobians, conductance times minor',scaled_minors[0])
    print('Shared empty-vertex edge current:',tuple(current[-1,:]))


if __name__=='__main__':
    print('Smaller biconnected high-cycle labeled graph/mark placements reduced:',check_smaller())
    check_obstruction()
