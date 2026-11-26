"""The fixed sl2 color channel, its alternating jet and its first exact moments."""
from fractions import Fraction
from itertools import permutations
from math import comb, factorial, prod
from pathlib import Path
from random import Random
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from matrix_unit_channel import color_function, parity, path
from rooted_trees import graphs
from wick_weights import wick_sign
from htk import Graph, SPReduction, brute
from htk.conventions import S2_kernel

H = {(0, 0): 1, (1, 1): -1}
E = {(0, 1): 1}
F = {(1, 0): 1}


def spoke_loads(graph, times, sources):
    currents = brute.currents(graph, times, sources)
    return [tuple((1 if target == graph.root else -1)*value for value in currents[k])
            for k, (origin, target) in enumerate(graph.edges)
            if graph.root in (origin, target)]


def color_check(n):
    """The only surviving graph has E and F at the ends of its interaction path."""
    survivors = []
    for edges in graphs(n):
        root = 2*n+1
        act = color_function(edges, root)
        count = 0
        for first, last in permutations(range(n+2), 2):
            outs = [H]*(n+2)
            outs[first], outs[last] = E, F
            count += bool(act([H]*(n+1)+outs))
        if count:
            survivors.append(count)
    assert survivors == [2]


def ordered_probe(n, rng):
    """Compare the physical spoke currents with all individual harmonic responses."""
    m = n+1
    graph = Graph(2*n+2, path(n))
    times = [Fraction(rng.randrange(1,24), rng.randrange(1,11)) for _ in graph.edges]
    vertices = list(range(1,2*n,2))
    endpoint = spoke_loads(graph, times, {0: (1,0), 2*n: (0,1)})
    reward = [a/c for a,c in endpoint]
    assert all(u > v for u,v in zip(reward, reward[1:]))
    response = [spoke_loads(graph, times, {v: (1,)}) for v in vertices]
    harmonic = [sum(weight*p[0] for weight,p in zip(reward, column)) for column in response]
    assert all(u > v for u,v in zip(harmonic, harmonic[1:]))
    assert all(sum(p[0] for p in column) == 1 for column in response)
    weights = tuple(range(n,0,-1))
    combined = spoke_loads(graph, times, {v: (t,) for v,t in zip(vertices,weights)})
    assert sum(p[0]*r for p,r in zip(combined,reward)) == sum(t*h for t,h in zip(weights,harmonic))
    d = comb(n,2)
    if n <= 5:
        alternated = sum(parity(order)*sum(weights[k]*h for k,h in zip(order,harmonic))**d
                         for order in permutations(range(n)))
        vandermonde = prod(u-v for j,u in enumerate(harmonic) for v in harmonic[j+1:])
        assert alternated == factorial(d)*vandermonde
    forward = color_function(graph.edges,graph.root)([H]*m+[E]+[H]*n+[F])
    reverse = color_function(graph.edges,graph.root)([H]*m+[F]+[H]*n+[E])
    assert -reverse == forward == -2**(2*n+1)
    assert forward*wick_sign(graph.edges,graph.root)*(-1)**comb(m,2) == 2**(2*n+1)
    assert (-1)**(n+1+comb(n+1,2)+comb(n+2,2)) == 1


def exact_moment(n, insertion):
    """Direct signed Wick sum including both endpoint assignments. Input jets act in slot order."""
    m = n+1
    graph = Graph(2*n+2, path(n))
    reduction = SPReduction(graph)
    evaluate = color_function(graph.edges,graph.root)
    colors = [E]+[H]*n+[F]
    momenta = [(1,0,0)]+[(0,k,0) for k in range(n,0,-1)]+[(0,0,1)]
    vertices = [0,*range(1,2*n,2),2*n]
    total = 0
    for middle in permutations(range(1,n+1)):
        for first,last in ((0,n+1),(n+1,0)):
            order = (first,)+middle+(last,)
            sources = {v:momenta[k] for v,k in zip(vertices,order)}
            color = evaluate([H]*m+[colors[k] for k in order])
            def alternating(currents):
                absorbed = [tuple((-1 if k<n else 1)*p for p in currents[3*k]) for k in range(m)]
                return sum(parity(perm)*insertion([absorbed[k] for k in perm])
                           for perm in permutations(range(m)))
            period = S2_kernel(reduction.kernel(sources,insertion=alternating,check=1),n)
            total += color*parity(order)*period
    return total*Fraction((-1)**comb(m,2)*wick_sign(graph.edges,graph.root),2)


if __name__ == '__main__':
    rng = Random(563)
    for n in range(1,9):
        if n <= 5:
            color_check(n)
        ordered_probe(n,rng)
    examples = [
        (1,lambda p:p[0][0],Fraction(2,3)),
        (1,lambda p:p[0][0]*p[1][2]-p[0][2]*p[1][0],Fraction(4,3)),
        (1,lambda p:(p[0][0]*p[1][2]-p[0][2]*p[1][0])
                    *p[0][0]*p[0][2]*p[1][0]*p[1][2],Fraction(1,45)),
        (2,lambda p:p[0][0]*p[1][1],Fraction(1,270)),
        (3,lambda p:p[0][0]*p[1][1]*p[2][1]**2,Fraction(499,29767500)),
    ]
    for n,insertion,expected in examples:
        result = exact_moment(n,insertion)
        assert result == expected
        print(n,result)
