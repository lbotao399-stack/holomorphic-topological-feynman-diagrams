"""Reconstruct the full K3,3 kernel for arbitrary conserved two-component sources."""
from math import prod

import networkx as nx
import sympy as sp

from replay import EP, INJECTIONS, K33, graph_terms, integrate, master_tables, run
from check import EDGE_INDEX, permutation_sign

# Each seed is a product of four brackets in the graph's vertex labels.
SEEDS = (
    ((0, 2), (0, 2), (0, 2), (1, 3)),
    ((0, 1), (0, 1), (0, 2), (0, 3)),
    ((0, 1), (0, 1), (0, 2), (2, 3)),
    ((0, 1), (0, 1), (0, 2), (3, 4)),
    ((0, 1), (0, 2), (1, 2), (3, 4)),
)
PROBES = (
    INJECTIONS,
    [(2, 0), (-1, -1), (0, 0), (0, 1), (-1, 1), (0, -1)],
    [(0, 1), (0, 0), (-1, -1), (2, 1), (-1, 0), (0, -1)],
    [(1, -1), (0, 1), (-1, -1), (1, 1), (-1, 1), (0, -1)],
    [(0, -1), (0, 1), (0, 1), (-1, 1), (-1, -1), (2, -1)],
    [(0, -2), (-1, 0), (1, -1), (1, 1), (1, 1), (-2, 1)],
    [(-3, -1), (-2, 2), (0, -1), (1, -1), (3, 1), (1, 0)],
)


def exact_coefficients():
    z3, z5 = sp.symbols('z3 z5')
    return [
        (193 + 3072*z3/7 - 4800*z5/7)/138240,
        (11189 + 28664*z3 - 44160*z5)/806400,
        (-1739 - sp.Rational(24168, 5)*z3 + 7296*z5)/96768,
        (-637 - sp.Rational(60696, 35)*z3 + sp.Rational(18432, 7)*z5)/13824,
        sp.Rational(1, 3456),
    ]


def character_and_evaluations():
    graph = nx.Graph(K33)
    rows = sp.zeros(len(PROBES), len(SEEDS))
    character = 0
    for g in nx.isomorphism.GraphMatcher(graph, graph).isomorphisms_iter():
        edge_sign = permutation_sign([EDGE_INDEX[frozenset((g[i], g[j]))] for i, j in K33])
        # h_k on the conserved source space. Multiply the complete-symmetric
        # generating function for the permutation representation by (1-t).
        h = [1, -1, 0, 0, 0, 0]
        seen = set()
        for v in range(6):
            if v in seen:
                continue
            length = 0
            while v not in seen:
                seen.add(v)
                length += 1
                v = g[v]
            for k in range(length, 6):
                h[k] += h[k - length]
        character += edge_sign * (h[4]**2 - h[3]*h[5])
        for i, source in enumerate(PROBES):
            assert tuple(map(sum, zip(*source))) == (0, 0)
            def bracket(pair):
                a,b = [source[g[v]] for v in pair]
                return a[0]*b[1] - a[1]*b[0]
            for j, seed in enumerate(SEEDS):
                rows[i,j] += edge_sign * prod(map(bracket, seed))
    assert character == 360 and rows[:5,:].det() == 2654208
    return rows


def solve(directory, executable):
    matrix = character_and_evaluations()
    tables = master_tables(directory, executable)
    values = []
    for i, source in enumerate(PROBES):
        terms = graph_terms(source, source_factor=1)
        pair = [integrate(directory, executable, tables, expression, f'{i}_{cut}').subs(EP, 0)
                for cut, expression in enumerate(terms)]
        assert sp.expand(pair[0] - pair[1]) == 0
        values.append(pair[0])
        print('probe', i, '=', values[-1], flush=True)
    coefficients = (matrix[:5,:].inv() * sp.Matrix(values[:5])).applyfunc(sp.expand)
    for i in (5, 6):
        assert sp.expand((matrix[i,:] * coefficients)[0] - values[i]) == 0
    assert sp.expand(values[6] + sp.Rational(789,640) - sp.Rational(59,80)*sp.Symbol('z3')) == 0
    assert all((matrix[5,:] * matrix[:5,:].inv()))
    expected = exact_coefficients()
    assert all(sp.expand(c-e) == 0 for c,e in zip(coefficients, expected))
    for seed, coefficient in zip(SEEDS, coefficients):
        print(seed, '=', coefficient, flush=True)


if __name__ == '__main__':
    run(solve)
