"""Enumerate admissible marked graphs and compare currents and their doubled cofactors."""

from fractions import Fraction as F
from math import factorial, prod
from pathlib import Path
from random import Random
import sys
import networkx as nx
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "code"))
from htk import Graph, sp_reduce
from htk.catalog import enumerate_tight
from htk.conventions import S2_kernel
from kernel import kernel


def marked_graph(m):
    sink = 2 * m + 1
    edges = [(0, sink)]
    for j in range(m):
        u, v = 2 * j + 1, 2 * j + 2
        edges.extend(((0, u), (u, v), (v, sink)))
    return Graph(sink + 1, edges, "scalar marked product")


def test_graphs():
    for m in (1, 2, 3):
        model = marked_graph(m).to_networkx()
        seen = False
        for g in enumerate_tight(3, m):
            for x, y in g.edges:
                core = g.copy()
                core.remove_nodes_from((x, y))
                if any(g.degree(v) > 3 for v in core) or not nx.is_forest(core):
                    continue
                attachments = [
                    sum(g.has_edge(v, x) + g.has_edge(v, y) for v in component)
                    for component in nx.connected_components(core)
                ]
                if max(attachments) > 2:
                    continue
                assert nx.is_isomorphic(g, model)
                assert all(
                    len(c) == 2 and sum(g.has_edge(v, x) for v in c) == 1
                    for c in nx.connected_components(core)
                )
                seen = True
        assert seen


def test_kernels():
    rng = Random(193)
    for m in (1, 2, 3):
        graph = marked_graph(m)
        reduction = sp_reduce(graph)
        for trial in range(6):
            source = tuple(rng.randrange(-3, 4) for j in range(3))
            pairs = [
                tuple(tuple(rng.randrange(-3, 4) for j in range(3)) for k in range(2))
                for i in range(m)
            ]
            waves = [source] + [k for pair in pairs for k in pair]
            powers = (
                {}
                if not trial
                else {
                    (rng.randrange(3 * m + 1), rng.randrange(3)): 1,
                    (rng.randrange(3 * m + 1), rng.randrange(3)): 1,
                }
            )
            value = kernel(source, pairs, powers)
            independent = S2_kernel(
                reduction.integrate(waves, powers=powers, check=1), m
            )
            assert value == independent, (m, powers, value, independent)
    for m in range(1, 13):
        sample = [((1, 0, 0), (0, 1, 0))] * m
        assert kernel((0, 0, 1), sample) == F((-1) ** m * 2, 2**m * factorial(m + 2))
    for m in range(1, 7):
        w = [F(i + 1, 2 * (m + 1) ** 2) for i in range(m)]
        pairs = [
            tuple(sp.Matrix([rng.randint(-4, 4) for j in range(3)]) for _ in range(2))
            for _ in range(m)
        ]
        L = sp.Matrix([0, 1, 2])
        jac = sp.zeros(3 * m)
        for i, (p, q) in enumerate(pairs):
            for j, (u, v) in enumerate(pairs):
                if i == j:
                    jac[3 * i : 3 * i + 3, 3 * j] = L
                jac[3 * i : 3 * i + 3, 3 * j + 1] = (w[i] - int(i == j)) * u
                jac[3 * i : 3 * i + 3, 3 * j + 2] = (w[i] - int(i == j)) * v
        expected = (1 - sum(w)) ** 2 * prod(L.dot(p.cross(q)) for p, q in pairs)
        assert jac.det() == expected


if __name__ == "__main__":
    test_graphs()
    test_kernels()
    print(
        "marked topologies through three loops, 18 raw jet integrals, six Jacobians and twelve exact all-order moments passed"
    )
