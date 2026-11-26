"""Independent topology enumeration and exact polarized graph checks."""

from collections import Counter
from fractions import Fraction as F

import networkx as nx

from kernels import Rooted, rooted_two, theta, typed_symbol
from rooted_graphs import canonical, graphs
from htk import polygon
from htk.catalog import enumerate_tight


def check_topologies():
    for n, expected in enumerate(((1, 0), (2, 3), (5, 18), (15, 95)), 1):
        generated = {key for _, _, key in graphs(n)}
        assert Counter(k[0] for k in generated) == Counter(T=expected[0], C=expected[1])
        if n > 3:
            continue
        independent = set()
        for graph in enumerate_tight(3, n):
            for root in graph.nodes():
                if graph.degree(root) not in (n, n + 1):
                    continue
                if any(graph.degree(v) > 3 for v in graph if v != root):
                    continue
                interior = graph.copy()
                interior.remove_node(root)
                assert nx.is_connected(interior)
                independent.add(canonical(interior, graph[root]))
        assert generated == independent


def check_color():
    basis = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    probes = ((1, 0, 1), (0, 1, 2), (1, 2, 0))
    box = Rooted(polygon(4), 3)
    shapes = rooted_two()
    expected = (F(-1, 180), F(7, 2160), 0)
    for axis in range(3):
        derivative = ((), (axis,))
        outgoing = list(zip(basis, (0, 1, 2), (True, False, False)))
        assert sum(box.gauge((i, i), derivative, outgoing) for i in range(3)) == 0
        outgoing = list(zip(basis, (0, 1, 2), (False,) * 3))
        assert (
            tuple(
                sum(g.gauge((i, i), derivative, outgoing) for i in range(3))
                for g in shapes
            )
            == expected
        )
    assert sum(expected) == F(-1, 432)
    for spec, want in zip(((1, 3, 3), (2, 2, 3)), (F(6835, 864), F(355, 54))):
        graph = Rooted(theta(*spec)[0], 0)
        outgoing = list(
            zip(probes + ((2, 1, -1),), (0, 0, 1, 1), (True, False, False, False))
        )
        assert graph.gauge((0, 1, 2), ((), (0,), ()), outgoing) == want

    def evaluate(root, kinds, **kwargs):
        return typed_symbol(
            box, root, ((), (0,)), [(p, *k) for p, k in zip(probes, kinds)], **kwargs
        )

    assert evaluate([("b", 0), ("b", 0)], [("b", 0), ("a", 1), ("a", 2)]) == F(-5, 12)
    assert evaluate(
        [("g", None)] * 2, [("g", None), ("h", None), ("h", None)], pure_poisson=True
    ) == F(1, 12)
    assert evaluate([("b", 0), ("g", None)], [("b", 1), ("a", 2), ("h", None)]) == F(
        2, 3
    )
    for graph, want in zip(
        shapes, (F(8411, 297675), F(110939, 14288400), F(72283, 14288400))
    ):
        insertion = [("g", None)] * 2
        outgoing = [(p, "h", None) for p in probes]
        assert (
            typed_symbol(graph, insertion, ((), (0,)), outgoing, pure_poisson=True)
            == want / 4
        )
        assert typed_symbol(graph, insertion, ((), (0,)), outgoing) == want
    assert F(8411, 297675) + F(110939, 14288400) + F(72283, 14288400) == F(559, 13608)


if __name__ == "__main__":
    check_topologies()
    check_color()
    print(
        "topologies through four loops; all five two-loop rooted classes; gauge, gravity and mixed contractions passed"
    )
