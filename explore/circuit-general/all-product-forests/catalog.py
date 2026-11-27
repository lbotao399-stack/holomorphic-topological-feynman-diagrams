"""Complete marked product-forest topologies, with arbitrary superpotential degree.

The two insertion vertices are ordered.  Interaction trees are unrooted;
their canonical forms retain the A/B spoke at every vertex.  Components
are unordered.  No source choice or current word enters this census.
"""
from functools import lru_cache
from itertools import product
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "code"))

from htk import Graph
from htk.forest_model import NotProductForest, recognize_product_forest


def size(tree):
    return 1 + sum(size(child) for child in tree[1])


def canonical_tree(adjacency, spokes):
    """AHU encoding of an unrooted tree with fixed, ordered boundary colors."""
    def rooted(vertex, parent):
        return (spokes[vertex], tuple(sorted(rooted(child, vertex)
                for child in adjacency[vertex] if child != parent)))
    # Taking all roots avoids choosing either vertex of a bicentered tree.
    return min(rooted(vertex, None) for vertex in range(len(adjacency)))


def forest_graph(forest):
    count = sum(size(tree) for tree in forest)
    A, B = count, count + 1
    edges = []
    next_vertex = 0

    def append(tree, parent=None):
        nonlocal next_vertex
        vertex = next_vertex
        next_vertex += 1
        if parent is not None:
            edges.append((parent, vertex))
        if tree[0] & 1:
            edges.append((vertex, A))
        if tree[0] & 2:
            edges.append((vertex, B))
        for child in tree[1]:
            append(child, vertex)
    for tree in forest:
        append(tree)
    edges.append((A, B))
    return Graph(count + 2, sorted(tuple(sorted(edge)) for edge in edges)), (A, B)


@lru_cache(maxsize=None)
def trees(count):
    """Every physical interaction-tree topology with exactly count vertices.

    At most one vertex may have full degree above three.  All admissible
    spoke assignments are generated, and the exact connected-subtree
    sparsity criterion is imposed before canonicalization.
    """
    if isinstance(count, bool) or not isinstance(count, int) or count < 1:
        raise ValueError("the number of interaction vertices must be positive")
    import networkx as nx
    underlying = [nx.empty_graph(1)] if count == 1 else nx.nonisomorphic_trees(count)
    result = set()
    for graph in underlying:
        adjacency = tuple(tuple(sorted(graph.neighbors(vertex))) for vertex in range(count))
        allowed = [tuple(mask for mask in range(4)
                         if len(adjacency[vertex]) + mask.bit_count() >= 2)
                   for vertex in range(count)]
        for spokes in product(*allowed):
            if sum(mask.bit_count() for mask in spokes) != count + 1:
                continue
            if sum(len(adjacency[vertex]) + mask.bit_count() > 3
                   for vertex, mask in enumerate(spokes)) > 1:
                continue
            edges = list(graph.edges())
            for vertex, mask in enumerate(spokes):
                if mask & 1:
                    edges.append((vertex, count))
                if mask & 2:
                    edges.append((vertex, count + 1))
            edges.append((count, count + 1))
            candidate = Graph(count + 2, edges)
            try:
                recognize_product_forest(candidate, (count, count + 1))
            except NotProductForest:
                continue
            result.add(canonical_tree(adjacency, spokes))
    return tuple(sorted(result))


@lru_cache(maxsize=None)
def forests(count):
    """Every unordered forest with ordered insertion vertices and count loops."""
    if isinstance(count, bool) or not isinstance(count, int) or count < 0:
        raise ValueError("the loop number must be a nonnegative integer")
    if count == 0:
        return ((),)
    choices = sorted((tree for order in range(1, count + 1) for tree in trees(order)))
    result = []

    def visit(remaining, start, prefix):
        if remaining == 0:
            result.append(tuple(prefix))
            return
        for index in range(start, len(choices)):
            tree = choices[index]
            order = size(tree)
            if order <= remaining:
                visit(remaining - order, index, prefix + [tree])
    visit(count, 0, [])
    return tuple(result)


def count_table(last):
    result = []
    for order in range(1, last + 1):
        cubic = sum(recognize_product_forest(*forest_graph(forest)).cubic
                    for forest in forests(order))
        result.append({"loops": order, "trees": len(trees(order)),
                       "forests": len(forests(order)), "cubic_forests": cubic,
                       "higher_valence_forests": len(forests(order)) - cubic})
    return result


if __name__ == "__main__":
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("last", type=int, nargs="?", default=6)
    arguments = parser.parse_args()
    print(json.dumps(count_table(arguments.last), indent=2))
