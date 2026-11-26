"""Marked product-forest graphs, exact canonical labels, and source validation.

The cubic encoding and enumerator adapt the existing implementations in
``explore/circuit-general/ht-d2/forests.py`` and
``explore/nonSP/verify/full_ward/forest_orders.py``.  An interaction is an
unordered pair of ports, each a child interaction or the boundary leaf
``'A'``/``'B'``; the pair ``('A', 'B')`` is the one-vertex tree.  The two
boundaries remain ordered.  See sources/LLC_main.tex, Theorems b-00249 and
b-00405.  This module evaluates no integral and discards no source.

For H=2,T=0, vertex relabeling and edge reversal preserve the bare kernel.
Edge permutation contributes its parity; reversed edge currents separately
contribute a sign to each inserted current component.  These are the
conventions of code/htk/CONVENTIONS.md, Section 2.2.
"""
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from math import factorial
from operator import index

from .graph import Graph

__all__ = [
    "CubicForest", "ForestNode", "ProductForestStructure", "NotProductForest",
    "NonCubicProductForest", "UnsupportedForestSources", "recognize_product_forest",
    "trees", "forests", "forest_counts", "forest_graph", "tree_size", "tree_label",
    "tree_automorphisms", "forest_automorphisms",
]


class NotProductForest(ValueError):
    """The marked graph does not satisfy the product-forest conditions."""


class NonCubicProductForest(NotImplementedError):
    """A valid higher-valence product forest needs the sourced reduction method."""

    def __init__(self, structure):
        self.structure = structure
        super().__init__("this product forest is not cubic: higher-valence interactions and "
                         "multiple surviving sources require sourced reductions")


class UnsupportedForestSources(NotImplementedError):
    """The cubic source-free-star method cannot eliminate a non-root source."""


def _integer(value, name, minimum=0):
    if isinstance(value, bool):
        raise TypeError(name + " must be an integer")
    try:
        value = index(value)
    except TypeError as error:
        raise TypeError(name + " must be an integer") from error
    if value < minimum:
        raise ValueError(name + " must be at least " + str(minimum))
    return value


def _parity(permutation):
    seen = set()
    result = 1
    for start in range(len(permutation)):
        if start in seen:
            continue
        vertex, length = start, 0
        while vertex not in seen:
            seen.add(vertex)
            length += 1
            vertex = permutation[vertex]
        if length % 2 == 0:
            result = -result
    return result


def _tree_key(tree):
    return (0, tree) if isinstance(tree, str) else (1, _tree_key(tree[0]), _tree_key(tree[1]))


def _canonical_tree(tree):
    if isinstance(tree, str):
        if tree not in ("A", "B"):
            raise ValueError("boundary leaves must be 'A' or 'B'")
        return tree
    if not isinstance(tree, (tuple, list)) or len(tree) != 2:
        raise ValueError("every interaction must be an unordered pair of ports")
    pair = tuple(sorted((_canonical_tree(port) for port in tree), key=_tree_key))
    if pair in (("A", "A"), ("B", "B")):
        raise ValueError("parallel spokes to one boundary are not allowed")
    return pair


def tree_size(tree):
    """Number of interaction vertices in a binary-port tree."""
    return 0 if isinstance(tree, str) else 1 + tree_size(tree[0]) + tree_size(tree[1])


def tree_label(tree):
    return tree if isinstance(tree, str) else "(" + ",".join(map(tree_label, tree)) + ")"


def trees(n):
    """Canonical cubic rooted product trees with exactly n interactions."""
    n = _integer(n, "tree size")
    return _trees(n)


@lru_cache(maxsize=None)
def _trees(n):
    if n == 0:
        return ()
    if n == 1:
        return (("A", "B"),)
    result = set()
    for child in trees(n - 1):
        result.add(("A", child))
        result.add(("B", child))
    for first_size in range(1, (n - 1) // 2 + 1):
        second_size = n - 1 - first_size
        for first in trees(first_size):
            for second in trees(second_size):
                result.add(tuple(sorted((first, second), key=_tree_key)))
    return tuple(sorted(result, key=_tree_key))


def forests(n):
    """Generate each multiset of cubic rooted trees of total size n once.

    Enumeration is intrinsically exponential in n; ``forest_counts`` obtains
    the counts without constructing any graphs or trees.
    """
    n = _integer(n, "loop number")
    available = tuple(sorted((tree for size in range(1, n + 1) for tree in trees(size)), key=_tree_key))
    sizes = tuple(map(tree_size, available))

    def generate(remaining, start, selected):
        if remaining == 0:
            yield tuple(selected)
            return
        for position in range(start, len(available)):
            size = sizes[position]
            if size <= remaining:
                selected.append(available[position])
                yield from generate(remaining - size, position, selected)
                selected.pop()

    yield from generate(n, 0, [])


def forest_counts(n):
    """Return exact tree and forest counts through n without enumeration.

    The tree recurrence counts an A/B unary extension or an unordered pair
    of nonempty children.  The forest recurrence is the Euler transform for
    multisets. Both lists include index zero: zero trees, one empty forest.
    """
    n = _integer(n, "loop number")
    tree_counts = [0] * (n + 1)
    if n:
        tree_counts[1] = 1
    for size in range(2, n + 1):
        pairs = sum(tree_counts[left] * tree_counts[size - 1 - left]
                    for left in range(1, size - 1))
        if (size - 1) % 2 == 0:
            pairs += tree_counts[(size - 1) // 2]
        tree_counts[size] = 2 * tree_counts[size - 1] + pairs // 2
    divisor_sums = [0] * (n + 1)
    for size in range(1, n + 1):
        for multiple in range(size, n + 1, size):
            divisor_sums[multiple] += size * tree_counts[size]
    counts = [1] + [0] * n
    for size in range(1, n + 1):
        numerator = sum(divisor_sums[k] * counts[size - k] for k in range(1, size + 1))
        if numerator % size:
            raise ArithmeticError("the exact forest recurrence did not divide integrally")
        counts[size] = numerator // size
    return tuple(tree_counts), tuple(counts)


def tree_automorphisms(tree):
    if isinstance(tree, str):
        return 1
    return (tree_automorphisms(tree[0]) * tree_automorphisms(tree[1])
            * (2 if tree[0] == tree[1] else 1))


def forest_automorphisms(forest):
    result = 1
    for tree, multiplicity in Counter(forest).items():
        result *= tree_automorphisms(tree) ** multiplicity * factorial(multiplicity)
    return result


def forest_graph(forest):
    """Build (Graph, roots) in canonical preorder, A=n, B=n+1, sorted edges."""
    forest = tuple(sorted((_canonical_tree(tree) for tree in forest), key=_tree_key))
    if any(isinstance(tree, str) for tree in forest):
        raise ValueError("a forest component must contain an interaction")
    n = sum(map(tree_size, forest))
    A, B = n, n + 1
    edges, roots = [(A, B)], []
    next_vertex = 0

    def visit(tree, parent=None):
        nonlocal next_vertex
        vertex = next_vertex
        next_vertex += 1
        if parent is None:
            roots.append(vertex)
        else:
            edges.append((parent, vertex))
        for port in tree:
            if isinstance(port, str):
                edges.append((vertex, A if port == "A" else B))
            else:
                visit(port, vertex)

    for tree in forest:
        visit(tree)
    return Graph(n + 2, sorted(edges)), tuple(roots)


@dataclass(frozen=True)
class ProductForestStructure:
    marked: tuple
    components: tuple
    degree_two_vertices: tuple
    higher_valence_vertices: tuple
    cubic: bool


def recognize_product_forest(G, marked):
    """Recognize the broader marked product topology without assuming cubic sources.

    On a connected interaction subtree, the H=2 sparsity condition involving
    both boundaries is ``sum(number_of_spokes - 1) <= 1``. A maximum connected
    subtree calculation verifies every such subset in linear time. Subsets
    containing at most one boundary already obey the bound because the
    interaction graph is a forest and no parallel spokes are allowed.
    """
    if not isinstance(G, Graph):
        raise TypeError("G must be an htk.Graph")
    try:
        marked = tuple(marked)
    except TypeError as error:
        raise ValueError("marked must be the ordered pair (A, B)") from error
    if len(marked) != 2:
        raise ValueError("marked must be the ordered pair (A, B)")
    A, B = (_integer(vertex, "marked vertex") for vertex in marked)
    if A == B or A >= G.nv or B >= G.nv:
        raise ValueError("marked vertices must be distinct vertices of G")
    edges = [frozenset(edge) for edge in G.edges]
    if len(set(edges)) != len(edges):
        raise NotProductForest("parallel edges are not product-forest edges")
    if frozenset((A, B)) not in edges:
        raise NotProductForest("the marked edge A--B is absent")
    if G.ne != 2 * G.nv - 3:
        raise NotProductForest("a product forest has n+2 vertices and 2n+1 edges")
    adjacent = [set() for _ in range(G.nv)]
    for u, v in G.edges:
        adjacent[u].add(v)
        adjacent[v].add(u)
    internal = set(range(G.nv)) - {A, B}
    remaining = set(internal)
    components, degree_two, high = [], [], []
    while remaining:
        start = min(remaining)
        parent, order = {start: None}, [start]
        for vertex in order:
            for neighbor in sorted(adjacent[vertex] & internal):
                if neighbor == parent[vertex]:
                    continue
                if neighbor in parent:
                    raise NotProductForest("deleting A and B does not leave a forest")
                parent[neighbor] = vertex
                order.append(neighbor)
        component = tuple(sorted(order))
        remaining.difference_update(component)
        if any(len(adjacent[vertex]) < 2 for vertex in component):
            raise NotProductForest("interaction vertices must have full graph degree at least two")
        spokes = {vertex: len(adjacent[vertex] & {A, B}) for vertex in component}
        if (not any(A in adjacent[vertex] for vertex in component)
                or not any(B in adjacent[vertex] for vertex in component)):
            raise NotProductForest("every interaction tree must meet both A and B")
        if sum(spokes.values()) != len(component) + 1:
            raise NotProductForest("each interaction tree must have one more spoke than vertices")
        best = {}
        for vertex in reversed(order):
            best[vertex] = spokes[vertex] - 1 + sum(
                max(0, best[child]) for child in adjacent[vertex] & internal if parent.get(child) == vertex)
            if best[vertex] > 1:
                raise NotProductForest("an induced interaction subtree with A,B violates the H=2 subgraph bound")
        roots = tuple(vertex for vertex in component if len(adjacent[vertex]) == 2)
        higher = tuple(vertex for vertex in component if len(adjacent[vertex]) > 3)
        if len(higher) > 1:
            raise NotProductForest("a product tree permits at most one higher-valence superpotential vertex")
        components.append(component)
        degree_two.append(roots)
        high.append(higher)
    cubic = all(len(roots) == 1 and not higher for roots, higher in zip(degree_two, high))
    return ProductForestStructure((A, B), tuple(components), tuple(degree_two), tuple(high), cubic)


@dataclass(frozen=True)
class ForestNode:
    vertex: int
    original_vertex: int
    ports: tuple
    encoding: tuple


class CubicForest:
    """Validated cubic marked forest with exact maps to canonical graph order."""

    def __init__(self, G, marked):
        structure = recognize_product_forest(G, marked)
        if not structure.cubic:
            raise NonCubicProductForest(structure)
        self.graph = self.G = Graph(G.nv, G.edges, G.name)
        self.structure = structure
        self.marked = structure.marked
        self.loops = G.nv - 2
        self.A, self.B = self.loops, self.loops + 1
        old_A, old_B = self.marked
        adjacent = [set() for _ in range(G.nv)]
        for u, v in G.edges:
            adjacent[u].add(v)
            adjacent[v].add(u)
        encodings, ports = {}, {}
        roots = [vertices[0] for vertices in structure.degree_two_vertices]
        for root in roots:
            parent, order = {root: None}, [root]
            for vertex in order:
                for child in sorted(adjacent[vertex] - {old_A, old_B}):
                    if child != parent[vertex]:
                        parent[child] = vertex
                        order.append(child)
            for vertex in reversed(order):
                items = [("A", "A")] if old_A in adjacent[vertex] else []
                if old_B in adjacent[vertex]:
                    items.append(("B", "B"))
                items += [(encodings[child], child) for child in adjacent[vertex] - {old_A, old_B}
                          if parent.get(child) == vertex]
                items.sort(key=lambda item: (_tree_key(item[0]), -1 if isinstance(item[1], str) else item[1]))
                if len(items) != 2:
                    raise ArithmeticError("a cubic rooted interaction does not have two outgoing ports")
                encodings[vertex] = tuple(item[0] for item in items)
                ports[vertex] = tuple(item[1] for item in items)
        roots.sort(key=lambda vertex: (_tree_key(encodings[vertex]), vertex))
        preorder, stack = [], list(reversed(roots))
        while stack:
            vertex = stack.pop()
            preorder.append(vertex)
            stack.extend(reversed([port for port in ports[vertex] if not isinstance(port, str)]))
        vertex_order = tuple(preorder) + (old_A, old_B)
        inverse_vertices = [0] * G.nv
        for canonical, original in enumerate(vertex_order):
            inverse_vertices[original] = canonical
        self.original_to_canonical_vertices = tuple(inverse_vertices)
        self.canonical_to_original_vertices = vertex_order
        self.forest = tuple(encodings[root] for root in roots)
        self.original_roots = tuple(roots)
        self.canonical_roots = tuple(inverse_vertices[root] for root in roots)
        self.canonical_nodes = tuple(ForestNode(
            inverse_vertices[vertex], vertex,
            tuple(port if isinstance(port, str) else inverse_vertices[port] for port in ports[vertex]),
            encodings[vertex]) for vertex in preorder)
        mapped = [(inverse_vertices[u], inverse_vertices[v]) for u, v in G.edges]
        canonical_edges = sorted(tuple(sorted(edge)) for edge in mapped)
        self.canonical_graph = Graph(G.nv, canonical_edges)
        positions = {edge: i for i, edge in enumerate(canonical_edges)}
        edge_map = tuple(positions[tuple(sorted(edge))] for edge in mapped)
        inverse_edges = [0] * G.ne
        for original, canonical in enumerate(edge_map):
            inverse_edges[canonical] = original
        self.original_to_canonical_edges = edge_map
        self.canonical_to_original_edges = tuple(inverse_edges)
        self.original_edge_signs = tuple(1 if u < v else -1 for u, v in mapped)
        self.kernel_sign = _parity(edge_map)

    def normalize_sources(self, lam):
        """Return full balanced Fraction pairs in canonical vertex order.

        Inputs follow htk's source convention: nv-1 rows omit the original
        grounded vertex, nv rows include its checked source, and dictionary
        entries omit zero vertices. Every non-root interaction source is
        checked; unsupported nonzero sources are never silently removed.
        """
        count = self.graph.nv

        def row(values):
            values = tuple(values)
            if len(values) != 2:
                raise ValueError("a C^2 source must have exactly two components")
            if any(isinstance(value, float) for value in values):
                raise TypeError("source components must be exact rational values")
            return tuple(Fraction(value) for value in values)

        if lam is None:
            full = [(Fraction(0), Fraction(0))] * count
        elif isinstance(lam, dict):
            if any(not isinstance(vertex, int) or isinstance(vertex, bool) or not 0 <= vertex < count for vertex in lam):
                raise ValueError("source dictionary contains an invalid vertex")
            full = [row(lam.get(vertex, (0, 0))) for vertex in range(count - 1)]
            ground = tuple(-sum(values[i] for values in full) for i in range(2))
            if count - 1 in lam and row(lam[count - 1]) != ground:
                raise ValueError("sources do not sum to zero")
            full.append(ground)
        else:
            full = [row(values) for values in lam]
            if len(full) == count - 1:
                full.append(tuple(-sum(values[i] for values in full) for i in range(2)))
            elif len(full) != count:
                raise ValueError("sources require nv-1 or nv two-component rows")
            if any(sum(values[i] for values in full) for i in range(2)):
                raise ValueError("sources do not sum to zero")
        allowed = set(self.original_roots) | set(self.marked)
        unsupported = tuple(vertex for vertex, values in enumerate(full) if vertex not in allowed and any(values))
        if unsupported:
            raise UnsupportedForestSources("nonzero sources at non-root cubic vertices %s require sourced star reductions" % (unsupported,))
        return tuple(full[vertex] for vertex in self.canonical_to_original_vertices)

    def map_powers(self, powers=None):
        """Map {(original edge, component): power} and return its complete sign."""
        mapped, sign = {}, self.kernel_sign
        for location, power in ({} if powers is None else powers).items():
            try:
                edge, component = location
            except (TypeError, ValueError) as error:
                raise ValueError("current powers are indexed by (edge, component)") from error
            edge = _integer(edge, "edge index")
            component = _integer(component, "current component")
            power = _integer(power, "current power")
            if edge >= self.graph.ne or component >= 2:
                raise ValueError("edge or current component lies outside this graph")
            if power:
                mapped[self.original_to_canonical_edges[edge], component] = power
                sign *= self.original_edge_signs[edge] ** power
        return mapped, sign
