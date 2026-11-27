"""Exact source-complete leaf elimination on every marked product forest.

The topology is compiled once.  Removing the marked vertices leaves trees;
postorder elimination consequently gives every removed vertex at most the
three retained neighbours parent,A,B.  Intermediate graphs need not be
Laman.  Positive conductances are combined by their exact Schur complement,
and EVERY original source is transferred with the same coefficients.

The result is an unexpanded source-linear response for every original edge,
the full current differential at fixed original sources, and the H=2,T=0
current-area density.  This is an exact pointwise network compiler, not an
assertion that every resulting resistance integral has been evaluated.

The local rule is the sourced circuit rule of sources/LLC_main.tex,
``local-B-sourced-current-map`` and ``local-B-sourced-operator``; its global
applicability follows from the product-forest theorem ``b-00249``.  The
density normalization is code/htk/CONVENTIONS.md Sections 2 and 2.3.
No spanning trees, expanded Kirchhoff polynomial, or period tables are
enumerated by this implementation.
"""
from collections.abc import Mapping
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import prod
from operator import index

from .forest_model import recognize_product_forest
from .graph import Graph
from .linalg import det

__all__ = ["ForestKronPlan", "ForestKronResponse", "KronStep", "KronTransfer",
           "KronBudgetExceeded"]


class KronBudgetExceeded(RuntimeError):
    """A complete coefficient buffer exceeds its explicit budget."""


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


def _exact(value, name):
    if isinstance(value, (float, complex, bool)):
        raise TypeError(name + " must be an exact rational value")
    return Fraction(value)


def _sources(graph, sources):
    def row(values):
        values = tuple(values)
        if len(values) != 2:
            raise ValueError("each H=2 source must retain both spinor entries")
        return tuple(_exact(value, "source") for value in values)

    if sources is None:
        return ((Fraction(0), Fraction(0)),) * graph.nv
    if isinstance(sources, Mapping):
        for vertex in sources:
            if _integer(vertex, "source vertex") >= graph.nv:
                raise ValueError("source vertex lies outside the original graph")
        rows = [row(sources.get(vertex, (0, 0))) for vertex in range(graph.nv - 1)]
        ground = tuple(-sum(values[a] for values in rows) for a in range(2))
        if graph.root in sources and row(sources[graph.root]) != ground:
            raise ValueError("original sources do not sum to zero")
        return tuple(rows) + (ground,)
    rows = [row(values) for values in sources]
    if len(rows) == graph.nv - 1:
        rows.append(tuple(-sum(values[a] for values in rows) for a in range(2)))
    elif len(rows) != graph.nv:
        raise ValueError("sources require nv-1 or nv complete spinor rows")
    if any(sum(values[a] for values in rows) for a in range(2)):
        raise ValueError("original sources do not sum to zero")
    return tuple(rows)


def _edge(a, b):
    return (a, b) if a < b else (b, a)


@dataclass(frozen=True)
class KronStep:
    vertex: int
    neighbors: tuple
    new_edges: tuple
    parallel_edges: tuple
    remaining_vertices: int
    remaining_edges: int


@dataclass(frozen=True)
class KronTransfer:
    """Full source transfer at one eliminated vertex, in original source labels.

    ``source_coefficients[v]`` multiplies the original source at v<nv-1.
    The original last source is balanced, never discarded.  Each retained
    neighbour receives ``weights[j] * source_coefficients``.  Backwards,
    phi_vertex = source_coefficients/pivot + sum_j weights[j]*phi_neighbor.
    """
    vertex: int
    neighbors: tuple
    pivot: Fraction
    weights: tuple
    source_coefficients: tuple


@dataclass(frozen=True)
class ForestKronResponse:
    graph: Graph
    marked: tuple
    resistances: tuple
    current_response: tuple
    potential_response: tuple
    transfers: tuple
    retained_source_response: tuple
    retained_conductance: Fraction
    laplacian_determinant: Fraction
    U: Fraction
    cotree: tuple
    max_entries: int
    stats: dict

    @property
    def pivots(self):
        return tuple(transfer.pivot for transfer in self.transfers)

    def _budget(self, count, label):
        if count > self.max_entries:
            raise KronBudgetExceeded(label + " requires %d exact entries, above max_entries=%d; no partial output is returned" %
                                     (count, self.max_entries))

    def _apply(self, response, sources):
        rows = _sources(self.graph, sources)
        return tuple(tuple(sum((coefficient * rows[vertex][a] for vertex, coefficient in enumerate(row)), Fraction(0))
                           for a in range(2)) for row in response)

    def currents(self, sources):
        """Every original oriented edge current, with both spinor entries."""
        return self._apply(self.current_response, sources)

    def potentials(self, sources):
        """Every original potential, grounded at the original last vertex."""
        return self._apply(self.potential_response, sources)

    def transferred_sources(self, sources):
        """Complete removed-vertex sources and final ordered A,B source rows."""
        removed = self._apply(tuple(transfer.source_coefficients for transfer in self.transfers), sources)
        retained = self._apply(self.retained_source_response, sources)
        return {"removed_vertices": tuple(transfer.vertex for transfer in self.transfers),
                "removed_sources": removed, "retained_vertices": self.marked,
                "retained_sources": retained}

    def _resistance_factor(self, edge, resistance_edge):
        tail, head = self.graph.edges[resistance_edge]
        response = self.current_response[edge]
        difference = ((response[tail] if tail != self.graph.root else 0)
                      - (response[head] if head != self.graph.root else 0))
        return (difference - int(edge == resistance_edge)) / self.resistances[resistance_edge]

    def source_jacobian(self):
        """All coefficient derivatives, ordered (edge,source vertex,time edge).

        This differentiates the complete original network at fixed ORIGINAL
        sources.  It therefore retains every differential of the transferred
        sources, without expanding a global numerator or a dual-number tensor.
        """
        edges, vertices = self.graph.ne, self.graph.nv - 1
        self._budget(edges * vertices * edges, "the complete source-response derivative tensor")
        result = []
        for edge in range(edges):
            factors = tuple(self._resistance_factor(edge, f) for f in range(edges))
            result.append(tuple(tuple(factors[f] * self.current_response[f][vertex] for f in range(edges))
                                for vertex in range(vertices)))
        return tuple(result)

    def current_jacobian(self, sources):
        """All current derivatives, ordered (original edge,spinor entry,time edge)."""
        edges = self.graph.ne
        self._budget(2 * edges * edges, "the complete current differential")
        currents = self.currents(sources)
        result = []
        for edge in range(edges):
            factors = tuple(self._resistance_factor(edge, f) for f in range(edges))
            result.append(tuple(tuple(factors[f] * currents[f][a] for f in range(edges)) for a in range(2)))
        return tuple(result)

    def density(self, sources, *, gauge_edge=None):
        """Exact htk H=2,T=0 scalar D(t), with every source-shift differential.

        The returned coefficient multiplies i_E(dt_0 wedge ... wedge dt_N-1).
        Gauge-edge independence and all original orientation signs follow
        from the full current-area determinant, not a local scalar ansatz.
        """
        if gauge_edge is None:
            gauge_edge = self.graph.ne - 1
        gauge_edge = _integer(gauge_edge, "gauge edge")
        if gauge_edge >= self.graph.ne:
            raise ValueError("gauge edge lies outside the original graph")
        dimension = 2 * self.graph.loops
        if dimension != self.graph.ne - 1 or len(self.cotree) != self.graph.loops:
            raise ArithmeticError("the original graph lost its H=2 density degree")
        self._budget(dimension * dimension, "the complete current-area matrix")
        currents = self.currents(sources)
        retained = tuple(edge for edge in range(self.graph.ne) if edge != gauge_edge)
        matrix = []
        for edge in self.cotree:
            factors = tuple(self._resistance_factor(edge, f) for f in retained)
            for a in range(2):
                matrix.append([factor * currents[f][a] for f, factor in zip(retained, factors)])
        loops = self.graph.loops
        sign = (-1) ** (gauge_edge + loops * (loops - 1) // 2)
        return sign * det(matrix) / self.resistances[gauge_edge]


class ForestKronPlan:
    """Reusable complete topology plan; every retained source is permitted.

    ``roots`` may specify one original interaction vertex per component.
    Otherwise a higher-valence vertex is preferred, then a degree-two one.
    The original marked labels, edge order, and orientations are preserved.
    ``max_entries`` bounds the main exact coefficient buffers and requested
    derivative outputs. It does not bound rational numerator bit lengths.
    """
    def __init__(self, graph, marked, *, roots=None, max_entries=1000000):
        self.structure = recognize_product_forest(graph, marked)
        self.graph = self.G = Graph(graph.nv, graph.edges, graph.name)
        self.marked = self.structure.marked
        self.max_entries = _integer(max_entries, "max_entries", 1)
        components = self.structure.components
        if roots is None:
            roots = tuple((higher or degree_two or component)[0]
                          for component, higher, degree_two in zip(components,
                              self.structure.higher_valence_vertices, self.structure.degree_two_vertices))
        else:
            roots = tuple(_integer(root, "component root") for root in roots)
            if len(roots) != len(components) or any(root not in component for root, component in zip(roots, components)):
                raise ValueError("roots must give one original vertex in each interaction component")
        self.roots = roots
        vertices, edges = self.graph.nv, self.graph.ne
        self.response_entry_bound = ((3 * vertices + edges) * (vertices - 1)
                                     + 12 * edges + 8 * vertices)
        if self.response_entry_bound > self.max_entries:
            raise KronBudgetExceeded("the complete source-linear response requires at most %d coefficient entries, above max_entries=%d" %
                                     (self.response_entry_bound, self.max_entries))
        adjacent = [set() for unused in range(vertices)]
        for a, b in self.graph.edges:
            adjacent[a].add(b)
            adjacent[b].add(a)
        original_internal = tuple(tuple(sorted(neighbors - set(self.marked))) for neighbors in adjacent)
        order = []
        for root in roots:
            parent, visit = {root: None}, [root]
            for vertex in visit:
                for neighbor in original_internal[vertex]:
                    if neighbor != parent[vertex]:
                        if neighbor in parent:
                            raise ArithmeticError("the validated interaction component gained a cycle")
                        parent[neighbor] = vertex
                        visit.append(neighbor)
            order.extend(reversed(visit))
        steps = []
        current_edges, current_vertices = edges, vertices
        for vertex in order:
            neighbors = tuple(sorted(adjacent[vertex]))
            if not 2 <= len(neighbors) <= 3:
                raise ArithmeticError("a product-forest leaf failed the two/three-neighbour elimination invariant")
            if len(set(neighbors) - set(self.marked)) > 1:
                raise ArithmeticError("a removed subtree has more than its one retained parent")
            pairs = tuple(combinations(neighbors, 2))
            new = tuple(pair for pair in pairs if pair[1] not in adjacent[pair[0]])
            parallel = tuple(pair for pair in pairs if pair not in new)
            for a, b in pairs:
                adjacent[a].add(b)
                adjacent[b].add(a)
            for neighbor in neighbors:
                adjacent[neighbor].remove(vertex)
            adjacent[vertex].clear()
            current_edges += len(new) - len(neighbors)
            current_vertices -= 1
            steps.append(KronStep(vertex, neighbors, new, parallel, current_vertices, current_edges))
        A, B = self.marked
        if (len(order) != vertices - 2 or current_vertices != 2 or current_edges != 1
                or adjacent[A] != {B} or adjacent[B] != {A}):
            raise ArithmeticError("the complete forest was not reduced to its marked edge")
        self.steps = tuple(steps)
        # This selects one basis, never enumerates spanning trees.
        self.cotree = tuple(self.graph.cotree())
        self.stats = {"original_vertices": vertices, "original_edges": edges,
                      "original_loops": self.graph.loops, "eliminated_vertices": len(steps),
                      "two_neighbor_steps": sum(len(step.neighbors) == 2 for step in steps),
                      "three_neighbor_steps": sum(len(step.neighbors) == 3 for step in steps),
                      "maximum_retained_neighbors": max((len(step.neighbors) for step in steps), default=0),
                      "non_Laman_intermediate_stages": sum(step.remaining_edges != 2 * step.remaining_vertices - 3 for step in steps),
                      "response_entry_bound": self.response_entry_bound,
                      "expanded_Kirchhoff_terms": 0, "spanning_trees_enumerated": 0}

    def response(self, resistances):
        resistances = tuple(_exact(value, "resistance") for value in resistances)
        if len(resistances) != self.graph.ne or any(value <= 0 for value in resistances):
            raise ValueError("one strictly positive exact resistance is required per original edge")
        conductances = {_edge(a, b): 1 / value for (a, b), value in zip(self.graph.edges, resistances)}
        size = self.graph.nv - 1
        sources = {vertex: tuple(Fraction(int(vertex == source)) for source in range(size))
                   for vertex in range(size)}
        sources[self.graph.root] = (Fraction(-1),) * size
        transfers = []
        pivot_product = Fraction(1)
        for step in self.steps:
            values = tuple(conductances[_edge(step.vertex, neighbor)] for neighbor in step.neighbors)
            pivot = sum(values, Fraction(0))
            if pivot <= 0:
                raise ArithmeticError("a positive network acquired a nonpositive pivot")
            weights = tuple(value / pivot for value in values)
            if sum(weights, Fraction(0)) != 1:
                raise ArithmeticError("a complete source transfer failed exact conservation")
            source = sources.pop(step.vertex)
            for neighbor, weight in zip(step.neighbors, weights):
                sources[neighbor] = tuple(old + weight * value for old, value in zip(sources[neighbor], source))
            for first, second in combinations(range(len(step.neighbors)), 2):
                edge = _edge(step.neighbors[first], step.neighbors[second])
                conductances[edge] = conductances.get(edge, Fraction(0)) + values[first] * values[second] / pivot
            for neighbor in step.neighbors:
                del conductances[_edge(step.vertex, neighbor)]
            transfers.append(KronTransfer(step.vertex, step.neighbors, pivot, weights, source))
            pivot_product *= pivot
        A, B = self.marked
        final = conductances[_edge(A, B)]
        if len(conductances) != 1 or set(sources) != {A, B} or any(a + b for a, b in zip(sources[A], sources[B])):
            raise ArithmeticError("the full source-transfer chain did not close on A,B")
        retained_sources = (sources[A], sources[B])
        potentials = {B: (Fraction(0),) * size, A: tuple(value / final for value in sources[A])}
        for transfer in reversed(transfers):
            potentials[transfer.vertex] = tuple(
                transfer.source_coefficients[source] / transfer.pivot
                + sum((weight * potentials[neighbor][source] for neighbor, weight in zip(transfer.neighbors, transfer.weights)), Fraction(0))
                for source in range(size))
        ground = potentials[self.graph.root]
        potential_response = tuple(tuple(value - reference for value, reference in zip(potentials[vertex], ground))
                                   for vertex in range(self.graph.nv))
        current_response = tuple(tuple((left - right) / value for left, right in zip(potentials[a], potentials[b]))
                                 for (a, b), value in zip(self.graph.edges, resistances))
        determinant = pivot_product * final
        U = determinant * prod(resistances)
        return ForestKronResponse(self.graph, self.marked, resistances, current_response, potential_response,
                                  tuple(transfers), retained_sources, final, determinant, U,
                                  self.cotree, self.max_entries, dict(self.stats, all_original_sources_retained=True))
