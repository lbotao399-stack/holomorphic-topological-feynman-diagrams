"""Bounded exact reduction of a whole linear combination of scalar periods.

Input rows are ``(edges, powers, coefficient)`` at a common integer q.
Every row denotes coefficient * I_q(edges; powers), with I_q as in engine.py.
The output always retains the whole expression: a finite zeta polynomial
plus explicit unresolved completed periods in the engine's Pi normalization.
An exhausted search budget is never an assertion that a residual is zero.

The weighted union-find records Pi(a) = factor * Pi(b) for the existing
uniqueness, inverse uniqueness, Fourier-duality and twist identities.  It
combines coefficients across graphs before requesting further identities.
This layer does not add a new integration formula or construct period tables.
"""
from dataclasses import dataclass
from fractions import Fraction
import heapq
from math import prod
from numbers import Integral
from types import MappingProxyType

import engine


class CombinationBudgetExceeded(RuntimeError):
    """The input itself exceeds a declared storage or graph-size budget."""


class IncompleteCombination(LookupError):
    """The exact expression still contains explicitly retained graph periods."""


@dataclass(frozen=True)
class CompletedPeriodTerm:
    """coefficient * Pi(vertices, weights), with weighted edges (i,j,beta)."""

    q: int
    vertices: int
    weights: tuple
    coefficient: Fraction


@dataclass(frozen=True)
class CombinationReduction:
    """An exact sum; finite_coefficients alone is the value only if complete."""

    finite_coefficients: object
    residuals: tuple
    statistics: object

    @property
    def complete(self):
        return not self.residuals

    @property
    def value(self):
        if self.residuals:
            raise IncompleteCombination('%d completed periods remain' % len(self.residuals))
        return self.finite_coefficients


def _integer(value, name, minimum):
    if isinstance(value, bool) or not isinstance(value, Integral) or value < minimum:
        raise ValueError('%s must be an integer >= %d' % (name, minimum))
    return int(value)


def _exact(value):
    if isinstance(value, (float, complex, bool)):
        raise ValueError('coefficients must be exact rational numbers')
    return Fraction(value)


class _WeightedOrbits:
    """weight[k] satisfies Pi(k)=weight[k]*Pi(parent[k])."""

    def __init__(self):
        self.parent = {}
        self.weight = {}
        self.size = {}
        self.known = {}
        self.coefficient = {}

    def add(self, key):
        if key not in self.parent:
            self.parent[key] = key
            self.weight[key] = Fraction(1)
            self.size[key] = 1
            self.coefficient[key] = Fraction(0)

    def find(self, key):
        parent = self.parent[key]
        if parent != key:
            root, weight = self.find(parent)
            self.weight[key] *= weight
            self.parent[key] = root
        return self.parent[key], self.weight[key]

    def set_known(self, key, value):
        root, weight = self.find(key)
        value = {label: coefficient / weight for label, coefficient in value.items() if coefficient}
        if root in self.known and self.known[root] != value:
            raise ArithmeticError('two finite routes disagree in their exact coefficient dictionaries')
        self.known[root] = value

    def join(self, first, second, factor):
        self.add(first)
        self.add(second)
        a, wa = self.find(first)
        b, wb = self.find(second)
        if a == b:
            if wa != factor * wb:
                raise ArithmeticError('inconsistent rational factor around a conformal identity cycle')
            return False
        factor = wa / (factor * wb)
        if self.size[a] < self.size[b]:
            a, b = b, a
            factor = 1 / factor
        self.parent[b] = a
        self.weight[b] = factor
        self.size[a] += self.size.pop(b)
        self.coefficient[a] += factor * self.coefficient.pop(b)
        value = self.known.pop(b, None)
        if value is not None:
            self.set_known(b, value)
        return True


def reduce_period_combination(terms, q, *, max_states=10000, max_nodes=50000,
                              max_input_terms=250000, max_edges=18,
                              max_vertices=None, use_inverse_uniqueness=True):
    """Preserve and reduce every exact coefficient within explicit budgets.

    max_states bounds expanded canonical graphs; max_nodes bounds all stored
    canonical graphs.  When either search limit is reached, the remaining
    completed periods stay in the result.  Inputs exceeding max_input_terms,
    max_edges, max_vertices or the initial max_nodes limit raise explicitly.
    No product factorization or numerator terminal is needed for the linear
    equivalence reduction here; engine.solve offers those additional routes.
    """
    q = _integer(q, 'q', 2)
    max_states = _integer(max_states, 'max_states', 0)
    max_nodes = _integer(max_nodes, 'max_nodes', 1)
    max_input_terms = _integer(max_input_terms, 'max_input_terms', 1)
    max_edges = _integer(max_edges, 'max_edges', 1)
    if max_vertices is not None:
        max_vertices = _integer(max_vertices, 'max_vertices', 3)
    if not isinstance(use_inverse_uniqueness, bool):
        raise ValueError('use_inverse_uniqueness must be boolean')
    coefficients = {}
    input_terms = 0
    input_vertex_maximum = 0
    for input_terms, (edges, powers, coefficient) in enumerate(terms, 1):
        if input_terms > max_input_terms:
            raise CombinationBudgetExceeded('max_input_terms exceeded')
        coefficient = _exact(coefficient)
        if not coefficient:
            continue
        edges = tuple(tuple(edge) for edge in edges)
        if not edges or len(edges) > max_edges:
            raise CombinationBudgetExceeded('input edge count is outside the declared max_edges budget')
        if any(len(edge) != 2 or edge[0] == edge[1] for edge in edges):
            raise ValueError('each edge must join two distinct vertices')
        powers = tuple(_integer(power, 'edge power', 1) for power in powers)
        if len(powers) != len(edges):
            raise ValueError('one positive integer power is required per edge')
        if any(power >= q for power in powers):
            raise ValueError('a convergent input must have every edge power < q')
        if not engine.nx.is_connected(engine.nx.Graph(edges)):
            raise ValueError('the input graph must be connected')
        if not engine.gutil.convergent(edges, powers, q):
            raise ValueError('every supplied scalar monomial must converge strictly')
        vertices, weighted, _ = engine.complete(edges, powers, q)
        input_vertex_maximum = max(input_vertex_maximum, vertices)
        if max_vertices is not None and vertices > max_vertices:
            raise CombinationBudgetExceeded('an input completion exceeds max_vertices')
        key = engine.canon(vertices, weighted)
        coefficient *= prod(engine.gam(q - power) for power in powers) / engine.gam(q)
        coefficients[key] = coefficients.get(key, Fraction(0)) + coefficient
        if len(coefficients) > max_nodes:
            raise CombinationBudgetExceeded('the input completions exceed max_nodes')
    coefficients = {key: value for key, value in coefficients.items() if value}
    if max_vertices is None:
        max_vertices = input_vertex_maximum + int(use_inverse_uniqueness)
    orbits = _WeightedOrbits()
    queue = []
    queued = set()

    def enqueue(key):
        if key not in queued:
            queued.add(key)
            heapq.heappush(queue, (key[0], key))

    for key in sorted(coefficients):
        orbits.add(key)
        orbits.coefficient[key] = coefficients[key]
        enqueue(key)
    expanded = merges = cycles = terminal_values = skipped_nodes = skipped_vertices = cancelled_skips = 0
    while queue and expanded < max_states:
        _, key = heapq.heappop(queue)
        root, _ = orbits.find(key)
        if root in orbits.known:
            continue
        if not orbits.coefficient[root]:
            cancelled_skips += 1
            continue
        expanded += 1
        vertices, weighted = engine.unkey(key)
        terminal = engine.terminal(vertices, weighted, q)
        if terminal is not None:
            orbits.set_known(key, terminal[0])
            terminal_values += 1
            continue
        for count, graph, factor, _ in engine.successors(
                vertices, weighted, q, use_m2=use_inverse_uniqueness):
            if count > max_vertices:
                skipped_vertices += 1
                continue
            other = engine.canon(count, graph)
            if other not in orbits.parent and len(orbits.parent) >= max_nodes:
                skipped_nodes += 1
                continue
            if orbits.join(key, other, factor):
                merges += 1
            else:
                cycles += 1
            if orbits.find(other)[0] not in orbits.known:
                enqueue(other)
    finite = {}
    residuals = {}
    finite_input_orbits = 0
    for key, coefficient in coefficients.items():
        root, weight = orbits.find(key)
        coefficient *= weight
        if root in orbits.known:
            finite_input_orbits += 1
            for label, value in orbits.known[root].items():
                finite[label] = finite.get(label, Fraction(0)) + coefficient * value
        else:
            residuals[root] = residuals.get(root, Fraction(0)) + coefficient
    cancelled = sum(not coefficient for coefficient in residuals.values())
    residual_terms = tuple(CompletedPeriodTerm(q, root[0], root[1], coefficient)
                           for root, coefficient in sorted(residuals.items()) if coefficient)
    labels = {('*'.join(label) or '1'): value for label, value in sorted(finite.items()) if value}
    statistics = dict(input_terms=input_terms, input_completed_orbits=len(coefficients),
                      expanded_states=expanded, stored_states=len(orbits.parent),
                      pending_states=len(queue), merged_states=merges,
                      checked_identity_cycles=cycles, terminal_values=terminal_values,
                      finite_input_orbits=finite_input_orbits,
                      cancelled_unresolved_components=cancelled,
                      cancelled_state_skips=cancelled_skips,
                      residual_components=len(residual_terms),
                      node_budget_skips=skipped_nodes, vertex_budget_skips=skipped_vertices,
                      state_budget_reached=bool(queue and expanded == max_states),
                      read_only_tables=True)
    return CombinationReduction(MappingProxyType(labels), residual_terms,
                                MappingProxyType(statistics))
