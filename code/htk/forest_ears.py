"""Exact sourced triangle-ear reductions in arbitrary tight H=2 cores.

An ear is a degree-two vertex v whose neighbours u,w are already joined.
Its source is retained by the affine transfer
    lambda_u -> lambda_u + (1-r)*lambda_v,
    lambda_w -> lambda_w + r*lambda_v.
The removed currents are also retained.  A bundle of k ears on the same
edge has a simplex of current-splitting fractions and only a two-component
total source transfer.  Exact Dirichlet moments remove the simplex; finite
source moments and total-degree Newton interpolation remove the remaining
source variables.  No polynomial in 2*k local ratios is expanded.

The residual core need not be series-parallel and its sources are arbitrary.
This module does not assert rationality of the residual periods.
"""
from collections.abc import Mapping
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import comb, factorial, prod
from numbers import Integral
import random

from . import brute
from .graph import Graph
from .linalg import perm_sign
from .forest_contractions import normalize_forest_sources
from .forest_periods import (CoreCurrentPeriods, ExactForestSeries, MissingForestPeriod,
                             ForestPeriodBudgetExceeded, _PeriodWorker, _finite)

try:
    from flint import fmpq, fmpq_poly, fmpq_mpoly_ctx
except ImportError as error:
    _FLINT_ERROR = error
else:
    _FLINT_ERROR = None


__all__ = ["EarBundle", "EarPushforward", "EarSample", "ForestEarMoments",
           "EarBudgetExceeded", "find_ear_bundles"]


class EarBudgetExceeded(RuntimeError):
    """An explicit exact source-sampling or polynomial budget was exceeded."""


def _integer(value, name, minimum=0):
    if isinstance(value, bool) or not isinstance(value, Integral) or value < minimum:
        raise ValueError("%s must be an integer at least %d" % (name, minimum))
    return int(value)


def _exact(value):
    if isinstance(value, float):
        raise TypeError("exact rational coefficients are required")
    return Fraction(value)


def _q(value):
    value = Fraction(value)
    return fmpq(value.numerator, value.denominator)


def _fraction(value):
    return Fraction(int(value.p), int(value.q))


def _bracket(left, right):
    return left[0] * right[1] - left[1] * right[0]


def _polynomial(graph, powers=None, polynomial=None):
    size = 2 * graph.ne
    if polynomial is not None:
        if powers is not None:
            raise ValueError("provide a current word or a polynomial, not both")
        if not isinstance(polynomial, Mapping):
            raise TypeError("current polynomials map full exponent tuples to rational coefficients")
        result = {}
        for exponent, coefficient in polynomial.items():
            exponent = tuple(_integer(power, "current exponent") for power in exponent)
            if len(exponent) != size:
                raise ValueError("one exponent is required for each original current component")
            coefficient = _exact(coefficient)
            if coefficient:
                result[exponent] = result.get(exponent, Fraction(0)) + coefficient
        return {exponent: coefficient for exponent, coefficient in result.items() if coefficient}
    if powers is not None and not isinstance(powers, Mapping):
        raise TypeError("current powers map (original edge, component) to integers")
    exponent = [0] * size
    for location, power in ({} if powers is None else powers).items():
        if not isinstance(location, tuple) or len(location) != 2:
            raise ValueError("a current location is (edge, component)")
        edge, component = (_integer(location[0], "edge"), _integer(location[1], "component"))
        if edge >= graph.ne or component >= 2:
            raise ValueError("current location lies outside the graph")
        exponent[2 * edge + component] = _integer(power, "current exponent")
    return {tuple(exponent): Fraction(1)}


def find_ear_bundles(graph, marked=()):
    """All sourced or unsourced ear bundles, largest first, without changing sources."""
    protected = set(marked)
    adjacent = [set() for _ in range(graph.nv)]
    for a, b in graph.edges:
        adjacent[a].add(b)
        adjacent[b].add(a)
    groups = {}
    for vertex in range(graph.nv):
        if vertex in protected or graph.degree(vertex) != 2:
            continue
        neighbours = tuple(sorted(adjacent[vertex]))
        if len(neighbours) == 2 and neighbours[1] in adjacent[neighbours[0]]:
            groups.setdefault(neighbours, []).append(vertex)
    return tuple((edge, tuple(vertices)) for edge, vertices in sorted(
        groups.items(), key=lambda item: (-len(item[1]), item[0])))


def _source_basis(vectors):
    nonzero = [vector for vector in vectors if any(vector)]
    if not nonzero:
        return (), tuple(() for _ in vectors)
    first = nonzero[0]
    second = next((vector for vector in nonzero if _bracket(first, vector)), None)
    if second is None:
        component = 0 if first[0] else 1
        return (first,), tuple((vector[component] / first[component],) for vector in vectors)
    denominator = _bracket(first, second)
    coordinates = tuple((_bracket(vector, second) / denominator,
                         _bracket(first, vector) / denominator) for vector in vectors)
    return (first, second), coordinates


@lru_cache(maxsize=None)
def _binomial_polynomials(degree):
    value, result = fmpq_poly([1]), []
    for order in range(degree + 1):
        result.append(tuple(value[index] for index in range(order + 1)))
        value = value * fmpq_poly([-order, 1]) / (order + 1)
    return tuple(result)


@lru_cache(maxsize=128)
def _source_weights(coordinates, beta_powers, degree):
    """Exact weighted box-source functional on a one/two-variable Newton grid.

    The exponential moment series of one ear has coefficient
      coordinate**exponent / exponent! *
      (|exponent|+alpha)! beta! / (|exponent|+alpha+beta+1)!.
    Its finite products contain every moment needed through ``degree``.
    Transposing the exact Newton finite-difference formula gives the sample
    weights. They may have either sign; this is polynomial integration.
    """
    dimension = len(coordinates[0])
    context = fmpq_mpoly_ctx.get(tuple("source%d" % i for i in range(dimension)))
    value = context.constant(1)
    indices = (tuple((i,) for i in range(degree + 1)) if dimension == 1
               else tuple((i, j) for i in range(degree + 1) for j in range(degree + 1 - i)))
    for vector, (alpha, beta) in zip(coordinates, beta_powers):
        terms = {}
        for exponent in indices:
            total = sum(exponent)
            coefficient = fmpq(factorial(total + alpha) * factorial(beta),
                               factorial(total + alpha + beta + 1))
            for coordinate, power in zip(vector, exponent):
                coefficient *= _q(coordinate) ** power / factorial(power)
            if coefficient:
                terms[exponent] = coefficient
        expanded = value * context.from_dict(terms)
        value = context.from_dict({exponent: coefficient for exponent, coefficient in expanded.terms()
                                   if sum(exponent) <= degree})
    moments = {}
    for exponent, coefficient in value.terms():
        for power in exponent:
            coefficient *= factorial(int(power))
        moments[tuple(map(int, exponent))] = coefficient
    falling = _binomial_polynomials(degree)
    weights = {index: fmpq(0) for index in indices}
    for index in indices:
        moment = fmpq(0)
        for exponent in product(*(range(order + 1) for order in index)):
            coefficient = moments.get(exponent, fmpq(0))
            for order, power in zip(index, exponent):
                coefficient *= falling[order][power]
            moment += coefficient
        if not moment:
            continue
        for sample in product(*(range(order + 1) for order in index)):
            coefficient = moment * (-1) ** (sum(index) - sum(sample))
            for order, at in zip(index, sample):
                coefficient *= comb(order, at)
            weights[sample] += coefficient
    return tuple((sample, _fraction(weight)) for sample, weight in weights.items() if weight)


@dataclass(frozen=True)
class EarSample:
    point: tuple
    sources: tuple
    polynomial: dict


@dataclass(frozen=True)
class EarPushforward:
    graph: Graph
    marked: tuple
    samples: tuple
    stats: dict


class EarBundle:
    """Remove any chosen bundle of degree-two vertices on one existing edge."""
    def __init__(self, graph, edge, vertices=None, marked=(), *, max_terms=1000000,
                 max_source_samples=10000):
        if _FLINT_ERROR is not None:
            raise ImportError("sourced ear reductions require python-flint") from _FLINT_ERROR
        self.G = self.graph = Graph(graph.nv, graph.edges, graph.name)
        if not self.graph.is_connected() or not self.graph.tight(2):
            raise ValueError("a connected tight H=2 graph is required")
        if len(set(map(frozenset, self.graph.edges))) != self.graph.ne:
            raise ValueError("the sourced-ear input must have no parallel edges")
        self.max_terms = _integer(max_terms, "max_terms", 1)
        self.max_source_samples = _integer(max_source_samples, "max_source_samples", 1)
        self.marked = tuple(_integer(v, "marked vertex") for v in marked)
        if len(set(self.marked)) != len(self.marked) or any(v >= graph.nv for v in self.marked):
            raise ValueError("protected vertices must be distinct graph vertices")
        self.left, self.right = tuple(_integer(v, "ear endpoint") for v in edge)
        if self.left == self.right or max(self.left, self.right) >= graph.nv:
            raise ValueError("ear endpoints must be distinct graph vertices")
        self.direct = next((i for i, ends in enumerate(graph.edges)
                            if set(ends) == {self.left, self.right}), None)
        if self.direct is None:
            raise ValueError("the retained ear edge is absent")
        available = tuple(v for ends, group in find_ear_bundles(graph, self.marked)
                          if set(ends) == {self.left, self.right} for v in group)
        self.vertices = available if vertices is None else tuple(_integer(v, "ear vertex") for v in vertices)
        if not self.vertices or len(set(self.vertices)) != len(self.vertices) or not set(self.vertices) <= set(available):
            raise ValueError("every removed vertex must be an unmarked degree-two ear on the retained edge")
        self.k = len(self.vertices)
        self.spokes = tuple(tuple(next(e for e, ends in enumerate(graph.edges) if set(ends) == set(pair))
                                   for pair in ((self.left, v), (v, self.right))) for v in self.vertices)
        local = (self.direct,) + tuple(e for pair in self.spokes for e in pair)
        self.others = tuple(e for e in range(graph.ne) if e not in local)
        self.kernel_sign = ((-1) ** (self.k * graph.loops - self.k * (self.k - 1) // 2)
                            * perm_sign(local + self.others))
        retained = tuple(v for v in range(graph.nv) if v not in self.vertices)
        mapping = {v: i for i, v in enumerate(retained)}
        self.vmap = tuple(mapping.get(v) for v in range(graph.nv))
        self.red = self.reduced_graph = Graph(len(retained),
            [(mapping[self.left], mapping[self.right])]
            + [(mapping[graph.edges[e][0]], mapping[graph.edges[e][1]]) for e in self.others])
        self.reduced_marked = tuple(mapping[v] for v in self.marked)
        self._locations = {self.direct: ("direct", 0, 1 if graph.edges[self.direct] == (self.left, self.right) else -1)}
        for i, vertex in enumerate(self.vertices):
            for side, (edge_index, canonical) in enumerate(zip(self.spokes[i],
                                          ((self.left, vertex), (vertex, self.right)))):
                self._locations[edge_index] = ("ear", (i, side), 1 if graph.edges[edge_index] == canonical else -1)
        self._locations.update({edge_index: ("outside", i + 1, 1) for i, edge_index in enumerate(self.others)})

    def shifted_sources(self, sources, ratios):
        sources = normalize_forest_sources(self.graph, sources)
        if len(ratios) != self.k:
            raise ValueError("one source ratio is required per ear")
        rows = [list(sources[v]) for v in range(self.graph.nv) if self.vmap[v] is not None]
        for vertex, ratio in zip(self.vertices, ratios):
            ratio = _exact(ratio)
            for a in range(2):
                rows[self.vmap[self.left]][a] += (1 - ratio) * sources[vertex][a]
                rows[self.vmap[self.right]][a] += ratio * sources[vertex][a]
        return tuple(map(tuple, rows))

    def time_map(self, ratios, weights, times):
        ratios = tuple(_exact(value) for value in ratios)
        weights = tuple(_exact(value) for value in weights)
        times = tuple(_exact(value) for value in times)
        if len(ratios) != self.k or len(weights) != self.k or len(times) != self.red.ne:
            raise ValueError("ear coordinates have the wrong length")
        direct_weight = 1 - sum(weights)
        if not 0 < direct_weight < 1 or any(not 0 < weight < 1 for weight in weights):
            raise ValueError("the splitting weights must lie inside their simplex")
        if any(not 0 < ratio < 1 for ratio in ratios) or any(time <= 0 for time in times):
            raise ValueError("source ratios and all times must be positive interior coordinates")
        result = [None] * self.graph.ne
        result[self.direct] = times[0] / direct_weight
        for (left, right), ratio, weight in zip(self.spokes, ratios, weights):
            result[left] = ratio * times[0] / weight
            result[right] = (1 - ratio) * times[0] / weight
        for i, edge_index in enumerate(self.others):
            result[edge_index] = times[i + 1]
        return tuple(result)

    def current_map(self, sources, ratios, weights, reduced_currents):
        sources = normalize_forest_sources(self.graph, sources)
        ratios = tuple(_exact(value) for value in ratios)
        weights = tuple(_exact(value) for value in weights)
        if len(ratios) != self.k or len(weights) != self.k or len(reduced_currents) != self.red.ne:
            raise ValueError("ear current coordinates have the wrong length")
        result, J = [], reduced_currents[0]
        for edge in range(self.graph.ne):
            kind, location, sign = self._locations[edge]
            if kind == "outside":
                vector = reduced_currents[location]
            elif kind == "direct":
                vector = tuple((1 - sum(weights)) * value for value in J)
            else:
                i, side = location
                source = sources[self.vertices[i]]
                vector = tuple(weights[i] * J[a] + (ratios[i] if side else ratios[i] - 1) * source[a]
                               for a in range(2))
            result.append(tuple(sign * value for value in vector))
        return tuple(result)

    def check(self, sources, count=1, seed=748192):
        """Exact original-current and density checks, including all removed sources."""
        sources = normalize_forest_sources(self.graph, sources)
        rng = random.Random(seed)
        for _ in range(_integer(count, "check count")):
            ratios = tuple(Fraction(rng.randint(1, 8), 9) for _ in range(self.k))
            raw = [rng.randint(1, 9) for _ in range(self.k + 1)]
            weights = tuple(Fraction(value, sum(raw)) for value in raw[1:])
            times = tuple(Fraction(rng.randint(1, 9), rng.randint(1, 5)) for _ in range(self.red.ne))
            original_times = self.time_map(ratios, weights, times)
            rows = self.shifted_sources(sources, ratios)
            currents = brute.currents(self.red, times, rows)
            if tuple(map(tuple, brute.currents(self.graph, original_times, sources))) != self.current_map(sources, ratios, weights, currents):
                raise ArithmeticError("a sourced-ear current map changed an original current")
            jacobian = times[0] ** (2 * self.k) / (1 - sum(weights)) ** 2
            factor = self.kernel_sign * (1 - sum(weights))
            for vertex, weight in zip(self.vertices, weights):
                jacobian /= weight ** 3
                factor *= _bracket(sources[vertex], currents[0])
            actual = brute.density(self.graph, original_times, sources, T=0).rational() * jacobian
            expected = factor * brute.density(self.red, times, rows, T=0).rational()
            if actual != expected:
                raise ArithmeticError("a sourced-ear density/sign identity failed")

    def pushforward(self, sources, powers=None, *, polynomial=None):
        """Exact finite lower-loop source samples for every original current word."""
        sources = normalize_forest_sources(self.graph, sources)
        original = _polynomial(self.graph, powers, polynomial)
        if len(original) > self.max_terms:
            raise EarBudgetExceeded("the input current polynomial exceeds max_terms")
        degree = max((sum(exponent) for exponent in original), default=0)
        # All transferred source variation is between two retained
        # endpoints. Every core current changes by h_e*mu, so the varying
        # part of the projective current determinant has rank at most the
        # number of retained loops. Its mu-degree is therefore <= loops,
        # although its unrestricted total source degree is twice that.
        source_degree = self.red.loops + self.k + degree
        ear_sources = tuple(sources[v] for v in self.vertices)
        basis, coordinates = _source_basis(ear_sources)
        stats = {"original_loops": self.graph.loops, "reduced_loops": self.red.loops,
                 "removed_vertices": self.vertices, "input_current_degree": degree,
                 "source_degree_bound": source_degree, "source_variables": len(basis),
                 "global_ratio_polynomial_variables": 0, "kernel_sign": self.kernel_sign,
                 "largest_polynomial_terms": len(original), "source_samples": 0}
        if not original or any(not any(vector) for vector in ear_sources):
            stats["identically_zero"] = True
            return EarPushforward(self.red, self.reduced_marked, (), stats)
        sample_bound = source_degree + 1 if len(basis) == 1 else (source_degree + 1) * (source_degree + 2) // 2
        if sample_bound > self.max_source_samples:
            raise EarBudgetExceeded("the exact source grid needs %d samples, exceeding max_source_samples=%d" %
                                    (sample_bound, self.max_source_samples))
        context = fmpq_mpoly_ctx.get(tuple("j%d" % i for i in range(2 * self.red.ne)))
        variables = context.gens()
        J = variables[:2]
        determinant_factor = context.constant(1)
        for source in ear_sources:
            determinant_factor *= _q(source[0]) * J[1] - _q(source[1]) * J[0]
        groups = {}
        for exponent, original_coefficient in original.items():
            # State = splitting exponents, beta exponents for source ratios,
            #         exponents of all retained current components.
            initial = ((0,) * (self.k + 1), (0,) * (2 * self.k), (0,) * (2 * self.red.ne))
            states = {initial: original_coefficient}
            for location, power in enumerate(exponent):
                if not power:
                    continue
                edge, component = divmod(location, 2)
                kind, where, orientation = self._locations[edge]
                new = {}
                for (splitting, beta, currents), coefficient in states.items():
                    choices = range(power + 1) if kind == "ear" else (power,)
                    for kept in choices:
                        split, source_powers, current_powers = list(splitting), list(beta), list(currents)
                        scalar = coefficient * orientation ** power
                        if kind == "outside":
                            current_powers[2 * where + component] += power
                        else:
                            current_powers[component] += kept
                            if kind == "direct":
                                split[0] += kept
                            else:
                                i, side = where
                                split[i + 1] += kept
                                removed = power - kept
                                scalar *= comb(power, kept) * ear_sources[i][component] ** removed
                                if not side:
                                    scalar *= (-1) ** removed
                                source_powers[2 * i + (0 if side else 1)] += removed
                        if scalar:
                            key = tuple(split), tuple(source_powers), tuple(current_powers)
                            new[key] = new.get(key, Fraction(0)) + scalar
                states = {key: value for key, value in new.items() if value}
                stats["largest_polynomial_terms"] = max(stats["largest_polynomial_terms"], len(states))
                if len(states) > self.max_terms:
                    raise EarBudgetExceeded("ear word expansion exceeds max_terms")
            for (splitting, beta, currents), coefficient in states.items():
                weight = Fraction(factorial(splitting[0] + 1), factorial(self.k + 1 + sum(splitting)))
                for power in splitting[1:]:
                    weight *= factorial(power)
                alpha_beta = tuple((beta[2 * i], beta[2 * i + 1]) for i in range(self.k))
                terms = groups.setdefault(alpha_beta, {})
                terms[currents] = terms.get(currents, Fraction(0)) + self.kernel_sign * coefficient * weight
            if sum(len(terms) for terms in groups.values()) > self.max_terms:
                raise EarBudgetExceeded("the complete ear word expansion exceeds max_terms")
        sample_polynomials = {}
        for alpha_beta, terms in groups.items():
            polynomial = context.from_dict({exponent: _q(coefficient) for exponent, coefficient in terms.items() if coefficient})
            for point, weight in _source_weights(coordinates, alpha_beta, source_degree):
                sample_polynomials[point] = sample_polynomials.get(point, context.constant(0)) + _q(weight) * polynomial
            if sum(len(value) for value in sample_polynomials.values()) > self.max_terms:
                raise EarBudgetExceeded("the complete source-sampled current expression exceeds max_terms")
        retained = [list(sources[v]) for v in range(self.graph.nv) if self.vmap[v] is not None]
        for source in ear_sources:
            for a in range(2):
                retained[self.vmap[self.left]][a] += source[a]
        samples, sample_term_count = [], 0
        for point, polynomial in sample_polynomials.items():
            polynomial *= determinant_factor
            stats["largest_polynomial_terms"] = max(stats["largest_polynomial_terms"], len(polynomial))
            if len(polynomial) > self.max_terms:
                raise EarBudgetExceeded("a sampled current polynomial exceeds max_terms")
            if not polynomial:
                continue
            sample_term_count += len(polynomial)
            if sample_term_count > self.max_terms:
                raise EarBudgetExceeded("the complete source-sampled current expression exceeds max_terms")
            rows = [list(row) for row in retained]
            displacement = tuple(sum(Fraction(point[j]) * basis[j][a] for j in range(len(basis))) for a in range(2))
            for a in range(2):
                rows[self.vmap[self.left]][a] -= displacement[a]
                rows[self.vmap[self.right]][a] += displacement[a]
            values = {tuple(map(int, exponent)): _fraction(coefficient)
                      for exponent, coefficient in polynomial.terms() if coefficient}
            samples.append(EarSample(point, tuple(map(tuple, rows)), values))
        stats["source_samples"] = len(samples)
        stats["sampled_current_terms"] = sample_term_count
        stats["source_moment_groups"] = len(groups)
        stats["source_basis"] = basis
        return EarPushforward(self.red, self.reduced_marked, tuple(samples), stats)


class ForestEarMoments:
    """Complete finite/series moments after source-preserving ear contractions.

    The geometric plan is fixed before any source interpolation: special
    source values never trigger node-specific graph simplifications. Every
    terminal node therefore uses the same reduced edge and period basis.
    Residual Schwinger/period polynomials are summed before the private
    read-only period worker is called, once per denominator power.

    ``moment`` accepts powers on the user's original edges. ``integrate``
    accepts an arbitrary finite current polynomial in edge-major order.
    ``series`` returns the complete exact convergent-series object, including
    a terminating constant when the reduced core is series-parallel.
    """
    def __init__(self, graph, sources, marked, *, check=1, max_terms=1000000,
                 max_source_samples=10000, max_tree_candidates=1000000,
                 max_period_terms=250000, max_series_edges=18,
                 period_timeout=120, max_matrix_entries=1000000):
        if _FLINT_ERROR is not None:
            raise ImportError("sourced ear reductions require python-flint") from _FLINT_ERROR
        self.graph = self.G = Graph(graph.nv, graph.edges, graph.name)
        if not self.graph.is_connected() or not self.graph.tight(2):
            raise ValueError("a connected tight H=2 graph is required")
        self.marked = tuple(_integer(v, "marked vertex") for v in marked)
        if len(self.marked) != 2 or len(set(self.marked)) != 2 or max(self.marked) >= graph.nv:
            raise ValueError("marked must contain two distinct graph vertices")
        self.sources = normalize_forest_sources(graph, sources)
        self.check = _integer(check, "check count")
        self.max_terms = _integer(max_terms, "max_terms", 1)
        self.max_source_samples = _integer(max_source_samples, "max_source_samples", 1)
        self.options = {"check": self.check, "max_terms": self.max_terms,
                        "max_tree_candidates": _integer(max_tree_candidates, "max_tree_candidates", 1),
                        "max_period_terms": _integer(max_period_terms, "max_period_terms", 1),
                        "max_series_edges": _integer(max_series_edges, "max_series_edges", 1),
                        "period_timeout": _integer(period_timeout, "period_timeout", 1),
                        "max_matrix_entries": _integer(max_matrix_entries, "max_matrix_entries", 1)}
        current, current_marked, current_sources = self.graph, self.marked, self.sources
        original_vertices = tuple(range(graph.nv))
        self._moves, records = [], []
        while True:
            bundles = find_ear_bundles(current, current_marked)
            if not bundles:
                break
            edge, vertices = bundles[0]
            move = EarBundle(current, edge, vertices, current_marked, max_terms=self.max_terms,
                             max_source_samples=self.max_source_samples)
            if self.check:
                move.check(current_sources, self.check)
            records.append({"original_vertices": tuple(original_vertices[v] for v in vertices),
                            "original_endpoints": tuple(original_vertices[v] for v in edge),
                            "before_edges": tuple(current.edges), "after_edges": tuple(move.red.edges),
                            "before_loops": current.loops, "after_loops": move.red.loops,
                            "kernel_sign": move.kernel_sign,
                            "vertex_map": move.vmap, "exact_point_checks": self.check})
            self._moves.append(move)
            current_sources = move.shifted_sources(current_sources, (Fraction(1, 2),) * move.k)
            original_vertices = tuple(original_vertices[v] for v in range(current.nv) if move.vmap[v] is not None)
            current, current_marked = move.red, move.reduced_marked
        self._moves = tuple(self._moves)
        self.steps = tuple(records)
        self.reduced_graph, self.reduced_marked = current, current_marked
        self._worker = None
        self._cache = {}
        self._period_geometry = None
        self._period_orbits = {}
        self.last_stats = {}

    def close(self):
        if self._worker is not None:
            self._worker.close()
            self._worker = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, traceback):
        self.close()

    def _aggregate(self, original):
        key = tuple(sorted(original.items()))
        if key in self._cache:
            return self._cache[key]
        schwinger, constant, template = {}, Fraction(0), None
        stats = {"original_loops": self.graph.loops, "reduced_loops": self.reduced_graph.loops,
                 "ear_bundles": len(self._moves), "removed_vertices": sum(move.k for move in self._moves),
                 "all_original_sources_retained": True, "original_current_terms": len(original),
                 "terminal_source_evaluations": 0, "source_samples_by_stage": [0] * len(self._moves),
                 "source_degree_bounds": [0] * len(self._moves),
                 "source_variables_by_stage": [0] * len(self._moves),
                 "largest_polynomial_terms": len(original), "largest_core_numerator_terms": 0,
                 "fixed_reduced_graph": tuple(self.reduced_graph.edges),
                 "global_ratio_polynomial_variables": 0,
                 "source_aggregation": "before_series_beta_and_graph_symmetry",
                 "exact_point_checks": self.check * len(self._moves)}

        def visit(stage, sources, polynomial):
            nonlocal constant, template
            if not polynomial:
                return
            if stage < len(self._moves):
                pushed = self._moves[stage].pushforward(sources, polynomial=polynomial)
                stats["source_samples_by_stage"][stage] += len(pushed.samples)
                stats["source_degree_bounds"][stage] = max(stats["source_degree_bounds"][stage],
                                                            pushed.stats["source_degree_bound"])
                stats["source_variables_by_stage"][stage] = max(stats["source_variables_by_stage"][stage],
                                                                pushed.stats["source_variables"])
                stats["largest_polynomial_terms"] = max(stats["largest_polynomial_terms"],
                                                          pushed.stats["largest_polynomial_terms"])
                for sample in pushed.samples:
                    visit(stage + 1, sample.sources, sample.polynomial)
                return
            stats["terminal_source_evaluations"] += 1
            if stats["terminal_source_evaluations"] > self.max_source_samples:
                raise EarBudgetExceeded("the complete fixed-core source evaluation exceeds max_source_samples=%d" %
                                        self.max_source_samples)
            engine = CoreCurrentPeriods(self.reduced_graph, sources, self.reduced_marked, **self.options)
            try:
                if engine._sp is not None:
                    constant += engine.integrate(polynomial, allow_series=False)
                    stats["largest_polynomial_terms"] = max(stats["largest_polynomial_terms"],
                                                              engine.last_stats.get("substituted_factor_terms", 0))
                    return
                # Every sample has the very same retained graph and time
                # variables. Sum its exact source polynomial first. Chain
                # beta factors and graph symmetries are then applied once,
                # instead of repeating identical combinatorics per sample.
                values = engine._schwinger(engine._polynomial(polynomial))
                if template is None:
                    template = engine
                stats["largest_core_numerator_terms"] = max(stats["largest_core_numerator_terms"],
                                                              sum(len(value) for value in values.values()))
                for q, numerator in values.items():
                    schwinger[q] = schwinger.get(q, engine._context.constant(0)) + numerator
                    if len(schwinger[q]) > self.max_terms:
                        raise EarBudgetExceeded("the source-integrated Schwinger numerator exceeds max_terms")
            finally:
                engine.close()

        if len(original) > self.max_terms:
            raise EarBudgetExceeded("the input current polynomial exceeds max_terms")
        visit(0, self.sources, original)
        schwinger = {q: numerator for q, numerator in schwinger.items() if numerator}
        total, period_edges = {}, ()
        if template is not None:
            if self._period_geometry is None:
                template._chain_data()
                self._period_geometry = (template._chains, tuple(template._core), template._permutations)
            chains, period_edges, permutations = self._period_geometry
            for q, numerator in schwinger.items():
                values = {}
                for exponent, coefficient in numerator.terms():
                    exponent = tuple(map(int, exponent))
                    alpha = tuple(sum(exponent[edge] + 1 for edge in chain.edges) for chain in chains)
                    if alpha not in self._period_orbits:
                        images = tuple(tuple(alpha[edge] for edge in permutation) for permutation in permutations)
                        orbit = min(images)
                        self._period_orbits.update({image: orbit for image in images})
                    orbit = self._period_orbits[alpha]
                    factor = _fraction(coefficient) * Fraction(prod(factorial(power) for power in exponent),
                                                               prod(factorial(power - 1) for power in alpha))
                    values[orbit] = values.get(orbit, Fraction(0)) + factor
                if len(values) > self.options["max_period_terms"]:
                    raise EarBudgetExceeded("aggregated weighted periods exceed max_period_terms")
                total[q] = {alpha: coefficient for alpha, coefficient in values.items() if coefficient}
        total = {q: terms for q, terms in total.items() if terms}
        stats["aggregated_schwinger_terms"] = sum(len(numerator) for numerator in schwinger.values())
        stats["aggregated_period_terms"] = sum(len(terms) for terms in total.values())
        stats["denominator_powers"] = sorted(total)
        result = period_edges, total, constant, stats
        self._cache[key] = result
        return result

    def _series(self, aggregate, reason=None):
        edges, total, constant, stats = aggregate
        self.last_stats = dict(stats, backend="exact_convergent_series" if total else "exact_terminating_series")
        if reason is not None:
            self.last_stats["finite_period_unavailable"] = reason
        components = tuple((q, {tuple(power - 1 for power in alpha): coefficient
                                for alpha, coefficient in terms.items()}) for q, terms in total.items())
        return ExactForestSeries(edges, components, constant=constant,
                                 max_edges=self.options["max_series_edges"], max_terms=self.max_terms,
                                 reason=reason)

    def integrate(self, polynomial, *, allow_series=True):
        original = _polynomial(self.graph, polynomial=polynomial)
        aggregate = self._aggregate(original)
        edges, total, constant, stats = aggregate
        if not total:
            self.last_stats = dict(stats, backend="exact_rational")
            return constant
        if self._worker is None:
            self._worker = _PeriodWorker(self.options["period_timeout"])
        result, records = {"1": constant}, []
        try:
            for q, terms in total.items():
                value, record = self._worker.integrate(edges, q, terms)
                records.append(record)
                for label, coefficient in value.items():
                    result[label] = result.get(label, Fraction(0)) + coefficient
        except MissingForestPeriod as error:
            if not allow_series:
                raise
            return self._series(aggregate, str(error))
        except (ForestPeriodBudgetExceeded, ArithmeticError, BrokenPipeError):
            self.close()
            raise
        self.last_stats = dict(stats, backend="finite_periods", period_records=records)
        return _finite(result)

    def moment(self, powers=None, *, allow_series=True):
        return self.integrate(_polynomial(self.graph, powers), allow_series=allow_series)

    def series(self, powers=None, *, polynomial=None):
        original = _polynomial(self.graph, powers, polynomial)
        return self._series(self._aggregate(original))
