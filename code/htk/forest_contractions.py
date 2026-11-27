"""Exact source-free local contractions inside a sourced product forest.

A retained parent may have any allowed valence and carry any source.  Only
the centre of a single star, or both centres of a paired star, must have
zero source.  The complete remaining source rows are retained.  Each step
introduces two or four local ratios, substitutes the original edge currents,
and integrates those ratios immediately.  No global ratio chart is built.

The result is a homogeneous polynomial in the remaining edge currents,
ready for a source-complete core-period evaluator.  Its coefficients include
the htk loop normalization and every edge-order/orientation sign.  See
``star.SingleMove`` and ``star.PairedMove`` for the exact circuit maps, and
``sources/LLC_main.tex`` for the source-free local identities.
"""

from collections.abc import Mapping
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from math import factorial
from operator import index
import random

from . import brute
from .forest_model import recognize_product_forest
from .graph import Graph
from .linalg import perm_sign
from .star import SingleMove, PairedMove, find_single, find_paired

try:
    from flint import fmpq, fmpq_mpoly_ctx
except ImportError as error:
    _FLINT_ERROR = error
else:
    _FLINT_ERROR = None

__all__ = ["ForestContractions", "ContractionResult", "ContractionBudgetExceeded",
           "contract_forest", "normalize_forest_sources"]


class ContractionBudgetExceeded(RuntimeError):
    """An explicit polynomial-term budget was exceeded; no partial result exists."""


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
    if isinstance(value, float):
        raise TypeError(name + " must be an exact rational value")
    return Fraction(value)


def normalize_forest_sources(G, lam):
    """Full conserved H=2 rows, without restrictions on source locations.

    Accept non-root rows, full balanced rows, or a vertex dictionary, with
    the same missing-root convention as ``htk.brute.sources``.  Floating
    input is rejected rather than silently converted to a binary rational.
    """
    def row(values):
        values = tuple(values)
        if len(values) != 2:
            raise ValueError("four-dimensional holomorphic sources have two components")
        return tuple(_exact(value, "source component") for value in values)

    if lam is None:
        return ((Fraction(0), Fraction(0)),) * G.nv
    if isinstance(lam, Mapping):
        for vertex in lam:
            if _integer(vertex, "source vertex") >= G.nv:
                raise ValueError("source vertex lies outside this graph")
        full = [row(lam.get(vertex, (0, 0))) for vertex in range(G.nv - 1)]
        ground = tuple(-sum(values[i] for values in full) for i in range(2))
        if G.root in lam and row(lam[G.root]) != ground:
            raise ValueError("sources do not sum to zero")
        return tuple(full) + (ground,)
    full = [row(values) for values in lam]
    if len(full) == G.nv - 1:
        full.append(tuple(-sum(values[i] for values in full) for i in range(2)))
    elif len(full) != G.nv:
        raise ValueError("sources require nv-1 or nv two-component rows")
    if any(sum(values[i] for values in full) for i in range(2)):
        raise ValueError("sources do not sum to zero")
    return tuple(full)


def _q(value):
    value = Fraction(value)
    return fmpq(value.numerator, value.denominator)


def _fraction(value):
    return Fraction(int(value.p), int(value.q))


def _bracket(left, right):
    return left[0] * right[1] - left[1] * right[0]


@lru_cache(maxsize=None)
def _ratio_weight(exponents):
    if len(exponents) == 2:
        first, second = exponents
        return fmpq(1, (first + 2) * (second + 2))
    first, second, third, fourth = exponents
    # Integral over the square of r0^a r1^b (r0-r1), times
    # the triangle integral r2^c r3^d (1-r2-r3).
    if first == second:
        return fmpq(0)
    return (fmpq(first - second,
                 (first + 1) * (first + 2) * (second + 1) * (second + 2))
            * fmpq(factorial(third) * factorial(fourth), factorial(third + fourth + 3)))


def _move_sign(move):
    """Kernel sign with the positive ratio-domain measure, not signed J_aff."""
    if isinstance(move, SingleMove):
        keys = ("ux", "wx", "vu", "vw", "vx")
        loop_sign = (-1) ** move.G.loops
    else:
        keys = ("1u", "1A", "1B", "2u", "2A", "2B", "AB")
        loop_sign = -1
    order = [move.e[key] for key in keys] + move.others
    return loop_sign * perm_sign(order)


def _current_map(move, reduced, ratios):
    """Current vectors in original stage edge order, valid only at zero centres."""
    result = [None] * move.G.ne
    if isinstance(move, SingleMove):
        marked, left, right = reduced[:3]
        first, second = ratios
        local = {
            "vu": tuple(-marked[i] - (1 - first) * left[i] for i in range(2)),
            "vw": tuple(marked[i] - (1 - second) * right[i] for i in range(2)),
            "vx": tuple((1 - first) * left[i] + (1 - second) * right[i] for i in range(2)),
            "ux": tuple(first * item for item in left),
            "wx": tuple(second * item for item in right),
        }
    else:
        left, right, marked = reduced[:3]
        first, second, third, fourth = ratios
        local = {
            "1u": tuple(-first * left[i] - second * right[i] for i in range(2)),
            "2u": tuple(-(1 - first) * left[i] - (1 - second) * right[i] for i in range(2)),
            "1A": tuple(first * left[i] - third * marked[i] for i in range(2)),
            "2A": tuple((1 - first) * left[i] - fourth * marked[i] for i in range(2)),
            "1B": tuple(second * right[i] + third * marked[i] for i in range(2)),
            "2B": tuple((1 - second) * right[i] + fourth * marked[i] for i in range(2)),
            "AB": tuple((1 - third - fourth) * item for item in marked),
        }
    for name, vector in local.items():
        result[move.e[name]] = tuple(move.s[name] * value for value in vector)
    for position, edge in enumerate(move.others):
        result[edge] = tuple(reduced[3 + position])
    return result


@dataclass(frozen=True)
class ContractionResult:
    """Complete local pushforward; integrate ``polynomial`` against the core kernel.

    ``polynomial`` maps length ``2*graph.ne`` nonnegative exponent tuples to
    Fraction coefficients, in edge-major, component-zero/component-one
    order.  The zero polynomial is an empty dictionary.  Vertex maps refer
    to the user's original labels; a removed vertex maps to ``None``.
    """
    graph: Graph
    sources: tuple
    marked: tuple
    polynomial: dict
    degree: int
    steps: tuple
    stats: dict
    original_to_reduced_vertices: tuple


class ForestContractions:
    """Plan and evaluate exact passive-star contractions for a marked forest.

    ``contract(powers)`` accepts ``{(original_edge, component): power}``.
    Alternatively, ``contract(polynomial=...)`` accepts an arbitrary
    homogeneous sparse polynomial in the original edge currents.  Sources
    can occupy every retained vertex, including higher-valence parents.

    The deterministic greedy plan removes zero-source single stars first,
    then paired stars, and never removes either marked insertion.  Failure
    to find another such move is a valid remainder, not a claim that the
    remainder is irreducible by every possible identity.  ``check`` exact
    rational circuit/density comparisons are performed per planned move.
    """

    def __init__(self, G, lam, marked, *, max_terms=250000, check=2):
        if _FLINT_ERROR is not None:
            raise ImportError("forest contractions require the optional python-flint package") from _FLINT_ERROR
        self.structure = recognize_product_forest(G, marked)
        self.graph = self.G = Graph(G.nv, G.edges, G.name)
        self.marked = tuple(self.structure.marked)
        self.sources = normalize_forest_sources(self.graph, lam)
        self.max_terms = _integer(max_terms, "max_terms", 1)
        self.check = _integer(check, "check count")
        self.last_stats = {}

        # The local star classes keep the grounded vertex.  Ground at B by
        # a vertex permutation only: original edge order and directions stay.
        vertex_order = tuple(v for v in range(G.nv) if v != self.marked[1]) + (self.marked[1],)
        initial_map = {old: new for new, old in enumerate(vertex_order)}
        current = Graph(G.nv, [(initial_map[a], initial_map[b]) for a, b in G.edges], G.name)
        current_sources = tuple(self.sources[v] for v in vertex_order)
        current_marked = tuple(initial_map[v] for v in self.marked)
        original_vertices = vertex_order
        planned, steps = [], []
        rng = random.Random(307984)
        while True:
            move = self._next_move(current, current_sources, current_marked)
            if move is None:
                break
            reduced_sources = tuple(current_sources[v] for v in range(current.nv)
                                    if move.vmap[v] is not None)
            reduced_marked = tuple(move.vmap[v] for v in current_marked)
            sign = _move_sign(move)
            if self.check:
                self._check_move(move, current_sources, reduced_sources, sign, rng)
            next_vertices = tuple(original_vertices[v] for v in range(current.nv)
                                  if move.vmap[v] is not None)
            steps.append({
                "kind": "single" if move.k == 2 else "paired",
                "original_centres": tuple(original_vertices[v] for v in move.centres),
                "stage_centres": tuple(move.centres),
                "before_loops": current.loops,
                "after_loops": move.red.loops,
                "before_edges": tuple(current.edges),
                "after_edges": tuple(move.red.edges),
                "local_edge_indices": dict(move.e),
                "local_edge_orientation_signs": dict(move.s),
                "kernel_sign": sign,
                "local_ratio_count": move.k,
                "remaining_original_vertices": next_vertices,
                "exact_point_checks": self.check,
            })
            planned.append((move, sign))
            current, current_sources, current_marked = move.red, reduced_sources, reduced_marked
            original_vertices = next_vertices

        self._planned = tuple(planned)
        self.steps = tuple(steps)
        self.reduced_graph = current
        self.reduced_sources = current_sources
        self.reduced_marked = current_marked
        remaining = {old: new for new, old in enumerate(original_vertices)}
        self.original_to_reduced_vertices = tuple(remaining.get(v) for v in range(G.nv))
        self.reduced_to_original_vertices = original_vertices

    @staticmethod
    def _next_move(G, sources, marked):
        protected = set(marked)
        for vertices in find_single(G):
            if vertices[0] not in protected and not any(sources[vertices[0]]):
                return SingleMove(G, *vertices)
        for vertices in find_paired(G):
            if not any(v in protected or any(sources[v]) for v in vertices[:2]):
                return PairedMove(G, *vertices)
        return None

    def _check_move(self, move, sources, reduced_sources, sign, rng):
        for _ in range(self.check):
            ratios, times = move.random_point(rng)
            original_times = move.tmap(ratios, times)
            original_currents = brute.currents(move.G, original_times, sources)
            reduced_currents = brute.currents(move.red, times, reduced_sources)
            lifted = _current_map(move, reduced_currents, ratios)
            if tuple(map(tuple, original_currents)) != tuple(lifted):
                raise ArithmeticError("source-free forest current map failed exact verification")
            if move.k == 2:
                factor = (sign * ratios[0] * ratios[1]
                          * _bracket(reduced_currents[1], reduced_currents[2]))
            else:
                factor = (sign * (ratios[0] - ratios[1]) * (1 - ratios[2] - ratios[3])
                          * _bracket(reduced_currents[0], reduced_currents[2])
                          * _bracket(reduced_currents[1], reduced_currents[2]))
            reference = (brute.density(move.G, original_times, sources, T=0).rational()
                         * abs(move.jacobian(ratios, times)))
            transformed = factor * brute.density(move.red, times, reduced_sources, T=0).rational()
            if reference != transformed:
                raise ArithmeticError("source-free forest density/sign failed exact verification")

    def _normalize_polynomial(self, powers, polynomial):
        size = 2 * self.graph.ne
        if polynomial is not None:
            if powers is not None:
                raise ValueError("provide either current powers or a current polynomial")
            if not isinstance(polynomial, Mapping):
                raise TypeError("current polynomial must map exponent tuples to exact coefficients")
            result, degrees = {}, set()
            for exponents, coefficient in polynomial.items():
                if not isinstance(exponents, tuple) or len(exponents) != size:
                    raise ValueError("current exponents require a tuple of length 2*graph.ne")
                exponents = tuple(_integer(value, "current power") for value in exponents)
                coefficient = _exact(coefficient, "polynomial coefficient")
                if coefficient:
                    result[exponents] = result.get(exponents, Fraction(0)) + coefficient
                    degrees.add(sum(exponents))
            if len(degrees) > 1:
                raise ValueError("the current polynomial must be homogeneous")
            return result, next(iter(degrees), 0)
        if powers is not None and not isinstance(powers, Mapping):
            raise TypeError("current powers must map (edge, component) to nonnegative integers")
        exponents = [0] * size
        for location, power in ({} if powers is None else powers).items():
            if not isinstance(location, tuple) or len(location) != 2:
                raise ValueError("current powers are indexed by (edge, component)")
            edge = _integer(location[0], "edge index")
            component = _integer(location[1], "current component")
            power = _integer(power, "current power")
            if edge >= self.graph.ne or component >= 2:
                raise ValueError("edge or current component lies outside this graph")
            exponents[2 * edge + component] = power
        return {tuple(exponents): Fraction(1)}, sum(exponents)

    def contract(self, powers=None, *, polynomial=None):
        """Return the exact reduced current polynomial and complete source data."""
        original, degree = self._normalize_polynomial(powers, polynomial)
        stats = {
            "original_loops": self.graph.loops,
            "reduced_loops": self.reduced_graph.loops,
            "input_current_degree": degree,
            "output_current_degree": degree + 2 * (self.graph.loops - self.reduced_graph.loops),
            "contraction_steps": len(self.steps),
            "single_steps": sum(move.k == 2 for move, _ in self._planned),
            "paired_steps": sum(move.k == 4 for move, _ in self._planned),
            "largest_polynomial_terms": len(original),
            "largest_local_ratio_count": max((move.k for move, _ in self._planned), default=0),
            "global_ratio_coordinates": 0,
            "exact_point_checks": self.check * len(self.steps),
            "step_term_counts": [],
        }

        def checked(poly):
            size = len(poly)
            stats["largest_polynomial_terms"] = max(stats["largest_polynomial_terms"], size)
            if size > self.max_terms:
                raise ContractionBudgetExceeded(
                    "forest contraction produced %d terms, exceeding max_terms=%d" % (size, self.max_terms))
            return poly

        context = fmpq_mpoly_ctx.get(tuple("j%d" % i for i in range(2 * self.graph.ne)))
        value = checked(context.from_dict({exponents: _q(coefficient) for exponents, coefficient in original.items()}))
        for move, sign in self._planned:
            size = 2 * move.red.ne
            names = tuple("j%d" % i for i in range(size))
            local = fmpq_mpoly_ctx.get(names + tuple("r%d" % i for i in range(move.k)))
            variables = local.gens()
            currents = [variables[2 * edge:2 * edge + 2] for edge in range(move.red.ne)]
            lifted = _current_map(move, currents, variables[size:])
            transformed = checked(value.compose(*(item for pair in lifted for item in pair)))
            coefficients = {}
            for exponents, coefficient in transformed.terms():
                weight = _ratio_weight(exponents[size:])
                if weight:
                    key = exponents[:size]
                    coefficients[key] = coefficients.get(key, fmpq(0)) + coefficient * weight
            context = fmpq_mpoly_ctx.get(names)
            value = checked(context.from_dict(coefficients))
            variables = context.gens()
            currents = [variables[2 * edge:2 * edge + 2] for edge in range(move.red.ne)]
            if move.k == 2:
                factor = _bracket(currents[1], currents[2])
            else:
                factor = _bracket(currents[0], currents[2]) * _bracket(currents[1], currents[2])
            value = checked(sign * value * factor)
            stats["step_term_counts"].append(len(value))
        result = {tuple(map(int, exponents)): _fraction(coefficient)
                  for exponents, coefficient in value.terms() if coefficient}
        final_degree = stats["output_current_degree"]
        if any(sum(exponents) != final_degree for exponents in result):
            raise ArithmeticError("forest contraction violated current homogeneity")
        stats["output_terms"] = len(result)
        self.last_stats = dict(stats)
        return ContractionResult(self.reduced_graph, self.reduced_sources, self.reduced_marked,
                                 result, final_degree, self.steps, stats,
                                 self.original_to_reduced_vertices)


def contract_forest(G, lam, marked, powers=None, *, polynomial=None, **options):
    """Convenience wrapper for ``ForestContractions(...).contract(...)``."""
    return ForestContractions(G, lam, marked, **options).contract(powers, polynomial=polynomial)
