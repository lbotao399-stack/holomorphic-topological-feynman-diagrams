"""Exact current moments on supported four-dimensional product forests.

Cubic forests with sources only at their component roots and the ordered
insertions use ``ForestMoments``. Higher-valence forests, or extra sources
on cubic topologies, use the existing exact series-parallel reduction when
available. That fallback retains every source and computes the generic
current-Jacobian density, checked against the original graph.

This adapter does not extend the source-free star theorem to sourced stars.
A sourced non-series-parallel remainder raises an explicit unsupported
reduction error. All returned values are complete finite current moments
at zero displacement, equivalently exact displacement Taylor data.
"""
from collections.abc import Mapping
from fractions import Fraction
from operator import index

from .forest_model import (CubicForest, UnsupportedForestSources,
                           recognize_product_forest)
from .forest_moments import ForestMoments
from .graph import Graph
from .spred import NotSeriesParallel, plan_generic, sp_reduce

__all__ = ["ProductForestMoments", "product_forest_moment",
           "UnsupportedProductForestReduction"]


class UnsupportedProductForestReduction(NotImplementedError):
    """The retained source pattern requires an unevaluated sourced non-SP reduction."""


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


def _balanced_sources(G, lam):
    """Exact H=2 source normalization without a restriction on source locations."""
    def row(values):
        values = tuple(values)
        if len(values) != 2:
            raise ValueError("four-dimensional holomorphic sources have two components")
        if any(isinstance(value, float) for value in values):
            raise TypeError("source components must be exact rational values")
        return tuple(Fraction(value) for value in values)

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


class ProductForestMoments:
    """Dispatch exact moments while preserving the complete original source data.

    ``marked=(A,B)`` keeps the two insertions ordered. Inputs and returned
    moments use the original htk graph's vertex labels, edge order, and
    orientations. ``method`` records the selected exact implementation:
    ``'cubic_boundary_moments'`` or ``'series_parallel_generic'``.

    The generic fallback reuses one density with ``check`` exact comparisons
    against the original Kirchhoff solution; default ``check=2``. The term
    and state budgets apply to the cubic boundary-moment implementation.
    """
    def __init__(self, G, lam, marked, *, check=2, max_terms=250000, max_states=100000):
        self.structure = recognize_product_forest(G, marked)
        self.graph = self.G = Graph(G.nv, G.edges, G.name)
        self.marked = self.structure.marked
        self.sources = _balanced_sources(self.graph, lam)
        self.check = _integer(check, "check count")
        max_terms = _integer(max_terms, "max_terms", 1)
        max_states = _integer(max_states, "max_states", 1)
        self._forest = self._reduction = self._density = None
        self._results = {}
        if self.structure.cubic:
            model = CubicForest(self.graph, self.marked)
            try:
                model.normalize_sources(self.sources)
            except UnsupportedForestSources:
                pass
            else:
                self._forest = ForestMoments(self.graph, self.sources, self.marked,
                                             max_terms=max_terms, max_states=max_states)
                self.method = "cubic_boundary_moments"
                return
        try:
            # A source-complete generic density does not require a tight
            # local plan. Avoid the backtracking tight-plan search, especially
            # when the graph has a non-series-parallel remainder.
            plan = plan_generic(self.graph)
            self._reduction = sp_reduce(self.graph, D=2, plan=plan)
        except NotSeriesParallel as error:
            raise UnsupportedProductForestReduction(
                "this product forest retains sources requiring a sourced non-series-parallel "
                "reduction; the cubic source-free-star formula does not apply") from error
        self.method = "series_parallel_generic"

    def _powers(self, powers):
        if powers is None:
            return ()
        if not isinstance(powers, Mapping):
            raise TypeError("current powers must be a mapping {(edge, component): power}")
        normalized = {}
        for location, power in powers.items():
            if not isinstance(location, tuple) or len(location) != 2:
                raise ValueError("current powers are indexed by (edge, component)")
            edge = _integer(location[0], "edge index")
            component = _integer(location[1], "current component")
            power = _integer(power, "current power")
            if edge >= self.graph.ne or component >= 2:
                raise ValueError("edge or current component lies outside this graph")
            if power:
                normalized[edge, component] = power
        return tuple(sorted(normalized.items()))

    def moment(self, powers=None):
        """Integrate the requested original-edge current monomial exactly."""
        key = self._powers(powers)
        if key not in self._results:
            if self._forest is not None:
                value = self._forest.moment(dict(key))
            else:
                if self._density is None:
                    self._density = self._reduction.density(
                        self.sources, T=0, method="generic", check=self.check)
                value = self._reduction.integrate(
                    self.sources, T=0, powers=dict(key), method="generic",
                    check=0, dens=self._density)
            self._results[key] = Fraction(value)
        return self._results[key]

    def moments(self, words):
        """Evaluate a finite sequence of words, reusing the same exact density."""
        return tuple(self.moment(word) for word in words)


def product_forest_moment(G, lam, marked, powers=None, **options):
    """Convenience wrapper for ``ProductForestMoments(...).moment(powers)``."""
    return ProductForestMoments(G, lam, marked, **options).moment(powers)
