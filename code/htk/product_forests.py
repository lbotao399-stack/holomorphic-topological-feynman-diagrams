"""Exact current moments on four-dimensional product forests.

Cubic forests with sources only at their component roots and the ordered
insertions use ``ForestMoments``. Higher-valence forests, or extra sources
on cubic topologies, use the existing exact series-parallel reduction when
available. That fallback retains every source and computes the generic
current-Jacobian density, checked against the original graph.

Remaining forests first contract source-free single and paired stars, with
all retained sources unchanged. The residual current polynomial is integrated
by the exact core-period engine. A finite result is a Fraction or ExactPeriod;
if a finite period is unavailable, the default returns an ExactForestSeries,
the full convergent rational series with certified remainder bounds. It never
returns an unresolved integral, a partial sum as the exact answer, or a zero
in place of an unavailable period. Explicit resource budgets still apply.
"""
from collections.abc import Mapping
from fractions import Fraction
from operator import index

from .forest_model import (CubicForest, UnsupportedForestSources,
                           recognize_product_forest)
from .forest_moments import ForestMoments
from .forest_contractions import ForestContractions, normalize_forest_sources
from .forest_periods import CoreCurrentPeriods, ExactForestSeries, MissingForestPeriod
from .graph import Graph
from .spred import NotSeriesParallel, plan_generic, sp_reduce

__all__ = ["ProductForestMoments", "product_forest_moment",
           "UnsupportedProductForestReduction"]


class UnsupportedProductForestReduction(NotImplementedError):
    """A requested finite-only evaluation needs the available exact series route."""


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


class ProductForestMoments:
    """Dispatch exact moments while preserving the complete original source data.

    ``marked=(A,B)`` keeps the two insertions ordered. Inputs and returned
    moments use the original htk graph's vertex labels, edge order, and
    orientations. ``method`` records the selected exact implementation:
    ``'cubic_boundary_moments'``, ``'series_parallel_generic'``,
    ``'sourcefree_contractions'``, or ``'sourced_core_periods'``.

    The generic fallback reuses one density with ``check`` exact comparisons
    against the original Kirchhoff solution; default ``check=2``. The term
    budget applies to the new contraction and core polynomial construction.
    The legacy direct series-parallel integrator retains its own allocation
    behavior. ``max_states`` is
    the cubic-message budget; further explicit budgets control tree enumeration,
    period requests, and the subset recursion of the exact series. A finite
    period expression is distinct from an exact convergent-series result.
    """
    def __init__(self, G, lam, marked, *, check=2, max_terms=250000, max_states=100000,
                 allow_series=True, max_tree_candidates=1000000,
                 max_period_terms=250000, max_series_edges=18, period_timeout=120,
                 max_matrix_entries=1000000):
        self.structure = recognize_product_forest(G, marked)
        self.graph = self.G = Graph(G.nv, G.edges, G.name)
        self.marked = self.structure.marked
        self.sources = normalize_forest_sources(self.graph, lam)
        self.check = _integer(check, "check count")
        max_terms = _integer(max_terms, "max_terms", 1)
        max_states = _integer(max_states, "max_states", 1)
        max_tree_candidates = _integer(max_tree_candidates, "max_tree_candidates", 1)
        max_period_terms = _integer(max_period_terms, "max_period_terms", 1)
        max_series_edges = _integer(max_series_edges, "max_series_edges", 1)
        period_timeout = _integer(period_timeout, "period_timeout", 1)
        max_matrix_entries = _integer(max_matrix_entries, "max_matrix_entries", 1)
        if not isinstance(allow_series, bool):
            raise TypeError("allow_series must be a boolean")
        self.allow_series = allow_series
        self._core_options = dict(check=self.check, max_terms=max_terms,
                                  max_tree_candidates=max_tree_candidates,
                                  max_period_terms=max_period_terms,
                                  max_series_edges=max_series_edges,
                                  period_timeout=period_timeout,
                                  max_matrix_entries=max_matrix_entries)
        self._forest = self._reduction = self._density = None
        self._contractions = self._core = None
        self._results = {}
        self.last_stats = {}
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
        except NotSeriesParallel:
            self._prepare_core()
            self.method = ("sourcefree_contractions" if self._contractions.steps
                           else "sourced_core_periods")
        else:
            self.method = "series_parallel_generic"

    def _prepare_core(self):
        if self._contractions is None:
            self._contractions = ForestContractions(
                self.graph, self.sources, self.marked,
                max_terms=self._core_options["max_terms"], check=self.check)

    def _core_for(self, result):
        if self._core is None:
            self._core = CoreCurrentPeriods(result.graph, result.sources,
                                            marked=result.marked, **self._core_options)
        return self._core

    def close(self):
        """Release the private period worker, if one was needed."""
        if self._core is not None:
            self._core.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, traceback):
        self.close()

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
        """Return a Fraction, finite ExactPeriod, or full ExactForestSeries."""
        key = self._powers(powers)
        if key not in self._results:
            if self._forest is not None:
                value = self._forest.moment(dict(key))
                stats = dict(self._forest.last_stats, backend="cubic_boundary_moments")
            elif self._reduction is not None:
                if self._density is None:
                    self._density = self._reduction.density(
                        self.sources, T=0, method="generic", check=self.check)
                value = self._reduction.integrate(
                    self.sources, T=0, powers=dict(key), method="generic",
                    check=0, dens=self._density)
                value = Fraction(value)
                stats = dict(backend="series_parallel_generic", original_loops=self.graph.loops,
                             all_original_sources_retained=True)
            else:
                result = self._contractions.contract(dict(key))
                core = self._core_for(result)
                try:
                    value = core.integrate(result.polynomial, allow_series=self.allow_series)
                except MissingForestPeriod as error:
                    raise UnsupportedProductForestReduction(
                        "finite-only product-forest evaluation has an unavailable core period; "
                        "use allow_series=True or .series(powers) for its exact convergent series") from error
                stats = dict(core.last_stats, contraction=result.stats,
                             all_original_sources_retained=True)
            self._results[key] = value, stats
        value, stats = self._results[key]
        self.last_stats = dict(stats)
        return value

    def moments(self, words):
        """Evaluate a finite sequence of words, reusing the same exact density."""
        return tuple(self.moment(word) for word in words)

    def series(self, powers=None):
        """Return an explicit convergent rational series, including finite termination.

        The returned object represents the complete series. Its ``terms`` and
        ``interval`` methods provide exact rational coefficients and certified
        finite-order bounds, without a decimal estimate or an equality claim
        for a partial sum. The cubic and series-parallel fast paths reuse their
        complete rational answer as a terminating series. A sourced non-SP
        core follows the full series route even if a finite period is known.
        """
        key = self._powers(powers)
        if self._forest is not None or self._reduction is not None:
            value = self.moment(dict(key))
            self.last_stats = dict(self.last_stats,
                                   finite_backend=self.last_stats["backend"],
                                   backend="exact_terminating_series",
                                   all_original_sources_retained=True)
            return ExactForestSeries((), (), constant=value,
                                     max_edges=self._core_options["max_series_edges"],
                                     max_terms=self._core_options["max_terms"])
        self._prepare_core()
        result = self._contractions.contract(dict(key))
        core = self._core_for(result)
        value = core.series(result.polynomial)
        self.last_stats = dict(core.last_stats, contraction=result.stats,
                               all_original_sources_retained=True)
        return value


def product_forest_moment(G, lam, marked, powers=None, **options):
    """Convenience wrapper for ``ProductForestMoments(...).moment(powers)``."""
    with ProductForestMoments(G, lam, marked, **options) as engine:
        return engine.moment(powers)
