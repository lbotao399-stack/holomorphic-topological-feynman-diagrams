"""Exact rational series for integer and half-integer bulk graph densities.

Every dominant-tree chart is retained.  Squaring its coordinates resolves
the half-integer boundary powers.  The normalized Kirchhoff denominator
lies between 1 and the number of spanning trees on the entire closed cube.
A binomial expansion about a positive integer square therefore has rational
coefficients and a finite rational bound for its complete omitted tail.
"""
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations, permutations
from math import ceil, comb, factorial, isqrt, prod
from types import SimpleNamespace

from flint import fmpq_mpoly_ctx

from .core_sector_periods import _cycle_count
from .forest_cubes import GreedyCube, _integer, _rational
from .forest_cube_expansion import _moment, _q
from .topological_polynomial import rational

__all__ = ["ExactBulkSeries", "BulkSeriesTerm", "BulkSeriesBudgetExceeded"]


class BulkSeriesBudgetExceeded(RuntimeError):
    """An explicit budget prevents a complete series calculation."""


@dataclass(frozen=True)
class BulkSeriesTerm:
    order: int
    term: Fraction
    partial_sum: Fraction
    error_bound: Fraction

    @property
    def lower(self):
        return self.partial_sum - self.error_bound

    @property
    def upper(self):
        return self.partial_sum + self.error_bound


class ExactBulkSeries:
    """The exact integral, represented by a complete convergent rational series.

    ``density`` is a BulkPolynomial with exact rational sources.  ``word``
    is a sequence of original (edge, spinor-index) current slots.  Resource
    failures raise; no incomplete atlas or finite partial sum is returned
    as the exact kernel.  This is a bulk construction at zero topological
    momentum, with no defect or inserted topological-current claim.
    """

    def __init__(self, density, word=(), *, max_tree_candidates=1_000_000,
                 max_charts=1_000_000, max_terms=1_000_000,
                 max_total_terms=3_000_000, max_order=10_000, reason=None):
        if density.symbolic_sources:
            raise ValueError("certified absolute bounds require rational source values")
        self.density = density
        self.graph = density.graph
        self.word = tuple(word)
        self.numerator, self.power = density.inserted(self.word)
        self.max_tree_candidates = _integer(max_tree_candidates, "max_tree_candidates", 1)
        self.max_charts = _integer(max_charts, "max_charts", 1)
        self.max_terms = _integer(max_terms, "max_terms", 1)
        self.max_total_terms = _integer(max_total_terms, "max_total_terms", 1)
        self.max_order = _integer(max_order, "max_order")
        self.reason = reason
        self.dimension = self.graph.ne - 1
        self.context = fmpq_mpoly_ctx.get(tuple("cube%d" % j for j in range(self.dimension)))
        self.components = None
        self.statistics = {}

    def _bounded(self, polynomial):
        if len(polynomial) > self.max_terms:
            raise BulkSeriesBudgetExceeded("complete polynomial exceeds max_terms")
        return polynomial

    def _buffer(self, values):
        if sum(len(p) for p in values) > self.max_total_terms:
            raise BulkSeriesBudgetExceeded("complete polynomial state exceeds max_total_terms")

    def _atlas(self):
        graph = self.graph
        if comb(graph.ne, graph.nv - 1) > self.max_tree_candidates:
            raise BulkSeriesBudgetExceeded("spanning-tree enumeration exceeds max_tree_candidates")
        trees = tuple(tree for tree in combinations(range(graph.ne), graph.nv - 1)
                      if _cycle_count(graph, tree) == 0)
        count = len(trees)*factorial(graph.loops)
        if count > self.max_charts:
            raise BulkSeriesBudgetExceeded("complete dominant-tree atlas exceeds max_charts")
        if not trees:
            raise ValueError("the graph is disconnected")
        self.trees = trees
        self.tree_count = len(trees)
        self.chart_count = count
        return SimpleNamespace(G=graph, dimension=self.dimension, loops=graph.loops, trees=trees)

    def _chart_polynomial(self, chart):
        context = self.context
        # A positive projective monomial chart has absolute Jacobian
        # product(t_e)/product(x_j).  Its orientation is already present
        # in the original density P; no new current-form sign is guessed.
        dilation = 2 if self.power.denominator == 2 else 1
        if self.power.denominator not in (1, 2):
            raise ValueError("only integer and half-integer bulk poles are supported")
        cotree_exponent = tuple(sum(chart.time_exponents[e][j] for e in chart.cotree_order)
                               for j in range(self.dimension))
        denominator_data = {}
        for exponent in chart.weight_exponents:
            key = tuple(dilation*x for x in exponent)
            denominator_data[key] = denominator_data.get(key, 0) + 1
        denominator = context.from_dict(denominator_data)
        numerator_data = {}
        for original, coefficient in self.numerator.to_dict().items():
            exponent = tuple(dilation*(sum((int(original[e])+1)*chart.time_exponents[e][j]
                                           for e in range(self.graph.ne))
                                        - self.power*cotree_exponent[j]) - 1
                             for j in range(self.dimension))
            if any(x.denominator != 1 or x < 0 for x in exponent):
                raise ValueError("a complete chart has a nonintegrable boundary power; "
                                 "no convergent-series value is asserted")
            key = tuple(int(x) for x in exponent)
            numerator_data[key] = numerator_data.get(key, Fraction(0)) + rational(coefficient)*dilation**self.dimension
        numerator = context.from_dict({a: _q(c) for a, c in numerator_data.items() if c})
        return self._bounded(denominator), self._bounded(numerator)

    def prepare(self):
        if self.components is not None:
            return self
        if not self.numerator:
            self.tree_count = self.chart_count = 0
            self.square_root_centre = self.centre = 1
            self.ratio = Fraction(0)
            self.components = ()
            self.statistics = {"complete_zero_density": True, "original_current_slots": self.word}
            return self
        atlas = self._atlas()
        centre_root = max(1, isqrt((self.tree_count+1)//2))
        while 2*centre_root**2 < self.tree_count+1:
            centre_root += 1
        self.square_root_centre = centre_root
        self.centre = centre_root**2
        self.ratio = max(Fraction(self.centre-1, self.centre),
                         Fraction(self.tree_count-self.centre, self.centre))
        if not 0 <= self.ratio < 1:
            raise ArithmeticError("the geometric ratio is not strictly below one")
        grouped = {}
        denominators = {}
        position = 0
        for tree_index, tree in enumerate(self.trees):
            cotree = tuple(e for e in range(self.graph.ne) if e not in tree)
            for order in permutations(cotree):
                chart = GreedyCube(atlas, position, tree_index, order)
                denominator, numerator = self._chart_polynomial(chart)
                position += 1
                key = tuple(sorted((tuple(map(int, a)), rational(c)) for a, c in denominator.to_dict().items()))
                denominators.setdefault(key, denominator)
                grouped[key] = self._bounded(grouped.get(key, self.context.constant(0)) + numerator)
                self._buffer(list(grouped.values()) + list(denominators.values()))
        if position != self.chart_count:
            raise ArithmeticError("the atlas is incomplete")
        self.components = tuple((denominators[k], p, _moment(p, absolute=True))
                                for k, p in grouped.items() if p)
        self.scale = Fraction(1, centre_root**int(2*self.power))
        self.statistics = {"spanning_trees": self.tree_count, "complete_charts": self.chart_count,
                           "combined_components": len(self.components),
                           "pole": str(self.power), "centre": self.centre,
                           "absolute_ratio": str(self.ratio),
                           "squared_coordinates": self.power.denominator == 2,
                           "all_original_sources_retained": True,
                           "all_original_current_slots_retained": True}
        return self

    def _order(self, order):
        order = _integer(order, "series order")
        if order > self.max_order:
            raise BulkSeriesBudgetExceeded("requested complete series exceeds max_order")
        return order

    def error_bound(self, order):
        order = self._order(order)
        self.prepare()
        if not self.components or not self.ratio:
            return Fraction(0)
        integer_power = ceil(self.power)
        # (q)_k/k! <= binomial(ceil(q)+k-1,k).  The following finite
        # rational expression is the exact tail of that majorant.
        majorant_tail = self.ratio**(order+1) * sum(
            (Fraction(comb(order+j, j), 1)/(1-self.ratio)**(integer_power-j)
             for j in range(integer_power)), Fraction(0))
        return self.scale*majorant_tail*sum((mass for _, _, mass in self.components), Fraction(0))

    def terms(self, last):
        last = self._order(last)
        self.prepare()
        if not self.components:
            for order in range(last+1):
                yield BulkSeriesTerm(order, Fraction(0), Fraction(0), Fraction(0))
            return
        states = [numerator for _, numerator, _ in self.components]
        factors = [self.context.constant(1) - denominator/_q(self.centre)
                   for denominator, _, _ in self.components]
        resident = [polynomial for denominator, numerator, _ in self.components
                    for polynomial in (denominator, numerator)] + factors
        self._buffer(resident + states)
        coefficient = self.scale
        partial = Fraction(0)
        for order in range(last+1):
            if order:
                states = [self._bounded(p*f) for p, f in zip(states, factors)]
                self._buffer(resident + states)
                coefficient *= (self.power + order - 1)/order
            term = coefficient*sum((_moment(p) for p in states), Fraction(0))
            partial += term
            yield BulkSeriesTerm(order, term, partial, self.error_bound(order))

    def interval(self, order):
        row = None
        for row in self.terms(order):
            pass
        return row.lower, row.upper

    def required_order(self, width):
        width = _rational(width)
        if width <= 0:
            raise ValueError("the enclosing width must be positive")
        if 2*self.error_bound(0) <= width:
            return 0
        low, high = 1, min(1, self.max_order)
        if not high:
            raise BulkSeriesBudgetExceeded("requested width exceeds max_order")
        while 2*self.error_bound(high) > width:
            if high == self.max_order:
                raise BulkSeriesBudgetExceeded("requested width exceeds max_order")
            low, high = high+1, min(2*high, self.max_order)
        while low < high:
            middle = (low+high)//2
            if 2*self.error_bound(middle) <= width:
                high = middle
            else:
                low = middle+1
        return low
