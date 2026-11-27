"""Explicit geometric rational series for every complete product-forest cube.

The dominant-tree atlas has normalized Kirchhoff polynomial 1 <= U <= N,
where N is the number of spanning trees.  On each compact cube, the full
current-area form times a finite original-current polynomial is P/U**s.
The augmented current determinant gives s=2*loops+1 before inserted
currents; exact polynomial divisions can lower this power further.

With c=(N+1)/2, expand U**(-s) in (c-U)/c.  Every coefficient is an EXACT
polynomial cube moment.  The absolute ratio is at most (N-1)/(N+1), and
a finite rational negative-binomial tail gives a certified remainder and
an a-priori sufficient order for every positive rational requested width.
This replaces spatial refinement by algebraic series order.  It proves
no universal finite expression in zeta values.

All sources, original currents, spinor entries, cube signs, and polynomial
degrees are retained.  Resource limits raise an exception; neither a
partial atlas nor a finite series sum is returned as the exact integral.
Only this module is new: no hashed historical backend is modified.
"""

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from math import comb, prod

from .forest_cubes import (ForestCubes, ForestCubeBudgetExceeded,
                           _integer, _rational)

try:
    from flint import fmpq, fmpq_mpoly_ctx
except ImportError as error:
    _FLINT_ERROR = error
else:
    _FLINT_ERROR = None


__all__ = ["ForestCubeExpansion", "ExactCubeExpansion", "CubeExpansionTerm",
           "ForestCubeExpansionBudgetExceeded"]


class ForestCubeExpansionBudgetExceeded(ForestCubeBudgetExceeded):
    """An explicit budget prevents a complete geometric-series computation."""


def _q(value):
    value = Fraction(value)
    return fmpq(value.numerator, value.denominator)


def _fraction(value):
    return Fraction(int(value.p), int(value.q))


def _moment(polynomial, absolute=False):
    """Integrate a complete polynomial, or its coefficientwise majorant."""
    return sum(((_fraction(abs(coefficient) if absolute else coefficient)
                 / prod(int(power)+1 for power in exponent))
                for exponent, coefficient in polynomial.terms()), Fraction(0))


def _polynomial_data(polynomial):
    return [[[int(power) for power in exponent], str(_fraction(coefficient))]
            for exponent, coefficient in sorted(polynomial.terms())]


@dataclass(frozen=True)
class CubeExpansionTerm:
    order: int
    term: Fraction
    partial_sum: Fraction
    error_bound: Fraction
    lower: Fraction
    upper: Fraction
    polynomial_terms: int
    charts: int


@dataclass(frozen=True)
class _Component:
    denominator: object
    power: int
    numerator: object
    absolute_moment: Fraction


class ForestCubeExpansion:
    """Reusable exact geometric integration of a full original-current input.

    The graph/source/marked arguments and remaining keyword budgets are
    those of ForestCubes.  max_terms bounds each sparse polynomial;
    max_total_terms bounds accumulated component or series-state support;
    max_order is an explicit requested-series resource budget.  These are
    computational limits, not restrictions of the mathematical formula.
    """

    def __init__(self, graph, sources, marked=None, *, max_terms=1_000_000,
                 max_total_terms=2_000_000, max_order=10_000, **cube_options):
        if _FLINT_ERROR is not None:
            raise ImportError("exact cube expansion requires python-flint") from _FLINT_ERROR
        self.max_terms = _integer(max_terms, "max_terms", 1)
        self.max_total_terms = _integer(max_total_terms, "max_total_terms", 1)
        self.max_order = _integer(max_order, "max_order")
        self.atlas = ForestCubes(graph, sources, marked, **cube_options)
        self.graph = self.G = self.atlas.G
        self.sources = self.atlas.sources
        self.marked = self.atlas.marked
        self.context = fmpq_mpoly_ctx.get(tuple("cube%d" % j
                                              for j in range(max(1, self.atlas.dimension))),
                                         ordering="lex")
        self.centre = Fraction(self.atlas.tree_count+1, 2)
        self.ratio = Fraction(self.atlas.tree_count-1, self.atlas.tree_count+1)

    def _bounded(self, polynomial):
        if len(polynomial) > self.max_terms:
            raise ForestCubeExpansionBudgetExceeded(
                "complete cube polynomial exceeds max_terms=%d" % self.max_terms)
        return polynomial

    def _buffer(self, polynomials):
        count = sum(len(polynomial) for polynomial in polynomials)
        if count > self.max_total_terms:
            raise ForestCubeExpansionBudgetExceeded(
                "complete cube coefficient buffer exceeds max_total_terms=%d" % self.max_total_terms)
        return count

    def _determinant(self, matrix):
        """Fraction-free polynomial determinant, with verified exact divisions."""
        context = self.context
        if not matrix:
            return context.constant(1)
        rows = [list(row) for row in matrix]
        previous = context.constant(1)
        sign = 1
        while len(rows) > 1:
            row_counts = [sum(bool(entry) for entry in row) for row in rows]
            col_counts = [sum(bool(entry) for entry in col) for col in zip(*rows)]
            chosen = min(((len(entry), row_counts[i]*col_counts[j], i, j)
                          for i, row in enumerate(rows)
                          for j, entry in enumerate(row) if entry), default=None)
            if chosen is None:
                return context.constant(0)
            unused, unused2, i, j = chosen
            if i:
                rows[i], rows[0] = rows[0], rows[i]
                sign = -sign
            if j:
                for row in rows:
                    row[j], row[0] = row[0], row[j]
                sign = -sign
            pivot = rows[0][0]
            following = []
            for row in rows[1:]:
                following_row = []
                for entry, top in zip(row[1:], rows[0][1:]):
                    numerator = self._bounded(pivot*entry-row[0]*top)
                    quotient, remainder = divmod(numerator, previous)
                    if remainder:
                        raise ArithmeticError("the augmented current determinant lost exact division")
                    following_row.append(self._bounded(quotient))
                following.append(following_row)
            self._buffer(entry for row in following for entry in row)
            rows, previous = following, pivot
        return sign*rows[0][0]

    def _chart(self, cube):
        context = self.context
        weights = [context.from_dict({exponent: 1}) for exponent in cube.weight_exponents]
        denominator = self._bounded(sum(weights, context.constant(0)))
        numerators = [context.constant(0) for unused in range(2*self.graph.ne)]
        for weight, flow in zip(weights, self.atlas._flows):
            for entry, value in enumerate(flow):
                if value:
                    numerators[entry] += _q(value)*weight
        self._buffer([denominator]+numerators)
        # det d(N/U) = det[[U,dU],[N,dN]] / U**(dimension+1).
        # Subtract each chosen reference-tree flow times the top row first.
        # This exact row operation removes constant terms without changing
        # the full determinant or discarding a source/current direction.
        chosen_flow = self.atlas._flows[self.atlas.trees.index(cube.tree)]
        selected = [self._bounded(numerators[2*edge+a]-_q(chosen_flow[2*edge+a])*denominator)
                    for edge in cube.cotree_order for a in range(2)]
        rows = [[polynomial]+[polynomial.derivative(j) for j in range(cube.dimension)]
                for polynomial in [denominator]+selected]
        numerator = self._bounded(cube.sign*self._determinant(rows))
        power = cube.dimension+1
        while numerator and power:
            quotient, remainder = divmod(numerator, denominator)
            if remainder:
                break
            numerator, power = self._bounded(quotient), power-1
        return denominator, tuple(numerators), numerator, power

    def series(self, powers=None, *, polynomial=None):
        return ExactCubeExpansion(self, self.atlas._polynomial(powers, polynomial))

    def spinor(self, edge, powers=None, *, polynomial=None):
        """Return both complete entries of one additional original edge current."""
        edge = _integer(edge, "original current edge")
        if edge >= self.graph.ne:
            raise ValueError("the current edge lies outside the original graph")
        base = self.atlas._polynomial(powers, polynomial)
        result = []
        for component in range(2):
            inserted = {}
            for exponent, coefficient in base.items():
                exponent = list(exponent)
                exponent[2*edge+component] += 1
                inserted[tuple(exponent)] = coefficient
            # Revalidate the increased complete degree before construction.
            result.append(self.series(polynomial=inserted))
        return tuple(result)


class ExactCubeExpansion:
    """An explicitly specified convergent rational series and exact tail bound."""

    def __init__(self, engine, polynomial):
        self.engine = engine
        self.atlas = engine.atlas
        self.polynomial = dict(polynomial)
        self.components = None
        self.constant = Fraction(0)
        self.last_stats = {}

    def __repr__(self):
        return "ExactCubeExpansion(cubes=%d, dimension=%d, ratio=%s, original_current_terms=%d)" % (
            self.atlas.chart_count, self.atlas.dimension, self.engine.ratio, len(self.polynomial))

    def _prepare(self):
        if self.components is not None:
            return
        if not self.polynomial:
            self.components = ()
            return
        if not self.atlas.dimension:
            self.constant = self.atlas.chart(0).density((), polynomial=self.polynomial)
            self.components = ()
            return
        engine, context = self.engine, self.engine.context
        groups, denominators = {}, {}
        constant = Fraction(0)
        powers_before, powers_after = Counter(), Counter()
        for cube in self.atlas.charts():
            denominator, currents, density, density_power = engine._chart(cube)
            if not density:
                continue
            denominator_key = tuple(sorted(denominator.terms()))
            denominators.setdefault(denominator_key, denominator)
            by_degree = {}
            for exponent, coefficient in self.polynomial.items():
                term = context.constant(_q(coefficient))
                for current, power in zip(currents, exponent):
                    if power:
                        term = engine._bounded(term*engine._bounded(current**power))
                degree = sum(exponent)
                by_degree[degree] = engine._bounded(by_degree.get(degree, context.constant(0))+term)
            for degree, insertion in by_degree.items():
                numerator = engine._bounded(density*insertion)
                if not numerator:
                    continue
                power = density_power+degree
                powers_before[power] += 1
                while power:
                    quotient, remainder = divmod(numerator, denominator)
                    if remainder:
                        break
                    numerator, power = engine._bounded(quotient), power-1
                powers_after[power] += 1
                if power == 0:
                    constant += _moment(numerator)
                    continue
                key = (denominator_key, power)
                groups[key] = engine._bounded(groups.get(key, context.constant(0))+numerator)
                if not groups[key]:
                    del groups[key]
            engine._buffer(list(groups.values())+list(denominators.values()))
        self.components = tuple(_Component(denominators[key], power, numerator, _moment(numerator, absolute=True))
                                for (key, power), numerator in groups.items() if numerator)
        self.constant = constant
        self.last_stats = {"complete_cubes": self.atlas.chart_count,
                           "combined_components": len(self.components),
                           "original_current_terms": len(self.polynomial),
                           "powers_before_insertion_cancellation": dict(sorted(powers_before.items())),
                           "powers_after_insertion_cancellation": dict(sorted(powers_after.items())),
                           "centre": str(engine.centre), "absolute_ratio": str(engine.ratio),
                           "all_original_sources_retained": True,
                           "complete_original_current_polynomial_retained": True}

    def _order(self, order):
        order = _integer(order, "series order")
        if order > self.engine.max_order:
            raise ForestCubeExpansionBudgetExceeded("complete series order exceeds max_order=%d" % self.engine.max_order)
        return order

    def error_bound(self, order):
        """Exact a-priori absolute tail after order, requiring no series terms."""
        order = self._order(order)
        self._prepare()
        centre, ratio = self.engine.centre, self.engine.ratio
        return ratio**(order+1)*sum((component.absolute_moment*
                    sum((Fraction(comb(order+j, j), 1)/centre**j
                         for j in range(component.power)), Fraction(0))
                    for component in self.components), Fraction(0))

    def required_order(self, width):
        """Sufficient order for a rational enclosing interval of at most width."""
        width = _rational(width)
        if width <= 0:
            raise ValueError("requested interval width must be positive")
        self._prepare()
        if 2*self.error_bound(0) <= width:
            return 0
        previous, high = 0, min(1, self.engine.max_order)
        if high == 0:
            raise ForestCubeExpansionBudgetExceeded("requested interval width requires an order above max_order")
        while 2*self.error_bound(high) > width:
            if high == self.engine.max_order:
                raise ForestCubeExpansionBudgetExceeded("requested interval width requires an order above max_order")
            previous, high = high, min(2*high, self.engine.max_order)
        low = previous+1
        while low < high:
            middle = (low+high)//2
            if 2*self.error_bound(middle) <= width:
                high = middle
            else:
                low = middle+1
        return low

    def terms(self, last):
        last = self._order(last)
        self._prepare()
        # Every state is the COMPLETE polynomial P*(c-U)**order.
        states = [component.numerator for component in self.components]
        factors = [self.engine.context.constant(_q(self.engine.centre))-component.denominator
                   for component in self.components]
        partial = Fraction(0)
        running = None
        for order in range(last+1):
            if order:
                states = [self.engine._bounded(polynomial*factor)
                          for polynomial, factor in zip(states, factors)]
            size = self.engine._buffer(states+factors)
            term = sum((Fraction(comb(component.power+order-1, order), 1)
                        *_moment(polynomial)/self.engine.centre**(component.power+order)
                        for component, polynomial in zip(self.components, states)), Fraction(0))
            if order == 0:
                term += self.constant
            partial += term
            error = self.error_bound(order)
            lower, upper = partial-error, partial+error
            if running is not None:
                lower, upper = max(lower, running[0]), min(upper, running[1])
                if lower > upper:
                    raise ArithmeticError("certified geometric-series intervals are disjoint")
            running = lower, upper
            self.last_stats.update({"complete_order": order, "polynomial_terms": size,
                                    "a_priori_error_bound": str(error)})
            yield CubeExpansionTerm(order, term, partial, error, lower, upper,
                                    size, self.atlas.chart_count)

    def interval(self, last):
        row = None
        for row in self.terms(last):
            pass
        return row.lower, row.upper

    def enclose(self, width):
        order = self.required_order(width)
        row = None
        for row in self.terms(order):
            pass
        if row.upper-row.lower > _rational(width):
            raise ArithmeticError("the explicit geometric remainder failed its requested width")
        return row

    def specification(self):
        """Complete JSON-safe coefficients defining the exact infinite answer."""
        self._prepare()
        return {"kind": "complete_geometric_rational_series", "version": 1,
                "original_graph": {"vertices": self.atlas.G.nv,
                                   "ordered_oriented_edges": [list(edge) for edge in self.atlas.G.edges],
                                   "marked": None if self.atlas.marked is None else list(self.atlas.marked)},
                "original_sources": [[str(value) for value in row] for row in self.atlas.sources],
                "original_current_polynomial": [[list(exponent), str(coefficient)]
                                                for exponent, coefficient in sorted(self.polynomial.items())],
                "dimension": self.atlas.dimension, "complete_cubes": self.atlas.chart_count,
                "constant": str(self.constant), "centre": str(self.engine.centre),
                "absolute_ratio": str(self.engine.ratio),
                "series_rule": "constant + sum_components sum_k>=0 binom(power+k-1,k) * cube_moment(numerator*(centre-denominator)^k) / centre^(power+k)",
                "monomial_cube_moment": "cube_moment(x^exponent) = product_j 1/(exponent_j+1)",
                "tail_after_order": "absolute_ratio^(order+1) * sum_components absolute_moment * sum_j=0..power-1 binom(order+j,j)/centre^j",
                "components": [{"power": component.power,
                                "denominator": _polynomial_data(component.denominator),
                                "numerator": _polynomial_data(component.numerator),
                                "absolute_moment": str(component.absolute_moment)}
                               for component in self.components]}
