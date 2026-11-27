"""Direct exact current integration on compact dominant-tree cubes.

For a bridgeless connected tight H=2 graph, choose a spanning tree and
order its cotree edges. Each tree resistance is an independent fraction
of the smallest cotree resistance crossing its fundamental cut. Cotree
resistances are successive ordered fractions, with the largest fixed to
one. These tau*n! cubes partition positive projective resistance space.

Every tree weight divided by the chosen cotree weight becomes a monomial
with nonnegative exponents. Its sum has constant coefficient one. Full
original currents are the corresponding weighted tree flows, so all
currents and their complete Jacobian are smooth on the CLOSED cube.
No globally expanded Schwinger current numerator, period table, or source
restriction is used. The normalized Kirchhoff tree weights remain explicit.

ExactCubeSeries integrates rational midpoint sums at successive dyadic
levels. The telescoping differences are exact rational series terms.
Interval arithmetic on each complete current form gives rigorous rational
remainder bounds. Uniform refinement converges: all denominators are at
least one and the finite arithmetic expressions are continuous on compact
cubes. Consequently every positive requested interval width is attained
after finitely many levels if the explicit resource budgets permit it.

This universal certified route is intentionally separate from finite zeta
reduction. Its chart count can still be large; no efficient all-loop bound
or finite transcendental closure is asserted.
"""

from collections.abc import Mapping
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations, permutations, product
from math import comb, factorial, prod
from operator import index

from . import brute
from .forest_contractions import normalize_forest_sources
from .forest_model import recognize_product_forest
from .graph import Graph
from .linalg import det
from .treeflow import tree_flow


__all__ = ["ForestCubes", "GreedyCube", "ExactCubeSeries", "CubeSeriesTerm",
           "ForestCubeBudgetExceeded"]


class ForestCubeBudgetExceeded(RuntimeError):
    """A complete computation exceeds an explicitly supplied resource budget."""


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


def _rational(value):
    if isinstance(value, (float, complex, bool)):
        raise TypeError("exact rational values are required")
    return Fraction(value)


@dataclass(frozen=True)
class _Interval:
    lower: Fraction
    upper: Fraction

    def __post_init__(self):
        if self.lower > self.upper:
            raise ValueError("interval endpoints are reversed")

    @staticmethod
    def exact(value):
        value = Fraction(value)
        return _Interval(value, value)

    def __add__(self, other):
        if not isinstance(other, _Interval):
            other = self.exact(other)
        return _Interval(self.lower+other.lower, self.upper+other.upper)

    __radd__ = __add__

    def __neg__(self):
        return _Interval(-self.upper, -self.lower)

    def __sub__(self, other):
        return self + (-other if isinstance(other, _Interval) else -Fraction(other))

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        if not isinstance(other, _Interval):
            other = self.exact(other)
        if self.lower >= 0 and other.lower >= 0:
            return _Interval(self.lower*other.lower, self.upper*other.upper)
        values = (self.lower*other.lower, self.lower*other.upper,
                  self.upper*other.lower, self.upper*other.upper)
        return _Interval(min(values), max(values))

    __rmul__ = __mul__

    def __truediv__(self, other):
        if not isinstance(other, _Interval):
            other = self.exact(other)
        if other.lower <= 0 <= other.upper:
            raise ZeroDivisionError("an interval denominator contains zero")
        return self * _Interval(1/other.upper, 1/other.lower)

    def __pow__(self, power):
        if power == 0:
            return self.exact(1)
        if power & 1:
            return _Interval(self.lower**power, self.upper**power)
        if self.lower >= 0:
            return _Interval(self.lower**power, self.upper**power)
        if self.upper <= 0:
            return _Interval(self.upper**power, self.lower**power)
        return _Interval(Fraction(0), max((-self.lower)**power, self.upper**power))

    def intersect(self, other):
        return _Interval(max(self.lower, other.lower), min(self.upper, other.upper))


def _det_interval(matrix, centre, centre_value, max_states):
    """A convergent enclosure, including matrices singular at the centre.

    Multilinearity in rows gives |det(C+R)-det(C)| <=
    product(row_l1(C)+row_l1(|R|))-product(row_l1(C)).
    This bound tends to zero with the entry radii, without dividing by a
    possibly vanishing determinant or pivot. Small matrices also use the
    exact subset determinant expression and intersect both valid ranges.
    """
    size = len(matrix)
    if not size:
        return _Interval.exact(1)
    base, enlarged, absolute = [], [], []
    for intervals, values in zip(matrix, centre):
        a = sum(abs(value) for value in values)
        radius = sum(max(abs(interval.lower-value), abs(interval.upper-value))
                     for interval, value in zip(intervals, values))
        base.append(a)
        enlarged.append(a+radius)
        absolute.append(sum(max(abs(entry.lower), abs(entry.upper)) for entry in intervals))
    error = prod(enlarged)-prod(base)
    bound = prod(absolute)
    answer = _Interval(centre_value-error, centre_value+error).intersect(_Interval(-bound, bound))
    if size >= max_states.bit_length():
        return answer
    states = {0: _Interval.exact(1)}
    for row in matrix:
        following = {}
        for used, value in states.items():
            for column, entry in enumerate(row):
                if used & (1 << column):
                    continue
                target = used | (1 << column)
                term = value * entry
                if (used >> (column+1)).bit_count() & 1:
                    term = -term
                following[target] = following.get(target, _Interval.exact(0)) + term
        states = following
    return answer.intersect(states[(1 << size)-1])


def _tree(graph, chosen):
    parent = list(range(graph.nv))
    def find(vertex):
        while parent[vertex] != vertex:
            vertex = parent[vertex]
        return vertex
    for edge in chosen:
        a, b = map(find, graph.edges[edge])
        if a == b:
            return False
        parent[b] = a
    return True


def _permutation_at(values, rank):
    values = list(values)
    result = []
    while values:
        size = factorial(len(values)-1)
        digit, rank = divmod(rank, size)
        result.append(values.pop(digit))
    return tuple(result)


def _word_value(currents, polynomial):
    value = 0
    for powers, coefficient in polynomial.items():
        term = coefficient
        for current, power in zip(currents, powers):
            if power:
                term *= current**power
        value += term
    return value


@dataclass(frozen=True)
class CubeSeriesTerm:
    order: int
    term: Fraction
    partial_sum: Fraction
    error_bound: Fraction
    lower: Fraction
    upper: Fraction
    evaluations: int


class ForestCubes:
    """A complete source-preserving dominant-tree cube atlas.

    ``marked`` optionally validates the original product-forest topology.
    Without it, connected tight bridgeless H=2 graphs are accepted as well.
    All resistance edge signs and ordering remain those of the input graph.
    ``max_charts`` bounds the complete atlas, even though charts are lazy.
    ``max_evaluations`` bounds cumulative midpoint cells of a series request.
    ``max_det_states`` controls an optional sharper interval determinant;
    larger matrices use a convergent row-perturbation bound instead.
    """

    def __init__(self, G, sources, marked=None, *, check=1,
                 max_tree_candidates=1_000_000, max_trees=100_000,
                 max_charts=1_000_000, max_chart_entries=1_000_000,
                 max_evaluations=1_000_000, max_det_states=4096,
                 max_input_degree=10_000):
        self.G = self.graph = Graph(G.nv, G.edges, G.name)
        if not self.G.is_connected() or self.G.ne != 2*self.G.loops+1:
            raise ValueError("a connected tight H=2 graph is required")
        self.marked = None if marked is None else tuple(marked)
        if self.marked is not None:
            recognize_product_forest(self.G, self.marked)
        self.sources = normalize_forest_sources(self.G, sources)
        self.check = _integer(check, "check count")
        self.max_tree_candidates = _integer(max_tree_candidates, "max_tree_candidates", 1)
        self.max_trees = _integer(max_trees, "max_trees", 1)
        self.max_charts = _integer(max_charts, "max_charts", 1)
        self.max_chart_entries = _integer(max_chart_entries, "max_chart_entries", 1)
        self.max_evaluations = _integer(max_evaluations, "max_evaluations", 1)
        self.max_det_states = _integer(max_det_states, "max_det_states", 1)
        self.max_input_degree = _integer(max_input_degree, "max_input_degree")
        self.dimension = self.G.ne-1
        self.loops = self.G.loops
        candidates = comb(self.G.ne, self.G.nv-1)
        if candidates > self.max_tree_candidates:
            raise ForestCubeBudgetExceeded("complete tree scan exceeds max_tree_candidates")
        if self.G.ne*self.dimension > self.max_chart_entries:
            raise ForestCubeBudgetExceeded("complete cube time matrix exceeds max_chart_entries")
        self.tree_count = int(brute.symanzik_U(self.G, [1]*self.G.ne))
        if self.tree_count > self.max_trees:
            raise ForestCubeBudgetExceeded("complete tree enumeration exceeds max_trees")
        self.orders_per_tree = factorial(self.loops)
        self.chart_count = self.tree_count*self.orders_per_tree
        if self.chart_count > self.max_charts:
            raise ForestCubeBudgetExceeded("complete dominant-tree atlas exceeds max_charts")
        if self.tree_count*max(2*self.G.ne, self.dimension) > self.max_chart_entries:
            raise ForestCubeBudgetExceeded("complete tree-flow or chart weights exceed max_chart_entries")
        trees, flows = [], []
        for selected in combinations(range(self.G.ne), self.G.nv-1):
            if _tree(self.G, selected):
                trees.append(selected)
                flows.append(tuple(value for row in tree_flow(self.G, self.sources, selected) for value in row))
        if len(trees) != self.tree_count:
            raise ArithmeticError("independent tree count and exact enumeration differ")
        if self.loops and any(all(edge in tree for tree in trees) for edge in range(self.G.ne)):
            raise ValueError("positive-loop cube atlas requires a bridgeless graph; every valid product forest is bridgeless")
        self.trees = tuple(trees)
        self._flows = tuple(flows)
        self._checked = set()
        self.stats = {"loops": self.loops, "dimension": self.dimension,
                      "spanning_trees": self.tree_count, "complete_cubes": self.chart_count,
                      "ordinary_Hepp_orders": factorial(self.G.ne),
                      "global_current_numerator_constructed": False,
                      "all_original_sources_retained": True}

    def chart(self, position):
        position = _integer(position, "chart position")
        if position >= self.chart_count:
            raise IndexError("chart position lies outside the complete atlas")
        tree_position, order_position = divmod(position, self.orders_per_tree)
        tree = self.trees[tree_position]
        cotree = tuple(edge for edge in range(self.G.ne) if edge not in tree)
        order = _permutation_at(cotree, order_position)
        cube = GreedyCube(self, position, tree_position, order)
        if position < self.check and position not in self._checked:
            cube.check_original()
            self._checked.add(position)
        return cube

    def charts(self):
        for position in range(self.chart_count):
            yield self.chart(position)

    def locate(self, times):
        """Locate an original positive resistance point in its greedy cube.

        Exact ties use original edge order; overlapping cube boundaries have
        zero integration measure. The returned fractions reproduce every
        original resistance up to the common positive projective scale.
        """
        times = tuple(_rational(value) for value in times)
        if len(times) != self.G.ne or any(value <= 0 for value in times):
            raise ValueError("one strictly positive exact resistance is required per original edge")
        parent = list(range(self.G.nv))
        def find(vertex):
            while parent[vertex] != vertex:
                vertex = parent[vertex]
            return vertex
        tree = []
        for edge in sorted(range(self.G.ne), key=lambda edge: (times[edge], edge)):
            a, b = map(find, self.G.edges[edge])
            if a != b:
                parent[b] = a
                tree.append(edge)
        tree = tuple(sorted(tree))
        natural = [edge for edge in range(self.G.ne) if edge not in tree]
        ordered = sorted(natural, key=lambda edge: (times[edge], edge))
        rank = 0
        remaining = list(natural)
        for edge in ordered:
            digit = remaining.index(edge)
            rank += digit*factorial(len(remaining)-1)
            remaining.pop(digit)
        cube = self.chart(self.trees.index(tree)*self.orders_per_tree+rank)
        return cube, cube.coordinates(times)

    def _polynomial(self, powers=None, polynomial=None):
        if powers is not None and polynomial is not None:
            raise ValueError("supply an original current word or polynomial, not both")
        if polynomial is None:
            if powers is not None and not isinstance(powers, Mapping):
                raise TypeError("current powers must be a mapping")
            exponent = [0]*(2*self.G.ne)
            for location, power in (powers or {}).items():
                if not isinstance(location, tuple) or len(location) != 2:
                    raise ValueError("current locations are (original edge, spinor entry)")
                edge = _integer(location[0], "current edge")
                component = _integer(location[1], "spinor entry")
                if edge >= self.G.ne or component >= 2:
                    raise ValueError("a current location lies outside the original graph")
                exponent[2*edge+component] = _integer(power, "current power")
            polynomial = {tuple(exponent): Fraction(1)}
        if not isinstance(polynomial, Mapping):
            raise TypeError("the full current polynomial must be a mapping")
        if len(polynomial)*max(1, 2*self.G.ne) > self.max_chart_entries:
            raise ForestCubeBudgetExceeded("the input current polynomial exceeds max_chart_entries")
        result = {}
        for exponent, coefficient in polynomial.items():
            if len(exponent) != 2*self.G.ne:
                raise ValueError("one exponent is required for every original current entry")
            exponent = tuple(_integer(power, "current exponent") for power in exponent)
            if sum(exponent) > self.max_input_degree:
                raise ForestCubeBudgetExceeded("the complete current word exceeds max_input_degree")
            coefficient = _rational(coefficient)
            if coefficient:
                result[exponent] = result.get(exponent, Fraction(0))+coefficient
        return {exponent: coefficient for exponent, coefficient in result.items() if coefficient}

    def series(self, powers=None, *, polynomial=None):
        return ExactCubeSeries(self, self._polynomial(powers, polynomial))


class GreedyCube:
    """One oriented compact cube; every map retains all original currents."""

    def __init__(self, atlas, position, tree_position, cotree_order):
        self.atlas = atlas
        self.position = position
        self.tree = atlas.trees[tree_position]
        self.cotree_order = cotree_order
        self.dimension = atlas.dimension
        edges, loops = atlas.G.ne, atlas.loops
        if not loops:
            self.time_exponents = ((),)
            self.weight_exponents = ((),)
            self.tree_reference_cotrees = ()
            self.orientation = 1
            self.sign = 1
            return
        cycles = atlas.G.fundamental_cycles(self.tree)
        natural_cotree = tuple(edge for edge in range(edges) if edge not in self.tree)
        cycle_of = {edge: cycles[j] for j, edge in enumerate(natural_cotree)}
        exponents = [[0]*self.dimension for _ in range(edges)]
        offset = len(self.tree)
        for order, edge in enumerate(cotree_order):
            for ratio in range(order, loops-1):
                exponents[edge][offset+ratio] = 1
        references = []
        for coordinate, edge in enumerate(self.tree):
            parent_cotree = next((candidate for candidate in cotree_order if cycle_of[candidate][edge]), None)
            if parent_cotree is None:
                raise ArithmeticError("a tree edge has no fundamental-cut cotree")
            exponents[edge] = list(exponents[parent_cotree])
            exponents[edge][coordinate] = 1
            references.append(parent_cotree)
        chosen = tuple(sum(exponents[edge][j] for edge in cotree_order) for j in range(self.dimension))
        weights = []
        for tree in atlas.trees:
            exponent = tuple(sum(exponents[edge][j] for edge in range(edges) if edge not in tree)-chosen[j]
                             for j in range(self.dimension))
            if any(power < 0 for power in exponent):
                raise ArithmeticError("a normalized tree weight is not a nonnegative monomial")
            weights.append(exponent)
        if sum(not any(exponent) for exponent in weights) != 1:
            raise ArithmeticError("the normalized Kirchhoff denominator lacks its unique constant one")
        self.time_exponents = tuple(map(tuple, exponents))
        self.weight_exponents = tuple(weights)
        self.tree_reference_cotrees = tuple(references)
        determinant = det([[1]+list(exponent) for exponent in exponents])
        if abs(determinant) != 1:
            raise ArithmeticError("the dominant-tree cube is not an invertible projective monomial chart")
        self.orientation = int(determinant)
        self.sign = self.orientation * (-1)**(loops*(loops-1)//2)

    def _point(self, point):
        point = tuple(_rational(value) for value in point)
        if len(point) != self.dimension or any(not 0 <= value <= 1 for value in point):
            raise ValueError("one exact closed-unit-interval coordinate is required per cube direction")
        return point

    def time_map(self, point):
        point = self._point(point)
        return tuple(prod(value**power for value, power in zip(point, exponent))
                     for exponent in self.time_exponents)

    def coordinates(self, times):
        """The inverse projective map on this complete dominant-tree cone."""
        times = tuple(_rational(value) for value in times)
        if len(times) != self.atlas.G.ne or any(value <= 0 for value in times):
            raise ValueError("one positive exact resistance is required per original edge")
        point = tuple(times[edge]/times[parent] for edge, parent in zip(self.tree, self.tree_reference_cotrees))
        point += tuple(times[a]/times[b] for a, b in zip(self.cotree_order, self.cotree_order[1:]))
        if any(value > 1 for value in point):
            raise ValueError("the resistance point lies outside this dominant-tree cube")
        return point

    def _currents(self, coordinates, derivative=True):
        zero = coordinates[0]*0 if coordinates else Fraction(0)
        denominator = zero
        numerator = [zero]*(2*self.atlas.G.ne)
        derivative_denominator = [zero]*self.dimension
        derivative_numerator = [[zero]*self.dimension for _ in numerator]
        for exponent, flow in zip(self.weight_exponents, self.atlas._flows):
            factors = [value**power for value, power in zip(coordinates, exponent)]
            weight = prod(factors)
            denominator += weight
            for entry, value in enumerate(flow):
                if value:
                    numerator[entry] += value*weight
            if derivative:
                for direction, power in enumerate(exponent):
                    if not power:
                        continue
                    changed = power*coordinates[direction]**(power-1)
                    changed *= prod(value for j, value in enumerate(factors) if j != direction)
                    derivative_denominator[direction] += changed
                    for entry, value in enumerate(flow):
                        if value:
                            derivative_numerator[entry][direction] += value*changed
        values = tuple(value/denominator for value in numerator)
        if not derivative:
            return denominator, values, ()
        jacobian = tuple(tuple((entry-value*change)/denominator for entry, change in zip(row, derivative_denominator))
                         for row, value in zip(derivative_numerator, values))
        return denominator, values, jacobian

    def currents(self, point):
        denominator, values, unused = self._currents(self._point(point), derivative=False)
        if denominator < 1:
            raise ArithmeticError("the normalized Kirchhoff denominator is below one")
        return tuple(tuple(values[2*edge:2*edge+2]) for edge in range(self.atlas.G.ne))

    def current_jacobian(self, point):
        """All original current derivatives, shaped (edges,2,cube dimension)."""
        unused, values, jacobian = self._currents(self._point(point))
        return tuple(tuple(jacobian[2*edge+component] for component in range(2))
                     for edge in range(self.atlas.G.ne))

    def density(self, point, powers=None, *, polynomial=None):
        polynomial = self.atlas._polynomial(powers, polynomial)
        unused, values, jacobian = self._currents(self._point(point))
        matrix = tuple(jacobian[2*edge+component] for edge in self.cotree_order for component in range(2))
        return self.sign*det(matrix)*_word_value(values, polynomial)

    def check_original(self):
        point = tuple(Fraction(j+2, 2*j+5) for j in range(self.dimension))
        times = self.time_map(point)
        currents = tuple(map(tuple, brute.currents(self.atlas.G, times, self.atlas.sources)))
        if self.currents(point) != currents:
            raise ArithmeticError("a complete cube current differs from the original Kirchhoff solve")
        # Projective monomial Jacobian: |det[time;dtime]| = prod(time)/prod(x).
        measure = Fraction(prod(times), prod(point))
        expected = brute.density(self.atlas.G, times, self.atlas.sources).rational()*measure
        if self.density(point) != expected:
            raise ArithmeticError("the complete signed cube density differs from the original graph")
        return {"point": point, "times": times, "density": expected,
                "all_original_currents_checked": True}

    def _cell(self, bounds, polynomial):
        centre = tuple((lower+upper)/2 for lower, upper in bounds)
        unused, centre_currents, centre_jacobian = self._currents(centre)
        centre_matrix = tuple(centre_jacobian[2*edge+component] for edge in self.cotree_order for component in range(2))
        centre_determinant = det(centre_matrix)
        value = self.sign*centre_determinant*_word_value(centre_currents, polynomial)
        coordinates = tuple(_Interval(lower, upper) for lower, upper in bounds)
        denominator, currents, jacobian = self._currents(coordinates)
        if isinstance(denominator, _Interval) and denominator.lower < 1:
            raise ArithmeticError("a complete cube interval lost the constant-one denominator")
        matrix = tuple(jacobian[2*edge+component] for edge in self.cotree_order for component in range(2))
        determinant = _det_interval(matrix, centre_matrix, centre_determinant, self.atlas.max_det_states)
        insertion = _word_value(currents, polynomial)
        enclosure = self.sign*determinant*insertion
        if not isinstance(enclosure, _Interval):
            enclosure = _Interval.exact(enclosure)
        if not enclosure.lower <= value <= enclosure.upper:
            raise ArithmeticError("the exact midpoint lies outside its complete current-form interval")
        return value, enclosure


class ExactCubeSeries:
    """Exact rational telescoping series with certified finite-level intervals.

    A partial sum is never returned as the exact integral. Every requested
    level includes every source, chart, cell, and original current term.
    ``enclose(width)`` terminates for every positive rational width when its
    explicit finite evaluation budget is sufficient.
    """

    def __init__(self, atlas, polynomial):
        self.atlas = atlas
        self.polynomial = dict(polynomial)
        self.last_stats = {}

    def __repr__(self):
        return "ExactCubeSeries(cubes=%d, dimension=%d, original_current_terms=%d)" % (
            self.atlas.chart_count, self.atlas.dimension, len(self.polynomial))

    def _cost(self, level):
        power = self.atlas.dimension*level
        if power >= self.atlas.max_evaluations.bit_length():
            raise ForestCubeBudgetExceeded("the complete dyadic level exceeds max_evaluations")
        cost = self.atlas.chart_count*(1 << power)
        if cost > self.atlas.max_evaluations:
            raise ForestCubeBudgetExceeded("the complete dyadic level exceeds max_evaluations")
        return cost

    def _level(self, level):
        cost = self._cost(level)
        subdivisions = 1 << level
        volume = Fraction(1, 1 << (self.atlas.dimension*level))
        partial = Fraction(0)
        enclosure = _Interval.exact(0)
        for cube in self.atlas.charts():
            for cell in product(range(subdivisions), repeat=self.atlas.dimension):
                bounds = tuple((Fraction(position, subdivisions), Fraction(position+1, subdivisions))
                               for position in cell)
                value, interval = cube._cell(bounds, self.polynomial)
                partial += volume*value
                enclosure += volume*interval
        return partial, enclosure, cost

    def terms(self, last):
        last = _integer(last, "dyadic series level")
        if last+1 > self.atlas.max_evaluations:
            raise ForestCubeBudgetExceeded("the complete series output exceeds max_evaluations")
        if not self.polynomial:
            for level in range(last+1):
                yield CubeSeriesTerm(level, Fraction(0), Fraction(0), Fraction(0), Fraction(0), Fraction(0), 0)
            return
        if not self.atlas.dimension:
            value = self.atlas.chart(0).density((), polynomial=self.polynomial)
            for level in range(last+1):
                yield CubeSeriesTerm(level, value if level == 0 else Fraction(0), value,
                                     Fraction(0), value, value, 1)
            return
        # Preflight the complete request before allocating any chart/cell arrays.
        total = sum(self._cost(level) for level in range(last+1))
        if total > self.atlas.max_evaluations:
            raise ForestCubeBudgetExceeded("the complete requested series exceeds cumulative max_evaluations")
        previous, running, used = Fraction(0), None, 0
        for level in range(last+1):
            partial, enclosure, count = self._level(level)
            running = enclosure if running is None else running.intersect(enclosure)
            used += count
            error = max(abs(partial-running.lower), abs(running.upper-partial))
            self.last_stats = {"complete_level": level, "evaluations": used,
                               "cubes_per_level": self.atlas.chart_count,
                               "all_original_sources_retained": True,
                               "global_current_numerator_constructed": False}
            yield CubeSeriesTerm(level, partial-previous, partial, error, running.lower, running.upper, used)
            previous = partial

    def interval(self, last):
        row = None
        for row in self.terms(last):
            pass
        return row.lower, row.upper

    def enclose(self, width):
        width = _rational(width)
        if width <= 0:
            raise ValueError("requested interval width must be positive")
        if not self.polynomial:
            return CubeSeriesTerm(0, Fraction(0), Fraction(0), Fraction(0), Fraction(0), Fraction(0), 0)
        previous, running, used, level = Fraction(0), None, 0, 0
        while True:
            count = self._cost(level)
            if used+count > self.atlas.max_evaluations:
                raise ForestCubeBudgetExceeded("certified enclosure requires more than max_evaluations; no partial value is returned")
            partial, enclosure, count = self._level(level)
            running = enclosure if running is None else running.intersect(enclosure)
            used += count
            error = max(abs(partial-running.lower), abs(running.upper-partial))
            row = CubeSeriesTerm(level, partial-previous, partial, error, running.lower, running.upper, used)
            self.last_stats = {"complete_level": level, "evaluations": used,
                               "cubes_per_level": self.atlas.chart_count,
                               "all_original_sources_retained": True,
                               "global_current_numerator_constructed": False}
            if running.upper-running.lower <= width:
                return row
            previous, level = partial, level+1
