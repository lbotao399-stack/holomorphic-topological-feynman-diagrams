"""Exact sparse relations for an entire source-supported period numerator.

The graph, all powers, and every rational coefficient remain explicit.  The
existing regulator-free diamond, star--triangle, inverse star--triangle, and
position-space vertex identities supply rows; this module closes only the
orbits reached from the supplied target and eliminates their whole linear
combination.  It never enumerates a complete level or builds a period table.

A stopped computation returns its finite contribution AND all unresolved
periods.  The finite contribution is exposed as a value only when the latter
is empty.  The method is additive to combination.py: its completed-graph
normalization and certificates are not changed.
"""
from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
import heapq
from numbers import Integral
from pathlib import Path
import sys
from types import MappingProxyType

from flint import fmpq

import engine

_FIVELOOP = Path(__file__).resolve().parents[1] / 'fiveloop'
if str(_FIVELOOP) not in sys.path:
    sys.path.insert(0, str(_FIVELOOP))
import msys5
import rel5
import startri


class FamilyTargetBudgetExceeded(RuntimeError):
    """Input or relation construction exceeded an explicit storage bound."""


class IncompleteFamilyTarget(LookupError):
    """The whole exact expression still contains scalar graph periods."""


@dataclass(frozen=True)
class ScalarPeriodTerm:
    q: int
    edges: tuple
    powers: tuple
    coefficient: Fraction


@dataclass(frozen=True)
class FamilyTargetResult:
    finite_coefficients: object
    residuals: tuple
    statistics: object

    @property
    def complete(self):
        return not self.residuals

    @property
    def value(self):
        if self.residuals:
            raise IncompleteFamilyTarget('%d scalar periods remain' % len(self.residuals))
        return self.finite_coefficients


def _integer(value, name, minimum):
    if isinstance(value, bool) or not isinstance(value, Integral) or value < minimum:
        raise ValueError('%s must be an integer >= %d' % (name, minimum))
    return int(value)


def _rational(value):
    if isinstance(value, (float, complex, bool)):
        raise ValueError('coefficients must be exact rational numbers')
    return Fraction(value)


def _labels(values):
    result = {}
    for label, coefficient in values.items():
        if isinstance(label, tuple):
            label = '*'.join(sorted(label)) or '1'
        else:
            label = '*'.join(sorted(label.split('*'))) if label != '1' else '1'
        result[label] = result.get(label, Fraction()) + _rational(coefficient)
    return {label: coefficient for label, coefficient in result.items() if coefficient}


def read_only_lower(edges, powers, q):
    """Existing finite terminals only; no automatic table construction."""
    result = engine.eval_I(edges, powers, q)
    return None if result is None else _labels(result[0])


class _Rows(list):
    def __init__(self, maximum, maximum_coefficients):
        super().__init__()
        self.maximum = maximum
        self.maximum_coefficients = maximum_coefficients
        self.coefficients = 0

    def append(self, value):
        count = len(value[0])
        if len(self) >= self.maximum:
            raise FamilyTargetBudgetExceeded('max_rows exceeded during relation construction')
        if self.coefficients + count > self.maximum_coefficients:
            raise FamilyTargetBudgetExceeded('max_row_coefficients exceeded during relation construction')
        super().append(value)
        self.coefficients += count


class _System(msys5.MultiSystem):
    def __init__(self, edges, q, maximum_variables, maximum_rows, maximum_coefficients):
        super().__init__([edges], q)
        self.maximum_variables = maximum_variables
        self.rows = _Rows(maximum_rows, maximum_coefficients)
        self.lower = read_only_lower
        self.counts.update(X=0, Y=0, finite_seed=0, dimension=0)

    def var(self, family, powers):
        representative = engine.gutil.orbit_rep(tuple(powers), self.fam[family]['autos'])
        key = (family, representative)
        if key not in self.index and len(self.reps) >= self.maximum_variables:
            raise FamilyTargetBudgetExceeded('max_variables exceeded during relation construction')
        return super().var(family, powers)


class FamilyTarget:
    """A reusable target on one simple core; powers use its supplied edge order.

    ``terms`` contains ``(powers, exact coefficient)``.  The scalar convention
    is I_q = integral product(t_e**(powers_e-1)) / U**q with the projective
    measure already used by engine.py and msys5.py.  Repeated inputs add.

    Construction bounds raise explicitly.  Elimination bounds instead return
    a sound incomplete result, because every remaining coefficient is kept.
    ``x_edge_limits`` only restricts which vertex identities are requested;
    no term from an accepted identity is omitted.
    """

    def __init__(self, edges, terms, q, *, max_input_terms=250000,
                 max_variables=100000, max_rows=300000,
                 max_row_coefficients=5000000):
        self.q = _integer(q, 'q', 2)
        max_input_terms = _integer(max_input_terms, 'max_input_terms', 1)
        max_variables = _integer(max_variables, 'max_variables', 1)
        max_rows = _integer(max_rows, 'max_rows', 1)
        max_row_coefficients = _integer(max_row_coefficients, 'max_row_coefficients', 1)
        self.edges = tuple(tuple(edge) for edge in edges)
        if not self.edges or any(len(edge) != 2 or edge[0] == edge[1] for edge in self.edges):
            raise ValueError('each edge must join two distinct vertices')
        if len({frozenset(edge) for edge in self.edges}) != len(self.edges):
            raise ValueError('the input core must be simple; reduce parallel edges before orbit canonicalization')
        if not engine.nx.is_connected(engine.nx.Graph(self.edges)):
            raise ValueError('the graph must be connected')
        self.system = _System(self.edges, self.q, max_variables, max_rows, max_row_coefficients)
        self.target = {}
        self.input_terms = 0
        for self.input_terms, (powers, coefficient) in enumerate(terms, 1):
            if self.input_terms > max_input_terms:
                raise FamilyTargetBudgetExceeded('max_input_terms exceeded')
            coefficient = _rational(coefficient)
            powers = tuple(_integer(power, 'edge power', 1) for power in powers)
            if len(powers) != len(self.edges):
                raise ValueError('one positive integer power is required per edge')
            if not engine.gutil.convergent(self.edges, powers, self.q):
                raise ValueError('every input scalar monomial must converge strictly')
            if coefficient:
                column = self.system.var(0, powers)
                self.target[column] = self.target.get(column, Fraction()) + coefficient
        self.target = {column: value for column, value in self.target.items() if value}
        self.input_columns = tuple(self.target)
        self.local_processed = 0
        self._initial_diamonds = False

    def seed_finite(self, *, search_states=40):
        """Read-only bounded conformal seeds for the input periods only."""
        search_states = _integer(search_states, 'search_states', 0)
        if not search_states:
            return
        for column in self.input_columns:
            powers = self.system.reps[column][1]
            result = engine.period(self.edges, powers, self.q, budget=search_states,
                                   use_m2=True)
            if result is not None:
                self.system.rows.append(({column: Fraction(1)}, _labels(result[0])))
                self.system.counts['finite_seed'] += 1

    def close_local(self, *, x_edge_limits=None):
        """Close reached orbits under complete regulator-free local rows.

        Diamond rows first expose the source-supported family.  Star--triangle
        boundaries and inverse moves are then closed, followed by the allowed
        vertex rows.  A zero-power vertex identity is rejected by rel5.x_items
        because it requires an additional boundary term; it is never truncated.
        """
        limits = {}
        for edge, limit in (x_edge_limits or {}).items():
            edge = _integer(edge, 'x edge index', 0)
            if edge >= len(self.edges):
                raise ValueError('x edge index is out of range')
            limits[edge] = _integer(limit, 'x edge power limit', 1)
        system = self.system
        triangles = rel5.triangles(self.edges)
        if not self._initial_diamonds:
            for column in self.input_columns:
                system.diamond_rows(0, list(system.reps[column][1]))
            self._initial_diamonds = True
        initial = set(self.input_columns)
        first = self.local_processed
        column = first
        while column < len(system.reps):
            powers = system.reps[column][1]
            for vertex in startri.unique_vertices(self.edges, powers, self.q):
                edges, transformed = startri.star_to_triangle(self.edges, list(powers), vertex, self.q)
                system.add_relation([(self.edges, powers, Fraction(1)),
                                     (edges, transformed, Fraction(-1))], 'S')
            for triangle in triangles:
                if sum(powers[index] for index in triangle[3:]) == 2 * self.q:
                    edges, transformed = rel5.tri_to_star(self.edges, list(powers), triangle, self.q)
                    system.add_relation([(self.edges, powers, Fraction(1)),
                                         (edges, transformed, Fraction(-1))], 'Y')
            if column not in initial:
                system.diamond_rows(0, list(powers))
            column += 1
        rules = rel5.x_rules(self.edges)
        processed = len(system.reps)
        for column in range(first, processed):
            powers = system.reps[column][1]
            if any(powers[edge] > limit for edge, limit in limits.items()):
                continue
            for rule in rules:
                items = rel5.x_items(self.edges, list(powers), rule, self.q)
                if items is not None:
                    system.add_relation(items, 'X')
        self.local_processed = processed
        # New orbits reached by X are retained and receive local rows too.
        if self.local_processed < len(system.reps):
            self.close_local(x_edge_limits=limits)

    def add_dimension_rows(self, lower_values, lower_q, *, max_polynomial_terms=100000):
        """Add I_lower_q(beta) = sum [t^m]U**shift I_q(beta+m).

        ``lower_values`` contains ``(positive powers, exact value dictionary)``
        from an independently known lower level on the SAME ordered graph.
        No gamma factor occurs: multiplying numerator and denominator by the
        same power of U preserves the original projective measure exactly.
        """
        lower_q = _integer(lower_q, 'lower_q', 2)
        if lower_q >= self.q:
            raise ValueError('lower_q must be smaller than q')
        maximum = _integer(max_polynomial_terms, 'max_polynomial_terms', 1)
        cotrees = [tuple(int(index in cotree) for index in range(len(self.edges)))
                   for cotree in self.system.fam[0]['cot']]
        polynomial = {tuple(0 for _ in self.edges): 1}
        for _ in range(self.q - lower_q):
            updated = Counter()
            for exponent, coefficient in polynomial.items():
                for cotree in cotrees:
                    updated[tuple(a + b for a, b in zip(exponent, cotree))] += coefficient
                    if len(updated) > maximum:
                        raise FamilyTargetBudgetExceeded('max_polynomial_terms exceeded in dimension shift')
            polynomial = updated
        for powers, value in lower_values:
            powers = tuple(_integer(power, 'lower edge power', 1) for power in powers)
            if len(powers) != len(self.edges) or not engine.gutil.convergent(self.edges, powers, lower_q):
                raise ValueError('every lower-level monomial must converge strictly')
            row = {}
            for exponent, coefficient in polynomial.items():
                shifted = tuple(a + b for a, b in zip(powers, exponent))
                column = self.system.var(0, shifted)
                row[column] = row.get(column, 0) + coefficient
            self.system.rows.append((row, _labels(value)))
            self.system.counts['dimension'] += 1

    def reduce(self, **budgets):
        result = eliminate_target(self.system.rows, self.target,
                                  [powers for _, powers in self.system.reps],
                                  self.edges, self.q, **budgets)
        statistics = dict(result.statistics)
        statistics.update(input_terms=self.input_terms,
                          input_orbits=len(self.input_columns),
                          relation_counts=dict(self.system.counts),
                          local_processed=self.local_processed,
                          read_only_tables=True)
        return FamilyTargetResult(result.finite_coefficients, result.residuals,
                                  MappingProxyType(statistics))


def eliminate_target(rows, target, representatives, edges, q, *,
                     max_operations=100000000, max_active_coefficients=5000000,
                     max_pivots=100000):
    """Sparse exact elimination with deterministic short-row pivoting.

    Rows mean sum(row[j] * period[j]) = constant.  Every pivot is substituted
    into all remaining equations and the entire target.  Incidence controls
    fill; a bound is checked before the next indivisible pivot.  The operation
    budget counts nonconstant period-coefficient products (including the
    target), excluding updates of finite constants.  It is not a wall-clock
    or total process-memory guarantee.
    """
    maximum_operations = _integer(max_operations, 'max_operations', 0)
    maximum_coefficients = _integer(max_active_coefficients, 'max_active_coefficients', 1)
    maximum_pivots = _integer(max_pivots, 'max_pivots', 0)
    incidence = defaultdict(set)
    queue, active, versions = [], [], []
    for identifier, (row, constant) in enumerate(rows):
        row = {int(column): fmpq(str(_rational(value))) for column, value in row.items() if value}
        constant = {label: fmpq(str(value)) for label, value in _labels(constant).items()}
        if any(column < 0 or column >= len(representatives) for column in row):
            raise ValueError('a row refers to an unregistered scalar period')
        active.append((row, constant))
        versions.append(0)
        for column in row:
            incidence[column].add(identifier)
        if row:
            heapq.heappush(queue, (len(row), identifier, 0))
        elif constant:
            raise ArithmeticError('an exact relation has a nonzero constant and no periods')
    count = sum(len(row) for row, _ in active)
    if count > maximum_coefficients:
        raise FamilyTargetBudgetExceeded('initial rows exceed max_active_coefficients')
    peak = count
    target = {int(column): fmpq(str(_rational(value))) for column, value in target.items() if value}
    if any(column < 0 or column >= len(representatives) for column in target):
        raise ValueError('the target refers to an unregistered scalar period')
    finite = {}
    complexity = {column: (-max(powers), -sum(power * power for power in powers),
                           tuple(-power for power in powers))
                  for column, powers in enumerate(representatives)}
    pivots = operations = 0
    stopped = None
    while queue and target:
        if pivots >= maximum_pivots:
            stopped = 'max_pivots'
            break
        _, identifier, version = heapq.heappop(queue)
        if version != versions[identifier] or active[identifier] is None:
            continue
        row, constant = active[identifier]
        if not row:
            continue
        column = min(row, key=lambda index: (len(incidence[index]), complexity[index]))
        denominator = row[column]
        pivot = {index: value / denominator for index, value in row.items() if index != column}
        pivot_constant = {label: value / denominator for label, value in constant.items()}
        affected = sorted(incidence[column] - {identifier})
        cost = len(pivot) * (len(affected) + int(column in target))
        # An upper bound before cancellations keeps each accepted pivot within
        # the requested active-coefficient budget; no partial row is applied.
        fill_bound = count - len(row) + len(affected) * (len(pivot) - 1)
        if operations + cost > maximum_operations:
            stopped = 'max_operations'
            break
        if fill_bound > maximum_coefficients:
            stopped = 'max_active_coefficients'
            break
        active[identifier] = None
        count -= len(row)
        for index in row:
            incidence[index].discard(identifier)
        if column in target:
            factor = target.pop(column)
            for index, value in pivot.items():
                updated = target.get(index, 0) - factor * value
                if updated:
                    target[index] = updated
                else:
                    target.pop(index, None)
            for label, value in pivot_constant.items():
                updated = finite.get(label, 0) + factor * value
                if updated:
                    finite[label] = updated
                else:
                    finite.pop(label, None)
        for other in affected:
            other_row, other_constant = active[other]
            factor = other_row.pop(column)
            incidence[column].discard(other)
            count -= 1
            for index, value in pivot.items():
                existed = index in other_row
                updated = other_row.get(index, 0) - factor * value
                if updated:
                    other_row[index] = updated
                    if not existed:
                        incidence[index].add(other)
                        count += 1
                elif existed:
                    other_row.pop(index)
                    incidence[index].discard(other)
                    count -= 1
            for label, value in pivot_constant.items():
                updated = other_constant.get(label, 0) - factor * value
                if updated:
                    other_constant[label] = updated
                else:
                    other_constant.pop(label, None)
            versions[other] += 1
            if other_row:
                heapq.heappush(queue, (len(other_row), other, versions[other]))
            else:
                if other_constant:
                    raise ArithmeticError('exact period relations have inconsistent finite constants')
                active[other] = None
        pivots += 1
        operations += cost
        peak = max(peak, count)
    residuals = tuple(ScalarPeriodTerm(q, tuple(edges), tuple(representatives[column]),
                                      Fraction(str(coefficient)))
                      for column, coefficient in sorted(target.items()) if coefficient)
    values = MappingProxyType({label: Fraction(str(value)) for label, value in sorted(finite.items())})
    statistics = dict(variables=len(representatives), rows=len(active), pivots=pivots,
                      rational_coefficient_products=operations,
                      peak_active_coefficients=peak, active_coefficients=count,
                      stopped_by=stopped, residual_periods=len(residuals))
    return FamilyTargetResult(values, residuals, MappingProxyType(statistics))
