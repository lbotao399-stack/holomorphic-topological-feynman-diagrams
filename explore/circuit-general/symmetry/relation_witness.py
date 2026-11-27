"""Sparse exact target elimination with an independently summable row witness.

Rows mean sum(row[j] * period[j]) = constant.  Besides retaining every
unresolved target coefficient, the solver returns rational weights on the
ORIGINAL rows.  Summing those rows verifies the result without repeating
elimination.  Reverse propagation through the elimination history avoids
expanding a separate dense proof vector for every active row.

The pivot policy follows family_target.py.  This additive implementation
does not change that module or any earlier certificate.
"""
from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction
import heapq
from types import MappingProxyType

from flint import fmpq

from family_target import (FamilyTargetBudgetExceeded, FamilyTargetResult,
                           ScalarPeriodTerm, _integer, _labels, _rational)


@dataclass(frozen=True)
class WitnessResult(FamilyTargetResult):
    row_weights: object


def eliminate_with_witness(rows, target, representatives, edges, q, *,
                           max_operations=10000000,
                           max_active_coefficients=3000000,
                           max_pivots=30000, max_history_entries=2000000):
    """Return a complete finite-plus-residual expression and its exact proof.

    Bounds are checked before an indivisible pivot.  ``max_operations``
    counts nonconstant period-coefficient products, as in family_target;
    ``max_history_entries`` additionally bounds recorded row dependencies.
    These are arithmetic/storage counters, not a wall-clock guarantee.
    """
    maximum_operations = _integer(max_operations, 'max_operations', 0)
    maximum_coefficients = _integer(max_active_coefficients, 'max_active_coefficients', 1)
    maximum_pivots = _integer(max_pivots, 'max_pivots', 0)
    maximum_history = _integer(max_history_entries, 'max_history_entries', 0)
    active, history, versions, queue = [], [], [], []
    incidence = defaultdict(set)
    for identifier, (row, constant) in enumerate(rows):
        row = {int(column): fmpq(str(_rational(value))) for column, value in row.items() if value}
        constant = {label: fmpq(str(value)) for label, value in _labels(constant).items()}
        if any(column < 0 or column >= len(representatives) for column in row):
            raise ValueError('a relation refers to an unregistered period')
        active.append((row, constant))
        history.append([])
        versions.append(0)
        for column in row:
            incidence[column].add(identifier)
        if row:
            heapq.heappush(queue, (len(row), identifier, 0))
        elif constant:
            raise ArithmeticError('a relation has a nonzero constant and no periods')
    count = sum(len(row) for row, _ in active)
    if count > maximum_coefficients:
        raise FamilyTargetBudgetExceeded('initial rows exceed max_active_coefficients')
    original_target = {int(column): _rational(value) for column, value in target.items() if value}
    if any(column < 0 or column >= len(representatives) for column in original_target):
        raise ValueError('the target refers to an unregistered period')
    target = {column: fmpq(str(value)) for column, value in original_target.items()}
    finite, target_weights, pivots = {}, {}, []
    complexity = {column: (-max(powers), -sum(power * power for power in powers),
                           tuple(-power for power in powers))
                  for column, powers in enumerate(representatives)}
    operations = history_entries = 0
    peak = count
    stopped = None
    while queue and target:
        if len(pivots) >= maximum_pivots:
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
        fill_bound = count - len(row) + len(affected) * (len(pivot) - 1)
        if operations + cost > maximum_operations:
            stopped = 'max_operations'
            break
        if fill_bound > maximum_coefficients:
            stopped = 'max_active_coefficients'
            break
        if history_entries + len(affected) > maximum_history:
            stopped = 'max_history_entries'
            break
        step = len(pivots)
        pivots.append((identifier, denominator, history[identifier]))
        history[identifier] = None
        active[identifier] = None
        count -= len(row)
        for index in row:
            incidence[index].discard(identifier)
        if column in target:
            factor = target.pop(column)
            target_weights[step] = factor
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
            history[other].append((step, factor))
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
                    raise ArithmeticError('exact relations have inconsistent finite constants')
                active[other] = None
                history[other] = None
        history_entries += len(affected)
        operations += cost
        peak = max(peak, count)

    # Each normalized pivot is its original row minus earlier normalized
    # pivots, divided by its pivot coefficient.  Reverse substitution gives
    # weights on original rows without a large proof vector at every step.
    proof = {}
    for step in reversed(range(len(pivots))):
        weight = target_weights.get(step, 0)
        if not weight:
            continue
        identifier, denominator, dependencies = pivots[step]
        weight /= denominator
        proof[identifier] = proof.get(identifier, 0) + weight
        for earlier, factor in dependencies:
            updated = target_weights.get(earlier, 0) - weight * factor
            if updated:
                target_weights[earlier] = updated
            else:
                target_weights.pop(earlier, None)
    proof = {identifier: Fraction(str(weight)) for identifier, weight in sorted(proof.items()) if weight}
    residuals = tuple(ScalarPeriodTerm(q, tuple(edges), tuple(representatives[column]),
                                      Fraction(str(coefficient)))
                      for column, coefficient in sorted(target.items()) if coefficient)
    finite = {label: Fraction(str(value)) for label, value in sorted(finite.items()) if value}
    statistics = dict(variables=len(representatives), rows=len(rows), pivots=len(pivots),
                      rational_coefficient_products=operations,
                      peak_active_coefficients=peak, active_coefficients=count,
                      history_entries=history_entries, witness_rows=len(proof),
                      stopped_by=stopped, residual_periods=len(residuals))
    result = WitnessResult(MappingProxyType(finite), residuals, MappingProxyType(statistics),
                           MappingProxyType(proof))
    # Independent summation checks every original target coefficient,
    # including the whole residual when a budget stops elimination.
    verify_row_witness(rows, original_target,
                       {column: Fraction(str(value)) for column, value in target.items()},
                       finite, proof)
    return result


def verify_row_witness(rows, target, residual, finite, weights):
    """Check target = finite + residual by direct exact original-row sums."""
    combined, constants = {}, {}
    for identifier, weight in weights.items():
        if not isinstance(identifier, int) or identifier < 0 or identifier >= len(rows):
            raise ValueError('a witness refers to an absent original row')
        weight = _rational(weight)
        row, constant = rows[identifier]
        for column, value in row.items():
            combined[column] = combined.get(column, Fraction()) + weight * _rational(value)
        for label, value in _labels(constant).items():
            constants[label] = constants.get(label, Fraction()) + weight * value
    expected = {column: _rational(value) for column, value in target.items() if value}
    for column, value in residual.items():
        expected[column] = expected.get(column, Fraction()) - _rational(value)
    combined = {column: value for column, value in combined.items() if value}
    expected = {column: value for column, value in expected.items() if value}
    constants = {label: value for label, value in constants.items() if value}
    if combined != expected or constants != _labels(finite):
        raise ArithmeticError('the original-row witness does not reproduce the complete target identity')
    return True
