"""Exact, boundary-aware downward dimension relation for completed K5.

This module constructs the Cayley--Menger determinant by its permutation
definition.  It never drops a term on the ground that a formal falling
factor vanishes: conversion to positive-weight momentum periods requires
all shifted position weights to remain positive.  The complete relation
in J normalization also accepts nonpositive position weights, provided
every shifted completed integral satisfies strict convergence.

The accompanying DERIVATION.md proves the dimension relation by a Gram
change of variables, without integration by parts across singular faces.
"""

from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations
from math import factorial


E5 = tuple(combinations(range(5), 2))
E6 = tuple(combinations(range(6), 2))
I5 = {edge: index for index, edge in enumerate(E5)}
I6 = {edge: index for index, edge in enumerate(E6)}
STARS6 = tuple(tuple(index for index, edge in enumerate(E6) if v in edge)
               for v in range(6))
TRIANGLES6 = tuple(tuple(I6[edge] for edge in combinations(vertices, 2))
                   for vertices in combinations(range(6), 3))
TRIANGLES5 = tuple(tuple(I5[edge] for edge in combinations(vertices, 2))
                   for vertices in combinations(range(5), 3))
CLIQUES5 = tuple(tuple(I5[edge] for edge in combinations(vertices, 2))
                 for vertices in combinations(range(5), 4))


class InvalidPeriod(ValueError):
    pass


class BoundaryContributionRequired(ValueError):
    """The proposed positive-weight scalar conversion would omit a boundary."""


class DivergentLowerTerm(ValueError):
    """A complete Gram term fails a strict lower-dimensional cut condition."""


def _integers(values, count):
    values = tuple(values)
    if len(values) != count or any(isinstance(value, bool) or not isinstance(value, int)
                                  for value in values):
        raise InvalidPeriod("weights must be an exactly sized tuple of integers")
    return values


def completed_convergent(beta, q):
    beta = _integers(beta, 15)
    return (isinstance(q, int) and not isinstance(q, bool) and q >= 2
            and all(sum(beta[index] for index in star) == 2 * q for star in STARS6)
            and max(beta) < q
            and all(sum(beta[index] for index in triangle) < 2 * q
                    for triangle in TRIANGLES6))


def momentum_convergent(alpha, q):
    alpha = _integers(alpha, 10)
    return (isinstance(q, int) and not isinstance(q, bool) and q >= 2
            and sum(alpha) == 6 * q and all(0 < a < q for a in alpha)
            and all(sum(alpha[index] for index in triangle) > q
                    for triangle in TRIANGLES5)
            and all(sum(alpha[index] for index in clique) > 3 * q
                    for clique in CLIQUES5))


def complete_momentum(alpha, q):
    alpha = _integers(alpha, 10)
    if not momentum_convergent(alpha, q):
        raise InvalidPeriod("the K5 scalar period is not strictly convergent")
    beta = [0] * 15
    for edge, a in zip(E5, alpha):
        beta[I6[edge]] = q - a
    for v in range(5):
        beta[I6[(v, 5)]] = 2 * q - sum(beta[index] for index in STARS6[v])
    beta = tuple(beta)
    if not completed_convergent(beta, q):
        raise AssertionError("scalar/completed convergence conventions disagree")
    return beta


@lru_cache(None)
def cayley_menger_terms():
    """All 130 nonzero terms of the symmetric 6x6 zero-diagonal determinant."""
    result = defaultdict(int)
    for permutation in permutations(range(6)):
        if any(permutation[i] == i for i in range(6)):
            continue
        powers = [0] * 15
        for i, j in enumerate(permutation):
            powers[I6[tuple(sorted((i, j)))]] += 1
        inversions = sum(permutation[i] > permutation[j]
                         for i in range(6) for j in range(i + 1, 6))
        result[tuple(powers)] += (-1) ** inversions
    result = {powers: value for powers, value in result.items() if value}
    if len(result) != 130:
        raise AssertionError("wrong complete Cayley--Menger support")
    for powers in result:
        if sum(powers) != 6 or any(sum(powers[i] for i in star) != 2 for star in STARS6):
            raise AssertionError("Cayley--Menger monomial is not balanced")
    return tuple(sorted(result.items()))


@lru_cache(None)
def gram_terms():
    """Q=-det(Cayley--Menger)/16 after the five infinity edges equal one."""
    result = defaultdict(Fraction)
    for powers, coefficient in cayley_menger_terms():
        reduced = tuple(powers[I6[edge]] for edge in E5)
        result[reduced] -= Fraction(coefficient, 16)
    result = {powers: value for powers, value in result.items() if value}
    if any(sum(powers) != 4 or max(powers) > 2 for powers in result):
        raise AssertionError("wrong four-dimensional Gram support")
    return tuple(sorted(result.items()))


def lower_completed(beta, q):
    """Return every exact lower J term, or refuse any divergent term.

    beta follows E6.  No quotient by relabelings and no table lookup occurs.
    The result is a tuple of (shifted_beta, coefficient), all at q-1.
    """
    beta = _integers(beta, 15)
    if q < 3 or not completed_convergent(beta, q):
        raise InvalidPeriod("a convergent completed period at q>=3 is required")
    prefactor = -Fraction(1, 4 * (q - 2) * (2 * q - 3) * (2 * q - 5))
    result = []
    for powers, coefficient in cayley_menger_terms():
        shifted = tuple(b - a for b, a in zip(beta, powers))
        if not completed_convergent(shifted, q - 1):
            raise DivergentLowerTerm((q - 1, shifted, powers))
        result.append((shifted, prefactor * coefficient))
    return tuple(result)


def lower_momentum(alpha, q):
    """Lower to positive-weight K5 I periods with all boundary checks explicit."""
    alpha = _integers(alpha, 10)
    if q < 3 or not momentum_convergent(alpha, q):
        raise InvalidPeriod("a convergent K5 scalar period at q>=3 is required")
    beta = tuple(q - a for a in alpha)
    prefactor = Fraction(4, (q - 1) * (q - 2) * (2 * q - 3) * (2 * q - 5))
    result = []
    for powers, coefficient in gram_terms():
        if any(b <= shift for b, shift in zip(beta, powers)):
            raise BoundaryContributionRequired((q - 1, alpha, powers))
        shifted = tuple(a - 1 + shift for a, shift in zip(alpha, powers))
        if not momentum_convergent(shifted, q - 1):
            raise DivergentLowerTerm((q - 1, shifted, powers))
        factorial_ratio = 1
        for b, shift in zip(beta, powers):
            factorial_ratio *= factorial(b - 1) // factorial(b - shift - 1)
        result.append((shifted, prefactor * coefficient * factorial_ratio))
    return tuple(result)


@lru_cache(None)
def _permutation_indices(vertices):
    edges = E5 if vertices == 5 else E6
    edge_index = I5 if vertices == 5 else I6
    return tuple(tuple(edge_index[tuple(sorted((permutation[u], permutation[v])))]
                       for u, v in edges)
                 for permutation in permutations(range(vertices)))


def canonical(weights, vertices):
    """Independent exhaustive vertex relabeling, used only on bounded controls."""
    weights = tuple(weights)
    return min(tuple(weights[index] for index in permutation)
               for permutation in _permutation_indices(vertices))


def scaled_sum(terms, lookup):
    """Whole exact result; missing input values raise, never become zero."""
    result = defaultdict(Fraction)
    for weights, coefficient in terms:
        for label, value in lookup(weights).items():
            result[label] += coefficient * value
    return {label: value for label, value in sorted(result.items()) if value}
