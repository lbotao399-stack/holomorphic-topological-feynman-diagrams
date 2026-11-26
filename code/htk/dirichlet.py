"""Exact rational cubature on products of standard simplices.

``product_simplex_rule((b_1, ..., b_k), degree)`` integrates every polynomial
of total degree at most ``degree`` on the product of the b_j-dimensional
standard simplices.  Points contain *all* barycentric coordinates of every
block; each coordinate is strictly positive.  Weights are signed Fractions,
and their sum is ``prod(1 / factorial(b_j))`` (Lebesgue, not probability,
normalization).  The empty product has one point and volume one.

The simplex rule is Grundmann--Moeller [GM78], of exact degree 2*s+1:

  sum_{i=0}^s (-1)^i (b+2*s-2*i+1)^(2*s+1)
                / (2^(2*s) i! (b+2*s-i+1)!)
       * sum_{|beta|=s-i} f((2*beta+1)/(b+2*s-2*i+1)).

For several blocks, write Delta_s = Q_s-Q_{s-1}, Q_{-1}=0, and use
``sum_{sum(s_j)<=s} tensor_j Delta_{s_j}``.  Delta_s annihilates block
polynomials of degree below 2*s, so this rule is exact through total degree
2*s+1.  The implementation uses the equivalent combination formula

  sum_{max(0,s-k+1)<=|a|<=s} (-1)^(s-|a|) binom(k-1,s-|a|)
                                    * tensor_j Q_{a_j}.

Zero-dimensional factors are restored after omitting them from this sum.
Duplicate points are merged and zero weights removed, all over Q.  This is
an exact integration rule for the stated polynomial space, not numerical
quadrature with an error tolerance.  No floating-point arithmetic is used.

[GM78] A. Grundmann and H. M. Moeller, Invariant Integration Formulas for
the n-Simplex by Combinatorial Methods, SIAM J. Numer. Anal. 15 (1978),
282--290.  https://doi.org/10.1137/0715019
"""

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import combinations_with_replacement, product
from math import comb, factorial
from numbers import Rational
from operator import index

__all__ = ["CubatureLimitError", "CubatureRule", "simplex_rule",
           "product_simplex_rule", "monomial_integral"]


class CubatureLimitError(ValueError):
    """The requested rule exceeds an explicit enumeration budget."""


def _nonnegative_integer(value, name):
    if isinstance(value, bool):
        raise TypeError("%s must be an integer, not a Boolean" % name)
    try:
        out = index(value)
    except TypeError as error:
        raise TypeError("%s must be an integer" % name) from error
    if out < 0:
        raise ValueError("%s must be nonnegative" % name)
    return out


def _positive_budget(value, name):
    value = _nonnegative_integer(value, name)
    if value == 0:
        raise ValueError("%s must be positive" % name)
    return value


def _compositions(total, parts):
    """Weak compositions without recursion in the number of parts."""
    if parts == 0:
        if total == 0:
            yield ()
        return
    for places in combinations_with_replacement(range(parts), total):
        counts = [0] * parts
        for place in places:
            counts[place] += 1
        yield tuple(counts)


@dataclass(frozen=True)
class CubatureRule:
    """Immutable exact rule; ``nodes`` contains ``(grouped_point, weight)``.

    A grouped point is a tuple of barycentric tuples, one per dimension in
    ``dimensions``. ``degree`` records the requested degree; ``exact_degree``
    records the odd degree guaranteed by this construction.  ``integrate``
    accepts integer or rational scalar values from its callback.
    """

    dimensions: tuple
    degree: int
    exact_degree: int
    nodes: tuple

    def __len__(self):
        return len(self.nodes)

    @property
    def volume(self):
        result = Fraction(1)
        for dimension in self.dimensions:
            result /= factorial(dimension)
        return result

    def integrate(self, function):
        result = Fraction(0)
        for point, weight in self.nodes:
            value = function(point)
            if not isinstance(value, Rational):
                raise TypeError("cubature callback must return an exact rational scalar")
            result += weight * Fraction(value)
        return result


@lru_cache(maxsize=64)
def _simplex_order_nodes(dimension, order):
    """Preflighted primitive GM rule, in ungrouped barycentric coordinates."""
    if dimension == 0:
        return (((Fraction(1),), Fraction(1)),)
    result = {}
    for i in range(order + 1):
        denominator = dimension + 2 * order - 2 * i + 1
        weight = Fraction((-1) ** i * denominator ** (2 * order + 1),
                          2 ** (2 * order) * factorial(i)
                          * factorial(dimension + 2 * order - i + 1))
        for beta in _compositions(order - i, dimension + 1):
            point = tuple(Fraction(2 * value + 1, denominator) for value in beta)
            result[point] = result.get(point, Fraction(0)) + weight
    return tuple(sorted((point, weight) for point, weight in result.items() if weight))


def _raw_term_bound(dimensions, order, max_terms):
    """Upper bound on tensor-node contributions before any duplicate merging."""
    coefficients = [1] + [0] * order
    cap = max_terms + 1
    for dimension in dimensions:
        counts = [min(comb(dimension + level + 1, dimension + 1), cap)
                  for level in range(order + 1)]
        updated = [0] * (order + 1)
        for total in range(order + 1):
            updated[total] = min(sum(coefficients[total - level] * counts[level]
                                     for level in range(total + 1)), cap)
        coefficients = updated
    first = max(0, order - len(dimensions) + 1)
    return min(sum(coefficients[first:]), cap)


@lru_cache(maxsize=16)
def _product_cached(dimensions, degree, max_nodes, max_terms):
    order = degree // 2
    active = tuple(dimension for dimension in dimensions if dimension)
    if not active:
        point = tuple((Fraction(1),) for _ in dimensions)
        return CubatureRule(dimensions, degree, 2 * order + 1,
                            ((point, Fraction(1)),))
    # Apart from the node count, bound the size of an individual stored point.
    if sum(dimensions) + len(dimensions) > max_terms:
        raise CubatureLimitError("barycentric point size exceeds max_terms")
    # Even the interval has (order+1)(order+2)/2 raw GM contributions.
    # This cheap check avoids allocating an unbounded polynomial in order.
    if (order + 1) * (order + 2) // 2 > max_terms:
        raise CubatureLimitError("simplex enumeration exceeds max_terms")
    bound = _raw_term_bound(active, order, max_terms)
    if bound > max_terms:
        raise CubatureLimitError("sparse cubature enumeration exceeds max_terms")
    result = {}
    blocks = len(active)
    for total in range(max(0, order - blocks + 1), order + 1):
        coefficient = (-1) ** (order - total) * comb(blocks - 1, order - total)
        for levels in _compositions(total, blocks):
            factors = [_simplex_order_nodes(dimension, level)
                       for dimension, level in zip(active, levels)]
            for entries in product(*factors):
                point = tuple(entry[0] for entry in entries)
                weight = Fraction(coefficient)
                for entry in entries:
                    weight *= entry[1]
                result[point] = result.get(point, Fraction(0)) + weight
                if len(result) > max_nodes:
                    raise CubatureLimitError("sparse cubature nodes exceed max_nodes")
    nodes = []
    for point, weight in result.items():
        if not weight:
            continue
        active_points = iter(point)
        full_point = tuple(next(active_points) if dimension else (Fraction(1),)
                           for dimension in dimensions)
        nodes.append((full_point, weight))
    return CubatureRule(dimensions, degree, 2 * order + 1, tuple(sorted(nodes)))


def product_simplex_rule(dimensions, degree, *, max_nodes=100000, max_terms=1000000):
    """Exact cubature for a polynomial of bounded *total* degree.

    ``dimensions`` is an iterable of nonnegative integers.  ``max_nodes``
    bounds distinct points during duplicate merging, including points whose
    weight may later cancel.  ``max_terms`` bounds a conservative count of
    primitive tensor-node contributions and the coordinate count of a point.
    Bounds are checked before enumeration where possible; exceeding either
    raises ``CubatureLimitError``.  Cached rules are immutable and keyed by
    the dimensions, degree and both budgets.
    """
    degree = _nonnegative_integer(degree, "degree")
    max_nodes = _positive_budget(max_nodes, "max_nodes")
    max_terms = _positive_budget(max_terms, "max_terms")
    dimensions = tuple(_nonnegative_integer(value, "simplex dimension")
                       for value in dimensions)
    # Count even zero-dimensional factors before constructing their point.
    if len(dimensions) > max_terms:
        raise CubatureLimitError("number of simplex factors exceeds max_terms")
    return _product_cached(dimensions, degree, max_nodes, max_terms)


def simplex_rule(dimension, degree, *, max_nodes=100000, max_terms=1000000):
    """One-block version of ``product_simplex_rule`` (points remain grouped)."""
    return product_simplex_rule((dimension,), degree,
                                max_nodes=max_nodes, max_terms=max_terms)


def monomial_integral(exponents):
    """Dirichlet integral of a barycentric monomial on a product of simplices.

    Each block gives one nonnegative integer exponent per barycentric
    coordinate.  Thus ``((2, 1, 0), (0, 3))`` means a triangle times an
    interval.  Empty blocks are invalid; the empty product has integral one.
    """
    result = Fraction(1)
    for block in exponents:
        powers = tuple(_nonnegative_integer(value, "monomial exponent") for value in block)
        if not powers:
            raise ValueError("a simplex monomial needs at least one barycentric exponent")
        numerator = 1
        for power in powers:
            numerator *= factorial(power)
        result *= Fraction(numerator, factorial(sum(powers) + len(powers) - 1))
    return result
