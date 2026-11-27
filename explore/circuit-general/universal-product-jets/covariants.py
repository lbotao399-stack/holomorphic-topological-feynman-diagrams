"""Exact complete symmetric spinor tensors of arbitrary finite rank.

Each returned term is (source multiplicities, three bracket powers,
coefficient). The three sources are gamma4, c, gamma5; the brackets are
<gamma4,c>, <gamma5,c>, <gamma4,gamma5>. Source multiplicities total the
tensor rank and bracket powers total three. Products of source spinors
mean their average over all permutations of the free lower spinor indices.

The calculation uses every symmetric tensor entry and every source-grid
coefficient. The two formal contraction variables below are only a way
to store the complete tensor; no input/output spinor is fixed or removed.
"""
from collections import defaultdict
from fractions import Fraction
from math import comb
from operator import index

from flint import fmpq, fmpq_mpoly_ctx


def integer(value, name, minimum=0):
    if isinstance(value, bool):
        raise TypeError(name + " must be an integer")
    try:
        value = index(value)
    except TypeError as error:
        raise TypeError(name + " must be an integer") from error
    if value < minimum:
        raise ValueError(name + " must be at least " + str(minimum))
    return value


def _fraction(value):
    return Fraction(int(value.p), int(value.q))


def _q(value):
    value = Fraction(value)
    return fmpq(value.numerator, value.denominator)


def restore(polynomials, order, *, max_terms=1000000):
    """Restore all free indices from the entire normalized source polynomial.

    Entry k has exactly k lower-minus indices and order-k lower-plus
    indices. Its binary coefficient therefore includes binomial(order,k).
    A spectator-homogeneous part of degree h>3 must be divisible by the
    spectator contraction to power h-3. Exact polynomial division checks
    that covariance constraint before restoring the three original sources.
    """
    order = integer(order, "tensor order")
    max_terms = integer(max_terms, "max_terms", 1)
    if len(polynomials) != order + 1:
        raise ValueError("every entry of the complete symmetric tensor is required")
    context = fmpq_mpoly_ctx.get(("spectator0", "spectator1", "index0", "index1"))
    first, second, index0, index1 = context.gens()
    homogeneous = defaultdict(dict)
    result = defaultdict(Fraction)
    for k, polynomial in enumerate(polynomials):
        for (a, b), value in polynomial.items():
            a, b = integer(a, "source exponent"), integer(b, "source exponent")
            value = Fraction(value) * comb(order, k)
            if not value:
                continue
            if a + b > order + 3:
                raise ArithmeticError("the complete source tensor exceeds the proven degree bound")
            if a + b <= 3:
                result[(order - k, k, 0), (3 - a - b, a, b)] += value
            else:
                homogeneous[a + b][a, b, order - k, k] = _q(value)
    for degree, terms in sorted(homogeneous.items()):
        if len(terms) > max_terms:
            raise RuntimeError("the full tensor exceeded max_terms before covariance division")
        power = degree - 3
        quotient, remainder = divmod(context.from_dict(terms), (first * index0 + second * index1) ** power)
        if remainder:
            raise ArithmeticError("the complete tensor violates polynomial GL(2) covariance")
        for (a, b, p, q), value in quotient.terms():
            a, b, p, q = map(int, (a, b, p, q))
            if a + b != 3 or p + q + power != order:
                raise ArithmeticError("covariance division changed a tensor or source degree")
            result[(p, q, power), (0, a, b)] += _fraction(value)
    if len(result) > max_terms:
        raise RuntimeError("the complete restored tensor exceeded max_terms")
    return tuple((counts, powers, coefficient) for (counts, powers), coefficient in sorted(result.items()) if coefficient)


def normalized_polynomials(terms, order):
    """Recover all normalized entries algebraically, independently of division."""
    order = integer(order, "tensor order")
    result = [defaultdict(Fraction) for k in range(order + 1)]
    for (p, q, r), (u, v, w), coefficient in terms:
        if p + q + r != order or u + v + w != 3:
            raise ValueError("invalid complete covariant tensor degree")
        for second in range(r + 1):
            k = q + second
            result[k][v + r - second, w + second] += Fraction(coefficient) * Fraction(comb(r, second), comb(order, k))
    return tuple({powers: value for powers, value in polynomial.items() if value} for polynomial in result)


def evaluate(terms, vectors, order):
    """Return every symmetric tensor entry at the complete original sources."""
    order = integer(order, "tensor order")
    if len(vectors) != 3 or any(len(vector) != 2 for vector in vectors):
        raise ValueError("all three complete two-component source spinors are required")
    if any(isinstance(value, float) for vector in vectors for value in vector):
        raise TypeError("source spinors must be exact")
    vectors = tuple(tuple(map(Fraction, vector)) for vector in vectors)
    bracket = lambda left, right: left[0] * right[1] - left[1] * right[0]
    brackets = (bracket(vectors[0], vectors[1]), bracket(vectors[2], vectors[1]), bracket(vectors[0], vectors[2]))
    result = [Fraction(0)] * (order + 1)
    for counts, powers, coefficient in terms:
        if sum(counts) != order or sum(powers) != 3:
            raise ValueError("invalid full covariant tensor degree")
        scalar = Fraction(coefficient)
        for value, power in zip(brackets, powers):
            scalar *= value ** power
        binary = [Fraction(1)]
        for vector, count in zip(vectors, counts):
            factor = [comb(count, k) * vector[0] ** (count - k) * vector[1] ** k for k in range(count + 1)]
            product = [Fraction(0)] * (len(binary) + count)
            for i, left in enumerate(binary):
                for j, right in enumerate(factor):
                    product[i + j] += left * right
            binary = product
        for k, value in enumerate(binary):
            result[k] += scalar * value / comb(order, k)
    return tuple(result)
