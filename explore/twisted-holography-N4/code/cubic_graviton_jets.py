"""Finite jet checks for both color sectors of cubic superfield cohomology on C^(2|3)."""

from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, combinations_with_replacement, permutations, product

from single_particle import bulk_basis_A2, div_A2


# Each jet variable has a site and a direction; the last three directions are odd.
def odd(variable):
    return variable % 5 >= 2


def times_variable(poly, variable):
    result = Counter()
    for term, value in poly.items():
        if odd(variable) and variable in term:
            continue
        inversions = odd(variable) and sum(odd(x) and x > variable for x in term)
        result[tuple(sorted((*term, variable)))] += (-1)**inversions * value
    return result


def total_derivative(poly, direction, sites):
    result = Counter()
    for site in sites:
        result.update(times_variable(poly, 5*site + direction))
    return clean(result)


def clean(poly):
    return {key: value for key, value in poly.items() if value}


@lru_cache(maxsize=8192)
def substituted_term(term, destinations):
    result = {(): 1}
    for variable in term:
        result = total_derivative(result, variable % 5, destinations[variable // 5])
    return result


def substitute(poly, destinations):
    result = Counter()
    for term, value in poly.items():
        for key, coefficient in substituted_term(term, destinations).items():
            result[key] += value * coefficient
    return clean(result)


def permute_slots(term, sites=3, alternating=True):
    result = Counter()
    for perm in permutations(range(sites)):
        inversions = alternating and sum(perm[i] > perm[j] for i, j in combinations(range(sites), 2))
        for key, value in substituted_term(term, tuple((i,) for i in perm)).items():
            result[key] += (-1)**inversions * value
    return clean(result)


def relative_symbols(weights, sites=3, alternating=True):
    choices = []
    for direction, degree in enumerate(weights):
        distribute = combinations_with_replacement if direction < 2 else combinations
        choices.append(list(distribute(range(sites), degree)))
    seen = set()
    for assignment in product(*choices):
        if len(set().union(*map(set, assignment))) != sites:
            continue
        term = tuple(sorted(5*site + direction for direction, sites in enumerate(assignment)
                            for site in sites))
        if term not in seen:
            polynomial = permute_slots(term, sites, alternating)
            seen.add(term)
            seen.update(polynomial)
            if polynomial:
                yield polynomial


def classical_variation(cubic):
    # Coefficient of the traceless Jordan tensor for the first pair. The three
    # tensors are independent; other coefficients are permutations of this one.
    variation = Counter(substitute(cubic, ((0, 1), (2,), (3,))))
    variation.subtract(substitute(cubic, ((0,), (2,), (1, 3))))
    variation.subtract(substitute(cubic, ((1,), (2,), (0, 3))))
    return clean(variation)


def commutator_variation(cubic):
    variation = Counter(substitute(cubic, ((0, 1), (2,), (3,))))
    variation.update(substitute(cubic, ((0,), (1,), (2, 3))))
    variation.subtract(substitute(cubic, ((0,), (1, 2), (3,))))
    variation.subtract(substitute(cubic, ((1,), (2,), (0, 3))))
    return clean(variation)


def quadratic_image(primitive):
    variation = Counter()
    for k in range(3):
        variation.update(substitute(primitive, ((k, (k+1) % 3), ((k+2) % 3,))))
    return clean(variation)


def graviton_components(weights):
    for slots in combinations_with_replacement(range(5), 3):
        remaining = [degree - slots.count(i) for i, degree in enumerate(weights)]
        if min(remaining) < 0 or max(remaining[2:]) > 1:
            continue
        polynomial = permute_slots(tuple(5*site+i for site, i in enumerate(slots)))
        for direction, count in enumerate(remaining):
            for _ in range(count):
                polynomial = total_derivative(polynomial, direction, range(3))
        if polynomial:
            yield polynomial


def rank(columns):
    pivots = {}
    for column in columns:
        reduced = {key: Fraction(value) for key, value in column.items()}
        while reduced:
            key = max(reduced)
            scale = reduced[key]
            if key not in pivots:
                pivots[key] = {i: value / scale for i, value in reduced.items()}
                break
            for i, value in pivots[key].items():
                reduced[i] = reduced.get(i, 0) - scale * value
                if not reduced[i]:
                    del reduced[i]
    return len(pivots)


def bulk_dimension(weights):
    # The two tangent polyvector slots map to the complementary tangent forms.
    # Divergence raises the form degree by one. Regular flux images have degree one.
    weight = tuple(x-1 for x in weights[:2]) + tuple(-x-1 for x in weights[2:])
    degree_three = [term for term in bulk_basis_A2(weight)
                    if 2 - sum(i < 2 for i in term[2]) + sum(term[1]) == 3]
    return len(degree_three) - rank(map(div_A2, degree_three))


def check():
    # Include the first several possible one-loop source weights of a two-ghost
    # graviton derivative, a repeated odd form direction, and mixed tangent jets.
    weights = [(k, 0, 2, 1, 1) for k in range(4)]
    weights += [(1, 1, 1, 0, 0), (2, 1, 1, 0, 0), (2, 2, 1, 0, 0),
                (0, 0, 3, 0, 0), (1, 1, 2, 1, 0), (1, 0, 2, 2, 1)]
    for degrees in weights:
        cubics = list(relative_symbols(degrees))
        components = list(graviton_components(degrees))
        assert all(not classical_variation(poly) for poly in components)
        kernel = len(cubics) - rank(map(classical_variation, cubics))
        assert kernel == rank(components)
        commutators = list(relative_symbols(degrees, alternating=False))
        closed = len(commutators) - rank(map(commutator_variation, commutators))
        boundaries = [quadratic_image(primitive) for primitive in relative_symbols(degrees, sites=2)]
        assert all(not commutator_variation(poly) for poly in boundaries)
        assert closed == rank(boundaries)
        print(degrees, "symmetric cochains", len(cubics), "closed gravitons", kernel,
              "commutator cochains", len(commutators), "closed and exact", closed)

    # Count all classes in each of the first three homogeneous jet degrees
    # without using the form-complex character.
    for total, expected in [(3, 25), (4, 84), (5, 156)]:
        count = 0
        for directions in combinations_with_replacement(range(5), total):
            degrees = tuple(directions.count(i) for i in range(5))
            cubics = list(relative_symbols(degrees))
            kernel = len(cubics) - rank(map(classical_variation, cubics))
            assert kernel == bulk_dimension(degrees)
            count += kernel
        assert count == expected
        print("Total degree", total, "independent closed cubic cochains", count)


if __name__ == "__main__":
    check()
