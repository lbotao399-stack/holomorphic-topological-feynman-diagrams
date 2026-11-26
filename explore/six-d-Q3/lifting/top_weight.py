"""Signed color interchange relations for five first modes at finite sl2 rank."""

from collections import Counter
from itertools import combinations, combinations_with_replacement, product

from quadratic import exact_rank


def canonical(flavors, colors):
    sign = 1
    result = []
    for flavor in range(5):
        block = tuple(c for f, c in zip(flavors, colors) if f == flavor)
        if flavor >= 2:
            if len(set(block)) != len(block):
                return 0, ()
            sign *= (-1) ** sum(first > second for first, second in combinations(block, 2))
        result.extend(sorted(block))
    return sign, tuple(result)


def quotient(flavors, charge):
    monomials = set()
    relations = []
    for colors in product(range(3), repeat=len(flavors)):
        if sum(1 - c for c in colors) != charge:
            continue
        sign, monomial = canonical(flavors, colors)
        if not sign or colors != monomial:
            continue
        monomials.add(monomial)
        for i, j in combinations(range(len(colors)), 2):
            if colors[i] == colors[j]:
                continue
            swapped = list(colors)
            swapped[i], swapped[j] = swapped[j], swapped[i]
            other_sign, other = canonical(flavors, swapped)
            relation = Counter({monomial: 1})
            if other_sign:
                relation[other] -= other_sign
            relations.append(relation)
    return len(monomials) - exact_rank(relations)


if __name__ == "__main__":
    count = 0
    total = 0
    for length in range(6):
        for flavors in combinations_with_replacement(range(5), length):
            singlets = quotient(flavors, 0) - quotient(flavors, 1)
            assert singlets == int(length % 2 == 0 and all(flavors.count(i) <= 1 for i in (2, 3, 4)))
            count += 1
            total += singlets
    print(f"{count} relative top multidegrees, {total} singlets, exact color-interchange ranks")
