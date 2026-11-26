"""Check the six component Cartan contractions and the source graph gradings."""
from collections import defaultdict


def act(mask, index, insert):
    bit = 1 << index
    if bool(mask & bit) == insert:
        return mask, 0
    sign = (-1) ** ((mask & (bit - 1)).bit_count())
    return mask ^ bit, sign


def component_product(mask):
    terms = {}
    left = mask
    while True:
        right = mask ^ left
        inversions = sum((right & ((1 << i) - 1)).bit_count() for i in range(6) if left >> i & 1)
        terms[left, right] = (-1) ** inversions
        if not left:
            return terms
        left = (left - 1) & mask


def verify_cartan():
    for mask in range(64):
        for i in range(6):
            added, sign = act(mask, i, True)
            lhs = {pair: sign * c for pair, c in component_product(added).items() if sign}
            rhs = defaultdict(int)
            for (left, right), c in component_product(mask).items():
                a, sa = act(left, i, True)
                b, sb = act(right, i, True)
                rhs[a, right] += c * sa
                rhs[left, b] += c * (-1) ** left.bit_count() * sb
            assert lhs == {pair: c for pair, c in rhs.items() if c}
            for j in range(6):
                raised, sign1 = act(mask, i, True)
                first, sign2 = act(raised, j, False)
                lowered, sign3 = act(mask, j, False)
                second, sign4 = act(lowered, i, True)
                total = defaultdict(int)
                total[first] += sign1 * sign2
                total[second] += sign3 * sign4
                assert {m: v for m, v in total.items() if v} == ({mask: 1} if i == j else {})


def verify_grading():
    for m in range(1, 6):
        vectors, strengths = 11 - 2 * m, 2 * m - 1
        taylor = 22 + m - strengths
        assert vectors + strengths == 10
        assert 2 * taylor - 2 * vectors - strengths == 25
    # A trivalent graph has three flags per vertex and two flags per edge.
    for n in range(1, 17):
        arity = 3 * (5 * n + 2) - 2 * (6 * n + 1)
        assert arity + 5 == 3 * n + 9
        assert arity + 5 - 3 * (n - 1) == 12


if __name__ == '__main__':
    verify_cartan()
    verify_grading()
    print('64 form components; 2304 Cartan anticommutators; 384 exterior-product identities')
    print('five strict input sectors at arity 10 and principal weight 25')
