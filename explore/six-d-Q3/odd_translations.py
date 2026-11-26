"""Component signs of the pure-ghost response at the smallest allowed insertion length."""
from functools import reduce
from itertools import combinations, product
from math import comb
from operator import or_


def wedge_variable(poly, bit, coefficient):
    result = {}
    for mon, value in poly.items():
        if not mon & bit:
            sign = (-1)**(mon // (2*bit)).bit_count()
            result[mon | bit] = coefficient * sign * value
    return result


def translation_volume(length):
    """Coefficients of the ordered product of odd-coordinate differences."""
    result = {0: 1}
    for slot in range(length-1):
        for flavor in range(2):
            direct = wedge_variable(result, 1 << (2*slot + flavor), 1)
            last = wedge_variable(result, 1 << (2*(length-1) + flavor), -1)
            for mon, value in last.items():
                direct[mon] = direct.get(mon, 0) + value
            result = {mon: value for mon, value in direct.items() if value}
    return result


def component_words(length):
    coefficients = {}
    for mon, value in translation_volume(length).items():
        degrees = [(mon >> (2*slot)) & 3 for slot in range(length)]
        # Collect all odd parameters on the left of all component fields.
        parity = 0
        crossings = 0
        for pair in degrees:
            crossings += parity * pair.bit_count()
            parity += (1+pair.bit_count()) % 2
        word = tuple(('c','g1','g2','b')[pair] for pair in degrees)
        coefficients[word] = (-1)**crossings * value
    return coefficients


def infinitesimal(poly, indices, raising=False):
    result = {}
    for mon, coefficient in poly.items():
        for j in indices:
            bit = 1 << j
            if not mon & bit:
                continue
            remainder = mon ^ bit
            sign = (-1)**(remainder & (bit-1)).bit_count()
            if raising:
                partner = bit >> 1
                if remainder & partner:
                    continue
                sign *= (-1)**(remainder & (partner-1)).bit_count()
                remainder |= partner
            result[remainder] = result.get(remainder, 0) + sign * coefficient
    return {mon: value for mon, value in result.items() if value}


def singlet_check(pairs, loops):
    from sympy import SparseMatrix
    subsets = list(combinations(range(pairs), loops))
    terms = {}
    rows = {}
    for col, (first, second) in enumerate(product(subsets, repeat=2)):
        mon = reduce(or_,(1 << (2*i) for i in first),0)
        mon |= reduce(or_,(1 << (2*i+1) for i in second),0)
        for target, sign in infinitesimal({mon:1}, range(1,2*pairs,2), raising=True).items():
            row = rows.setdefault(target,len(rows))
            terms[row,col] = sign
    rank = SparseMatrix(len(rows),len(subsets)**2,terms).rank()
    return len(subsets)**2-rank


if __name__ == '__main__':
    for length in range(1,10):
        poly = translation_volume(length)
        for flavor in range(2):
            assert not infinitesimal(poly, range(flavor,2*length,2))
        assert not infinitesimal(poly, range(1,2*length,2), raising=True)
        coeffs = component_words(length)
        assert len(coeffs) == length**2 and set(map(abs,coeffs.values())) == {1}
        assert coeffs[('b',)*(length-1)+('c',)] == 1
    coeffs = component_words(4)
    assert coeffs['g1','g2','b','b'] == -coeffs['b','b','b','c'] == -1
    for pairs in range(6):
        for n in range(pairs+1):
            expected = comb(pairs,n)*comb(pairs+1,n)//(n+1)
            assert singlet_check(pairs,n) == expected
    print('C10 signs:', coeffs['g1','g2','b','b'], coeffs['b','b','b','c'])
    print('Translation identities through nine slots and all flavor ranks through six slots verified.')
