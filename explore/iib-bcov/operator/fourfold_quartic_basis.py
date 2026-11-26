"""Seven conserved fourfold waves: symmetric genus-four primary symbols and their collinear image."""
from collections import Counter
from fractions import Fraction
from itertools import combinations, combinations_with_replacement, permutations
from math import factorial, prod
from random import Random

import numpy as np
from flint import fmpz_mat, nmod_mat


def slot_series(values, cycle, minimum):
    """Multiply by the single-slot torus character with its low degrees removed."""
    result = values.copy()
    for axis, size in enumerate(values.shape):
        for residue in range(min(cycle, size)):
            selected = [slice(None)] * values.ndim
            selected[axis] = slice(residue, None, cycle)
            section = result[tuple(selected)]
            np.cumsum(section, axis=axis, out=section)
    if minimum:
        result -= values
    if minimum == 2:
        for axis in range(values.ndim):
            lower, upper = [slice(None)] * values.ndim, [slice(None)] * values.ndim
            lower[axis], upper[axis] = slice(None, -cycle), slice(cycle, None)
            result[tuple(upper)] -= values[tuple(lower)]
    return result


def weight_eight_multiplicity(array):
    result = array.copy()
    for axis in range(4):
        upper, lower = [slice(None)] * 4, [slice(None)] * 4
        upper[axis], lower[axis] = slice(1, None), slice(None, -1)
        result[tuple(upper)] -= result[tuple(lower)].copy()
    answer = 0
    for perm in permutations(range(4)):
        powers = tuple(8-i+j for i, j in enumerate(perm))
        sign = (-1)**sum(a>b for a, b in combinations(perm, 2))
        answer += sign * result[powers]
    return answer


def partitions(n, cap):
    if not n:
        yield ()
    for part in range(1, min(n, cap)+1):
        for rest in partitions(n-part, part):
            yield (part, *rest)


def dimension(minimum):
    # The total-momentum regular sequence contributes the four factors (1-x_i).
    # The conserved single-soft ideals have product resolution length <=21
    # in dimension 24. The quotient's length <=22 excludes origin torsion,
    # so the product equals the intersection of the required slotwise powers.
    unit = np.zeros((12, 11, 10, 9), dtype=object)
    unit[0, 0, 0, 0] = 1
    sym = [unit]
    for m in range(1, 8):
        total = sum(slot_series(sym[m-j], j, minimum) for j in range(1, m+1))
        assert not np.any(total % m)
        sym.append(total // m)
    answer = weight_eight_multiplicity(sym[-1])
    by_cycle = Fraction()
    for cycle in partitions(7, 7):
        result = unit
        for length in cycle:
            result = slot_series(result, length, minimum)
        order = prod(k**m * factorial(m) for k, m in Counter(cycle).items())
        by_cycle += Fraction(weight_eight_multiplicity(result), order)
    assert by_cycle == answer
    return answer


SUBSETS = tuple(combinations(range(7), 4))
SUBSET_IDS = {subset: i for i, subset in enumerate(SUBSETS)}
ORDINARY_WORDS = ((13, 14, 15, 16, 18, 24, 27, 29), (16, 17, 17, 17, 19, 23, 23, 28))


def permuted_minors():
    indices, signs = [], []
    for perm in permutations(range(7)):
        row = [tuple(perm[i] for i in s) for s in SUBSETS]
        indices.append([SUBSET_IDS[tuple(sorted(t))] for t in row])
        signs.append([(-1)**sum(a>b for a, b in combinations(t, 2)) for t in row])
    return np.asarray(indices, dtype=np.int16), np.asarray(signs, dtype=np.int16)


def square_orbits(indices):
    masks = [sum(1<<k for k in s) for s in SUBSETS]
    words = {word for word in combinations_with_replacement(range(35), 4)
             if masks[word[0]] | masks[word[1]] | masks[word[2]] | masks[word[3]] == 127}
    groups = {}
    for word in sorted(words):
        if word in words:
            orbit = set(map(tuple, np.sort(indices[:, word], axis=1).tolist()))
            words -= orbit
            groups[word] = np.asarray(sorted(orbit), dtype=np.int16)
    return groups


def frames():
    rng = Random(47851)
    for _ in range(28):
        rows = [[int(i==j) for j in range(4)] for i in range(4)]
        rows += [[rng.randrange(-6, 7) for _ in range(4)] for _ in range(2)]
        yield rows + [[-sum(c) for c in zip(*rows)]]


def collinear(strengths):
    rows = [[t, 0, 0, 0] for t in strengths]
    rows += [[int(i==j) for j in range(4)] for i in (1, 2, 3)]
    return rows + [[-sum(c) for c in zip(*rows)]]


def minors(frame):
    return [int(fmpz_mat([frame[i] for i in s]).det()) for s in SUBSETS]


def square_values(groups, values):
    squares = np.asarray([m*m for m in values], dtype=object)
    return [int(np.prod(squares[orbit], axis=1).sum()) for orbit in groups.values()]


def pivots(rows):
    mat, rank = nmod_mat(rows, 1009).rref()
    return [next(j for j in range(mat.ncols()) if mat[i,j]) for i in range(rank)]


def verify_basis(groups, indices, signs):
    values = [minors(p) for p in frames()]
    rows = [square_values(groups, p) for p in values]
    selected = pivots(rows)
    assert len(selected) == fmpz_mat(rows).rank() == 26
    cols = [[row[j] % 1009 for row in rows] for j in selected]
    # Full S7 sums for the two ordinary words include signs from permuting rows.
    table = np.asarray(values, dtype=np.int32).T % 1009
    tables = np.transpose(table[indices] * signs[:,:,None], (1, 0, 2)) % 1009
    for word in ORDINARY_WORDS:
        assert all(sum(i in SUBSETS[m] for m in word) >= 2 for i in range(7))
        product = tables[word[0]].copy()
        for k in word[1:]:
            product *= tables[k]
            product %= 1009
        cols.append((product.sum(axis=0) % 1009).tolist())
    assert int(nmod_mat(cols, 1009).det()) == 870
    return [list(groups)[j] for j in selected]


def verify_collinear(groups):
    keys = ((0, 0, 1, 9), (0, 1, 2, 3))
    selected = {word: groups[word] for word in keys}
    for strengths in ((1, 1, 1), (1, 1, -2), (2, 3, 4)):
        q = collinear(strengths)
        result = square_values(selected, minors(q))
        assert result == [24*prod(strengths)**2*sum(t*t for t in strengths),
                          4*prod(strengths)**2*sum(strengths)**2]
    # Quadratic softness divides out (abc)^2; symmetry and weight eight
    # leave only a symmetric homogeneous quadratic in (a,b,c). Thus the
    # first two exact evaluations also identify each formula polynomially.
    assert int(fmpz_mat([[72, 36], [576, 0]]).det()) == -20736


if __name__ == '__main__':
    dims = [dimension(minimum) for minimum in (0, 1, 2)]
    assert dims == [40, 34, 28]
    indices, signs = permuted_minors()
    groups = square_orbits(indices)
    assert len(groups) == 41
    selected = verify_basis(groups, indices, signs)
    verify_collinear(groups)
    print('Unrestricted, soft and doubly soft dimensions:', *dims)
    print('26 squared-minor orbit representatives:', selected)
    print('Two signed representatives:', ORDINARY_WORDS)
    print('Full determinant modulo 1009: 870. Collinear image rank: 2.')
