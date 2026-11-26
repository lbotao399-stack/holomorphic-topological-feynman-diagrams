"""Fundamental single-trace color in the LLC fan with three beta letters."""
from collections import Counter
from functools import cache
from itertools import permutations
from math import comb, prod


def cycles(successor):
    seen = set()
    for i in range(len(successor)):
        if i in seen:
            continue
        orbit = []
        while i not in seen:
            seen.add(i)
            orbit.append(i)
            i = successor[i]
        yield orbit


def face_orders(size):
    """Fix the second vertex order. The m first-vertex gaps are counted separately."""
    for middle in permutations(range(1, size)):
        first = (0,) + middle
        following = [0] * size
        for i, j in zip(first, first[1:] + (0,)):
            following[i] = j
        yield following[1:] + following[:1]


def handle_coefficients(size):
    """Count the cyclic trace DEF, retaining empty fundamental index loops."""
    result = Counter()
    for successor in face_orders(size):
        faces = list(cycles(successor))
        containing_end = next(face for face in faces if size - 1 in face)
        pos = containing_end.index(size - 1)
        face = containing_end[pos:] + containing_end[:pos]
        weight = len(face)
        for distance, gap in enumerate(face[1:], 1):
            k = size - 1 - gap
            weight += (-1)**k * comb(size, k) * distance
        handle = (size - len(faces)) // 2
        result[handle] += weight
    return dict(result)


@cache
def moving_gap_counts(size, handle):
    """Count only permutations moving every gap, with the given cycle defect."""
    count = excess = 0
    distances = [0] * (size - 1)
    for successor in face_orders(size):
        if any(i == j for i, j in enumerate(successor)):
            continue
        faces = list(cycles(successor))
        if size - len(faces) != 2 * handle:
            continue
        count += 1
        gap = successor[-1]
        distance = 1
        while gap != size - 1:
            distances[size - 2 - gap] += distance
            excess += 1
            gap = successor[gap]
            distance += 1
    return count, excess, distances


def fixed_handle(size, handle):
    """The elementary binomial sum only needs permutations on at most 4h gaps."""
    if handle == 0:
        return 1
    total = 0
    for moved in range(2 * handle + 1, min(size, 4 * handle) + 1):
        count, excess, distances = moving_gap_counts(moved, handle)
        total += count * comb(size, moved) + excess * comb(size - 1, moved - 1)
        for j, distance in enumerate(distances):
            other = moved - 2 - j
            total -= distance * ((-1)**j * comb(size - 1, other)
                                 + (-1)**(size + other) * comb(size - 1, j))
    return total


def one_handle(n):
    return comb(n + 2, 4) + n * (n - 1) + comb(n, 3) + 2 + (-1)**n * (n - 1)


def check_unmarked(size, rank):
    colors = sum(rank ** len(list(cycles(order))) for order in face_orders(size))
    difference = prod(rank + j for j in range(size + 1)) - prod(
        rank - j for j in range(size + 1))
    assert colors * size * (size + 1) == difference


if __name__ == "__main__":
    for n in range(2, 10):
        direct = handle_coefficients(n + 1)
        reduced = {h: fixed_handle(n + 1, h) for h in range(n // 2 + 1)}
        assert direct == reduced and direct[0] == 1 and direct[1] == one_handle(n)
        print(n, [direct[h] for h in sorted(direct)])
        if n <= 7:
            check_unmarked(n + 1, 3)
    assert handle_coefficients(3) == {0: 1, 1: 6}
