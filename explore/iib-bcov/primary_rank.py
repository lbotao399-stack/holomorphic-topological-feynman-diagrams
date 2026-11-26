"""Exact rational-momentum nonvanishing tests for the polarized worldsheet primary."""
from collections import deque
from fractions import Fraction
from functools import cache
from itertools import combinations
from math import lcm
from random import Random

from flint import fmpz_mat


def integer_rows(waves):
    rows = [tuple(map(Fraction, p)) for p in waves]
    if not rows or not rows[0] or any(len(p) != len(rows[0]) for p in rows):
        raise ValueError("Expected a nonempty matrix of rational waves")
    if any(map(sum, zip(*rows))):
        raise ValueError("The sum of the momenta must vanish")
    # Rescaling a row here preserves its linear independence relations.
    integer = []
    for row in rows:
        denominator = lcm(*(x.denominator for x in row))
        integer.append(tuple(int(x*denominator) for x in row))
    return integer


def rank_of(rows, slots):
    return fmpz_mat([list(rows[i]) for i in slots]).rank() if slots else 0


def partition_witness(rows, genus):
    """Return independent classes or a subset of size greater than genus times its rank."""
    @cache
    def independent(slots):
        return rank_of(rows, slots) == len(slots)

    classes = [frozenset() for _ in range(genus)]
    for new_slot in range(len(rows)):
        predecessor = {new_slot: None}
        queue = deque([new_slot])
        while queue:
            slot = queue.popleft()
            for color, group in enumerate(classes):
                if slot in group:
                    continue
                extended = group | {slot}
                if independent(extended):
                    while True:
                        classes[color] |= {slot}
                        previous = predecessor[slot]
                        if previous is None:
                            break
                        parent, color = previous
                        classes[color] -= {slot}
                        slot = parent
                    break
                for displaced in sorted(group):
                    if displaced not in predecessor and independent(extended - {displaced}):
                        predecessor[displaced] = slot, color
                        queue.append(displaced)
            else:
                continue
            break
        else:
            obstruction = tuple(sorted(predecessor))
            assert len(obstruction) > genus*rank_of(rows, obstruction)
            return None, obstruction
        assert all(independent(group) for group in classes)
    assert sorted(i for group in classes for i in group) == list(range(len(rows)))
    return [tuple(sorted(group)) for group in classes], None


def primary_positive(genus, momenta):
    """Li–Wang worldsheet normalization, g >= 2, primary arity (N-3)g+3."""
    rows = integer_rows(momenta)
    dimension = len(rows[0])
    if genus < 2 or len(rows) != (dimension-3)*genus + 3:
        raise ValueError("Expected stable primary arity at genus at least two")
    if rank_of(rows, range(len(rows))) < dimension:
        return False
    classes, _ = partition_witness(rows, genus)
    return classes is not None


def check():
    rng = Random(192404)
    count = 0
    for dimension, genus in ((4, 2), (4, 4), (5, 2), (5, 3), (5, 4), (6, 3)):
        basis = [tuple(int(i == j) for i in range(dimension)) for j in range(dimension)]
        for _ in range(14):
            size = (dimension-3)*genus+3
            waves = [rng.choice(basis) if rng.randrange(4) else
                     tuple(rng.randrange(-1, 2) for _ in range(dimension)) for _ in range(size-1)]
            waves.append(tuple(-sum(c) for c in zip(*waves)))
            rows = integer_rows(waves)
            ok = primary_positive(genus, waves)
            full = rank_of(rows, range(size)) == dimension
            brute = full and all(len(subset) <= genus*rank_of(rows, subset)
                                 for k in range(1, size+1) for subset in combinations(range(size), k))
            assert ok == brute
            if dimension in (4, 5):
                simple = full and all(any(p) for p in rows)
                if dimension == 5:
                    simple = simple and all(
                        sum(rank_of(rows, (i, j)) == 1 for j in range(size)) <= genus
                        for i in range(size))
                assert ok == simple
            count += 1
    basis = [tuple(int(i == j) for i in range(5)) for j in range(5)]
    # Genus four: two full-rank, conserved examples on opposite sides of the line cutoff.
    for length in (4, 5):
        waves = [basis[0]]*length + basis[1:] + [basis[1]]*(6-length)
        waves.append(tuple(-sum(c) for c in zip(*waves)))
        assert primary_positive(4, waves) == (length == 4)
    print("Exact partition, violating-subset and exhaustive checks:", count)


if __name__ == '__main__':
    check()
