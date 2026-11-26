"""All connected two-loop cubic color contractions with three external legs."""
from collections import Counter
from fractions import Fraction
from itertools import combinations, combinations_with_replacement, permutations


def graphs():
    pairs = list(combinations(range(5), 2))
    relabelings = list(permutations(range(5)))
    seen = set()
    for edges in combinations_with_replacement(pairs, 6):
        degrees = Counter(i for pair in edges for i in pair)
        if max(degrees.values()) > 3 or len(degrees) != 5:
            continue
        reached = {0}
        for _ in range(4):
            reached.update(j for i, j in edges if i in reached)
            reached.update(i for i, j in edges if j in reached)
        if len(reached) != 5:
            continue
        canonical = min(tuple(sorted(tuple(sorted((p[i], p[j]))) for i, j in edges))
                        for p in relabelings)
        if canonical not in seen:
            seen.add(canonical)
            yield canonical


def structure_constants(rank):
    # Orthogonal SU(2) and Gell-Mann SU(3) bases. Entries of 2f are
    # integer * sqrt(3)^power, so the full contraction uses integers only.
    entries = {(1, 2, 3): (2, 0)}
    if rank == 3:
        for i, j, k, value in ((1, 4, 7, 1), (2, 4, 6, 1), (2, 5, 7, 1),
                              (3, 4, 5, 1), (1, 5, 6, -1), (3, 6, 7, -1)):
            entries[i, j, k] = (value, 0)
        entries[4, 5, 8] = entries[6, 7, 8] = (1, 1)
    result = {}
    for triple, (value, power) in entries.items():
        for p in permutations(triple):
            inversions = sum(p[i] > p[j] for i, j in combinations(range(3), 2))
            result[p] = ((-1)**inversions * value, power)
    return result


def contract(edges, f):
    slots = [[] for _ in range(5)]
    for edge, (i, j) in enumerate(edges):
        slots[i].append(edge)
        slots[j].append(edge)
    external = 6
    for vertex in slots:
        while len(vertex) < 3:
            vertex.append(external)
            external += 1
    # Visit adjacent vertices first and index nonzero entries by all
    # previously assigned incident colors at that vertex.
    known = set()
    order = []
    while len(order) < 5:
        current = max((i for i in range(5) if i not in order),
                      key=lambda i: len(known.intersection(slots[i])))
        order.append(current)
        known.update(slots[current])
    lookups = []
    known.clear()
    for i in order:
        previous = [j for j, edge in enumerate(slots[i]) if edge in known]
        groups = {}
        for triple, value in f.items():
            key = tuple(triple[j] for j in previous)
            groups.setdefault(key, []).append((triple, value))
        lookups.append(([slots[i][j] for j in previous], groups))
        known.update(slots[i])

    labels = [0] * 9
    result = Counter()
    def visit(depth, coefficient, power):
        if depth == 5:
            result[tuple(labels[6:]), power] += coefficient
            return
        vertex = slots[order[depth]]
        previous, groups = lookups[depth]
        for triple, (factor, root_power) in groups.get(tuple(labels[k] for k in previous), ()):
            for edge, color in zip(vertex, triple):
                labels[edge] = color
            combined = power + root_power
            visit(depth + 1, coefficient * factor * (3 if combined == 2 else 1), combined % 2)
    visit(0, 1, 0)
    return {entry: value for entry, value in result.items() if value}


def main():
    count = Counter()
    for edges in graphs():
        universal = None
        for rank in (2, 3):
            f = structure_constants(rank)
            values = contract(edges, f)
            scale = Fraction(values.get(((1, 2, 3), 0), 0), 16 * rank**2 * 2)
            expected = {(triple, power): scale * 16 * rank**2 * value
                        for triple, (value, power) in f.items() if scale}
            assert values == expected
            if universal is not None:
                assert scale == universal
            universal = scale
        count[abs(universal)] += 1
    assert count == {0: 2, Fraction(1, 4): 1, Fraction(1, 2): 3, 1: 3}
    print('Unlabeled connected cores per absolute coefficient of C_A^2 f:', dict(sorted(count.items())))
    return count


if __name__ == '__main__':
    main()
