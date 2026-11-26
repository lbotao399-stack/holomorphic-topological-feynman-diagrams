"""Ranks of the transverse linear-field tensors allowed by each theta shape."""

from collections import Counter
from itertools import combinations, combinations_with_replacement, permutations, product

from sympy import Matrix, zeros


def traceless_basis():
    basis = [{(i, j): 1} for i in range(3) for j in range(3) if i != j]
    basis.extend([{(i, i): 1, (2, 2): -1} for i in range(2)])
    return basis


def sign(indices):
    return (-1) ** sum(y > x for i, x in enumerate(indices) for y in indices[:i])


def column(matrices):
    result = Counter()
    # Each coordinate is a weakly ordered monomial of M and a strictly
    # ordered exterior monomial in the lower, linear-field input slots.
    for entries in product(*(list(a.items()) for a in matrices)):
        ix = tuple(key[0] for key, _ in entries)
        jx = tuple(key[1] for key, _ in entries)
        if len(set(jx)) < len(jx):
            continue
        coefficient = sign(jx)
        for _, value in entries:
            coefficient *= value
        result[tuple(sorted(ix)), tuple(sorted(jx))] += coefficient
    return result


def lift_all_target_monomials(n):
    ranks = []
    for k in range(2, n + 1):
        rank = 0
        for upper in combinations_with_replacement(range(n), k):
            for lower in combinations(range(n), k):
                order = next(
                    (
                        u
                        for u in permutations(upper)
                        if all(i != j for i, j in zip(u, lower))
                    ),
                    None,
                )
                if order is not None:
                    entries = [{(i, j): 1} for i, j in zip(order, lower)]
                else:
                    assert len(set(upper)) == 1 and upper[0] in lower
                    i = upper[0]
                    other = next(j for j in lower if j != i)
                    entries = [
                        {(i, j): 1} if j != i else {(i, i): 1, (other, other): -1}
                        for j in lower
                    ]
                image = {entry: v for entry, v in column(entries).items() if v}
                assert image == {(upper, lower): 1}
                rank += 1
        ranks.append(rank)
    return ranks


def cartan_three_cycle():
    basis = []
    raw = traceless_basis()
    for entries in raw:
        m = zeros(3)
        for (i, j), v in entries.items():
            m[i, j] = v
        basis.append(m)
    inverse = Matrix([[(a * b).trace() for b in basis] for a in basis]).inv()
    dual = [
        sum((inverse[i, j] * b for j, b in enumerate(basis)), zeros(3))
        for i in range(8)
    ]

    def bracket(a, b):
        return a * b - b * a

    def three(a, b, c):
        return (a * bracket(b, c)).trace()

    # Jacobi and trace invariance make the trace three-form closed.
    for four in combinations(basis, 4):
        cocycle = 0
        for i, j in combinations(range(4), 2):
            other = [v for k, v in enumerate(four) if k not in (i, j)]
            cocycle += (-1) ** (i + j) * three(bracket(four[i], four[j]), *other)
        assert cocycle == 0
    boundary, image = Counter(), Counter()
    pairing = 0
    for slots in combinations(range(8), 3):
        coeff = three(*(dual[k] for k in slots))
        if not coeff:
            continue
        pairing += coeff * three(*(basis[k] for k in slots))
        for key, val in column([raw[k] for k in slots]).items():
            image[key] += coeff * val
        for i, j in combinations(range(3), 2):
            remaining = slots[3 - i - j]
            comm = bracket(basis[slots[i]], basis[slots[j]])
            for n, v in enumerate(dual):
                s = sign((n, remaining)) if n != remaining else 0
                boundary[tuple(sorted((n, remaining)))] += (
                    coeff * (-1) ** (i + j) * s * (v * comm).trace()
                )
    assert not any(boundary.values()) and not any(image.values())
    assert pairing == -8
    print("Cartan three-cycle, trace pairing", pairing, "universal graph tensor zero")


def run():
    basis = traceless_basis()
    expected = [8, 18, 10, 0]
    for k in range(1, 5):
        data = [column(selected) for selected in combinations(basis, k)]
        keys = sorted(set().union(*data))
        mat = Matrix([[row.get(key, 0) for row in data] for key in keys])
        rank = mat.rank()
        assert rank == expected[k - 1]
        print(k, len(data), rank, len(data) - rank)
    # The traceless basis is spanned by u tensor q with q(u) = 0.
    for i, j in permutations(range(3), 2):
        u = tuple(int(t == i) + int(t == j) for t in range(3))
        q = tuple(int(t == i) - int(t == j) for t in range(3))
        assert sum(x * y for x, y in zip(u, q)) == 0
        assert u[i] * q[i] == 1 and u[j] * q[j] == -1


if __name__ == "__main__":
    run()
    cartan_three_cycle()
    ranks = lift_all_target_monomials(5)
    assert ranks == [150, 350, 350, 126]
    print("five-dimensional ranks for two through five inputs", ranks)
