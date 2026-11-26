"""Exact signed one-linear-field Ward identity for the sum of five theta shapes."""

from fractions import Fraction
from itertools import permutations, product
from math import factorial, prod

from theta_denominators import SHAPES


def inversions(elements):
    return (-1) ** sum(x > y for j, y in enumerate(elements) for x in elements[:j])


def arrows(paths, mask, out):
    r = []
    for h, group in enumerate(paths):
        f = bool(out >> h & 1)
        a = [f]
        for j in group:
            if j >= 0 and mask >> j & 1:
                f = not f
            a.append(f)
        r.append(a)
    return r


def sign(paths, arr, mask):
    words = sum(paths, ())
    ext = inversions(words) * inversions(
        tuple(j for j in words if j < 0 or not mask >> j & 1)
    )
    total = len(words)
    n = 0
    tails = []
    for row, ar in zip(paths, arr):
        vertices = [0, *range(n + 1, n + len(row) + 1), total + 1]
        n += len(row)
        for k, f in enumerate(ar):
            tails.append(vertices[k if f else k + 1])
            if k < len(row) and row[k] == -1 and f:
                ext = -ext
    return ext * inversions(tails)


def groups(s):
    return tuple(tuple(range(sum(s[:i]), sum(s[: i + 1]))) for i in range(3))


def endpoint_permutation(old, to_new):
    sign = 1
    for right in (False, True):
        selected = [i for i in range(3) if old[i][-1 if right else 0] == right]
        if len(selected) == 3:
            # The fiber epsilon and the endpoint momentum area both change
            # by the path permutation sign, since the three momenta sum to zero.
            continue
        if len(selected) != 1:
            raise ArithmeticError("endpoint orientation")
        sign *= inversions([to_new[i] for i in range(3) if i not in selected])
    return sign


def verify():
    shapes = {shape: (den, multiplicity) for shape, den, multiplicity in SHAPES}
    for shape, (_, automorphisms) in shapes.items():
        assert automorphisms == 2 * prod(factorial(shape.count(i)) for i in set(shape))
    c = (1, -1, 1)
    lifted = ((1, 0), (1, 1), (0, 1))
    assert all(sum(s * v[i] for s, v in zip(c, lifted)) == 0 for i in range(2))
    rates_checked = 0
    for sizes in product(range(6), repeat=3):
        if sum(sizes) != 9:
            continue
        rates = []
        for h in range(3):
            shape = tuple(m + 1 + int(i == h) for i, m in enumerate(sizes))
            if max(shape) <= 6:
                den = shapes[tuple(sorted(shape))][0]
                rates.append(Fraction((sizes[h] + 1) * (13 - 2 * shape[h]), 13 * den))
        assert rates and len(set(rates)) == 1
        rates_checked += len(rates)
    print("scalar weights", rates_checked)
    counts = [0] * 5
    for s in product(range(6), repeat=3):
        if sum(s) != 9:
            continue
        g = groups(s)
        for mask in range(1, 512):
            k = mask.bit_count()
            if k % 2 == 0:
                continue
            for out in (0, 3, 5, 6):
                seen = set()
                for h in range(3):
                    for pos in range(s[h] + 1):
                        w = list(g)
                        w[h] = w[h][:pos] + (-1,) + w[h][pos:]
                        a = arrows(w, mask, out)
                        seen.add(((-1) ** (h == 1)) * sign(w, a, mask))
                        counts[k // 2] += 1
                assert len(seen) == 1, (s, mask, out, seen)
    print("insertion", counts, flush=True)
    counts = [0] * 5
    for shape, _, _ in SHAPES:
        g = groups(tuple(n - 1 for n in shape))
        for mask in range(1, 1024):
            k = mask.bit_count()
            if k % 2 == 0:
                continue
            for out in (0, 3, 5, 6):
                a = arrows(g, mask, out)
                for perm in permutations(range(3)):
                    inv = [perm.index(i) for i in range(3)]
                    g2 = tuple(g[i] for i in perm)
                    a2 = [a[i] for i in perm]
                    d = inversions(perm) ** 5
                    for h, row in enumerate(g):
                        d *= (-1) ** (len(row) * (int(h == 1) + int(inv[h] == 1)))
                    assert sign(g, a, mask) == d * endpoint_permutation(a, inv) * sign(
                        g2, a2, mask
                    )
                    counts[k // 2] += 1
    print("permutation", counts)


if __name__ == "__main__":
    verify()
