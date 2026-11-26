"""Exact character bounds for the extreme current profile and its scalar predecessor."""

from collections import Counter
from itertools import permutations

import numpy as np
from relative_dimensions import nodes, torus_data

PRIME = 998244353
WEIGHTS = [(5, -3, -3, -3, -4), (4, -2, -3, -3, -4), (4, -3, -3, -3, -3)]


def exterior(characters, degree, count):
    e = [1]
    for k in range(1, count + 1):
        v = (
            sum(
                (-1) ** (r - 1) * characters[r][degree] * e[k - r] % PRIME
                for r in range(1, k + 1)
            )
            % PRIME
        )
        e.append(v * pow(k, -1, PRIME) % PRIME)
    return e[count]


def main():
    factors, current, characters = torus_data(5, 32, 9, 15, -3, PRIME, 3)
    for profile, expected in [
        ((2,) * 8 + (7,), 175),
        ((2,) * 6 + (3, 7), 321),
        ((2,) * 7 + (8,), 8),
    ]:
        v = factors * current % PRIME
        for order, k in Counter(profile).items():
            v = v * exterior(characters, order - 1, k) % PRIME
        assert int(v.sum() % PRIME) == expected
        print(profile, expected)
    # Weyl's character formula for the dual of each specified highest weight.
    factors, _, characters = torus_data(5, 32, 8, 8, 0, PRIME, 3)
    points = np.array([n for n, _ in nodes(5, 32)])
    roots = np.array(
        [pow(3, ((PRIME - 1) // 32) * i, PRIME) for i in range(32)], dtype=np.int64
    )

    def alternating(w):
        v = np.zeros(len(points), dtype=np.int64)
        for p in permutations(range(5)):
            sign = (-1) ** sum(i > j for k, i in enumerate(p) for j in p[k + 1 :])
            v += sign * roots[(-points @ np.array([w[i] for i in p])) % 32]
        return v % PRIME

    rho = [4, 3, 2, 1, 0]
    inv = np.array([pow(int(x), -1, PRIME) for x in alternating(rho)], dtype=np.int64)
    for n, shift, expected in [(8, 0, [14, 145, 16]), (7, 1, [0, 8, 0])]:
        char = exterior(characters, 1, n)
        result = []
        for w in WEIGHTS:
            w = (w[0] + shift,) + w[1:]
            dual = alternating([x + y for x, y in zip(w, rho)]) * inv % PRIME
            result.append(int((factors * char % PRIME * dual % PRIME).sum() % PRIME))
        assert result == expected
        print("quadratic exterior degree", n, "multiplicities", result)


if __name__ == "__main__":
    main()
