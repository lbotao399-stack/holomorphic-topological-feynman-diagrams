"""Scalar boundary constraints and quadratic-bracket rank certificates by Taylor profile."""

from collections import Counter
from functools import cache
from itertools import combinations, combinations_with_replacement, permutations
from math import comb
from pathlib import Path
import shlex
import subprocess
import os
from tempfile import TemporaryDirectory
from extreme_fusion import Jet

HERE = Path(__file__).resolve().parent


def pack(weight):
    return sum(x * 64**i for i, x in enumerate(weight))


RHO = (4, 3, 2, 1, 0)
WEYL = [
    (
        (-1) ** sum(i < j for k, i in enumerate(p) for j in p[k + 1 :]),
        pack([x - y for x, y in zip(RHO, p)]),
    )
    for p in permutations(RHO)
]


@cache
def exteriors(order, count, dual=False):
    weights = Jet(order).weight
    # A five-coordinate weight is encoded linearly; each coordinate difference
    # in these homogeneous character products is strictly smaller than 64.
    d = [Counter() for _ in range(count + 1)]
    d[0][0] = 1
    for k, m in Counter(pack(w) * (-1 if dual else 1) for w in weights).items():
        for j in range(count, 0, -1):
            for n in range(1, min(j, m) + 1):
                f = comb(m, n)
                for w, c in d[j - n].items():
                    d[j][w + n * k] += f * c
    return d[count]


def multiply(a, b):
    c = Counter()
    for x, i in a.items():
        for y, j in b.items():
            c[x + y] += i * j
    return c


def multiplicity(character, weight):
    return sum(sign * character[pack(weight) + k] for sign, k in WEYL)


def partners(profile, current):
    q = profile.count(2)
    high = Counter({pack([-3] * 5): 1})
    if current:
        high = multiply(high, Counter({64**i: 1 for i in range(5)}))
    for order, count in Counter(x for x in profile if x != 2).items():
        high = multiply(high, exteriors(order, count, True))
    low = exteriors(2, q)
    output = []
    # Each high jet contributes at most one negative coordinate to its dual.
    # The low exterior power contributes at most q to any positive coordinate.
    for rev in combinations_with_replacement(range(-3 - len(profile) + q, q + 1), 5):
        w = rev[::-1]
        if sum(w) != -q:
            continue
        a = multiplicity(low, w)
        b = multiplicity(high, w)
        if a and b:
            output.append((q, a, b, low[pack(w)], *w))
    return output


def brackets(profile):
    return {
        tuple(
            sorted(
                [d for k, d in enumerate(profile) if k not in (i, j)]
                + [profile[i] + profile[j] - 1]
            )
        )
        for i, j in combinations(range(len(profile)), 2)
    }


def run():
    small = (2,) * 6 + (4, 7)
    third = (2,) * 7 + (4, 5)
    scalar_ranks = {}
    with TemporaryDirectory(prefix="twist-profile-ranks-") as tmp:
        exe = Path(tmp) / "quadratic-rank"
        subprocess.run(
            [
                *shlex.split(os.environ.get("CXX", "c++")),
                "-O3",
                "-std=c++17",
                "-fsanitize=signed-integer-overflow",
                "-fno-sanitize-recover=all",
                HERE / "quadratic_rank.cpp",
                "-o",
                exe,
            ],
            check=True,
        )
        for profile, current, threshold, expected in [
            (small, False, -7, (76, 53)),
            ((2,) * 6 + (5, 6), False, -7, (135, 98)),
            ((2,) * 5 + (4, 4, 5), False, -7, (1729, 1076)),
            (third, True, 4, (998, 867)),
        ]:
            rows = [r for r in partners(profile, current) if r[4] >= threshold]
            data = "".join(" ".join(map(str, r)) + "\n" for r in rows)
            output = subprocess.run(
                [exe], input=data, text=True, capture_output=True, check=True
            ).stdout.splitlines()
            out = [tuple(map(int, row.split())) for row in output]
            sizes = (sum(r[-3] * r[-2] for r in out), sum(r[-2] * r[-1] for r in out))
            assert sizes == expected, (profile, sizes)
            if not current:
                scalar_ranks[profile] = sizes[1]
            print("profile", profile, "tested directions, rank", sizes, flush=True)
    profiles = [
        (2,) * 8 + (8,),
        (2,) * 7 + (3, 7),
        (2,) * 7 + (4, 6),
        (2,) * 6 + (3, 3, 6),
        (2,) * 6 + (4, 4, 4),
        (2,) * 5 + (3, 3, 3, 5),
    ]
    # Large profiles are integrated via the torus character implementation.
    from extreme_character import PRIME, exterior, torus_data

    factors, _, chars = torus_data(5, 32, 9, 15, -3, PRIME, 3)

    @cache
    def dim(profile):
        a = factors.copy()
        for d, k in Counter(profile).items():
            a = a * exterior(chars, d - 1, k) % PRIME
        return int(a.sum() % PRIME)

    for count in range(1, len(profiles) + 1):
        best = 1
        for size in range(1, count + 1):
            for subset in combinations(profiles[:count], size):
                targets = set().union(*(brackets(p) for p in subset))
                outgoing = sum(r for q, r in scalar_ranks.items() if q in targets)
                best = max(
                    best, sum(map(dim, subset)) - sum(map(dim, targets)) + outgoing
                )
        print(
            "scalar profiles through",
            profiles[count - 1],
            "independent cycles at least",
            best,
        )
        if count == 2:
            assert best == 120


if __name__ == "__main__":
    run()
