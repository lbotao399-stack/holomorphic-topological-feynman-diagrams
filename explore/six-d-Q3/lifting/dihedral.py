"""Exact polarization ranks for the two Cartan scalars and one Cartan antifield."""

from itertools import combinations_with_replacement, product
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "s2code"))
from salg import add, const, gen, mul
from quadratic import exact_rank


def word(letters):
    out = const()
    for letter in letters:
        out = mul(out, letter)
    return out


def check(m, limit):
    colors = [[gen((1 if a == 2 else 3, a, i, 0, 0)) for i in range(2)] for a in range(3)]
    generators = []
    for k in (2, m):
        for flavors in combinations_with_replacement(range(3), k):
            if k == 2:
                signs = (0, 1)
            else:
                signs = (0,) * k
            positive = word([colors[a][s] for a, s in zip(flavors, signs)])
            negative = word([colors[a][1-s] for a, s in zip(flavors, signs)])
            p = add(positive, negative)
            if p:
                generators.append((tuple(flavors.count(i) for i in range(3)), p))
    pieces = {(0, 0, 0): [const()]}
    for weight, p in generators:
        for total in range(limit + 1):
            for (x, y, b), basis in list(pieces.items()):
                if x + y + b != total or total + sum(weight) > limit:
                    continue
                key = tuple(s + t for s, t in zip((x, y, b), weight))
                if key[2] > 2:
                    continue
                rows = pieces.setdefault(key, [])
                old_rank = len(rows)
                for candidate in (mul(p, f) for f in basis):
                    if candidate and exact_rank([*rows, candidate]) > old_rank:
                        rows.append(candidate)
                        old_rank += 1
    sectors = 0
    for a, d, b in product(range(limit + 1), range(limit + 1), range(3)):
        if a + d + b > limit:
            continue
        dims = 0
        seen = set()
        for xp, yp in product(range(a + 1), range(d + 1)):
            for plus in (0, 1):
                minus = b - plus
                if minus not in (0, 1) or (2*(xp + yp + plus) - a - d - b) % m:
                    continue
                mon = (xp, yp, plus, minus)
                reflected = (a - xp, d - yp, minus, plus)
                if mon in seen:
                    continue
                seen.update((mon, reflected))
                if reflected != mon or b < 2:
                    dims += 1
        assert len(pieces.get((a, d, b), [])) == dims, (m, a, d, b, dims)
        sectors += 1
    return sectors


if __name__ == "__main__":
    print(f"G2 dihedral polarizations generate the exact Cartan answer in {check(6, 9)} graded sectors")
