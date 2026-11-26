"""Component and matrix replays for the mixed quartic and its one-loop completion."""

from fractions import Fraction
from functools import reduce
from itertools import permutations
from pathlib import Path
import sys
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "s2code"))
from ops import Theory, L
from salg import EPS, add, mul, scale, const


def word(*factors):
    return reduce(mul, factors, const())


def jet(i):
    return tuple(int(a == i) for a in range(3))


def contraction(first, second):
    result = {}
    for i in range(3):
        result = add(result, mul(L(first, i), L(second, i)))
    return result


class GhostBox(Theory):
    """For c,g1,g2-only inputs the two scalar derivatives saturate the Berezin measure."""

    def T1(self, left, right, momenta):
        key = (left, right, momenta)
        if key in self._T1:
            return self._T1[key]
        result = {}
        for d1, b, c1 in self.fby[left]:
            for d2, c, c2 in self.fby[b]:
                for d3, a, c3 in self.fby[c]:
                    if a != right:
                        continue
                    for axes, sign in EPS.items():
                        indices = [
                            tuple(i + int(j == axis) for j, i in enumerate(p))
                            for axis, p in zip(axes, momenta)
                        ]
                        mon = word(
                            *(L("c", d, k) for d, k in zip((d1, d2, d3), indices))
                        )
                        result = add(result, scale(mon, sign * c1 * c2 * c3))
        self._T1[key] = result
        return result


def check_components():
    theory = GhostBox(3, {p: Fraction(s) for p, s in EPS.items()})
    xx = contraction("g1", "g1")
    yy = contraction("g2", "g2")
    xy = contraction("g1", "g2")
    quartic = add(mul(xx, yy), mul(xy, xy), -1)
    primitive = {}
    for (i, j, k), sign in EPS.items():
        for a in range(3):
            for d, b, f1 in theory.fby[a]:
                for e, c, f2 in theory.fby[b]:
                    term = word(
                        L("g1", a),
                        L("c", d, jet(i)),
                        L("c", e, jet(j)),
                        L("g2", c, jet(k)),
                    )
                    primitive = add(primitive, scale(term, sign * f1 * f2))
    assert not theory.Q0(quartic)
    obstruction = theory.Q1(quartic)
    assert obstruction and obstruction == scale(theory.Q0(primitive), 2)
    assert not theory.Q1(primitive)
    return len(obstruction), len(primitive)


def check_matrices():
    x = sp.diag(1, -1, 0)
    y = sp.diag(1, 1, -2)
    ghosts = []
    for i, j in ((0, 1), (1, 2), (2, 0)):
        c = sp.zeros(3)
        c[i, j] = 1
        ghosts.append(c)
    pieces = []
    for seq in permutations(ghosts):
        result = y
        for c in reversed(seq):
            result = c * result - result * c
        pieces.append(sp.trace(x * result))
    assert pieces == [6, 6, 3, 0, 3, 0]
    assert 2 * sum(pieces) == 36
    assert sp.trace(x * x) * sp.trace(y * y) - sp.trace(x * y) ** 2 == 12
    return pieces


if __name__ == "__main__":
    print("nonzero terms in reduced box and primitive", *check_components())
    print("sl3 probe nested brackets", check_matrices())
