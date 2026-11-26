"""Native contractions on differential-form components before color traces."""

from itertools import product
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "s2code"))
from ops import L
from salg import TH1, TH2, add, dpar, gen, lderiv, mul, scale, sderiv

ETA = ((0, 3, 0, 0, 0), (0, 4, 0, 0, 0))
DZ = tuple((0, 5 + i, 0, 0, 0) for i in range(3))
DTHETA = ((5, 1, 0, 0, 0), (5, 2, 0, 0, 0))
JETS = ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))


def field(color):
    eta1, eta2 = map(gen, ETA)
    return add(add(L("c", color), mul(eta1, L("g1", color))),
               add(mul(eta2, L("g2", color)), mul(mul(eta1, eta2), L("b", color))))


def form(color):
    C = field(color)
    result = {}
    # dpar acts only on C, which contains source eta and component fields.
    for i, dz in enumerate(DZ):
        result = add(result, mul(gen(dz), dpar(C, i)))
    for dtheta, eta in zip(DTHETA, ETA):
        result = add(result, mul(gen(dtheta), lderiv(C, eta)))
    return result


def root(poly, color, jet):
    result = {}
    for t in (1, 2, 3, 4):
        result = add(result, sderiv(poly, (t, color) + jet))
    return result


def check():
    delta = mul(add(gen(TH1), gen(ETA[0]), -1), add(gen(TH2), gen(ETA[1]), -1))
    contracted = [root(form(0), 0, jet) for jet in JETS]
    odd_slot = {}
    for dtheta, eta in zip(DTHETA, ETA):
        odd_slot = add(odd_slot, mul(gen(dtheta), lderiv(delta, eta)), -1)
    assert contracted == [odd_slot] + [scale(mul(gen(dz), delta), -1) for dz in DZ]
    assert all(contracted)
    assert all(not mul(left, right) for left, right in product(contracted, repeat=2))
    # Differentiate the actual three-factor product. The unhit middle factor
    # exercises all graded spectator signs for each ordered root assignment.
    cubic = mul(mul(form(0), form(1)), form(2))
    assert cubic
    count = 0
    for i, j in ((0, 1), (0, 2), (1, 0), (1, 2), (2, 0), (2, 1)):
        for first, second in product(JETS, repeat=2):
            assert not root(root(cubic, i, first), j, second)
            count += 1
    return count


if __name__ == "__main__":
    print(f"{check()} ordered cubic root contractions vanish with independent form variables")
