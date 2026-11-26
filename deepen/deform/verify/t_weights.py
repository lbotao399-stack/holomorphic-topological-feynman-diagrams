"""Check own W_133, W_223 against run-1 c4 W_tot formulas (gl_N, sl_N) and c4 sl_N individual forms."""
from fractions import Fraction as Fr
from vce import perm_from_cycles
from vw import W133, W223


def P(*cycles):
    return perm_from_cycles([tuple(x - 1 for x in c) for c in cycles], 4)


def lin(*pairs):
    out = {}
    for s, w in pairs:
        for k, v in w.items():
            out[k] = out.get(k, 0) + s * v
    return {k: v for k, v in out.items() if v != 0}


def T(coef_by_npow, *cycles):
    sig = P(*cycles)
    return {(sig, n): Fr(c) for n, c in coef_by_npow.items() if c != 0}


def add_all(*ws):
    return lin(*[(1, w) for w in ws])


# c4 gl_N W_tot
c4gl = add_all(
    T({0: -10}, (1, 2, 3, 4)), T({0: -10, 2: -3}, (1, 2, 4, 3)), T({0: 20}, (1, 3, 2, 4)),
    T({0: -10, 2: -3}, (1, 3, 4, 2)), T({0: 20}, (1, 4, 2, 3)), T({0: -10}, (1, 4, 3, 2)),
    T({1: -26}, (1, 2), (3, 4)), T({1: 4}, (1, 3), (2, 4)), T({1: 4}, (1, 4), (2, 3)),
    T({1: 3}, (1,), (2, 3, 4)), T({1: 3}, (1,), (2, 4, 3)), T({1: 3}, (1, 2, 3), (4,)), T({1: 3}, (1, 2, 4), (3,)),
    T({1: 3}, (1, 3, 2), (4,)), T({1: 3}, (1, 3, 4), (2,)), T({1: 3}, (1, 4, 2), (3,)), T({1: 3}, (1, 4, 3), (2,)),
    T({0: 20}, (1,), (2,), (3, 4)), T({0: 20}, (1, 2), (3,), (4,)),
    T({0: -10}, (1,), (2, 3), (4,)), T({0: -10}, (1,), (2, 4), (3,)), T({0: -10}, (1, 3), (2,), (4,)),
    T({0: -10}, (1, 4), (2,), (3,)),
)
c4sl = add_all(
    T({0: -10}, (1, 2, 3, 4)), T({0: -10, 2: -3}, (1, 2, 4, 3)), T({0: 20}, (1, 3, 2, 4)),
    T({0: -10, 2: -3}, (1, 3, 4, 2)), T({0: 20}, (1, 4, 2, 3)), T({0: -10}, (1, 4, 3, 2)),
    T({1: -26}, (1, 2), (3, 4)), T({1: 4}, (1, 3), (2, 4)), T({1: 4}, (1, 4), (2, 3)),
)
c4sl133 = add_all(T({0: 2}, (1, 2, 3, 4)), T({0: 2, 2: 1}, (1, 2, 4, 3)), T({0: -4}, (1, 3, 2, 4)),
                  T({0: 2, 2: 1}, (1, 3, 4, 2)), T({0: -4}, (1, 4, 2, 3)), T({0: 2}, (1, 4, 3, 2)),
                  T({1: 6}, (1, 2), (3, 4)))
c4sl223 = add_all(T({0: -2}, (1, 2, 3, 4)), T({0: -2}, (1, 2, 4, 3)), T({0: 4}, (1, 3, 2, 4)),
                  T({0: -2}, (1, 3, 4, 2)), T({0: 4}, (1, 4, 2, 3)), T({0: -2}, (1, 4, 3, 2)),
                  T({1: -4}, (1, 2), (3, 4)), T({1: 2}, (1, 3), (2, 4)), T({1: 2}, (1, 4), (2, 3)))

for sl, ref in ((False, c4gl), (True, c4sl)):
    w133 = W133(sl)
    w223 = W223(sl)
    tot = lin((-3, w133), (2, w223))
    print("sl" if sl else "gl", "own W_tot == c4 W_tot:", lin((1, tot), (-1, ref)) == {})
    if sl:
        print("sl own W_133 == c4:", lin((1, w133), (-1, c4sl133)) == {})
        print("sl own W_223 == c4:", lin((1, w223), (-1, c4sl223)) == {})
    print("  N-powers W_133:", sorted({n for (_, n) in w133}), " W_223:", sorted({n for (_, n) in w223}))
