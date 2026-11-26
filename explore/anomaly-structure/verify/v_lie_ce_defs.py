"""Shared definitions (same as in v_lie_ce.py)."""
import lie_eval as le

THETA133 = (6, [(0, 5), (0, 1), (1, 2), (2, 5), (0, 3), (3, 4), (4, 5)], [1, 2, 3, 4])
THETA223 = (6, [(0, 1), (1, 5), (0, 2), (2, 5), (0, 3), (3, 4), (4, 5)], [1, 2, 3, 4])
BOX = (4, [(0, 1), (1, 2), (2, 3), (3, 0)], [0, 1, 2, 3])


def br2(p, q):
    return p[0] * q[1] - p[1] * q[0]


def det3(a, b, c):
    return (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]))


P133 = lambda ps: -3 * br2(ps[0], ps[1]) * br2(ps[2], ps[3])
P223 = lambda ps: 2 * br2(ps[0], ps[1]) * br2(ps[2], ps[3])
PF = lambda ps: br2(ps[0], ps[1]) * br2(ps[2], ps[3])
PBOX = lambda ps: det3(ps[0], ps[1], ps[2])


def colour(alg):
    W1, s1, _ = le.orientation_rule(alg, *THETA133)
    W2, s2, _ = le.orientation_rule(alg, *THETA223)
    Wb, sb, _ = le.orientation_rule(alg, *BOX)
    return W1, W2, Wb
