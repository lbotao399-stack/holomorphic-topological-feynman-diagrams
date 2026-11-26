"""Explicit-tensor check (lie.py, cochain.py; no ce.py) of the Moyal identities on gl_N:
j_M = (1/2)[p1p2] tr(1{23}),  i_M = (1/8)[p1p2]^2 tr(1[23]):   d j_M = 0,  (1/2){j_M,j_M} = d i_M,
and of the sl_N identity (1/2){D,D} = F_223/(12 N^2) + d eta_N, eta_N = ((N^2+4)/(2N^2)) [p1p2]^2 tr(1[23])."""
import sys
import numpy as np
from fractions import Fraction as Fr
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/deform/code')
import verify_explicit as V
import lie, cochain as cc
from c_obstruction import kernel_poly

P = V.P
mom = V.mom
u = P[0][0] * P[1][1] - P[0][1] * P[1][0]


def frac_cochain_eq(lhs_terms, rhs_terms, g):
    """lhs, rhs: lists of (Fraction, cochain dict mono->int tensor); exact comparison."""
    keys = set()
    for _, c in lhs_terms + rhs_terms:
        keys |= set(c)
    for m in keys:
        L = sum((x * c[m].astype(object) for x, c in lhs_terms if m in c), np.zeros((g.dim,) * 4, dtype=object))
        R = sum((x * c[m].astype(object) for x, c in rhs_terms if m in c), np.zeros((g.dim,) * 4, dtype=object))
        if not all(a == b for a, b in zip(L.reshape(-1), R.reshape(-1))):
            return False
    return True


for g in [lie.gl(2), lie.gl(3)]:
    D = V.tr_tensor(g, [[(0, 1, 2)], [(0, 2, 1)]])
    f = V.tr_tensor(g, [[(0, 1, 2)]]) - V.tr_tensor(g, [[(0, 2, 1)]])
    JD = V.J_cochain(D)                                  # [p1p2] D  (j_M = JD / 2)
    iM8 = cc.alternate(mom, u * u, f, 3)                 # alt([p1p2]^2 f) = 6 [p1p2]^2 f
    dJ = cc.d_CE(mom, g, JD, 3)
    br, den = V.bracket33(g, JD, JD)                     # {JD,JD} = br / den
    di = cc.d_CE(mom, g, iM8, 3)                         # d(6 [p1p2]^2 f) = 48 d i_M
    # (1/2){j_M,j_M} = (1/8){JD,JD} = br/(8 den);  d i_M = di/48
    ok = frac_cochain_eq([(Fr(1, 8 * den), br)], [(Fr(1, 48), di)], g)
    print("%s: d j_M = 0: %s;  (1/2){j_M,j_M} = d i_M (explicit tensors, exact): %s;  {j_M,j_M} != 0: %s" % (
        g.name, cc.is_zero(dJ), ok, bool(br)))

for g in [lie.sl(3)]:
    N = 3
    D = V.tr_tensor(g, [[(0, 1, 2)], [(0, 2, 1)]])
    f = V.tr_tensor(g, [[(0, 1, 2)]]) - V.tr_tensor(g, [[(0, 2, 1)]])
    JD = V.J_cochain(D)
    br, den = V.bracket33(g, JD, JD)
    G223, ext, K = kernel_poly((2, 2, 3))
    W, sc = cc.colour_weight(g, G223.nv, G223.edges, ext)
    br4 = (P[0][0] * P[1][1] - P[0][1] * P[1][0]) * (P[2][0] * P[3][1] - P[2][1] * P[3][0])
    F = cc.alternate(mom, br4, W, 4)                     # sc * F_223
    eta6 = cc.alternate(mom, u * u, f, 3)                # 6 [p1p2]^2 f
    deta = cc.d_CE(mom, g, eta6, 3)
    c_eta = Fr(N * N + 4, 2 * N * N) / 6
    ok = frac_cochain_eq([(Fr(1, 2 * den), br)], [(Fr(1, 12 * N * N * sc), F), (c_eta, deta)], g)
    print("%s: (1/2){D,D} = F_223/(12N^2) + d[((N^2+4)/(2N^2)) [p1p2]^2 tr(1[23])] (explicit tensors, exact): %s" % (g.name, ok))
