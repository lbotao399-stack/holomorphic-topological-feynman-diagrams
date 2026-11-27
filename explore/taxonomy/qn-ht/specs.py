"""Theory registry with gradings for the sector decomposition of H(Q_0) and d_2 = [Q_2^hat].
spec fields: build() -> Theory; spin(n, i) -> H-tuple of Fractions (letter weight = spin + alpha);
ghost(n, i) -> int (Z-grading, Q raises it by gstep) or None (Z/2 by parity);
charge(n, i) -> tuple (gauge Cartan weights and flavour charges; conserved by Q_0, Q_2)."""
from fractions import Fraction as Fr
import hteng
import theories as T

F2 = Fr(1, 2)
HADJ = [2, -2, 0]  # sl2 Cartan (h) weights of e, f, h
HF = [1, -1]


def sl2_killing():
    """trace form in the fundamental: k(e,f) = k(f,e) = 1, k(h,h) = 2"""
    return {(0, 1): 1, (1, 0): 1, (2, 2): 2}


def W_N8_sl2():
    """W = k(Phi, [X, Y]) = Phi^a X^b Y^c f_bc^d k_da, multiplets 0,1,2 = Phi, X, Y"""
    k = sl2_killing()
    W = {}
    for (b, c, d), v in T.F_SL2.items():
        for a in range(3):
            kv = k.get((d, a), 0)
            if kv:
                key = ((0, a), (1, b), (2, c))
                W[key] = W.get(key, 0) + v * kv
    return {kk: vv for kk, vv in W.items() if vv}


def su2_N8_sl2():
    adj = T.adj_rho(T.F_SL2, 3)
    th = hteng.Theory(1, 3, T.F_SL2, bf=[(3, adj)] * 3, W=W_N8_sl2(), name="3d N=8 SU(2) (sl2)")
    return th


def su2_N4_pure_sl2():
    return hteng.Theory(1, 3, T.F_SL2, bf=[(3, T.adj_rho(T.F_SL2, 3))], name="3d N=4 SU(2) pure (sl2)")


# ------------------------------------------------------------------ 3d
def spec3d(build, rch, cartan):
    """rch[k]: R-charge of phi_k; cartan(k, i) -> tuple of conserved charges of phi_k^i (psi opposite)."""
    def spin(n, i):
        if n == 'c':
            return (Fr(0),)
        if n == 'B':
            return (Fr(1),)
        k = int(n[1:])
        return (rch[k] / 2,) if n[0] == 'p' else (1 - rch[k] / 2,)
    S = 1
    from math import lcm
    for r in rch:
        S = lcm(S, (r / 2).denominator, r.denominator)
    def ghost(n, i):
        if n == 'c':
            return S
        if n == 'B':
            return 0
        k = int(n[1:])
        return int(S * (rch[k] if n[0] == 'p' else 1 - rch[k]))
    return dict(build=build, H=1, spin=spin, ghost=ghost, gstep=S, charge=None, cartan=cartan)


def _gauge_charge_sl2(n, i):
    if n == 'c':
        return HADJ[i]
    if n == 'B':
        return -HADJ[i]
    return None


def make_charge(gauge_ch, mat_ch):
    """gauge_ch(n, i) for c, B; mat_ch(k, i) -> tuple for phi_k^i (psi_k^i: negative)"""
    def charge(n, i):
        if n in ('c', 'B'):
            g = gauge_ch(n, i)
            return (g,) + mat_ch(None, None, zero=True)[1:]
        k = int(n[1:])
        t = mat_ch(k, i)
        return t if n[0] in 'pu' else tuple(-x for x in t)
    return charge


SPECS = {}

# 3d N=4 SQED with Nf hypers (charge q): multiplets 0 = Phi, 1..Nf = Q_i, Nf+1..2Nf = Qt_i
def _sqed(Nf, q=1):
    rch = [Fr(1)] + [F2] * (2 * Nf)
    def mat(k, i, zero=False):
        if zero:
            return (0,) * (Nf + 3)
        if k == 0:
            return (0, -2) + (0,) * (Nf + 1)
        if k <= Nf:
            return (q, 1) + tuple(1 if j == k - 1 else 0 for j in range(Nf)) + (0,)
        return (-q, 1) + tuple(-1 if j == k - 1 - Nf else 0 for j in range(Nf)) + (0,)
    ch = make_charge(lambda n, i: 0, mat)
    sp = spec3d(lambda: T.sqed_3d(Nf, q), rch, None)
    sp['charge'] = ch
    return sp

SPECS["3d_sqed1"] = _sqed(1)
SPECS["3d_sqed2"] = _sqed(2)
SPECS["3d_sqed1_q2"] = _sqed(1, 2)


def _su2N4(Nf):
    rch = [Fr(1)] + [F2] * (2 * Nf)
    def mat(k, i, zero=False):
        if zero:
            return (0, 0) + (0,) * Nf
        if k == 0:
            return (HADJ[i], -2) + (0,) * Nf
        if k <= Nf:
            return (HF[i], 1) + tuple(1 if j == k - 1 else 0 for j in range(Nf))
        return (-HF[i], 1) + tuple(-1 if j == k - 1 - Nf else 0 for j in range(Nf))
    sp = spec3d(lambda: T.su2_N4_3d(Nf), rch, None)
    sp['charge'] = make_charge(_gauge_charge_sl2, mat)
    return sp

SPECS["3d_su2N4_1"] = _su2N4(1)
SPECS["3d_su2N4_2"] = _su2N4(2)


def _su2pure4():
    rch = [Fr(1)]
    def mat(k, i, zero=False):
        if zero:
            return (0, 0)
        return (HADJ[i], 1)
    sp = spec3d(su2_N4_pure_sl2, rch, None)
    sp['charge'] = make_charge(_gauge_charge_sl2, mat)
    return sp

SPECS["3d_su2N4_pure"] = _su2pure4()


def _su2N8():
    rch = [Fr(1), F2, F2]  # W = Phi [X,Y]: R(Phi) + R(X) + R(Y) = 2  (ESW 11.8: Phi in K, X,Y in K^0 is r = (2,0,0))
    def mat(k, i, zero=False):
        if zero:
            return (0, 0, 0)
        # U(1)^2 flavour preserving W: Phi (-1,-1), X (1,0), Y (0,1)
        f = [(-1, -1), (1, 0), (0, 1)][k]
        return (HADJ[i],) + f
    sp = spec3d(su2_N8_sl2, rch, None)
    sp['charge'] = make_charge(_gauge_charge_sl2, mat)
    return sp

SPECS["3d_su2N8"] = _su2N8()


# ------------------------------------------------------------------ 5d, 7d
def spec_odd(build, H, gauge, mats, zgrade=True):
    """gauge in {'sl2','u1'}; mats: list over BF/sym multiplets of (kind, rep) with kind 'fund','adj',
    ('u1', q) for BF phi, or 'half' for the symplectic half-hyper.  Spins (9.2 form): phi, psi, u: 1/2."""
    D = H + 1
    one = (Fr(1),) * H
    half = (F2,) * H
    zero = (Fr(0),) * H
    nm = len(mats)
    def spin(n, i):
        if n == 'c':
            return zero
        if n == 'B':
            return one
        return half
    def ghost(n, i):
        if not zgrade:
            return None
        if n == 'c':
            return 1
        if n == 'B':
            return D - 2
        if n[0] == 'p':
            return 0
        if n[0] == 's':
            return D - 1
        return None
    def charge(n, i):
        if n in ('c', 'B'):
            g = (_gauge_charge_sl2(n, i),) if gauge == 'sl2' else (0,)
            return g + (0,) * nm
        k = int(n[1:])
        kind = mats[k]
        if kind == 'fund':
            g = HF[i]
        elif kind == 'adj':
            g = HADJ[i]
        elif kind == 'half':
            g = HF[i]
        else:
            g = kind[1]
        fl = tuple(1 if j == k else 0 for j in range(nm)) if kind != 'half' else (0,) * nm
        t = (g,) + fl
        if n[0] == 's':
            t = tuple(-x for x in t)
        return t
    return dict(build=build, H=H, spin=spin, ghost=ghost, gstep=1, charge=charge)


def _bf(H, gauge, reps, name):
    def build():
        if gauge == 'sl2':
            return hteng.Theory(H, 3, T.F_SL2, bf=reps, name=name)
        return hteng.Theory(H, 1, {}, bf=reps, name=name)
    return build

SPECS["5d_pure"] = spec_odd(_bf(2, 'sl2', [], "5d N=1 SU(2) pure (BF Bg)"), 2, 'sl2', [])
SPECS["5d_fund1"] = spec_odd(_bf(2, 'sl2', [(2, T.FUND_SL2)], "5d N=1 SU(2) + 1 hyper (BF R/g)"), 2, 'sl2', ['fund'])
SPECS["5d_N2"] = spec_odd(_bf(2, 'sl2', [(3, T.adj_rho(T.F_SL2, 3))], "5d N=2 SU(2) (BF g/g)"), 2, 'sl2', ['adj'])
SPECS["5d_u1q1"] = spec_odd(_bf(2, 'u1', T.u1_reps([1]), "5d N=1 U(1) + 1 hyper q=1 (BF R/g)"), 2, 'u1', [('u1', 1)])
SPECS["5d_u1q2"] = spec_odd(_bf(2, 'u1', T.u1_reps([2]), "5d N=1 U(1) + 1 hyper q=2 (BF R/g)"), 2, 'u1', [('u1', 2)])
SPECS["5d_half"] = spec_odd(lambda: T.sym_theory(2, 'sl2', [T.half_hyper_sl2()], "5d N=1 SU(2) + half-hyper (U//g)"),
                            2, 'sl2', ['half'], zgrade=False)
SPECS["7d_pure"] = spec_odd(_bf(3, 'sl2', [], "7d N=1 SU(2) (BF Bg)"), 3, 'sl2', [])


# positive control (3dN2 REPORT Sec. 4): U(1)_0 + one chiral q = 1, r = 1/2; lifted classes at J = 3, 4
def _u1chiral():
    def mat(k, i, zero=False):
        if zero:
            return (0, 0)
        return (1, 1)
    sp = spec3d(lambda: T.u1_3d([1]), [F2], None)
    sp['charge'] = make_charge(lambda n, i: 0, mat)
    return sp

SPECS["3d_u1chiral"] = _u1chiral()
