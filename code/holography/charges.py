"""LLC charge map and the 1/16-BPS index of Kinney-Maldacena-Minwalla-Raju (hep-th/0510251).

LLC (sources/LLC_main.tex, eq-E-sector, eq-E-charges): sector (theta1,theta2,theta3; n+, n-)_Y,
    Theta = theta1+theta2+theta3, n = n+ + n-,
    E = (Theta+Y)/2 + n,  R_I = theta_I - (Theta-Y)/2,  J_L = (Theta-Y+n)/2,  J_R = (n+ - n-)/2.
KMMR index (their (4.2)):  I = Tr (-1)^F t^{2(E+j1)} y^{2 j2} v^{R2} w^{R3}, with SU(4) Dynkin labels
    [R1,R2,R3] = [q2+q3, q1-q2, q2-q3] in terms of SO(6) Cartans q (their Table 2), j1 = J_L, j2 = J_R,
    q = (R_1, R_2, R_3) of LLC (gamma^I = phibar_I has q = e_I).  Hence the monomial of a sector is
    (-1)^F t^{2 Theta + 3n} y^{n+ - n-} v^{theta1-theta2} w^{theta2-theta3},  (-1)^F = (-1)^{Y+Theta}
    (letter parity: gamma even, beta odd, b even, dc odd; = form parity on C^{2|3}).
Series are dicts {(t, y, v, w) exponents: Fraction}, truncated at t-degree <= tmax.
"""
from fractions import Fraction

from .linalg import vadd


def llc_charges(theta, npm, Y):
    th1, th2, th3 = theta
    Th = th1 + th2 + th3
    n = npm[0] + npm[1]
    E = Fraction(Th + Y, 2) + n
    R = tuple(Fraction(x) - Fraction(Th - Y, 2) for x in theta)
    JL = Fraction(Th - Y + n, 2)
    JR = Fraction(npm[0] - npm[1], 2)
    return E, R, JL, JR


def kmmr_monomial(theta, npm, Y):
    """(sign, (t, y, v, w) exponents) of a state in sector (theta; n+, n-)_Y."""
    E, R, JL, JR = llc_charges(theta, npm, Y)
    assert E == R[0] + R[1] + R[2] + 2 * JL      # BPS relation
    q1, q2, q3 = R
    D1, D2, D3 = q2 + q3, q1 - q2, q2 - q3
    tex = 2 * (E + JL)
    for x in (tex, 2 * JR, D2, D3):
        assert x.denominator == 1
    Th = sum(theta)
    sign = (-1) ** ((Y + Th) % 2)
    return sign, (int(tex), int(2 * JR), int(D2), int(D3))


# ------------------------------------------------------------ truncated series in (t; y, v, w)
def smul(a, b, tmax):
    out = {}
    for ka, ca in a.items():
        for kb, cb in b.items():
            if ka[0] + kb[0] > tmax:
                continue
            k = tuple(x + y for x, y in zip(ka, kb))
            vadd(out, {k: ca * cb})
    return out


def geom(mono, tmax):
    """mono/(1-mono) for a monomial key with positive t-degree."""
    out = {}
    j = 1
    while mono[0] * j <= tmax:
        out[tuple(x * j for x in mono)] = Fraction(1)
        j += 1
    return out


def inv_one_minus(mono, tmax):
    out = geom(mono, tmax)
    out[(0, 0, 0, 0)] = Fraction(1)
    return out


SCALARS = [(2, 0, 1, 0), (2, 0, 0, -1), (2, 0, -1, 1)]     # t^2 v, t^2/w, t^2 w/v
DERIVS = [(3, 1, 0, 0), (3, -1, 0, 0)]                     # t^3 y, t^3/y


def kmmr_Isp(tmax):
    """KMMR (4.14) = (4.8): sum_scalars u/(1-u) - sum_derivs u/(1-u)."""
    out = {}
    for m in SCALARS:
        vadd(out, geom(m, tmax))
    for m in DERIVS:
        vadd(out, geom(m, tmax), -1)
    return out


def kmmr_f(tmax):
    """KMMR (4.2): single-letter index."""
    num = {}
    for m in SCALARS:
        vadd(num, {m: Fraction(1)})
    for m in DERIVS:
        vadd(num, {m: Fraction(1)}, -1)
    for m in [(4, 0, 0, 1), (4, 0, -1, 0), (4, 0, 1, -1)]:  # t^4 (w + 1/v + v/w)
        vadd(num, {m: Fraction(1)}, -1)
    vadd(num, {(6, 0, 0, 0): Fraction(2)})
    den = smul(inv_one_minus(DERIVS[0], tmax), inv_one_minus(DERIVS[1], tmax), tmax)
    return smul(num, den, tmax)


def sector_list(tmax):
    """All weights beta = (n+, n-; theta) with 2 Theta + 3 n <= tmax."""
    out = []
    for n in range(tmax // 3 + 1):
        for npl in range(n + 1):
            rest = tmax - 3 * n
            for Th in range(rest // 2 + 1):
                for t1 in range(Th + 1):
                    for t2 in range(Th - t1 + 1):
                        out.append((npl, n - npl, t1, t2, Th - t1 - t2))
    return out
