"""Normalization of the C^4 graph integrals (H = 4, T = 0) and of the twisted action.

(a) Omega = dx^1 ^ dx^2 ^ dx^3 ^ dx^4:  Omega ^ dxbar^1 ^ .. ^ dxbar^4 = 16 d^8 x_E  (exterior algebra, exact).
(b) Edge form on C^d (S2 eq-y generalized):  [K_t] - dt ^ dbar^*[K_t] = (2 pi)^{-d} e^{-y.x} d^d y,
    [K_t] = e^{-x.xbar/2t} dxbar^{1..d} / (2^d (2 pi t)^d),  dbar^* = -2 sum_a iota_{d/dxbar^a} d_a,  y = xbar/(2t);
    checked as an identity of forms in (t, xbar^a) with x^a as parameters, d = 2, 3, 4.
(c) Sign s_d(n) of  G = s_d(n) a^m b^{k-1} K_Gamma  (derivation chain in the docstring of sign_total), compared
    with S1 (d = 2) and S2 (d = 3); varsigma^ext for the C^4 graphs by an own Jacobian determinant and by htk.
(d) Superfield identity  int dtheta (1/3!) f_{ABC} C^A C^B C^C = (1/2) f_{ABC} gamma^A c^B c^C  (C = c + theta gamma).
"""
import sys, random, itertools
from fractions import Fraction as Fr
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/msym-D8/code')
import sympy as sp
from htk import conventions, brute as hbrute
from htk.graph import Graph
from graphs_c4 import pentagon, theta_O, THETAS
import brute_c4 as bf

LOG = open((__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/msym-D8/code/norm_c4.log', 'w')


def log(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    LOG.write(s + '\n')


# ------------------------------------------------------------------ minimal exterior algebra
class Form:
    """sum coeff * e_{i1} ^ ... ^ e_{ik} (indices sorted), coefficients sympy."""
    def __init__(self, d=None):
        self.d = dict(d or {})

    @staticmethod
    def one(i):
        return Form({(i,): sp.Integer(1)})

    @staticmethod
    def scalar(c):
        return Form({(): sp.sympify(c)})

    def __add__(self, o):
        r = dict(self.d)
        for k, v in o.d.items():
            r[k] = sp.expand(r.get(k, 0) + v)
        return Form({k: v for k, v in r.items() if v != 0})

    def __sub__(self, o):
        return self + o * Form.scalar(-1)

    def __mul__(self, o):
        r = {}
        for k1, v1 in self.d.items():
            for k2, v2 in o.d.items():
                if set(k1) & set(k2):
                    continue
                seq = list(k1) + list(k2)
                inv = sum(1 for i in range(len(seq)) for j in range(i + 1, len(seq)) if seq[i] > seq[j])
                key = tuple(sorted(seq))
                r[key] = sp.expand(r.get(key, 0) + (-1) ** inv * v1 * v2)
        return Form({k: v for k, v in r.items() if v != 0})

    def simplify(self):
        return Form({k: sp.simplify(v) for k, v in self.d.items() if sp.simplify(v) != 0})


def wedge(*fs):
    r = Form.scalar(1)
    for f in fs:
        r = r * f
    return r


def check_a():
    # basis: dxE^1..dxE^8 = indices 0..7 ; x^a = xE^{2a-1} + i xE^{2a}
    dx = [Form.one(2 * a) + Form.one(2 * a + 1) * Form.scalar(sp.I) for a in range(4)]
    dxb = [Form.one(2 * a) - Form.one(2 * a + 1) * Form.scalar(sp.I) for a in range(4)]
    vol = wedge(*dx, *dxb)
    log('(a) dx^1^..^dx^4^dxbar^1^..^dxbar^4 =', vol.d, ' -> Omega = dx^{1234}: Omega ^ dxbar^{1234} = 16 d^8x_E:',
        vol.d == {tuple(range(8)): 16})
    for d in (2, 3):
        dx = [Form.one(2 * a) + Form.one(2 * a + 1) * Form.scalar(sp.I) for a in range(d)]
        dxb = [Form.one(2 * a) - Form.one(2 * a + 1) * Form.scalar(sp.I) for a in range(d)]
        log('    d = %d: dx^{1..d} ^ dxbar^{1..d} =' % d, wedge(*dx, *dxb).d)


def check_b(d):
    # coordinates: t (index 0), xbar^a (indices 1..d); x^a are parameters
    t = sp.Symbol('t', positive=True)
    xb = sp.symbols('xb1:%d' % (d + 1))
    x = sp.symbols('x1:%d' % (d + 1))
    E = sp.exp(-sum(x[a] * xb[a] for a in range(d)) / (2 * t))
    dt = Form.one(0)
    dxb = [Form.one(a + 1) for a in range(d)]
    K = Form.scalar(E / (2 ** d * (2 * sp.pi * t) ** d)) * wedge(*dxb)
    # dbar^* K = -2 sum_a (d/dx^a coefficient) iota_{d/dxbar^a} dxbar^{1..d}
    coefK = E / (2 ** d * (2 * sp.pi * t) ** d)
    dstar = Form()
    for a in range(d):
        rest = wedge(*[dxb[b] for b in range(d) if b != a])
        dstar = dstar + Form.scalar(-2 * sp.diff(coefK, x[a]) * (-1) ** a) * rest
    lhs = K - dt * dstar
    y = [xb[a] / (2 * t) for a in range(d)]
    dy = [Form.scalar(1 / (2 * t)) * dxb[a] + Form.scalar(-xb[a] / (2 * t ** 2)) * dt for a in range(d)]
    rhs = Form.scalar(E / (2 * sp.pi) ** d) * wedge(*dy)
    ok = (lhs - rhs).simplify().d == {}
    log('(b) d = %d: [K_t] - dt ^ dbar^*[K_t] == (2 pi)^{-%d} e^{-y.x} d^%dy:' % (d, d, d), ok)
    return ok


def sign_total(d, n):
    """G = s_d(n) * a^m * b^{k-1} * K_Gamma,  m = d n + 1 edges, k = (d-1) n + 2 vertices,
    a = edge-form factor, b = vertex factor  (int Omega_v ^ dxbar_v^{1..d} e^{x mu} = b delta^d(mu) d^dxbar_v).
    Chain:  (i)  wedge_e P_e = (-1)^{m + (d-1) C(m,2)} int_F wedge_e F_e,  P_e = int dt_e dbar^* K;
            (ii) dbar wedge_e P_e |_UV = (-1)^{(d-1) C(m,2)} int_{RP, omega_m} Phi_{m-1}  (Stokes; facets t_j = eps);
            (iii) wedge_e d^d y_e = prod (2t_e)^{-d} varsigma^ext 2^{dn} vol_X ^ prod psi + ...  (S3 eq-extraction,
                  CONVENTIONS Sec. 4, X = xbar, c_a = 2);  vol_X ^ prod psi = (-1)^{d(k-1)(m-1)} prod psi ^ vol_X;
            (iv) (wedge_v Omega_v) ^ prod_a prod_v dxbar^a_v = (-1)^{C(d,2) C(k-1,2) + d C(k-1,2)} prod_v (Omega_v ^ dxbar_v);
            (v)  vertex integrals b^{k-1}, delta functions 2^{d(k-1)} (prod t / U)^d;  prod (2t)^{-d} 2^{dn} 2^{d(k-1)} (prod t)^d = 1."""
    m, k = d * n + 1, (d - 1) * n + 2
    C2 = lambda x: x * (x - 1) // 2
    ext = (C2(m) * C2(d) + n * d + (k - 1) * n * C2(d)) % 2
    e = ((d - 1) * C2(m) + d * (k - 1) * (m - 1) + C2(d) * C2(k - 1) + d * C2(k - 1) + ext) % 2
    return -1 if e else 1


def check_c():
    log('(c) s_d(n):')
    for d in (2, 3, 4, 5):
        log('    d = %d: s_d(n), n = 1..6:' % d, [sign_total(d, n) for n in range(1, 7)])
    # S1 (d = 2): a = -1/(8 pi^2), b = 2 (2 pi i)^2 = -8 pi^2  (compendium table, S1 Sec. 3.7)
    #   a^m b^{k-1} = (-1)^n (8 pi^2)^{-n};  S1: G = (-8 pi^2)^{-n} I^{S1}, I^{S1} = (-1)^{n(n+1)/2} K
    ok1 = all(sign_total(2, n) * (-1) ** n == (-1) ** n * conventions.sign_n(n) for n in range(1, 9))
    log('    d = 2 against S1 eq-1-7 (G = (-8 pi^2)^{-n} I^{S1}, I^{S1} = (-1)^{n(n+1)/2} K): agree for n = 1..8:', ok1)
    # S2 (d = 3): a = 1/(8 pi^3), b = 8 pi^3;  S2 eq-B-norm: G = (-1/(8 pi^3))^n I^{S2}, I^{S2} = (-1)^{n(n+1)/2} K
    rel = [sign_total(3, n) * (-1) ** n * conventions.sign_n(n) for n in range(1, 9)]
    log('    d = 3 against S2 eq-B-norm: ratio (this derivation)/(S2), n = 1..8:', rel)
    log('    d = 4: a = 1/(16 pi^4), b = (2 pi)^4 = 16 pi^4 (c_Omega i^4 = 1):  G = s_4(n) (16 pi^4)^{-n} K, s_4(n) =',
        [sign_total(4, n) for n in range(1, 9)])
    # varsigma^ext for the C^4 graphs: formula, htk, own determinant
    graphs = [('pent', pentagon()[0])] + [(pqr, theta_O(*pqr)[0]) for pqr in THETAS]
    for name, G in graphs:
        f = conventions.extraction_sign_formula(G, 4)
        e = conventions.extraction_sign(G, 4, c=[2] * 4, k=[2] * 4)
        own = own_extraction_sign(G, 4, G.cycle_matrix())
        log('    varsigma^ext(%s, D = 4): formula %+d, htk %+d, own determinant %+d' % (name, f, e, own))
    # K (positive cycle matrix) versus int wedge_{c in cotree} d^4 y_c (own brute force)
    rng = random.Random(4)
    for name, G in graphs:
        cot = bf.cotree_of(G.nv, G.edges)
        for _ in range(2):
            t = [Fr(rng.randint(1, 30), rng.randint(1, 7)) for _ in range(G.ne)]
            lam = {v: [rng.randint(-6, 6) for _ in range(4)] for v in range(G.nv - 1)}
            Down, _ = bf.density_D(G.nv, G.edges, t, lam, cot)
            Dh = hbrute.density(G, t, [lam[v] for v in range(G.nv - 1)], T=0)
            assert Dh.rational() == Down, (name, Dh, Down)
        log('    %s: own D(t) for wedge_{c in %s} d^4y_c == htk D(t) of K (positive C), 2 random points: True' % (name, cot))


def own_extraction_sign(G, D, C, c=2, kk=2, seed=3):
    """varsigma from  dg ^ wedge_e wedge_a theta^a_e = varsigma c^{nD} (-1)^{dim X} D_psi g vol_X ^ vol_t,
    theta^a_e = dX^a_e - c J^a_e dt_e, J = X_e/(kk t_e), basis (X^a_v a-major, t_e)."""
    rng = random.Random(seed)
    N, g, n = G.ne, G.nv - 1, len(C)
    while True:
        X = [[Fr(rng.randint(-9, 9), rng.randint(1, 5)) for _ in range(g)] for _ in range(D)]
        t = [Fr(rng.randint(1, 20), rng.randint(1, 7)) for _ in range(N)]
        B = [[(1 if v == a else 0) - (1 if v == b else 0) for v in range(g)] for (a, b) in G.edges]
        Xe = [[sum(B[e][v] * X[a][v] for v in range(g)) for e in range(N)] for a in range(D)]
        J = [[Xe[a][e] / (kk * t[e]) for e in range(N)] for a in range(D)]
        rows = [[Fr(0)] * (D * g) + [Fr(1)] * N]
        for e in range(N):
            for a in range(D):
                row = [Fr(0)] * (D * g + N)
                for v in range(g):
                    row[a * g + v] = Fr(B[e][v])
                row[D * g + e] = -c * J[a][e]
                rows.append(row)
        detM = bf.fdet(rows)
        j0 = N - 1
        prow = [[C[r][e] * J[a][e] for e in range(N) if e != j0] for a in range(D) for r in range(n)]
        Dpsi = Fr((-1) ** j0) / t[j0] * bf.fdet(prow)
        if Dpsi != 0:
            break
    val = detM * (-1) ** (D * g) / (Fr(c) ** (n * D) * Dpsi * sum(t))
    assert val in (1, -1), val
    return int(val)


def check_d():
    import ops_c41 as O
    n, f = O.structure_constants('so5')
    lhs, rhs = {}, {}
    for (A, B, Cc), v in f.items():
        lhs = O.add(lhs, O.scale(O.mul(O.mul(O.csuper(A), O.csuper(B)), O.csuper(Cc)), v / 6))
        rhs = O.add(rhs, O.scale(O.mul(O.mul(O.L('g', A), O.L('c', B)), O.L('c', Cc)), v / 2))
    ok = O.theta_coeff(lhs, '1') == rhs
    log('(d) so(5): int dtheta (1/3!) f C C C == (1/2) f gamma c c :', ok)


if __name__ == '__main__':
    check_a()
    for d in (2, 3, 4):
        check_b(d)
    check_c()
    check_d()
