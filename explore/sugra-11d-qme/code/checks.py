"""Exact checks of the g^1 mu^9 gamma machinery (flows.py) against brute force (htk.brute, no reduction)."""
import sys, random, time
from fractions import Fraction as Fr
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/sugra-11d-qme/code')
import flows as F
from htk import brute
from htk.graph import theta
from htk.linalg import perm_sign, det

rng = random.Random(20261120)
def rq(a=1, b=9):
    return Fr(rng.randint(-b, b), rng.randint(a, b))
def rvec():
    return [rq() for _ in range(F.H)]
def rsimplex(n):
    a = [Fr(rng.randint(1, 9)) for _ in range(n)]
    s = sum(a); return [x / s for x in a]

def check_vertex_factor():
    """S3 sec-assembly: dL/dmu^i dL/dgamma_j (mu^a mu^b d_a gamma_b) = delta^j_b d_i mu^b - mu^a d_a delta^j_i."""
    L = F.TERMS['b'][1]
    ok = True
    u = rvec(); y = rvec()
    for i in range(5):
        for j in range(5):
            val = Fr(0)
            for a in range(5):
                for b in range(5):
                    letters = [('mu', a, 0), ('mu', b, 1), ('ga', b, 2)]
                    for k in range(2):           # remove mu^i at position k (sign: odd letters before)
                        if letters[k][1] != i: continue
                        s = (-1) ** k
                        if b != j: continue      # remove gamma_j
                        ext = letters[1 - k]
                        # derivative d_a on gamma (contracted) -> y_a; external mu -> u
                        val += s * y[a] * u[ext[1]]
            ref = (u[j] * y[i] if True else 0) - (sum(u[a] * y[a] for a in range(5)) if i == j else 0)
            ok &= (val == ref)
    return ok

def check_currents(pqr, ntest=3):
    G, interior = theta(*pqr)
    ok = True
    for _ in range(ntest):
        lam = {v: rvec() for v in range(G.nv - 1)}
        x = [rsimplex(n) for n in pqr]; w = rsimplex(3)
        yp = F.poly_currents(pqr, interior, x, w, lam)
        t = F.tmap(pqr, x, w)
        yb = brute.currents(G, t, [lam[v] for v in range(G.nv - 1)])
        ok &= (yp == yb)
    return ok

def jac_det(pqr, x, w):
    """det[t; dt/dxi], xi = (x^A_2.., x^B_2.., x^C_2.., w_B, w_C), x^P_1 = 1 - sum, w_A = 1 - w_B - w_C."""
    import sympy as sp
    syms = []; X = []
    for P, n in enumerate(pqr):
        s = sp.symbols('x%d_2:%d' % (P, n + 1)) if n > 1 else ()
        syms += list(s); X.append([1 - sum(s)] + list(s))
    wB, wC = sp.symbols('wB wC'); syms += [wB, wC]
    W = [1 - wB - wC, wB, wC]
    t = [X[P][j] / W[P] for P in range(3) for j in range(pqr[P])]
    sub = {}
    for P, n in enumerate(pqr):
        for j in range(1, n):
            sub[sp.Symbol('x%d_%d' % (P, j + 1))] = sp.Rational(x[P][j].numerator, x[P][j].denominator)
    sub[wB] = sp.Rational(w[1].numerator, w[1].denominator); sub[wC] = sp.Rational(w[2].numerator, w[2].denominator)
    rows = [[ti.subs(sub) for ti in t]] + [[sp.diff(ti, s).subs(sub) for ti in t] for s in syms]
    d = sp.Matrix(rows).det()
    return Fr(int(sp.numer(d)), int(sp.denom(d)))

def check_density(pqr, ntest=2):
    """D(t) |det[t; dt/dxi]| = -(D_theta/2) prod_P w_P^{11/2-n_P}  (K^S3 = -K^htk, S3 eq-Ktheta)."""
    G, interior = theta(*pqr)
    ok = True
    for _ in range(ntest):
        lam = [rvec() for v in range(G.nv - 1)]
        x = [rsimplex(n) for n in pqr]; w = rsimplex(3)
        t = F.tmap(pqr, x, w)
        Dn = brute.density(G, t, lam, T=1)
        J = abs(jac_det(pqr, x, w))
        lhs = Dn.coeff * J                       # times radicand^{-1/2}
        Dth = brute.theta_D(lam, interior)
        # rhs = -(Dth/2) prod w^{11/2-n} ; compare squares and signs: lhs^2/U == rhs^2
        rhs_sq = (Dth / 2) ** 2 * w[0] ** (11 - 2 * pqr[0]) * w[1] ** (11 - 2 * pqr[1]) * w[2] ** (11 - 2 * pqr[2])
        ok &= (lhs ** 2 / Dn.radicand == rhs_sq) and ((lhs > 0) == (-Dth / 2 > 0)) and Dn.expo == 0
    return ok

def ext_data(G, typ, fl, lam_all):
    """External letter data at every vertex for the flow (one letter per interior vertex)."""
    ext = {}
    for v in range(G.nv):
        ext[v] = dict(u=rvec(), vg=rvec(), lam=lam_all[v])
    return ext

def check_flow_pointwise(pqr, which, ntest=1):
    """Flow integrand from polynomial currents == from Kirchhoff currents at t(x,w); edge-order and
    vertex-order sign rules of W_Gamma."""
    G, interior = theta(*pqr)
    fls = F.m1_flows(G)
    fl, typ = fls[which]
    lam = {v: rvec() for v in range(G.nv - 1)}; lam[0] = [Fr(0)] * F.H     # lambda_N = 0 (no letter at N)
    lam_root = [-sum(lam[v][i] for v in range(G.nv - 1)) for i in range(F.H)]
    lam_all = dict(lam); lam_all[G.root] = lam_root
    ext = ext_data(G, typ, fl, lam_all)
    x = [rsimplex(n) for n in pqr]; w = rsimplex(3)
    t0 = time.time()
    pats = F.patterns(G, fl, typ)
    yp = F.poly_currents(pqr, interior, x, w, lam)
    yb = brute.currents(G, F.tmap(pqr, x, w), [lam[v] for v in range(G.nv - 1)])
    v1 = F.flow_integrand(G, interior, fl, typ, pats, yp, ext)
    v2 = F.flow_integrand(G, interior, fl, typ, pats, yb, ext)
    # edge transposition (0 1): W_Gamma -> -W_Gamma
    order = list(range(G.ne)); order[0], order[1] = order[1], order[0]
    pats_t = F.patterns(G, fl, typ, order=order)
    v3 = F.flow_integrand(G, interior, fl, typ, pats_t, yp, ext)
    # vertex order reversed: unchanged (W even)
    pats_v = F.patterns(G, fl, typ, vorder=list(reversed(range(G.nv))))
    v4 = F.flow_integrand(G, interior, fl, typ, pats_v, yp, ext)
    return dict(npat=len(pats[0]), same_currents=(yp == yb), poly_eq_kirchhoff=(v1 == v2),
                edge_swap=(v3 == -v1), vertex_order=(v4 == v1), nonzero=(v1 != 0), sec=round(time.time() - t0, 1))

if __name__ == '__main__':
    print('vertex factor of S3 sec-assembly reproduced:', check_vertex_factor())
    for pqr in F.PQR:
        print(pqr, 'polynomial currents == Kirchhoff (exact, 3 points):', check_currents(pqr))
    for pqr in F.PQR:
        print(pqr, 'reduced density == brute density x Jacobian (sign, square; 2 points):', check_density(pqr))
    for pqr in F.PQR:
        G, interior = theta(*pqr)
        fls = F.m1_flows(G)
        roles = []
        for fl, typ in fls:
            c = [v for v in typ if typ[v] == 'c'][0]
            roles.append('c@N' if c == 0 else ("c@N'" if c == G.root else 'c@interior'))
        from collections import Counter
        print(pqr, 'm=1 flows:', len(fls), dict(Counter(roles)))
