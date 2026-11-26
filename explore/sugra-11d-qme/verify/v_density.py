"""Claim b2 (and the current formula used in b1): own Kirchhoff solver, own exterior algebra, own Jacobian.

(1) Polynomial currents of S3 eq-series / eq-parallel == Kirchhoff currents at t^P_j = x^P_j / w_P.
(2) Omega = U^{-5} rho_2(M) ^ prod_{i=1..5} prod_{r=1,2} psi^i_r = D(t) iota_E dt_0...dt_12,
    rho_2 = U^{-1/2} Pf(-1/2 dM M^{-1} dM);  Pf also compared with S3 eq-rho2.
    g = D(t(xi)) |det[t; dt/dxi]|  compared with  -(D_theta/2) prod_P w_P^{11/2-n_P}  (sign and square).
(3) varsigma^ext of S3 eq-extraction from the 78 x 78 Jacobian (S3 Section 5.3), hence
    varsigma_theta = (-1)^n (-1)^{(5n+1)5n/2} varsigma^ext (S3 eq-master).
Exact rationals throughout.
"""
import sys, os, random, time
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vcommon import theta_graph, PQR, H, gauss_det, kirchhoff, D_theta, perm_parity_sign, rq, rsimplex

rng = random.Random(777)


def poly_currents(pqr, paths, x, w, lam):
    """S3 eq-series, eq-parallel; direction N -> N' along every path."""
    Qs, mh = [], []
    for P, (pv, pe) in enumerate(paths):
        n_ = pqr[P]
        Q = [[Fr(0)] * H]
        for k in range(1, n_):
            Q.append([Q[-1][i] + lam[pv[k - 1]][i] for i in range(H)])
        Qs.append(Q)
        mh.append([sum(x[P][j] * Q[j][i] for j in range(n_)) for i in range(H)])
    S = [lam[0][i] + sum(m[i] for m in mh) for i in range(H)]
    y = {}
    for P, (pv, pe) in enumerate(paths):
        ell = [w[P] * S[i] - mh[P][i] for i in range(H)]
        for j, e in enumerate(pe):
            y[e] = [ell[i] + Qs[P][j][i] for i in range(H)]
    return [y[e] for e in range(len(y))]


def cycle_matrix(ne, paths):
    A, B, C = [pe for pv, pe in paths]
    c1 = [0] * ne; c2 = [0] * ne
    for e in A: c1[e] = 1
    for e in B: c1[e] = -1; c2[e] = 1
    for e in C: c2[e] = -1
    return [c1, c2]


def density_parts(ne, t, y, Cm, e0):
    """Returns (r, U): Omega coefficient on dt_{S}, S = all but e0, equals r * U^{-1/2};
    and D(t) = (-1)^e0 r U^{-1/2} / t_e0.  Also the Pf two-form and the eq-rho2 two-form."""
    n = 2
    M = [[sum(Fr(Cm[a][e] * Cm[b][e]) * t[e] for e in range(ne)) for b in range(n)] for a in range(n)]
    U = M[0][0] * M[1][1] - M[0][1] * M[1][0]
    Minv = [[M[1][1] / U, -M[0][1] / U], [-M[1][0] / U, M[0][0] / U]]
    dM = [[[Cm[a][e] * Cm[b][e] for e in range(ne)] for b in range(n)] for a in range(n)]
    # Pf(-1/2 dM M^{-1} dM) = A_12 = -1/2 sum_ab Minv_ab dM_{1a} ^ dM_{b2}; coefficient on dt_e ^ dt_f (e<f)
    pf = {}
    for e in range(ne):
        for f in range(e + 1, ne):
            s = Fr(0)
            for a in range(n):
                for b in range(n):
                    s += Minv[a][b] * (dM[0][a][e] * dM[b][1][f] - dM[0][a][f] * dM[b][1][e])
            pf[(e, f)] = -s / 2
    # eq-rho2: rho_2 = -(M11 dM12 dM22 - M12 dM11 dM22 + M22 dM11 dM12) / (2 U^{3/2})  => Pf = U^{1/2} rho_2
    def w2(u, v):
        return {(e, f): Fr(u[e] * v[f] - u[f] * v[e]) for e in range(ne) for f in range(e + 1, ne)}
    T1 = w2(dM[0][1], dM[1][1]); T2 = w2(dM[0][0], dM[1][1]); T3 = w2(dM[0][0], dM[0][1])
    rho2 = {k: -(M[0][0] * T1[k] - M[0][1] * T2[k] + M[1][1] * T3[k]) / (2 * U) for k in pf}
    pf_ok = all(pf[k] == rho2[k] for k in pf)
    # psi rows (i outer, r inner)
    Psi = [[Cm[r][e] * y[e][i] for e in range(ne)] for i in range(H) for r in range(n)]
    S = [e for e in range(ne) if e != e0]
    tot = Fr(0)
    for a_ in range(len(S)):
        for b_ in range(a_ + 1, len(S)):
            e, f = S[a_], S[b_]
            if pf[(e, f)] == 0:
                continue
            rest = [g for g in S if g not in (e, f)]
            sg = perm_parity_sign([e, f] + rest)
            tot += pf[(e, f)] * sg * gauss_det([[row[g] for g in rest] for row in Psi])
    r = tot / U ** 5
    return r, U, pf_ok


def jacobian(pqr, x, w):
    """det[t; dt/dxi], xi = (x^A_2.., x^B_2.., x^C_2.., w_B, w_C), t^P_j = x^P_j / w_P."""
    cols = []  # t index list (P, j)
    for P, n_ in enumerate(pqr):
        for j in range(n_):
            cols.append((P, j))
    rows = [[x[P][j] / w[P] for (P, j) in cols]]
    for Pk, n_ in enumerate(pqr):
        for k in range(1, n_):
            rows.append([(Fr(int(j == k) - int(j == 0)) / w[P]) if P == Pk else Fr(0) for (P, j) in cols])
    for Wv in (1, 2):   # w_B, w_C ; w_A = 1 - w_B - w_C
        dw = [Fr(-1), Fr(int(Wv == 1)), Fr(int(Wv == 2))]
        rows.append([-x[P][j] / w[P] ** 2 * dw[P] for (P, j) in cols])
    return gauss_det(rows)


def varsigma_ext(nv, edges, Cm, tvals, tau, zb, e0):
    """Coefficient of [prod_v dtau_v][prod_v dzb^1_v]..[prod_v dzb^5_v] ^ dt_S in
    wedge_e (theta^tau_e theta^1_e .. theta^5_e), theta^tau_e = dtau_e - w_e dt_e, theta^i_e = dzb^i_e - 4 y^i_e dt_e,
    divided by 4^{5n} (prod_r psi^tau_r prod_i prod_r psi^i_r)_S."""
    ne = len(edges); g = nv - 1; n = ne - nv + 1
    S = [e for e in range(ne) if e != e0]
    ncol = 6 * g + len(S)
    rows = []
    wv = []; yv = []
    for e, (a, b) in enumerate(edges):
        ta = tau[a] - tau[b]
        we = ta / (2 * tvals[e])
        ye = [(zb[a][i] - zb[b][i]) / (4 * tvals[e]) for i in range(H)]
        wv.append(we); yv.append(ye)
        for comp in range(6):
            row = [Fr(0)] * ncol
            for v, sgn in ((a, 1), (b, -1)):
                if v < g:
                    row[comp * g + v] += sgn
            coef = we if comp == 0 else 4 * ye[comp - 1]
            if e in S:
                row[6 * g + S.index(e)] = -coef
            rows.append(row)
    lhs = gauss_det(rows)
    Psi = [[Cm[r][e] * wv[e] for e in range(ne)] for r in range(n)]
    Psi += [[Cm[r][e] * yv[e][i] for e in range(ne)] for i in range(H) for r in range(n)]
    rhs = Fr(4) ** (5 * n) * gauss_det([[row[e] for e in S] for row in Psi])
    return lhs / rhs


if __name__ == '__main__':
    t0 = time.time()
    for pqr in PQR:
        nv, edges, paths = theta_graph(*pqr)
        ne = len(edges)
        Cm = cycle_matrix(ne, paths)
        ok_cur = True; ok_dens = True; ok_pf = True; ok_basic = True
        for trial in range(3):
            lam = {v: [rq(rng) for _ in range(H)] for v in range(nv - 1)}   # lambda_N != 0 included
            x = [rsimplex(rng, n_) for n_ in pqr]
            w = rsimplex(rng, 3)
            t = [x[P][j] / w[P] for P in range(3) for j in range(pqr[P])]
            yk = kirchhoff(nv, edges, t, lam)
            yp = poly_currents(pqr, paths, x, w, lam)
            ok_cur &= (yk == yp)
            # KCL / KVL of yk
            for v in range(nv - 1):
                for i in range(H):
                    ok_cur &= (sum(yk[e][i] for e, (a, b) in enumerate(edges) if a == v)
                               - sum(yk[e][i] for e, (a, b) in enumerate(edges) if b == v) == lam[v][i])
            for row in Cm:
                for i in range(H):
                    ok_cur &= (sum(row[e] * t[e] * yk[e][i] for e in range(ne)) == 0)
            r0, U, pfok = density_parts(ne, t, yk, Cm, 0)
            r1, U1, _ = density_parts(ne, t, yk, Cm, ne - 1)
            D0 = (-1) ** 0 * r0 / t[0]            # times U^{-1/2}
            D1 = (-1) ** (ne - 1) * r1 / t[ne - 1]
            ok_basic &= (D0 == D1)
            ok_pf &= pfok
            J = jacobian(pqr, x, w)
            g = D0 * abs(J)                         # times U^{-1/2}
            Dth = D_theta(lam, paths)
            rhs_sq = (Dth / 2) ** 2
            for P in range(3):
                rhs_sq *= w[P] ** (11 - 2 * pqr[P])
            ok_dens &= (g ** 2 / U == rhs_sq) and ((g > 0) == (-Dth > 0)) and Dth != 0
        print(pqr, 'poly currents == own Kirchhoff (+KCL,KVL), 3 pts:', ok_cur,
              '| Pf == eq-rho2:', ok_pf, '| D(t) same for e0=0,12:', ok_basic,
              '| D|J| = -(Dtheta/2) prod w^(11/2-n) (square & sign):', ok_dens)
    print('varsigma^ext (S3 eq-extraction), 2 random integer points x 2 choices of S per graph:')
    for pqr in PQR:
        nv, edges, paths = theta_graph(*pqr)
        Cm = cycle_matrix(len(edges), paths)
        vals = []
        for trial in range(2):
            tv = [Fr(rng.randint(1, 9)) for _ in edges]
            tau = [Fr(rng.randint(-9, 9)) for _ in range(nv - 1)] + [Fr(0)]
            zb = [[Fr(rng.randint(-9, 9)) for _ in range(H)] for _ in range(nv - 1)] + [[Fr(0)] * H]
            for e0 in (0, len(edges) - 1):
                vals.append(varsigma_ext(nv, edges, Cm, tv, tau, zb, e0))
        n = 2
        ext = vals[0]
        vs = (-1) ** n * (-1) ** ((5 * n + 1) * 5 * n // 2) * ext
        print(pqr, 'ratios', [str(v) for v in vals], '-> varsigma^ext =', ext, ', varsigma_theta =', vs)
    print('sec', round(time.time() - t0, 1))


def control_heptagon():
    """Control: S3 states varsigma^ext = -1 for C_7 (edges 1->2->..->7->1, root 7)."""
    nv = 7
    edges = [(i, i + 1) for i in range(6)] + [(6, 0)]
    Cm = [[1] * 7]
    vals = []
    for trial in range(2):
        tv = [Fr(rng.randint(1, 9)) for _ in edges]
        tau = [Fr(rng.randint(-9, 9)) for _ in range(nv - 1)] + [Fr(0)]
        zb = [[Fr(rng.randint(-9, 9)) for _ in range(H)] for _ in range(nv - 1)] + [[Fr(0)] * H]
        for e0 in (0, 6):
            vals.append(varsigma_ext(nv, edges, Cm, tv, tau, zb, e0))
    return vals


if __name__ == '__main__':
    print('control C_7: varsigma^ext ratios', [str(v) for v in control_heptagon()])
