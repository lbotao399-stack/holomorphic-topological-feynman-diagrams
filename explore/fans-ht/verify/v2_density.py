"""Pointwise exact checks of the densities (claims a4, a5, a7, d1 pointwise part).
Brute force: vlib.brute_density with the path-tree cycle basis, mesh-analysis currents, own Pfaffian;
g_brute(r) = D(t(r)) |det[t; dt/dr]|."""
import random
import sys
import time
from fractions import Fraction as Fr

from vlib import (fan_cycles, t_of_r, rand_r, jac_det, brute_density, claim_density, claim_Pi, claim_LM,
                  claim_gh, claim_R, half_power, same_halfint, laplacian, finv, period, fdet, ia, ib)

rng = random.Random(7)
NS = [int(x) for x in sys.argv[1:]] or [2, 4, 6]


def det2(x, y):
    return x[0] * y[1] - x[1] * y[0]


def s1_eq_1_40_integrand(n, u, v, lam2):
    """S1 eq-1-40 at z = 0 with r_{2k-1} = u_k, r_{2k} = v_k; lam2[i] (i = 1..n+1) 2-vectors;
    lambda_{+a} mu_+^a = det(lambda, mu) (htk CONVENTIONS Section 3); currents from eq-1-38."""
    r = {}
    for k in range(1, n + 1):
        r[2 * k - 1], r[2 * k] = u[k], v[k]

    def lred(i):
        out = [Fr(0), Fr(0)]
        for j in range(1, i + 1):
            p = Fr(1)
            for k in range(j, i):
                p *= r[2 * k - 1]
            out = [out[0] + p * lam2[j][0], out[1] + p * lam2[j][1]]
        return out
    y = {}
    ln = lred(n)
    y[n] = [r[2 * n - 1] * ln[c] - r[2 * n] * lam2[n + 1][c] for c in range(2)]
    for i in range(n - 1, 0, -1):
        li = lred(i)
        y[i] = [r[2 * i - 1] * r[2 * i] * li[c] - (1 - r[2 * i]) * lam2[i + 1][c] + (1 - r[2 * i]) * y[i + 1][c]
                for c in range(2)]
    val = Fr(1)
    for i in range(1, n):
        li = lred(i)
        val *= r[2 * i] * det2(li, [lam2[i + 1][c] - y[i + 1][c] for c in range(2)])
    val *= det2(lred(n), lam2[n + 1])
    return val, y


def potential_form(n, u, v, t, Gm, lam2):
    """J(r) prod_i det(Phi_i, Phi_{i+1}),  J = prod_{i<n} v_i R_i/(a_i B_i B_{i+1}) / (a_n B_n b_{n+1})."""
    R = claim_R(n, u, v)
    B = {i: u[i] * R[i] for i in range(1, n + 1)}
    a = {i: t[ia(i)] for i in range(1, n + 1)}
    J = Fr(1)
    for i in range(1, n):
        J *= v[i] * R[i] / (a[i] * B[i] * B[i + 1])
    J /= a[n] * B[n] * v[n]
    Phi = {w: [sum(Gm[w - 1][x - 1] * lam2[x][c] for x in range(1, n + 2)) for c in range(2)] for w in range(1, n + 2)}
    prod = Fr(1)
    for i in range(1, n + 1):
        prod *= det2(Phi[i], Phi[i + 1])
    return J * prod, J


if __name__ == "__main__":
    for n in NS:
        t0 = time.time()
        G, C = fan_cycles(n, "pathtree")
        sgn = (-1) ** (n * (n + 1) // 2)
        npts = 3 if n <= 4 else (2 if n == 6 else 1)
        for p in range(npts):
            u, v = rand_r(n, rng)
            t = t_of_r(n, u, v)
            Jr = jac_det(n, u, v)
            aJ = abs(Jr)
            # ---------- (2,0): general 2-vector sources
            lam2 = {i: [Fr(rng.randint(-5, 5)), Fr(rng.randint(-5, 5))] for i in range(1, n + 2)}
            D20, _ = brute_density(G, C, t, [[lam2[i][0] for i in range(1, n + 2)], [lam2[i][1] for i in range(1, n + 2)]], 0)
            g20 = D20 * aJ
            s1val, _ = s1_eq_1_40_integrand(n, u, v, lam2)
            Gm = finv(laplacian(G, t))
            pf, Jpot = potential_form(n, u, v, t, Gm, lam2)
            assert g20 == sgn * s1val == sgn * pf, (g20, sgn * s1val, sgn * pf)
            # ---------- (1,1): general scalar sources
            lam = {i: Fr(rng.randint(-6, 6)) for i in range(1, n + 2)}
            lam[1] = lam[1] or Fr(1)
            lam[n + 1] = lam[n + 1] or Fr(2)
            D11, U = brute_density(G, C, t, [[lam[i] for i in range(1, n + 2)]], 1)
            gb = (D11 * aJ, 1 / U)          # g = D11 |J| U^{-1/2}
            gc = claim_density(n, u, v, lam)
            assert same_halfint(gb[0], gb[1], gc[0], gc[1]), (gb, gc)
            # Pi = U^{1/2} J prod a_i
            pr, prad = claim_Pi(n, u, v)
            A = Fr(1)
            for i in range(1, n + 1):
                A *= t[ia(i)]
            assert same_halfint(pr, prad, Jpot * A, U)
            # M_j = sum_{w>j} h_w lambda_w
            g_, h_, P_, R_, c_ = claim_gh(n, u, v)
            L, M = claim_LM(n, u, v, lam)
            for j in range(1, n + 1):
                assert M[j] == sum(h_[w] * lam[w] for w in range(j + 1, n + 2))
                lred = Fr(0)
                for w in range(1, j + 1):
                    pw = Fr(1)
                    for k in range(w, j):
                        pw *= u[k]
                    lred += pw * lam[w]
                assert L[j] * P_[j] == lred          # L_j = lambda^red_j / P_j (S1 eq-B-fan-first-cycle)
            # ---------- (1,1) endpoint sources
            l1, ln1 = Fr(rng.randint(1, 5)), Fr(-rng.randint(1, 5))
            lamE = [l1] + [Fr(0)] * (n - 1) + [ln1]
            DE, U = brute_density(G, C, t, [lamE], 1)
            pairs = [(v[n], n - 2)]
            for k in range(1, n):
                pairs += [(v[k], 1), (1 - v[k], k - 2)]
            for m in range(1, n + 1):
                pairs.append((u[m], n - 1 - m))
            er, erad = half_power(pairs)
            dfact = 1
            for k in range(n - 1, 0, -2):
                dfact *= k
            er = er * (l1 * ln1) ** (n // 2) * dfact / 2 ** (n // 2)
            assert same_halfint(DE * aJ, 1 / U, er, erad)
            # ---------- (2,0) endpoint: (-1)^{n(n+1)/2} S1 eq-1-41
            x1 = [Fr(rng.randint(-4, 4)), Fr(rng.randint(1, 4))]
            xn = [Fr(rng.randint(1, 4)), Fr(rng.randint(-4, 4))]
            lamE2 = {i: [Fr(0), Fr(0)] for i in range(1, n + 2)}
            lamE2[1], lamE2[n + 1] = x1, xn
            D20E, _ = brute_density(G, C, t, [[lamE2[i][0] for i in range(1, n + 2)], [lamE2[i][1] for i in range(1, n + 2)]], 0)
            w41 = det2(x1, xn) ** n * v[n] ** (n - 1)
            for j in range(1, n):
                w41 *= u[j] ** (n - j) * v[j] * (1 - v[j]) ** (j - 1)
            assert D20E * aJ == sgn * w41
            print("n=%d point %d: sign det[t;dt/dr]=%+d; (2,0) brute = (-1)^{n(n+1)/2} eq-1-40 = (-1)^{n(n+1)/2} J prod det(Phi_i,Phi_i+1);"
                  " (1,1) brute = 2^{-n/2} Pi Haf[L,M] (sign and square); Pi = U^{1/2} J prod a_i; M_j = sum_{w>j} h_w lam_w;"
                  " endpoint (1,1) closed density; (2,0) endpoint = (-1)^{n(n+1)/2} eq-1-41;  g11 = %s * sqrt(%s)"
                  % (n, p, 1 if Jr > 0 else -1, gc[0], gc[1]))
        print("   n=%d: %.1f s" % (n, time.time() - t0))
