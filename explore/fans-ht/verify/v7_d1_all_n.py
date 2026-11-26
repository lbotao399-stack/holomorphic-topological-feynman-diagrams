"""V15 (claim d1, chain steps; claim a4 for odd n): (2,0) brute force at odd and even n, and the
endpoint chain  det(Phi_i,Phi_{i+1}) = det(lam_1,lam_{n+1}) P_n u_n v_n a_i,
J prod a_i (P_n u_n v_n)^n = v_n^{n-1} prod_{j<n} u_j^{n-j} v_j (1-v_j)^{j-1}."""
import random
from fractions import Fraction as Fr

from vlib import currents_mesh, fan_cycles, t_of_r, rand_r, jac_det, brute_density, claim_gh, laplacian, finv, ia
from v2_density import s1_eq_1_40_integrand, potential_form, det2

rng = random.Random(77)


def prod_u(u, j, i):
    p = Fr(1)
    for k in range(j, i):
        p *= u[k]
    return p


for n in range(1, 8):
    G, C = fan_cycles(n, "pathtree")
    sgn = (-1) ** (n * (n + 1) // 2)
    for p in range(2):
        u, v = rand_r(n, rng)
        t = t_of_r(n, u, v)
        aJ = abs(jac_det(n, u, v))
        lam2 = {i: [Fr(rng.randint(-5, 5)), Fr(rng.randint(-5, 5))] for i in range(1, n + 2)}
        D20, _ = brute_density(G, C, t, [[lam2[i][0] for i in range(1, n + 2)], [lam2[i][1] for i in range(1, n + 2)]], 0)
        s1val, _ = s1_eq_1_40_integrand(n, u, v, lam2)
        Gm = finv(laplacian(G, t))
        pf, J = potential_form(n, u, v, t, Gm, lam2)
        assert D20 * aJ == sgn * s1val == sgn * pf
        # potential-form intermediates of claim a4 (one component suffices; check both)
        _, ycur = s1_eq_1_40_integrand(n, u, v, lam2)
        Rr = claim_gh(n, u, v)[3]
        Bsp = {i: u[i] * Rr[i] for i in range(1, n + 1)}
        Phi2 = {w: [sum(Gm[w - 1][x - 1] * lam2[x][k] for x in range(1, n + 2)) for k in range(2)] for w in range(1, n + 2)}
        for k in range(2):
            lam_c = [lam2[i][k] for i in range(1, n + 2)]
            ymesh = currents_mesh(G, t, lam_c)
            for i in range(1, n + 1):
                lred = sum(lam2[j][k] * prod_u(u, j, i) for j in range(1, i + 1))
                assert lred == (Phi2[i][k] - Phi2[i + 1][k]) / t[ia(i)] + Phi2[i][k] / Bsp[i]
                assert ymesh[ia(i)] == ycur[i][k]          # eq-1-38 currents = Kirchhoff currents
            for i in range(1, n):
                lred1 = sum(lam2[j][k] * prod_u(u, j, i + 1) for j in range(1, i + 2))
                assert lred1 - ycur[i + 1][k] == Phi2[i + 1][k] / Bsp[i + 1]
        x1 = [Fr(rng.randint(-4, 4)), Fr(rng.randint(1, 4))]
        xn = [Fr(rng.randint(1, 4)), Fr(rng.randint(-4, 4))]
        lamE = {i: [Fr(0), Fr(0)] for i in range(1, n + 2)}
        lamE[1], lamE[n + 1] = x1, xn
        D20E, _ = brute_density(G, C, t, [[lamE[i][0] for i in range(1, n + 2)], [lamE[i][1] for i in range(1, n + 2)]], 0)
        w41 = det2(x1, xn) ** n * v[n] ** (n - 1)
        for j in range(1, n):
            w41 *= u[j] ** (n - j) * v[j] * (1 - v[j]) ** (j - 1)
        assert D20E * aJ == sgn * w41
        g, h, P, R, c = claim_gh(n, u, v)
        Phi = {w: [sum(Gm[w - 1][x - 1] * lamE[x][k] for x in range(1, n + 2)) for k in range(2)] for w in range(1, n + 2)}
        prodA = Fr(1)
        for i in range(1, n + 1):
            a_i = t[ia(i)]
            prodA *= a_i
            assert det2(Phi[i], Phi[i + 1]) == det2(x1, xn) * P[n] * u[n] * v[n] * a_i
        mono = v[n] ** (n - 1)
        for j in range(1, n):
            mono *= u[j] ** (n - j) * v[j] * (1 - v[j]) ** (j - 1)
        assert J * prodA * (P[n] * u[n] * v[n]) ** n == mono
    print("n=%d: (2,0) brute = (-1)^{n(n+1)/2} eq-1-40 = (-1)^{n(n+1)/2} J prod det(Phi_i,Phi_i+1) (general sources);"
          " eq-1-38 currents = Kirchhoff currents; lam^red_i = (Phi_i-Phi_i+1)/a_i + Phi_i/B_i; y^red_{i+1,O} = Phi_{i+1}/B_{i+1};"
          " endpoint brute = (-1)^{n(n+1)/2} eq-1-41; det(Phi_i,Phi_i+1) = det(l1,l_n+1) P_n u_n v_n a_i;"
          " J prod a_i (P_n u_n v_n)^n = eq-1-41 monomial: exact at 2 points" % n)
