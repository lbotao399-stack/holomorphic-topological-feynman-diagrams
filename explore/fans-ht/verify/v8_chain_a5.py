"""a5 chain steps: with phi = Lambda^{-1} lam (Kirchhoff potentials, own Laplacian),
alphahat_i = g_i phi_{i+1} - g_{i+1} phi_i = -a_i L_i,  betahat_j = h_j phi_{j+1} - h_{j+1} phi_j = a_j M_j,
and <l_i l_j> = 2 alphahat_i betahat_j (i < j) for l_i = tau_i phi_{i+1} - phi_i tau_{i+1}, <tau tau^T> = 2 Lambda^{-1}."""
import random
from fractions import Fraction as Fr
from vlib import fan_cycles, t_of_r, rand_r, laplacian, finv, claim_gh, claim_LM, ia

rng = random.Random(99)
for n in (2, 4, 6, 8):
    G, C = fan_cycles(n, "pathtree")
    for p in range(3):
        u, v = rand_r(n, rng)
        t = t_of_r(n, u, v)
        Gm = finv(laplacian(G, t))
        lam = {i: Fr(rng.randint(-6, 6)) for i in range(1, n + 2)}
        phi = {w: sum(Gm[w - 1][x - 1] * lam[x] for x in range(1, n + 2)) for w in range(1, n + 2)}
        g, h, P, R, c = claim_gh(n, u, v)
        L, M = claim_LM(n, u, v, lam)
        for i in range(1, n + 1):
            assert g[i] * phi[i + 1] - g[i + 1] * phi[i] == -t[ia(i)] * L[i]
            assert h[i] * phi[i + 1] - h[i + 1] * phi[i] == t[ia(i)] * M[i]
        # <l_i l_j> from the covariance 2 G (exact bilinear form)
        for i in range(1, n + 1):
            for j in range(i + 1, n + 1):
                # l_i = tau_i phi_{i+1} - tau_{i+1} phi_i
                ci = {i: phi[i + 1], i + 1: -phi[i]}
                cj = {j: phi[j + 1], j + 1: -phi[j]}
                cov = sum(2 * Gm[a - 1][b - 1] * x * y for a, x in ci.items() for b, y in cj.items())
                assert cov == 2 * (g[i] * phi[i + 1] - g[i + 1] * phi[i]) * (h[j] * phi[j + 1] - h[j + 1] * phi[j])
    print("n=%d: alphahat_i = -a_i L_i, betahat_j = a_j M_j, <l_i l_j> = 2 alphahat_i betahat_j (i<j): exact at 3 points" % n)
