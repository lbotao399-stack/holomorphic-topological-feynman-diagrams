"""E1: M_j / (P_n u_n v_n) = lam_{n+1} + v_n^{-1} sum_{k=j+1}^n (1-u_k) lam^{(j)}_k / prod_{m=k}^{n-1} (1-v_m),
       lam^{(j)}_k = sum_{w=j+1}^k (P_k/P_w) lam_w;   L_i = lam^red_i / P_i,  lam^red_i = sum_{w<=i} (P_i/P_w) lam_w
       (exact at random rational r, random integer sources, n = 2, 4, 6, 8)."""
import sys, random, math
from fractions import Fraction as Fr
sys.path.insert(0, ".")
import fanlib as F
from integrate_general import beta2
rng = random.Random(41)
for n in (2, 4, 6, 8):
    for _ in range(3):
        u, v = F.random_r(n, rng)
        lam = [Fr(rng.randint(-9, 9)) for _ in range(n + 1)]
        L, M = F.LM(n, u, v, lam)
        P = {1: Fr(1)}
        for k in range(2, n + 2):
            P[k] = P[k - 1] * u[k - 1]
        for i in range(1, n + 1):
            assert L[i] == sum(P[i] / P[w] * lam[w - 1] for w in range(1, i + 1)) / P[i]
        for j in range(1, n + 1):
            s = Fr(0)
            for k in range(j + 1, n + 1):
                lk = sum(P[k] / P[w] * lam[w - 1] for w in range(j + 1, k + 1))
                den = Fr(1)
                for m in range(k, n):
                    den *= 1 - v[m]
                s += (1 - u[k]) * lk / den
            assert M[j] / (P[n] * u[n] * v[n]) == lam[n] + s / v[n]
print("E1 holds, n = 2, 4, 6, 8")
