"""Taylor coefficients of the (1,1) fan kernel at endpoint sources (lam_2 = .. = lam_n = 0), in the
format of S1 eq-sum-N4-coeff / eq-fan-full-taylor.

  -sum_e y_e z_e = -lam_1 z_{1,O} - lam_{n+1} z_{n+1,O} - sum_i y_{i,i+1} Z_i ,  Z_i = z_{i,i+1} + z_{i+1,O} - z_{i,O},
  y_{i,i+1} = sum_{k=i}^n mu_{ik} lam_1 - nu_i lam_{n+1},
  mu_{ik} = prod_{m=i}^{k-1}(1-v_m) prod_{l<=k} u_l * (v_k if k < n else 1),  nu_i = v_n prod_{m=i}^{n-1} (1-v_m).

  K(lam; z) = (lam_1 lam_{n+1})^{n/2} sum  [(-lam_1 z_{1,O})^{r0}/r0!] [(-lam_{n+1} z_{n+1,O})^{s0}/s0!]
              prod_{i<=k} [(-lam_1 Z_i)^{r_ik}/r_ik!] prod_i [(lam_{n+1} Z_i)^{s_i}/s_i!]  c^{(n)}_{(1,1)}[r, s],

  c^{(n)}_{(1,1)} = (n-1)!!/2^{n/2} * Gamma(1/2 + sum_i r_in) Gamma(n/2 + sum_i s_i) / Gamma((n+3)/2 + sum_i r_in + sum_i s_i)
                  * prod_{j=1}^{n-1} 2/(n+1-j+2 X_j)
                  * prod_{j=1}^{n-1} B(3/2 + sum_{i<=j} r_ij,  j/2 + Y_j),
  X_j = sum_{k>=j} sum_{i<=k} r_ik,   Y_j = sum_{i<=j} [s_i + sum_{k>j} r_ik].
"""
import sys
import math
import itertools
from fractions import Fraction as Fr

sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/fans-ht/code")
from integrate_general import gamma_half, beta2


def dfact(k):
    out = 1
    while k > 1:
        out *= k
        k -= 2
    return out


def coeff_11(n, r, s):
    """c^{(n)}_{(1,1)}; r: dict {(i,k): r_ik} (1 <= i <= k <= n), s: dict {i: s_i} (1 <= i <= n).
    Returns (q, kpi): value = q * sqrt(pi)^kpi.  (r0, s0 do not enter.)"""
    R = lambda i, k: r.get((i, k), 0)
    S = lambda i: s.get(i, 0)
    X = {j: sum(R(i, k) for k in range(j, n + 1) for i in range(1, k + 1)) for j in range(1, n + 1)}
    Y = {j: sum(S(i) + sum(R(i, k) for k in range(j + 1, n + 1)) for i in range(1, j + 1)) for j in range(1, n)}
    q = Fr(dfact(n - 1), 2 ** (n // 2))
    kpi = 0
    rn = sum(R(i, n) for i in range(1, n + 1))
    sn = sum(S(i) for i in range(1, n + 1))
    qa, ka = gamma_half(1 + 2 * rn)
    qb, kb = gamma_half(n + 2 * sn)
    qc, kc = gamma_half(n + 3 + 2 * rn + 2 * sn)
    q *= qa * qb / qc
    kpi += ka + kb - kc
    for j in range(1, n):
        q *= Fr(2, n + 1 - j + 2 * X[j])
        dq, dk = beta2(3 + 2 * sum(R(i, j) for i in range(1, j + 1)), j + 2 * Y[j])
        q *= dq
        kpi += dk
    return q, kpi


def coeff_20(n, r, s):
    """S1 eq-sum-N4-coeff (the (2,0) coefficient), for comparison."""
    R = lambda i, k: r.get((i, k), 0)
    S = lambda i: s.get(i, 0)
    f = math.factorial
    rn = sum(R(i, n) for i in range(1, n + 1))
    sn = sum(S(i) for i in range(1, n + 1))
    c = Fr(f(rn) * f(n - 1 + sn), f(n + 1 + rn + sn))
    for j in range(1, n):
        c /= n - j + 1 + sum(R(i, k) for k in range(j, n + 1) for i in range(1, k + 1))
    for j in range(1, n):
        a = sum(R(i, j) for i in range(1, j + 1))
        b = sum(S(i) + sum(R(i, k) for k in range(j + 1, n + 1)) for i in range(1, j + 1))
        bb = sum(S(i) + sum(R(i, k) for k in range(j, n + 1)) for i in range(1, j + 1))
        c *= Fr(f(1 + a) * f(j - 1 + b), f(j + 1 + bb))
    return c


def multi_indices(n, m):
    """All (r0, s0, {r_ik}, {s_i}) with total order m."""
    slots = [("r0",), ("s0",)] + [("r", i, k) for k in range(1, n + 1) for i in range(1, k + 1)] + \
            [("s", i) for i in range(1, n + 1)]
    L = len(slots)
    for comp in itertools.combinations_with_replacement(range(L), m):
        cnt = [0] * L
        for c in comp:
            cnt[c] += 1
        r0 = s0 = 0
        r, s = {}, {}
        for sl, c in zip(slots, cnt):
            if not c:
                continue
            if sl[0] == "r0":
                r0 = c
            elif sl[0] == "s0":
                s0 = c
            elif sl[0] == "r":
                r[(sl[1], sl[2])] = c
            else:
                s[sl[1]] = c
        yield r0, s0, r, s
