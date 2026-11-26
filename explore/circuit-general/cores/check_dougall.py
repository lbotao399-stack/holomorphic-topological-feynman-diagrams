from fractions import Fraction as F
from math import comb, factorial
# Gegenbauer C_m^lam as exact polynomial coefficients (lam integer): C_m(x)=sum_k (-1)^k (lam)_{m-k} (2x)^{m-2k}/(k!(m-2k)!)
def rf(a, n):
    r = F(1)
    for i in range(n): r *= (a + i)
    return r
def geg(m, lam):
    c = [F(0)] * (m + 1)
    for k in range(m // 2 + 1):
        c[m - 2 * k] += F((-1) ** k) * rf(lam, m - k) * F(2) ** (m - 2 * k) / (factorial(k) * factorial(m - 2 * k))
    return c
def pmul(a, b):
    r = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            r[i + j] += x * y
    return r
def moment(k, lam):
    # int_{-1}^{1} x^k (1-x^2)^(lam-1/2) dx / same with k=0 : = 0 for odd k, else (1/2)_j/(lam+1)_j with j=k/2
    if k % 2: return F(0)
    j = k // 2
    return rf(F(1, 2), j) / rf(F(lam + 1), j)
def norm_int(m1, m2, m3, lam):
    p = pmul(pmul(geg(m1, lam), geg(m2, lam)), geg(m3, lam))
    return sum(c * moment(k, lam) for k, c in enumerate(p))
def dougall(m1, m2, m3, lam):
    if (m1 + m2 + m3) % 2: return F(0)
    s = (m1 + m2 + m3) // 2; a, b, c = s - m1, s - m2, s - m3
    if min(a, b, c) < 0: return F(0)
    r = F(lam, s + lam) * rf(2 * lam, s) / rf(lam, s)
    for x in (a, b, c): r *= rf(lam, x) / factorial(x)
    return r
ok = True
for lam in (1, 2, 3, 4):
    for m1 in range(0, 6):
        for m2 in range(0, 6):
            for m3 in range(0, 6):
                if norm_int(m1, m2, m3, lam) != dougall(m1, m2, m3, lam):
                    ok = False; print('fail', lam, m1, m2, m3, norm_int(m1, m2, m3, lam), dougall(m1, m2, m3, lam))
print('Dougall normalized triple-product formula exact for lam=1..4, m_i<=5:', ok)
