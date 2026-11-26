"""n4 (Killing), m4 - n4, dim (g*^(x)3)^g, g in S^2 g for classical families from root data (vbk.analyse),
compared with the explorer's polynomials (tr_V normalization, Killing = r tr_V):
  so_N: n4_tr = -3N(N-1)(N-2)(N^2-14N+32)/16, r = N-2;  sp_N: 3N(N+1)(N+2)(N^2+14N+32)/16, r = N+2;
  sl_N: 24N^2(N^2-1), r = 2N.   n4_K = n4_tr / r^4.
"""
import sys, time
from fractions import Fraction as Fr
import vbk

def pred(fam, N):
    if fam == 'so':
        return Fr(-3 * N * (N - 1) * (N - 2) * (N * N - 14 * N + 32), 16) / Fr(N - 2) ** 4
    if fam == 'sp':
        return Fr(3 * N * (N + 1) * (N + 2) * (N * N + 14 * N + 32), 16) / Fr(N + 2) ** 4
    return Fr(24 * N * N * (N * N - 1)) / Fr(2 * N) ** 4

def pred83(fam, N):
    if fam == 'so':
        return Fr(-3 * N * (N - 1) * (N - 3) * (N - 4) * (N + 2), 2 * (N - 2) ** 3)
    if fam == 'sp':
        return Fr(-3 * N * (N - 2) * (N + 1) * (N + 3) * (N + 4), 2 * (N + 2) ** 3)
    return Fr(-3 * (N * N - 1) * (N * N - 4), N * N)

cases = []
for N in range(5, int(sys.argv[1]) + 1):
    if N % 2:
        cases.append(('so', N, 'B', (N - 1) // 2))
    else:
        cases.append(('so', N, 'D', N // 2))
for N in range(4, int(sys.argv[2]) + 1, 2):
    cases.append(('sp', N, 'C', N // 2))
for N in range(2, int(sys.argv[3]) + 1):
    cases.append(('sl', N, 'A', N - 1))
allok = True
for fam, N, typ, n in cases:
    t0 = time.time()
    r = vbk.analyse(typ, n)
    d = r['d']
    ok1 = r['n'][4] == pred(fam, N)
    ok2 = 8 * r['n'][4] - 3 * d == pred83(fam, N)
    ok3 = r['m'][4] - r['n'][4] == Fr(d, 8)
    allok &= ok1 and ok2 and ok3
    print('%s%-3d (%s%d) d=%4d  n4 = %-14s poly %s  8n4-3d = %-16s formula %s  m4-n4 = d/8 %s  dim(g*^3)^g = %d  g in S^2 g %d  g in L^2 g %d  w(L^2 g - g) %s  (%.1fs)'
          % (fam, N, typ, n, d, r['n'][4], ok1, 8 * r['n'][4] - 3 * d, ok2, ok3, r['adj_gg'], r['adj_S2'], r['adj_L2'],
             sorted(set(str(x[1]) for x in r['X2_w'])), time.time() - t0), flush=True)
print('all family checks True:', allok)
