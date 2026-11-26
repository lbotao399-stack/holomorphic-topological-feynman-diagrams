"""Per-variable minimal exponents of the (1,1) general-source density g = 2^{-n/2} Pi Haf[L M] over all
monomials of its expansion in (u_l, 1-u_l, v_k, 1-v_k, lam) (exact expansion of integrate_general.build)."""
import sys
sys.path.insert(0, ".")
from fractions import Fraction as Fr
from integrate_general import build, hafnian
for n in (2, 4, 6, 8):
    S, L, M = build(n)
    H = hafnian(n, L, M)
    off = [0] * S.nv
    off[S.iv(n)] -= 2
    for k in range(1, n):
        off[S.iv(k)] += 1
        off[S.iV(k)] += k - 2
    for k in range(1, n + 1):
        off[S.iu(k)] -= k + 1
    mins = None
    for key in H:
        e = [x + y for x, y in zip(key, off)]
        mins = e if mins is None else [min(a, b) for a, b in zip(mins, e)]
    names = ["u%d" % k for k in range(1, n + 1)] + ["(1-u%d)" % k for k in range(1, n + 1)] + \
            ["v%d" % k for k in range(1, n + 1)] + ["(1-v%d)" % k for k in range(1, n + 1)]
    mono = " ".join("%s^%s" % (nm, Fr(m, 2)) for nm, m in zip(names, mins[:4 * n]) if m != 0)
    print("n=%d: %d terms; minimal monomial m_0 = %s" % (n, len(H), mono))
