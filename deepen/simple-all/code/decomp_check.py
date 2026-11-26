"""Cochain-level check of F_223 = x beta + a G_kk (Killing normalization), G_kk = alt([p0p1][p2p3] kappa_01 kappa_23),
with (x, a) read off from the probe solution (logs/b_analysis.log), on all entries and all 6 momentum coefficients.
Form used kappa, Killing = r kappa:  F_kappa = x r^2 beta_kappa + a r^3 Gkk_kappa  (F_kappa = s^-7 F_int).
Also: the same with F_133 (x', a') from F_133 = beta/2 - F_223 (Killing)."""
import sys, json, re, time
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/simple-all/code')
import numpy as np
from fractions import Fraction as Fr
import lie, t4

LOG = (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/simple-all/logs/'
sol = {}
for line in open(LOG + 'b_analysis.log'):
    m = re.match(r'(.{16}) d\s+\d+.*x = (\S+), a = (\S+),', line)
    if m:
        sol[m.group(1).strip()] = (Fr(m.group(2)), Fr(m.group(3)))
for key in sys.argv[1:]:
    t0 = time.time()
    L = lie.build(key)
    r, _ = L.killing_ratio()
    E = t4.Engine(L)
    W133, W223 = E.weights()
    F223 = t4.alt_bb(W223)
    B = t4.beta_t(E.ff)
    G = t4.alt_bb(E.kk)
    x, a = sol[L.name]
    keys, vals = t4.merge(F223 + B + G)
    Fv, Bv, Gv = vals[0:6], vals[6:12], vals[12:18]
    s7 = L.s ** 7
    cx, ca = x * r ** 2 * s7, a * r ** 3 * s7
    ok = all(all(Fr(int(f)) == cx * int(b) + ca * int(g) for f, b, g in zip(Fv[mu], Bv[mu], Gv[mu])) for mu in range(6))
    G_in_Cbeta = all(np.count_nonzero(Gv[mu] * int(Bv[3][np.nonzero(Bv[3])[0][0]]) - Bv[mu] * int(Gv[3][np.nonzero(Bv[3])[0][0]])) == 0 for mu in range(6))
    print("%-16s x = %s, a = %s : F_223 = x beta + a G_kk on all %d entries x 6 coefficients: %s ; G_kk in C beta: %s (%.1fs)"
          % (L.name, x, a, len(keys), ok, G_in_Cbeta, time.time() - t0), flush=True)
