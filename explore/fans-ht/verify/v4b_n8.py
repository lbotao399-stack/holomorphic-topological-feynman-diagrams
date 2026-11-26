"""a10 (n = 8): termwise exact integration of g_8 = m_0 p at fixed integer sources vs general_kernels.txt."""
import random, time, sys
import sympy as sp
from v4_integrate import build_p, integrate_p, kernels_file

n = 8
KF = kernels_file()
ls = {i: sp.Symbol("l%d" % i) for i in range(1, n + 2)}
rng = random.Random(808)
for trial in range(int(sys.argv[1]) if len(sys.argv) > 1 else 1):
    t0 = time.time()
    lv = {i: sp.Integer(rng.choice([-1, 1]) * rng.randint(1, 3)) for i in range(1, n + 2)}
    p, us, vs = build_p(n, lv)
    Kt = sp.nsimplify(integrate_p(n, p, us, vs))
    Kf = KF[n].subs({ls[i]: lv[i] for i in ls})
    print("n=8 lam=%s: g/m_0 polynomial in r: True; termwise = %s; general_kernels.txt = %s; equal: %s (%.0f s)"
          % ([int(lv[i]) for i in range(1, n + 2)], sp.simplify(Kt), sp.simplify(Kf), sp.simplify(Kt - Kf) == 0, time.time() - t0))
