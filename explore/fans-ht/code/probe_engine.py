import sys, time
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
import sympy as sp
from htk import fan, sp_reduce
for n in (2, 4):
    G = fan(n)
    lam = [[sp.Symbol("l%d" % (v + 1))] for v in range(G.nv - 1)]
    for pref in ("small", "large"):
        t0 = time.time()
        red = sp_reduce(G, D=2, prefer=pref)
        print(red.describe())
        dens = red.density(lam, T=1, check=2)
        print("exps", {k: str(v) for k, v in dens.exponents().items()})
        print("num", sp.factor(dens.num.as_expr()))
        K = dens.integrate()
        print("n=%d prefer=%s K =" % (n, pref), sp.factor(K), "time %.1f" % (time.time() - t0))
