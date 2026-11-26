import sys, time
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
import sympy as sp
from htk import fan, sp_reduce
n = 6
G = fan(n)
lam = [[sp.Symbol("l%d" % (v + 1))] for v in range(G.nv - 1)]
t0 = time.time()
red = sp_reduce(G, D=2, prefer="large")
K = red.kernel(lam, T=1, check=1)
print("engine n=6 general sources, prefer=large, time %.1f s" % (time.time() - t0))
print("K =", sp.factor(K))
