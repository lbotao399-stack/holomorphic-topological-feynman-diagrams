"""htk engine, n = 6, (H,T) = (1,1), integer sources (endpoint and two general vectors), compared with
general_kernels.txt (K_6 from (3.4)/(5.1))."""
import sys, time, pickle
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
import sympy as sp
from htk import fan, sp_reduce
n = 6
G = fan(n)
lam_s = sp.symbols("l1:%d" % (n + 2))
K = sp.sympify(pickle.load(open("../general_kernels.pkl", "rb"))[n], locals={str(x): x for x in lam_s})
red = sp_reduce(G, D=2, prefer="large")
for lam in ([1, 0, 0, 0, 0, 0, 1], [2, -1, 3, 1, -2, 1, 1], [1, 2, -1, 1, 3, -1, 2]):
    t0 = time.time()
    Ke = red.kernel([[x] for x in lam], T=1, check=1)
    Ke = sp.sympify(str(Ke)) if not isinstance(Ke, sp.Basic) else Ke
    mine = sp.expand(K.subs(dict(zip(lam_s, lam))))
    print("lam=%s engine %r, from (3.4): %s (%.0f s)" % (lam, Ke, mine, time.time() - t0), flush=True)
