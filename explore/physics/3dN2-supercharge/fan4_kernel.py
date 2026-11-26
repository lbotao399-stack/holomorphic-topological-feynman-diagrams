"""htk exact T=1 kernels at z = 0: fan(2), fan(4) (endpoint sources), tadpole(3,0); log fan4_kernel.log"""
import sys, time
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/circuit-general/ht-d2")
import sympy as sp
from htk import graph, sp_reduce
a, b = sp.symbols("a b")
for n in (2, 4):
    t = time.time()
    G = graph.fan(n)
    lam = [[a]] + [[0]] * (n - 1) + [[b]]
    K = sp.factor(sp_reduce(G, prefer="large").kernel(lam, T=1))
    Cn = (4 * sp.pi) ** sp.Rational(n, 2) / ((n + 1) * sp.factorial(n) ** 2)
    print("fan(%d) edges %s: K = %s ; C_n = %s ; K/(C_n (ab)^{n/2}) = %s  [%.1fs]" % (n, G.edges, K, Cn, sp.simplify(K / (Cn * (a * b) ** (n // 2))), time.time() - t))
from forests import tadpole_graph
for c in (3,):
    G, q = tadpole_graph(c, 0)
    lam = [[0]] * (G.nv - 1); lam[q] = [a]
    print("tadpole(%d,0) edges %s q=%d: K = %s" % (c, G.edges, q, sp.factor(sp_reduce(G, prefer="large").kernel(lam, T=1))))
G, q = tadpole_graph(5, 0)
print("tadpole(5,0) edges", G.edges, "q =", q, "(not series-parallel; exact value not available)")
