import sys, time
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
import sympy as sp
import htk
from htk import catalog, spred
import enumerate_graphs as eg

gs = [G for G in eg.enumerate_tight(2, 4) if eg.is_sp(G)]
for k, Gx in enumerate(gs):
    G, labels = htk.from_networkx(Gx, root=None) if False else htk.from_networkx(Gx, None)
    lam = [[sp.Symbol("l%d" % v)] for v in range(G.nv - 1)]
    t0 = time.time()
    red = spred.sp_reduce(G, prefer="large")
    print("graph", k, G.edges, "root", G.root)
    print(red.describe())
    dens = red.density(lam, T=1, check=3)
    K = dens.integrate()
    print("  exps", dens.exponents())
    print("  K =", sp.factor(K), " (%.1fs)" % (time.time() - t0))
