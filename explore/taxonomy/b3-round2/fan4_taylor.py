"""htk exact Taylor moments of fan(4) (H=1, T=1) in the edge order (-v,w) of qn-ht/hteng.h1_q4_shapes,
sources l0 = a, l4 = b (l1 = l2 = l3 = 0: no survivors there in cubic theories), z only on the root edges
(v,5).  K_m = int g (-sum y.z)^m/m!.  Output fan4_taylor.pkl, fan4_taylor.log.  Usage: python3 fan4_taylor.py M"""
import sys, os, time, pickle
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "code"))
import sympy as sp
from htk import graph, sp_reduce
M = int(sys.argv[1]) if len(sys.argv) > 1 else 2
log = open(os.path.join(HERE, "fan4_taylor.log"), "w")
def out(*x):
    s = " ".join(str(y) for y in x); print(s); log.write(s + "\n"); log.flush()
R = 5
fe = []
for i in range(4):
    fe += [(i, i + 1), (i, R)]
fe.append((4, R))
mine = sorted(fe, key=lambda e: (-e[0], e[1]))
a, b = sp.symbols("a b")
G = graph.Graph(6, mine, "fan4-mine")
zs = {}
for e, (v, w) in enumerate(mine):
    zs[e] = [sp.Symbol("z%d" % v)] if w == R else [0]
t0 = time.time()
red = sp_reduce(G, prefer="large")
Ks = red.taylor([[a], [0], [0], [0], [b]], zs, M, T=1, check=2)
Ks = [sp.expand(K) for K in Ks]
out("fan(4) edges (order -v,w):", mine)
for m, K in enumerate(Ks):
    out("K_%d/pi^2 =" % m, sp.factor(sp.expand(K / sp.pi ** 2)))
out("time %.1fs" % (time.time() - t0))
pickle.dump({"edges": mine, "K": Ks}, open(os.path.join(HERE, "fan4_taylor.pkl"), "wb"))
