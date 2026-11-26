"""Exact T=1 (H=1) Taylor tables of the six labelled bitriangles on vertices {0,1,2,O=3}
(the unique 2-Laman graph with 4 vertices = all n=2 graphs).  Edges: the 5 pairs (i<j),
oriented i->j, in lexicographic order.  Sources l0,l1,l2 symbolic, displacements z_e on all edges.
K_m = int g (-sum_e y_e z_e)^m/m!  (htk CONVENTIONS sec 2, 5).  Output: kernels_M.pkl, kernels.log"""
import sys, time, pickle, itertools
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
import sympy as sp
from htk import graph, sp_reduce

M = int(sys.argv[1]) if len(sys.argv) > 1 else 2
l = sp.symbols("l0 l1 l2")
pairs = [(i, j) for i in range(4) for j in range(i + 1, 4)]
out = {}
for missing in pairs:
    edges = [p for p in pairs if p != missing]
    G = graph.Graph(4, edges, "bitri-missing%d%d" % missing)
    z = {e: [sp.Symbol("z%d" % e)] for e in range(5)}
    t0 = time.time()
    red = sp_reduce(G, prefer="large")
    Ks = red.taylor([[l[0]], [l[1]], [l[2]]], z, M, T=1, check=2)
    Ks = [sp.expand(K) for K in Ks]
    out[tuple(edges)] = Ks
    print("edges", edges, "missing", missing, "time %.1fs" % (time.time() - t0))
    print("  K_0 =", sp.factor(Ks[0]))
    if M >= 1:
        print("  K_1 =", sp.factor(Ks[1]))
    sys.stdout.flush()
pickle.dump(out, open("kernels_%d.pkl" % M, "wb"))
