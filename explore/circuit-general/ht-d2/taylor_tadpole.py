"""Exact Taylor coefficients of the n = 2 single-source Q_2 graph (tadpole(3,0) = bitriangle with O
at a degree-two vertex) at (H,T) = (1,1), symbolic source and displacements on all edges.
Output: taylor_tadpole.log"""
import sympy as sp
from forests import tadpole_graph
from htk import sp_reduce
G, q = tadpole_graph(3, 0)
l = sp.Symbol("l")
lam = [[0] for _ in range(G.nv - 1)]
lam[q] = [l]
z = {e: [sp.Symbol("z_%d%d" % G.edges[e])] for e in range(G.ne)}
red = sp_reduce(G, prefer="large")
Ks = red.taylor(lam, z, 3, T=1, check=2)
print(G, "q =", q, "root O =", G.root)
for m, K in enumerate(Ks):
    print("K_%d =" % m, sp.factor(K))
# fan(2) with endpoint sources for comparison (same graph, root at a degree-3 vertex)
from htk import graph
F = graph.fan(2)
a, b = sp.symbols("a b")
print("fan(2), endpoint sources a, b:", sp_reduce(F, prefer="large").kernel([[a], [0], [b]], T=1))
