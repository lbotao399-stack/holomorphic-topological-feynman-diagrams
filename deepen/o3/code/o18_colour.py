"""(d) Colour factors W_G^(v) = sgn(L -> L') prod_{u != v} f(L'_u) d(L'_v) prod_e kappa^{-1} (kappa = tr, double line),
formal gl_N, trace basis over the external labels 1..4 (= ext vertices 1,2,3,4); sl_N: restriction (single-trace
cycles of length 1 dropped).  Coefficients are polynomials in N."""
from common import *

for pqr in [(1, 3, 3), (2, 2, 3)]:
    G, ext, lam = graph_data(pqr)
    for v in range(G.nv):
        W = colour_weight(G.nv, G.edges, ext, {v: 'd'})
        Wsl = {t: d for t, d in W.items() if not ce.has_singleton(t)}
        terms = []
        for t in sorted(Wsl, key=lambda t: (len(t), t)):
            d = Wsl[t]
            poly = " + ".join("%s N^%d" % (x, p) if p else "%s" % x for p, x in sorted(d.items()))
            terms.append("(%s) %s" % (poly, ce.tau_str(tuple(tuple(i for i in c) for c in t))))
        print("theta%s v=%d  W^(v)|_sl = %s" % (pqr, v, " + ".join(terms)))
        print("      gl_N extra (terms with tr of a single X): %d" % (len(W) - len(Wsl)))
