"""10d, one loop: [O_1] is proportional to the local class of Tr_ad(X^6) (hexagon, type-I README).  On a Cartan
subalgebra Tr_ad(ad_h^6) = sum_roots alpha(h)^6 =: P_6, kappa(h,h) ~ P_2 = sum alpha(h)^2.
(a) P_6 != 0 for every simple g (one-loop class nonzero: Sym^6(g*)^g -> H^7 injective [src GW 1810.06534]).
(b) GS-type factorization Tr_ad X^6 = kappa(X,X) Y_4(X), Y_4 in Sym^4(g*)^g  <=>  P_2 | P_6 in Q[h]
    (Chevalley restriction; the quotient is W-invariant).  Exact polynomial division (sympy)."""
import sys
import numpy as np
import sympy as sp
from common import Log
import lie

log = Log(__file__)
keys = sys.argv[1:] or ['A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7', 'B2', 'B3', 'B4', 'B5', 'B6', 'C3', 'C4', 'C5', 'C6',
                        'D4', 'D5', 'D6', 'D7', 'D8', 'G2', 'F4', 'E6', 'E7', 'E8']
for key in keys:
    L = lie.build(key)
    r = len(L.cartan)
    x = sp.symbols('x0:%d' % r)
    roots = [tuple(int(v) for v in L.wt[a]) for a in range(L.dim) if any(L.wt[a])]
    lin = [sum(c * xi for c, xi in zip(w, x)) for w in roots]
    P2 = sp.expand(sum(l ** 2 for l in lin))
    P6 = sp.expand(sum(l ** 6 for l in lin))
    q, rem = sp.div(sp.Poly(P6, *x), sp.Poly(P2, *x))
    fact = rem.is_zero
    log('%s %-16s rank %d roots %d: P_6 != 0: %s;  P_2 | P_6 (GS-type factorization): %s'
        % (key, L.name, r, len(roots), P6 != 0, fact))
# so(N) family from the vector representation: Tr_ad X^6 = (N-32) tr X^6 + 15 tr X^2 tr X^4 [src type-I README]
log('so(N), V = C^N: Tr_ad X^6 = (N-32) tr_V X^6 + 15 tr_V X^2 tr_V X^4  [src explore/type-i/README.md, closed exchange]')
log.close()
