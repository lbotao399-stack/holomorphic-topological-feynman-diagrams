"""hCS on C^3 (H = 3, T = 0): is the one-loop (box) obstruction exact?
omega_box = alt(K_box W_box), K_box = -Delta/6 (engine orientation of polygon(4)); weight (1,1,1).
Spanning set of invariant 3-cochains: alt(T (x) mu), T over trace monomials and f, mu over the
canonical monomials of degree (1,1,1) in (p_1, p_2).  Ranks exact (cochain.exact_rank)."""
import sys, time
import numpy as np
import sympy as sp
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/anomaly-structure/code')
import htk
from htk import polygon
import lie, cochain as cc
from c_obstruction import trace_tensors

H = 3
mom = cc.Momenta(6, H)
P = mom.p
box = polygon(4)
K = sp.expand(htk.kernel_exact(box, {0: list(P[0]), 1: list(P[1]), 2: list(P[2])}, T=0) * 6)
print("6 K_box =", sp.factor(K))
for g in [lie.sl(2), lie.gl(2), lie.sl(3), lie.gl(3)]:
    t0 = time.time()
    W, sc = cc.colour_weight(g, box.nv, box.edges, [0, 1, 2, 3])
    w = cc.alternate(mom, K, W, 4)
    monos4 = cc.all_monos(mom, 4, (1, 1, 1))
    monos3 = cc.all_monos(mom, 3, (1, 1, 1))
    cols = []
    for T in [T for s, T in trace_tensors(g, 3)] + [g.fl]:
        for mu in monos3:
            eta = cc.alternate(mom, mom.mono_expr(mu, 3, [P[0], P[1]]), T, 3)
            cols.append(cc.flatten(cc.d_CE(mom, g, eta, 3), monos4, g.dim, 4) if eta else
                        np.zeros(len(monos4) * g.dim ** 4, dtype=np.int64))
    A = np.stack(cols, axis=1)
    v = cc.flatten(w, monos4, g.dim, 4)
    r1 = cc.exact_rank(A)
    r2 = cc.exact_rank(np.concatenate([A, v[:, None]], axis=1))
    print("%s: omega_box zero %s, closed %s; rank d(C3) = %d, rank [d(C3) | omega_box] = %d -> exact: %s (%.1fs)" % (
        g.name, cc.is_zero(w), cc.is_zero(cc.d_CE(mom, g, w, 4)), r1, r2, r1 == r2, time.time() - t0), flush=True)
