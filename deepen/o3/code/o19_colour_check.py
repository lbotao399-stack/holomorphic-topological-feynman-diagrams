"""(d) Independent check of the colour factors W_G^(v) (d at vertex v): explicit structure tensors of lie.py
(fl = kappa([X_a,X_b],X_c), dl = tr(X_a{X_b,X_c}), kappa^{-1} exact), einsum over the graph with the local orders L'
and the sign sgn(L -> L'), compared entrywise with the double-line formal result evaluated at N
(formal_gl.tensor_from_formal), gl_3, sl_3, gl_2."""
import numpy as np
from common import *
import lie
import formal_gl as fg

LET = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
for g in (lie.gl(3), lie.sl(3), lie.gl(2)):
    d = g.dim
    M = np.array(g.mats, dtype=np.int64)
    # dl[a,b,c] = tr(X_a X_b X_c) + tr(X_a X_c X_b)
    t1 = np.einsum('aij,bjk,cki->abc', M, M, M)
    dl = t1 + np.transpose(t1, (0, 2, 1))
    allok = True
    for pqr in [(1, 3, 3), (2, 2, 3)]:
        G, ext, lam = graph_data(pqr)
        local = local_orders(G.nv, G.edges, ext)
        L = [('x', v) for v in ext] + [('e', e, s) for e in range(G.ne) for s in (0, 1)]
        let = {h: LET[i] for i, h in enumerate(L)}
        sgn = orientation_sign(G.nv, G.edges, ext)
        for v in range(G.nv):
            ops, terms = [], []
            for u in range(G.nv):
                terms.append(''.join(let[h] for h in local[u]))
                ops.append(dl if u == v else g.fl)
            for e in range(G.ne):
                terms.append(let[('e', e, 0)] + let[('e', e, 1)])
                ops.append(g.kinv)
            out = ''.join(let[('x', w)] for w in ext)
            T = np.einsum(','.join(terms) + '->' + out, *ops, optimize='greedy').astype(object) * sgn
            W = colour_weight(G.nv, G.edges, ext, {v: 'd'})
            Wf = fg.tensor_from_formal({t: {p: int(x) for p, x in dd.items()} for t, dd in W.items()}, g, 4).astype(object)
            ok = np.array_equal(T, Wf * g.kinv_scale ** G.ne)
            allok &= ok
            print("  %s theta%s v=%d: einsum tensor == double-line W^(v): %s ; nonzero %s" % (g.name, pqr, v, ok, bool(np.any(T != 0))))
    print("%s: all equal: %s" % (g.name, allok), flush=True)
