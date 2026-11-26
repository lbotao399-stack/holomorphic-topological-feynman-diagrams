"""Check of the sparse engine against simple-g/code/core.py (read-only import): for the same structure constants
(simple-g/code/algebras.build), s^7 W_133, s^7 W_223 (slot 0 = b_theta), s^7 F_133, s^7 F_223 and beta (slot 0 = b_theta,
all 6 momentum coefficients) are compared entry by entry; full tensors W are compared for the smallest algebras."""
import sys, time
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/simple-all/code')
import numpy as np
import lie, t4
sys.path.insert(1, lie.SG)
import algebras as al, core


def dense_slot0(T, th, d):
    m = T.idx[:, 0] == th
    out = np.zeros((d, d, d), dtype=object)
    for (a, b, c, e), v in zip(T.idx[m], T.val[m]):
        out[b, c, e] = int(v)
    return out


for key in sys.argv[1:]:
    t0 = time.time()
    Ld = al.build(key)
    L = lie.from_dense(Ld)
    E = t4.Engine(L)
    W133, W223 = E.weights()
    T = core.Tensors(Ld)
    th = Ld.theta
    d = L.dim
    ok = []
    for name, W in (('133', W133), ('223', W223)):
        ref = T.W(name, [th, None, None, None]).astype(object)
        ok.append(np.array_equal(dense_slot0(W, th, d), ref) and W.sp == 7)
    F = {g: t4.alt_bb(W) for g, W in (('133', W133), ('223', W223))}
    okF = []
    for g in ('133', '223'):
        ref = core.F_eval(T, g, th)
        okF.append(all(np.array_equal(dense_slot0(F[g][mu], th, d), ref[mu].astype(object)) for mu in range(6)))
    B = t4.beta_t(E.ff)
    refb = core.beta_eval(Ld, th)
    okb = all(np.array_equal(dense_slot0(B[mu], th, d), refb[mu].astype(object)) for mu in range(6))
    full = ''
    if d <= 15:
        for name, W in (('133', W133), ('223', W223)):
            ref = T.W(name, [None, None, None, None]).astype(object)
            D = np.zeros((d,) * 4, dtype=object)
            for (a, b, c, e), v in zip(W.idx, W.val):
                D[a, b, c, e] = int(v)
            full += ' full W%s equal: %s ;' % (name, np.array_equal(D, ref))
    print("%-14s dim %3d s %d : W133(b_theta,...) equal %s, W223 equal %s ; F133, F223 equal %s ; beta equal %s ;%s nnz W133 %d W223 %d (%.1fs)"
          % (L.name, d, L.s, ok[0], ok[1], okF, okb, full, W133.nnz(), W223.nnz(), time.time() - t0), flush=True)
