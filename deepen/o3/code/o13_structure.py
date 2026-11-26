"""(e), (f) Structure of the SL2-invariant subcomplex at weight (3,3) under charge conjugation C (formal N):
dimensions of C^3_inv,+-, C^4_inv,+-, Z^4_inv,+-, B^4_inv,+-; every single piece alt(K_G^(v) W_G^(v)) (one graph,
one vertex) closed / exact; Jacobi {D,{D,D}} = 0; d{D,iota} = 0;  primitives:
    sl_N: theta' = d(-(N^2/60)[p1p2]^3 d),  {D,[p1p2]^2 f} = d((1/3)[p1p2]^3 d);
    gl_N: theta'(d) = d(x1 [p1p2]^3 d + x2 [p1p2]^3 s1 + x3 [p1p2]^3 s0), {J_T, [p1p2]^2 f} = d(...)."""
import time
from common import *

t0 = time.time()
KD = load('kernels.pkl')
th = load('theta.pkl')
AUT = {(1, 3, 3): 4, (2, 2, 3): 4}
deg = (3, 3)


def npow_parts(c):
    pw = sorted({p for d in c.values() for p in d})
    return {p: {kk: Fr(d[p]) for kk, d in c.items() if d.get(p)} for p in pw}


for sl in (True, False):
    nm = 'sl_N' if sl else 'gl_N'
    k3, k4, k5 = spaces.keys_of(3, deg, sl), spaces.keys_of(4, deg, sl), spaces.keys_of(5, deg, sl)
    out = {}
    for k, keys in ((3, k3), (4, k4)):
        S = inv_span(k, sl, deg)
        for sgn in (+1, -1):
            P = [spaces.const_part(ce.clean(ce.scale(ce.lin((1, spaces.to_np(c)), (sgn, Cop(spaces.to_np(c), k))), Fr(1, 2))))
                 for c in S]
            P = [x for x in P if x]
            if P:
                _, M = la.matrix_from(P, keys)
                cc = la.certified_rank(M)
                B = [P[j] for j in cc['pivcols']]
            else:
                B = []
            out[(k, sgn)] = B
    for sgn in (+1, -1):
        C3s, C4s = out[(3, sgn)], out[(4, sgn)]
        img = [spaces.const_part(ce.d_ce(spaces.to_np(c), 3, sl)) for c in C3s]
        rB = la.certified_rank(la.matrix_from(img, k4)[1])['rank'] if img else 0
        d4 = [spaces.const_part(ce.d_ce(spaces.to_np(c), 4, sl)) for c in C4s]
        d4 = [x for x in d4]
        rZ = la.certified_rank(la.matrix_from(d4, k5)[1])['rank'] if any(d4) else 0
        print("%s C-parity %+d: dim C^3_inv = %d, dim C^4_inv = %d, rank d3 = %d, rank d4 = %d, dim Z^4_inv = %d, dim H^4_inv = %d"
              % (nm, sgn, len(C3s), len(C4s), rB, rZ, len(C4s) - rZ, len(C4s) - rZ - rB), flush=True)
    # single pieces
    cnt = [0, 0]
    for pqr in AUT:
        G, ext, lam = graph_data(pqr)
        for v in range(G.nv):
            Wv = colour_weight(G.nv, G.edges, ext, {v: 'd'})
            piece = alt_cochain(KD[pqr]['Kv'][v], Wv, sl=sl)
            if is_zero(piece):
                continue
            cnt[0] += 1
            cnt[1] += is_zero(ce.d_ce(piece, 4, sl))
    print("%s: single pieces alt(K^(v) W^(v)) nonzero: %d, closed: %d" % (nm, cnt[0], cnt[1]), flush=True)
    # primitives in span{[p1p2]^3 T}
    u3 = ce.poly_mul(ce.poly_mul(ce.bracket_poly(0, 1, 3), ce.bracket_poly(0, 1, 3)), ce.bracket_poly(0, 1, 3))
    Ts = {'d': {((0, 1, 2),): {0: Fr(1)}, ((0, 2, 1),): {0: Fr(1)}},
          's1': {((0,), (1, 2)): {0: Fr(1)}, ((0, 2), (1,)): {0: Fr(1)}, ((0, 1), (2,)): {0: Fr(1)}},
          's0': {((0,), (1,), (2,)): {0: Fr(1)}}}
    names = ['d'] if sl else ['d', 's1', 's0']
    Ps = [ce.canon_raw(ce.raw_from(u3, Ts[T]), 3, sl=sl) for T in names]
    dPs = [spaces.const_part(ce.d_ce(P_, 3, sl)) for P_ in Ps]
    _, M = la.matrix_from(dPs, k4)
    cc = la.certified_rank(M)
    print("%s: rank d{[p1p2]^3 T, T in %s} = %d" % (nm, names, cc['rank']))
    targets = [("theta'(d)", ce.restrict_sl(th['th']['d']) if sl else th['th']['d'])]
    u2f = f_bracket_sq(sl)
    for T in names:
        targets.append(("{J_%s, [p1p2]^2 f}" % T, ce.bracket(J_cochain(T, sl), 3, u2f, 3, sl=sl)))
    for tn, tc in targets:
        res = {}
        okall = True
        for p, x in npow_parts(tc).items():
            ok, sol = la.in_span(M, la.column(x, k4), cc)
            okall &= ok
            res[p] = [str(s) for s in sol] if ok else None
        print("  %s = d( sum_T x_T [p1p2]^3 T ), x over %s by N-power: %s (solvable: %s)" % (tn, names, res, okall), flush=True)
    if sl:
        D = J_cochain('d', sl)
        DD = ce.bracket(D, 3, D, 3, sl=sl)
        print("sl_N: {D,{D,D}} == 0: %s" % is_zero(ce.bracket(D, 3, DD, 4, sl=sl)))
        i2 = load('i2_sl.pkl')
        print("sl_N: d{D, iota} == 0: %s" % is_zero(ce.d_ce(ce.bracket(D, 3, i2['iota'], 3, sl=sl), 4, sl)))
    else:
        jr = ce.lin((2, J_cochain('d', sl)), (3, J_cochain('s1', sl)), (5, J_cochain('s0', sl)))
        X = ce.bracket(jr, 3, jr, 3, sl=sl)
        print("gl_N: j = 2D + 3S1 + 5S0: {j,{j,j}} == 0: %s ; {j,j} != 0: %s" % (is_zero(ce.bracket(jr, 3, X, 4, sl=sl)), not is_zero(X)))
print("time %.0fs" % (time.time() - t0))
