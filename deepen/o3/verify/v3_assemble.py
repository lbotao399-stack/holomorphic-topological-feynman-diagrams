"""Assemble o' (order hbar^2 graphs) and theta'(T) (order hbar^3, one J_1-type vertex T in {d, s1, s0}) as formal
trace-basis 4-cochains (own colour calculus + own Wick sign + own kernels), and test them:
 (I)  sl_N: o' - (N^2/8){D,D} = -((3N^2+10)/20) d(u2f)     [defines iota = -((3N^2+10)/20) u2f]
 (II) sl_N: theta'(d) + {D, iota} = d(xi), xi = -((2N^2+5)/30) D3
plus the pieces.  Colour factors cross-checked against explicit einsum contraction at N = 3 (gl_3, sl_3).
Covariance: theta', o' recomputed with local order choice C (rotation at every vertex) -- must be identical."""
import pickle, time, random, itertools
from fractions import Fraction as Fr
import numpy as np
from vlib import *
from vgraph import *
import vexp

t0 = time.time()
KD = pickle.load(open('data/kernels.pkl', 'rb'))
AUT = 4


def cochain_from(Kp, col, sl):
    poly = poly_to_full(Kp)
    terms = [(c, cycles_of(perm), poly) for perm, c in col.items()]
    return from_raw(4, terms, sl)


def graph_sums(sl, choice='A', kernels=None):
    """returns o', theta'(T) for T in d, s1, s0, and per-piece data."""
    op = {}
    th = {'d': {}, 's1': {}, 's0': {}}
    pieces = {}
    F = {}
    for pqr in [(1, 3, 3), (2, 2, 3)]:
        G, ext, lam = theta_data(pqr)
        loc = local_orders(G, ext, choice)
        sg = wick_sign(G, ext, loc)
        if choice == 'A':
            K0, Kv = KD[pqr]['K0'], KD[pqr]['Kv']
        else:
            K0 = kernel_over_pi(G, lam)
            Kv = {v: kernel_over_pi(G, lam, insertion(G, loc, v, (0, 1))) for v in range(G.nv)}
        col0 = {p: c * sg for p, c in colour(G, ext, loc, {}, sl).items()}
        w0 = alt(cochain_from(K0, col0, sl), 4, sl)
        op = lin((1, op), (Fr(1, AUT), w0))
        # F_G = alt([p1p2][p3p4] W_G)
        bb = fpmul(bracket_poly(0, 1, 4), bracket_poly(2, 3, 4))
        F[pqr] = alt(from_raw(4, [(c, cycles_of(p), bb) for p, c in col0.items()], sl), 4, sl)
        for T in ('d', 's1', 's0'):
            if sl and T != 'd':
                continue
            for v in range(G.nv):
                colv = {p: c * sg for p, c in colour(G, ext, loc, {v: T}, sl).items()}
                wv = alt(cochain_from(Kv[v], colv, sl), 4, sl)
                pieces[(pqr, v, T)] = (colv, wv)
                th[T] = lin((1, th[T]), (Fr(1, AUT), wv))
    return op, th, pieces, F


def show_col(col):
    return " + ".join("(%s)%s" % (c, perm_str(p)) for p, c in sorted(col.items()))


for sl in (True, False):
    tag = 'sl_N' if sl else 'gl_N'
    print("=" * 30, tag, "=" * 30, flush=True)
    op, th, pieces, F = graph_sums(sl)
    D = cubic('d', 1, sl)
    u2f = cubic('f', 2, sl)
    D3 = cubic('d', 3, sl)
    N2 = LP.N(2)
    print("o' closed: %s;  number of terms %d" % (is_zero(d(op, 4, sl)), len(op)))
    print("o' = (1/240)(-3 F_133 + 2 F_223): %s" % is_zero(lin((1, op), (Fr(3, 240), F[(1, 3, 3)]), (Fr(-2, 240), F[(2, 2, 3)]))))
    du2f = d(u2f, 3, sl)
    print("F_133 + F_223 = 2N^2 d(u2f): %s" % is_zero(lin((1, F[(1, 3, 3)]), (1, F[(2, 2, 3)]), (-2 * N2, du2f))))
    DD = bracket(D, 3, D, 3, sl)
    if sl:
        eta = scal((N2 + 4) * LP.N(-2) * Fr(1, 2), u2f)
        print("(1/2){D,D} = F_223/(12N^2) + d eta_N, eta_N = ((N^2+4)/(2N^2)) u2f: %s" %
              is_zero(lin((Fr(1, 2), DD), (-LP.N(-2) * Fr(1, 12), F[(2, 2, 3)]), (-1, d(eta, 3, sl)))))
        iota = scal(-(3 * N2 + 10) * Fr(1, 20), u2f)
        I_ = lin((1, op), (-N2 * Fr(1, 8), DD), (-1, d(iota, 3, sl)))
        print("(I)  o' - (N^2/8){D,D} - d(iota) = 0, iota = -((3N^2+10)/20) u2f: %s" % is_zero(I_))
        for v in range(6):
            for pqr in [(1, 3, 3), (2, 2, 3)]:
                pass
        for pqr in [(1, 3, 3), (2, 2, 3)]:
            for v in range(6):
                colv, wv = pieces[(pqr, v, 'd')]
                print("  W^(v,d) theta%s v=%d: %s" % (pqr, v, show_col(colv) if colv else "0"))
        thd = th['d']
        print("theta'(d) closed: %s; nonzero: %s" % (is_zero(d(thd, 4, sl)), not is_zero(thd)))
        print("theta'(d) = -(N^2/60) d(D3): %s" % is_zero(lin((1, thd), (N2 * Fr(1, 60), d(D3, 3, sl)))))
        DI = bracket(D, 3, iota, 3, sl)
        Du = bracket(D, 3, u2f, 3, sl)
        print("{D,u2f} = (1/3) d(D3): %s" % is_zero(lin((1, Du), (Fr(-1, 3), d(D3, 3, sl)))))
        xi = scal(-(2 * N2 + 5) * Fr(1, 30), D3)
        II = lin((1, thd), (1, DI), (-1, d(xi, 3, sl)))
        print("(II) theta'(d) + {D,iota} = d(xi), xi = -((2N^2+5)/30) D3: %s" % is_zero(II))
        # opposite relative sign
        II2 = lin((1, thd), (-1, DI))
        print("     opposite relative sign: theta'(d) - {D,iota} = d(((N^2+5)/30) D3):",
              is_zero(lin((1, II2), (-1, d(scal((N2 + 5) * Fr(1, 30), D3), 3, sl)))))
        print("C-parity: C theta' = -theta': %s;  C{D,iota} = -{D,iota}: %s;  C D = -D: %s;  C u2f = u2f: %s" %
              (is_zero(lin((1, Cop(thd, 4)), (1, thd))), is_zero(lin((1, Cop(DI, 4)), (1, DI))),
               is_zero(lin((1, Cop(D, 3)), (1, D))), is_zero(lin((1, Cop(u2f, 3)), (-1, u2f)))))
        print("SL2: E theta' = F theta' = 0: %s %s" % (is_zero(sl2_E(thd, 4)), is_zero(sl2_F(thd, 4))))
        pickle.dump(dict(op=op, th=th, F=F, DD=DD, DI=DI, iota=iota), open('data/assembled_sl.pkl', 'wb'))
    else:
        S1 = cubic('s1', 1, sl)
        S0 = cubic('s0', 1, sl)
        S1_3 = cubic('s1', 3, sl)
        S0_3 = cubic('s0', 3, sl)
        print("gl: (1/2){D,D} = (1/2) d(u2f): %s" % is_zero(lin((Fr(1, 2), DD), (Fr(-1, 2), du2f))))
        # F_223 = F_0 + N F_1
        F223 = F[(2, 2, 3)]
        F0 = {k: LP({0: v.c.get(0, Fr(0))}) for k, v in F223.items() if v.c.get(0)}
        F1 = {k: LP({0: v.c.get(1, Fr(0))}) for k, v in F223.items() if v.c.get(1)}
        print("F_223 = F_0 + N F_1 (no other powers): %s" % all(set(v.c) <= {0, 1} for v in F223.values()))
        beta = du2f
        checks = [("(1/2){D,D}", scal(Fr(1, 2), DD), (0, 0, Fr(1, 2))),
                  ("{D,S1}", bracket(D, 3, S1, 3, sl), (0, Fr(-1, 12), 0)),
                  ("{D,S0}", bracket(D, 3, S0, 3, sl), (Fr(1, 24), 0, 1)),
                  ("(1/2){S1,S1}", scal(Fr(1, 2), bracket(S1, 3, S1, 3, sl)), (Fr(1, 48), -LP.N(1) * Fr(1, 48), Fr(1, 2))),
                  ("{S1,S0}", bracket(S1, 3, S0, 3, sl), (LP.N(1) * Fr(1, 48), 0, LP.N(1) * Fr(1, 2))),
                  ("{S0,S0}", bracket(S0, 3, S0, 3, sl), (0, 0, 0))]
        for name, val, (x0, x1, xb) in checks:
            ok = is_zero(lin((1, val), (-1 * (x0 if isinstance(x0, LP) else LP.const(x0)), F0),
                             (-1 * (x1 if isinstance(x1, LP) else LP.const(x1)), F1),
                             (-1 * (xb if isinstance(xb, LP) else LP.const(xb)), beta)))
            print("  %s = %s F_0 + %s F_1 + %s beta: %s" % (name, x0, x1, xb, ok))
        print("o' = (1/48) F_223 - (N^2/40) beta: %s" % is_zero(lin((1, op), (Fr(-1, 48), F223), (N2 * Fr(1, 40), beta))))
        for T in ('d', 's1', 's0'):
            print("theta'(%s): zero: %s; closed: %s" % (T, is_zero(th[T]), is_zero(d(th[T], 4, sl))))
        print("theta'(d) = d(-(N^2/60) D3 + (N/30) S1_3): %s" %
              is_zero(lin((1, th['d']), (-1, d(lin((-N2 * Fr(1, 60), D3), (LP.N(1) * Fr(1, 30), S1_3)), 3, sl)))))
        print("{D,u2f} = (1/3) d D3: %s;  {S1,u2f} = (1/3) d S1_3: %s;  {S0,u2f} = 0: %s;  d S0_3 = 0: %s" %
              (is_zero(lin((1, bracket(D, 3, u2f, 3, sl)), (Fr(-1, 3), d(D3, 3, sl)))),
               is_zero(lin((1, bracket(S1, 3, u2f, 3, sl)), (Fr(-1, 3), d(S1_3, 3, sl)))),
               is_zero(bracket(S0, 3, u2f, 3, sl)), is_zero(d(S0_3, 3, sl))))
        pickle.dump(dict(op=op, th=th, F=F, F0=F0, F1=F1, beta=beta), open('data/assembled_gl.pkl', 'wb'))
    print("time %.1fs" % (time.time() - t0), flush=True)

# covariance: local order choice C
for sl in (True,):
    op, th, _, _ = graph_sums(sl)
    opC, thC, _, _ = graph_sums(sl, 'C')
    print("local order C: o' identical: %s;  theta'(d) identical: %s" %
          (is_zero(lin((1, op), (-1, opC))), is_zero(lin((1, th['d']), (-1, thC['d'])))))

# explicit colour check at N = 3
random.seed(7)
for sl in (False, True):
    alg = vexp.Algebra(3, sl)
    ok, n = True, 0
    for pqr in [(1, 3, 3), (2, 2, 3)]:
        G, ext, lam = theta_data(pqr)
        loc = local_orders(G, ext, 'A')
        for vt in [{}, {0: 'd'}, {1: 'd'}, {3: 'd'}, {5: 'd'}] + ([] if sl else [{2: 's1'}, {0: 's0'}, {4: 's1'}]):
            col = colour(G, ext, loc, vt, sl)
            Xs = [alg.random_element(2) for _ in range(4)]
            formal = sum(c.at(3) * vexp.tau_eval(p, Xs) for p, c in col.items())
            expl = explicit_colour = None
            # explicit einsum over basis indices
            B = alg.basis
            dim = alg.dim
            def T_eval(T, mats):
                if T == 'f':
                    return np.trace(mats[0].dot(mats[1].dot(mats[2]) - mats[2].dot(mats[1])))
                if T == 'd':
                    return np.trace(mats[0].dot(mats[1].dot(mats[2]) + mats[2].dot(mats[1])))
                if T == 's1':
                    tr = np.trace
                    return tr(mats[0]) * tr(mats[1].dot(mats[2])) + tr(mats[1]) * tr(mats[0].dot(mats[2])) + tr(mats[2]) * tr(mats[0].dot(mats[1]))
                if T == 's0':
                    return np.trace(mats[0]) * np.trace(mats[1]) * np.trace(mats[2])
            # vertex tensors as arrays over internal half-edge indices
            nE = len(G.edges)
            tot = Fr(0)
            # brute force contraction: loop over assignments of basis indices to edges' tail ends and use kappa^{-1}
            # sum_{a_e, b_e} prod_e kinv[a_e][b_e] prod_v T_v(...) ; dim^14 too large -> contract vertex by vertex
            vt_arrays = []
            for v in range(G.nv):
                hs = loc[v]
                internal = [h for h in hs if h[0] == 'e']
                shape = [dim] * len(internal)
                arr = np.zeros(shape, dtype=object)
                for idx in itertools.product(range(dim), repeat=len(internal)):
                    mats, it = [], iter(idx)
                    for h in hs:
                        mats.append(Xs[ext.index(h[1])] if h[0] == 'x' else B[next(it)])
                    arr[idx] = T_eval(vt.get(v, 'f'), mats)
                vt_arrays.append((arr, internal))
            Ki = np.array([[alg.kinv[a][b] for b in range(dim)] for a in range(dim)], dtype=object)
            letters = {}
            for e in range(nE):
                letters[('e', e, 0)] = chr(97 + 2 * e)
                letters[('e', e, 1)] = chr(97 + 2 * e + 1)
            ops, subs = [], []
            for arr, internal in vt_arrays:
                ops.append(arr)
                subs.append("".join(letters[h] for h in internal))
            for e in range(nE):
                ops.append(Ki)
                subs.append(letters[('e', e, 0)] + letters[('e', e, 1)])
            expl = np.asarray(np.einsum(",".join(subs) + "->", *ops, optimize='greedy')).item()
            n += 1
            ok = ok and (Fr(expl) == Fr(formal)) and (n == 0 or True)
            if vt == {1: 'd'} or vt == {}:
                print('   example %s %s %s: explicit %s formal %s' % ('sl' if sl else 'gl', pqr, vt, expl, formal))
    print("%s_3: formal colour == explicit einsum (%d cases): %s" % ('sl' if sl else 'gl', n, ok), flush=True)
print("time %.1fs" % (time.time() - t0))
