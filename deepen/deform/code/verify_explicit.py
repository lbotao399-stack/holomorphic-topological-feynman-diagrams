"""Independent check of the decisive classes with explicit structure constants (no ce.py):
lie.py bases, cochain.alternate / cochain.d_CE / cochain.colour_weight (anomaly-structure, explicit
tensors), bracket by explicit kappa^{-1} contraction (numpy einsum, written here).
sl3, gl3, sl2, gl2: for J = a D + b S1 + c S0 (explicit tensors) solve
    (1/2){J,J} = lam F_223 + d(eta)   exactly (eta in the span of alt(mu T), T invariant 3-tensors),
or certify that no lam works."""
import sys, time, itertools
import numpy as np
import sympy as sp
from fractions import Fraction as Fr
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/anomaly-structure/code')
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/deform/code')
import lie, cochain as cc
from c_obstruction import kernel_poly, trace_tensors
import la

mom = cc.Momenta(6, 2)
P = mom.p
monos3 = cc.all_monos(mom, 3, (2, 2))
monos4 = cc.all_monos(mom, 4, (2, 2))


def tr_tensor(g, word_sets):
    """sum over listed cycle structures of prod tr; word_sets: list of lists of cycles (0-based slots)."""
    d = g.dim
    mats = g.mats
    T = np.zeros((d, d, d), dtype=np.int64)
    for idx in itertools.product(range(d), repeat=3):
        val = 0
        for cycles in word_sets:
            v = 1
            for c in cycles:
                m = np.eye(mats[0].shape[0], dtype=np.int64)
                for j in c:
                    m = m @ mats[idx[j]]
                v *= int(np.trace(m))
            val += v
        T[idx] = val
    return T


def J_cochain(T):
    """canonical [p1p2] T: monomials over gens(3) = (p0_0, p0_1, p1_0, p1_1)."""
    return {(1, 0, 0, 1): T.copy(), (0, 1, 1, 0): -T.copy()}


def bracket33(g, A, B):
    """{A,B} for canonical 3-cochains (dict mono -> tensor): alt(gamma_0)/4, gamma_0 = A(x1,x2,e) kinv B(y1,y2,e').
    returns (canonical 4-cochain with integer tensors, denominator)."""
    out = {}
    for ma, TA in A.items():
        for mb, TB in B.items():
            G = np.einsum('abe,ef,cdf->abcd', TA, g.kinv, TB)
            expr = mom.mono_expr(ma, 3, [P[0], P[1]]) * mom.mono_expr(mb, 3, [P[2], P[3]])
            alt = cc.alternate(mom, expr, G, 4)
            for m, T in alt.items():
                out[m] = out.get(m, 0) + T
    return {m: T for m, T in out.items() if np.any(T)}, 4 * g.kinv_scale


def run(g):
    t0 = time.time()
    print("==== %s (dim %d) ====" % (g.name, g.dim), flush=True)
    # F_223 explicit
    G223, ext, K = kernel_poly((2, 2, 3))
    W, sc = cc.colour_weight(g, G223.nv, G223.edges, ext)
    br4 = (P[0][0] * P[1][1] - P[0][1] * P[1][0]) * (P[2][0] * P[3][1] - P[2][1] * P[3][0])
    F = cc.alternate(mom, br4, W, 4)       # = sc * F_223
    # B^4 spanning set
    cols = []
    for T in [T for s, T in trace_tensors(g, 3)] + [g.fl]:
        for mu in monos3:
            eta = cc.alternate(mom, mom.mono_expr(mu, 3, [P[0], P[1]]), T, 3)
            if eta:
                cols.append(cc.flatten(cc.d_CE(mom, g, eta, 3), monos4, g.dim, 4))
    Bm = np.stack(cols, axis=1)
    Fv = cc.flatten(F, monos4, g.dim, 4)
    keep = np.any(Bm != 0, axis=1) | (Fv != 0)
    # J tensors
    Ts = {"D": tr_tensor(g, [[(0, 1, 2)], [(0, 2, 1)]])}
    if g.name.startswith("gl"):
        Ts["S1"] = tr_tensor(g, [[(0,), (1, 2)], [(1,), (0, 2)], [(2,), (0, 1)]])
        Ts["S0"] = tr_tensor(g, [[(0,), (1,), (2,)]])
    nonzeroT = {k: bool(np.any(T)) for k, T in Ts.items()}
    print("  invariant tensors nonzero:", nonzeroT)
    # closedness of J
    for k_, T in Ts.items():
        J = J_cochain(T)
        alt6 = cc.alternate(mom, P[0][0] * P[1][1] - P[0][1] * P[1][0], T, 3)
        same = all(np.array_equal(alt6.get(m, 0) * 1, 6 * J.get(m, 0)) for m in set(alt6) | set(J))
        dJ = cc.d_CE(mom, g, J, 3)
        print("  J_%s = [p1p2] T: alt([p1p2]T) = 6 J: %s, d J = 0: %s" % (k_, same, cc.is_zero(dJ)), flush=True)
    # brackets
    names = list(Ts)
    Qv = {}
    for i, X in enumerate(names):
        for Y in names[i:]:
            Bc, den = bracket33(g, J_cochain(Ts[X]), J_cochain(Ts[Y]))
            v = cc.flatten(Bc, monos4, g.dim, 4) if Bc else np.zeros(len(monos4) * g.dim ** 4, dtype=np.int64)
            closed = cc.is_zero(cc.d_CE(mom, g, Bc, 4)) if Bc else True
            Qv[(X, Y)] = (v, den)
            print("  {J_%s, J_%s}: closed %s" % (X, Y, closed), flush=True)
    keep = keep | np.any(np.stack([v for v, d in Qv.values()], axis=1) != 0, axis=1)
    Bk, Fk = Bm[keep], Fv[keep]
    Mbase = np.concatenate([Bk, Fk[:, None]], axis=1)
    cert = la.certified_rank(Mbase)
    certB = la.certified_rank(Bk)
    print("  rank B^4 = %d, rank [B^4 | F_223] = %d  -> [F_223] != 0: %s" % (certB['rank'], cert['rank'],
          cert['rank'] > certB['rank']), flush=True)
    # test combinations (a,b,c)
    if g.name.startswith("gl"):
        N = g.mats[0].shape[0]
        tests = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, Fr(-2, N), 0), (1, Fr(-2, N), Fr(4, N * N)), (2, 1, 1), (0, 1, Fr(-2, N))]
    else:
        tests = [(1,)]
    for abc in tests:
        coef = dict(zip(names, abc))
        # (1/2){J,J} = sum_{X<=Y} w_XY {J_X,J_Y}, w = c_X c_Y (1/2 if X == Y)
        tot = None
        for (X, Y), (v, den) in Qv.items():
            w = Fr(coef[X]) * Fr(coef[Y]) * (Fr(1, 2) if X == Y else 1) / den
            if w == 0:
                continue
            term = [w * int(x) for x in v[keep]]
            tot = term if tot is None else [a_ + b_ for a_, b_ in zip(tot, term)]
        if tot is None:
            tot = [Fr(0)] * int(keep.sum())
        exact, _ = la.in_span(Bk, tot, certB) if certB['rank'] else (all(x == 0 for x in tot), None)
        ok, x = la.in_span(Mbase, tot, cert)
        lam = x[-1] * sc if ok else None      # F column is sc * F_223
        print("  J = %s: (1/2){J,J} exact: %s; = lam F_223 + exact: %s, lam = %s" % (
            dict((k_, str(v_)) for k_, v_ in coef.items()), exact, ok, lam), flush=True)
    print("  time %.0fs" % (time.time() - t0))


if __name__ == '__main__':
    for g in [lie.sl(2), lie.gl(2), lie.sl(3), lie.gl(3)]:
        run(g)
