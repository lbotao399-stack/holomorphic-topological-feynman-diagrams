"""VERIFY (c6), (c7) for all N.  Own code.

Formal argument.  Invariant cochains of gl_N (sl_N) are restrictions of formal trace-monomial
cochains alt(mu tr_sigma) (sl_N: sigma without fixed points).  The formal d is N-independent
(inserting [X_i, X_j] into a trace word produces no N).  The cochains are polynomial in N:
omega(N) = sum_k N^k omega_k, the omega_k formal (trace formulas of (c4), verified in v_lie_ce.py at
N = 2, 3, 4 on top of the degree bound <= 2).  At an N0 where the 4-argument (and 3-argument)
monomials are linearly independent functions (checked below), evaluation is injective on formal
cochains, so for every N >= 4:  omega(N) exact  <=>  (1, N, N^2) lies in
R = {(a,b,c) : a omega_0 + b omega_1 + c omega_2 in d(C^3)} computed at N0.
gl_N vs sl_N: W_G vanishes on the centre, so the gl_N cochains are pull-backs along the Lie
homomorphism pi: gl_N -> sl_N, X -> X - tr(X)/N; restriction along sl_N -> gl_N and pull-back along
pi both commute with d, hence a gl_N cochain of this form is exact iff its sl_N restriction is.
"""
import sys, time, random, itertools
from fractions import Fraction as Fr
import numpy as np
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/anomaly-structure/verify')
import lie_eval as le
import ce_eval as ce

rng = random.Random(2718)


def br2(p, q):
    return p[0] * q[1] - p[1] * q[0]


PF = lambda ps: br2(ps[0], ps[1]) * br2(ps[2], ps[3])


def T(ws, X):
    v = 1
    for w in ws:
        m = X[w[0]]
        for j in w[1:]:
            m = m @ X[j]
        v *= int(np.trace(m))
    return v


class FuncCochain:
    def __init__(self, P, Wf):
        self.P, self.Wf = P, Wf

    def __call__(self, args):
        tot = 0
        for s, sg in ce.PERMS[4]:
            pv = self.P([args[s[i]][1] for i in range(4)])
            if pv:
                tot += sg * pv * self.Wf([args[s[i]][2] for i in range(4)])
        return tot


t = lambda *ws: (lambda X: T(list(ws), X))


def lin(*terms):
    return lambda X: sum(c * f(X) for c, f in terms)


# W_tot = -3 W_133 + 2 W_223 (c4), split by powers of N
Wtot = {
    0: lin((-10, t((0, 1, 2, 3))), (-10, t((0, 1, 3, 2))), (20, t((0, 2, 1, 3))), (-10, t((0, 2, 3, 1))), (20, t((0, 3, 1, 2))), (-10, t((0, 3, 2, 1))),
           (20, t((0,), (1,), (2, 3))), (20, t((0, 1), (2,), (3,))),
           (-10, t((0,), (1, 2), (3,))), (-10, t((0,), (1, 3), (2,))), (-10, t((0, 2), (1,), (3,))), (-10, t((0, 3), (1,), (2,)))),
    1: lin((-26, t((0, 1), (2, 3))), (4, t((0, 2), (1, 3))), (4, t((0, 3), (1, 2))),
           (3, t((0,), (1, 2, 3))), (3, t((0,), (1, 3, 2))), (3, t((0, 1, 2), (3,))), (3, t((0, 1, 3), (2,))), (3, t((0, 2, 1), (3,))),
           (3, t((0, 2, 3), (1,))), (3, t((0, 3, 1), (2,))), (3, t((0, 3, 2), (1,)))),
    2: lin((-3, t((0, 1, 3, 2))), (-3, t((0, 2, 3, 1)))),
}
# sl_N restrictions (c4)
W223sl = {
    0: lin((-2, t((0, 1, 2, 3))), (-2, t((0, 1, 3, 2))), (4, t((0, 2, 1, 3))), (-2, t((0, 2, 3, 1))), (4, t((0, 3, 1, 2))), (-2, t((0, 3, 2, 1)))),
    1: lin((-4, t((0, 1), (2, 3))), (2, t((0, 2), (1, 3))), (2, t((0, 3), (1, 2)))),
}
W133sl = {
    0: lin((2, t((0, 1, 2, 3))), (2, t((0, 1, 3, 2))), (-4, t((0, 2, 1, 3))), (2, t((0, 2, 3, 1))), (-4, t((0, 3, 1, 2))), (2, t((0, 3, 2, 1)))),
    1: lin((6, t((0, 1), (2, 3))),),
    2: lin((1, t((0, 1, 3, 2))), (1, t((0, 2, 3, 1)))),
}
WtotSl = {k: (lambda X, k=k: Wtot[k](X)) for k in Wtot}  # on traceless matrices the tr(X_i) terms vanish
Wsum_sl = {k: (lambda X, k=k: W133sl.get(k, lambda X: 0)(X) + W223sl.get(k, lambda X: 0)(X)) for k in (0, 1, 2)}


def faithful(alg, k, fixed_point_free):
    cyc = [c for c in ce.trace_types(k) if (not fixed_point_free or all(len(x) > 1 for x in c))]
    rows = []
    for _ in range(len(cyc) + 15):
        X = [alg.mat([rng.randint(-3, 3) for _ in range(alg.dim)]) for _ in range(k)]
        rows.append([ce.tr_eval(c, X) for c in cyc])
    return len(cyc), ce.rank_Q(rows)


def relation_space(D, cols):
    """R = {a : sum_k a_k cols_k in colspan(D)}; returns a basis (list of vectors)."""
    m = len(D[0])
    K = len(cols)
    A = [d + [c[i] for c in cols] for i, d in enumerate(D)]
    # null space of A (exact), project on the last K coordinates
    n = len(A)
    M = [[Fr(x) for x in r] for r in A]
    ncol = m + K
    piv = []
    r = 0
    for c in range(ncol):
        p = next((i for i in range(r, n) if M[i][c] != 0), None)
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        pv = M[r][c]
        M[r] = [x / pv for x in M[r]]
        for i in range(n):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [x - f * y for x, y in zip(M[i], M[r])]
        piv.append(c)
        r += 1
    free = [c for c in range(ncol) if c not in piv]
    null = []
    for fc in free:
        v = [Fr(0)] * ncol
        v[fc] = Fr(1)
        for i, pc in enumerate(piv):
            v[pc] = -M[i][fc]
        null.append(v)
    proj = [v[m:] for v in null]
    # basis of the projection
    basis = []
    for v in proj:
        if ce.rank_Q(basis + [v]) > len(basis):
            basis.append(v)
    return basis


def analyse(alg, kind, nsamp=110):
    t0 = time.time()
    k4 = faithful(alg, 4, kind == "sl")
    k3 = faithful(alg, 3, kind == "sl")
    print(" %s: 4-argument trace monomials: %d, rank as functions %d; 3-argument: %d, rank %d" % (alg.name, k4[0], k4[1], k3[0], k3[1]), flush=True)
    mono22 = ce.monomials(2, (2, 2))
    cyc3 = [c for c in ce.trace_types(3) if (kind == "gl" or all(len(x) > 1 for x in c))]
    etas = [ce.TraceCochain(alg, mu, c, 3) for mu in mono22 for c in cyc3]
    if kind == "gl":
        targets = {"omega": [FuncCochain(PF, Wtot[k]) for k in (0, 1, 2)]}
    else:
        targets = {"omega": [FuncCochain(PF, WtotSl[k]) for k in (0, 1, 2)],
                   "F_223": [FuncCochain(PF, W223sl[k]) for k in (0, 1)],
                   "F_133+F_223": [FuncCochain(PF, Wsum_sl[k]) for k in (0, 1, 2)]}
    D, vals = [], {nm: [] for nm in targets}
    for _ in range(nsamp):
        args = ce.random_args(alg, rng, 4, 2)
        D.append([ce.d_eval(alg, f, args) for f in etas])
        for nm, fs in targets.items():
            vals[nm].append([f(args) for f in fs])
    rD = ce.rank_Q(D)
    print("   rank d(C^3) = %d (%d spanning cochains, %d samples)" % (rD, len(etas), nsamp))
    for nm, v in vals.items():
        cols = [[row[k] for row in v] for k in range(len(v[0]))]
        R = relation_space(D, cols)
        rk = ce.rank_Q([d + [c[i] for c in cols] for i, d in enumerate(D)])
        print("   %s = sum_k N^k %s_k: rank [D | %s_0..] = %d; relation space R (coefficient vectors a with sum a_k %s_k exact): %s" % (
            nm, nm, nm, rk, nm, [[str(x) for x in b] for b in R]), flush=True)
    print("   (%.0fs)" % (time.time() - t0))


print("=== formal-in-N analysis at N0 = 5 (sl and gl), cross-check at N0 = 6 (sl) ===", flush=True)
analyse(le.sl(5), "sl")
analyse(le.gl(5), "gl")
analyse(le.sl(6), "sl", nsamp=80)
