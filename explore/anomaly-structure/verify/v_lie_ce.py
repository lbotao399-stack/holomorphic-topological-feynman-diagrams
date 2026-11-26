"""VERIFY (b2), (b3 colour part), (b4), (c4), (c5), (c6), (c7), (c8).  Own code: lie_eval.py, ce_eval.py.

Method: cochains modulo conservation are evaluated as functions at random integer arguments
(c_i, p_i), sum p_i = 0; d by its defining formula on these values; invariant 3-cochains spanned
by alt(mu tr_sigma) (first fundamental theorem for GL_N; Reynolds); ranks of value matrices exact
over Q.  A rank increase rank[D|w] > rank[D] of the value matrix is a certificate that w is not a
combination of the columns of D as functions.  Exactness: exact solve on samples, then the identity
is re-evaluated at fresh random samples.
"""
import sys, time, random, itertools
from fractions import Fraction as Fr
import numpy as np
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/anomaly-structure/verify')
import lie_eval as le
import ce_eval as ce

rng = random.Random(31415)
T00 = time.time()
THETA133 = (6, [(0, 5), (0, 1), (1, 2), (2, 5), (0, 3), (3, 4), (4, 5)], [1, 2, 3, 4])
THETA223 = (6, [(0, 1), (1, 5), (0, 2), (2, 5), (0, 3), (3, 4), (4, 5)], [1, 2, 3, 4])
BOX = (4, [(0, 1), (1, 2), (2, 3), (3, 0)], [0, 1, 2, 3])


def br2(p, q):
    return p[0] * q[1] - p[1] * q[0]


def det3(a, b, c):
    return (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]))


P133 = lambda ps: -3 * br2(ps[0], ps[1]) * br2(ps[2], ps[3])   # (60/pi) K_133
P223 = lambda ps: 2 * br2(ps[0], ps[1]) * br2(ps[2], ps[3])    # (60/pi) K_223
PF = lambda ps: br2(ps[0], ps[1]) * br2(ps[2], ps[3])          # F_G kernel
PBOX = lambda ps: det3(ps[0], ps[1], ps[2])                   # 6 K_box up to sign


def colour(alg):
    W1, s1, _ = le.orientation_rule(alg, *THETA133)
    W2, s2, _ = le.orientation_rule(alg, *THETA223)
    Wb, sb, _ = le.orientation_rule(alg, *BOX)
    return W1, W2, Wb


# ------------------------------------------------------------------ part 1: (b2)
print("=== 1. (b2) antisymmetric parts of the colour weights; invariant 4-forms ===", flush=True)


def gl1gl2_form(X, Y):
    # gl2 + gl1 block diag (3x3): kappa = tr(AB) + 2 tr A tr B + 3 b b' + (tr A) b' + b (tr B)
    A, B = X[:2, :2], Y[:2, :2]
    a, b = int(X[2, 2]), int(Y[2, 2])
    return int(np.trace(A @ B)) + 2 * int(np.trace(A)) * int(np.trace(B)) + 3 * a * b + int(np.trace(A)) * b + a * int(np.trace(B))


def so4_form(X, Y):
    # sl2 + sl2 block diag with levels 1 and 3
    return int(np.trace(X[:2, :2] @ Y[:2, :2])) + 3 * int(np.trace(X[2:, 2:] @ Y[2:, 2:]))


def gl3_form(X, Y):
    return int(np.trace(X @ Y)) + 2 * int(np.trace(X)) * int(np.trace(Y))


sl2m = [le.E(2, 0, 1), le.E(2, 1, 0), le.E(2, 0, 0) - le.E(2, 1, 1)]
gl2m = [le.E(2, i, j) for i in range(2) for j in range(2)]
algs1 = [
    ("sl2", le.sl(2), 0),
    ("gl2", le.gl(2), 1),
    ("sl3", le.sl(3), 0),
    ("gl3 (twisted basis, kappa = tr XY + 2 trX trY)", le.Alg("gl3t", [sum(int(u) * m for u, m in zip(row, [le.E(3, i, j) for i in range(3) for j in range(3)])) for row in le.unimodular(9, rng)], gl3_form), 1),
    ("so4 = sl2+sl2 (levels 1, 3)", le.block_sum([sl2m, sl2m], [2, 2], "sl2+sl2", so4_form), 0),
    ("so5", le.so(5), 0),
    ("gl2+gl1 (non-diagonal invariant form)", le.block_sum([gl2m, [np.array([[1]])]], [2, 1], "gl2+gl1", gl1gl2_form), 2),
]


def rank_mod_p(rows, p=2147483647):
    M = [[x % p for x in r] for r in rows]
    n = len(M)
    if n == 0:
        return 0
    m = len(M[0])
    r = 0
    for c in range(m):
        piv = next((i for i in range(r, n) if M[i][c]), None)
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        inv = pow(M[r][c], p - 2, p)
        M[r] = [(x * inv) % p for x in M[r]]
        for i in range(n):
            if i != r and M[i][c]:
                f = M[i][c]
                M[i] = [(x - f * y) % p for x, y in zip(M[i], M[r])]
        r += 1
    return r


def invariant_alt4_nullity_modp(alg):
    d = alg.dim
    basis = list(itertools.combinations(range(d), 4))
    index = {b: i for i, b in enumerate(basis)}
    rows = []
    fu = alg.fu_int
    for e in range(d):
        # (e . w)(a,b,c,d) = - sum_slots w(.., [T_e, T_slot], ..)
        for (a, b, c, dd) in basis:
            row = [0] * len(basis)
            tup = [a, b, c, dd]
            for slot in range(4):
                for x in np.nonzero(fu[e, tup[slot]])[0]:
                    new = tup[:]
                    new[slot] = int(x)
                    if len(set(new)) < 4:
                        continue
                    row[index[tuple(sorted(new))]] += ce.perm_sign(sorted(range(4), key=lambda i: new[i])) * int(fu[e, tup[slot], x])
            if any(row):
                rows.append(row)
    from sympy.polys.matrices import DomainMatrix
    from sympy import QQ
    nul_p = len(basis) - rank_mod_p(rows)
    DM = DomainMatrix([[QQ(int(x)) for x in r] for r in rows], (len(rows), len(basis)), QQ)
    return nul_p, len(basis) - DM.rank()


for name, g, zdim in algs1:
    t0 = time.time()
    W1, W2, Wb = colour(g)
    a1 = not le.alt_part(W1).any()
    a2 = not le.alt_part(W2).any()
    # centrality: contract one external slot with a central element (if any)
    cz = g.center_dim()
    nul = invariant_alt4_nullity_modp(g) if g.dim <= 10 else None
    print("%-48s dim %2d centre %d: Alt(W_133) = 0: %s, Alt(W_223) = 0: %s; dim(Lambda^4 g*)^g (mod p, exact) = %s (predicted from H^4(s) = 0: dim z*(x)H^3(s) + dim Lambda^4 z* = %d) (%.1fs)" % (
        name, g.dim, cz, a1, a2, nul, zdim, time.time() - t0), flush=True)
# exhibit the invariant 4-form tr(1) tr([2,3]4) alternated on gl2+gl1 / gl3 (lower bound for the nullity)


# ------------------------------------------------------------------ part 2: (b4)
def value_matrix(alg, funcs, target_funcs, npts_args, nsamp, H, d_of=True):
    """rows: samples; columns: d(func) for func in funcs (if d_of) then target funcs (plain)."""
    rows = []
    for _ in range(nsamp):
        args = ce.random_args(alg, rng, npts_args, H)
        row = [ce.d_eval(alg, f, args) for f in funcs] if d_of else [f(args) for f in funcs]
        row += [t(args) for t in target_funcs]
        rows.append(row)
    return rows


def exactness(alg, eta_funcs, targets, H, nsamp, label, verify=25):
    t0 = time.time()
    rows = value_matrix(alg, eta_funcs, [tf for _, tf in targets], 4, nsamp, H)
    m = len(eta_funcs)
    D = [r[:m] for r in rows]
    rD = ce.rank_Q(D)
    out = {}
    for k, (tn, tf) in enumerate(targets):
        col = [r[m + k] for r in rows]
        rDw = ce.rank_Q([d + [c] for d, c in zip(D, col)])
        res = {"rank D": rD, "rank [D|w]": rDw, "exact": rDw == rD, "w zero on samples": not any(col)}
        if rDw == rD and any(col):
            x = ce.solve_Q(D, col)
            ok = True
            for _ in range(verify):
                args = ce.random_args(alg, rng, 4, H)
                lhs = sum(xj * ce.d_eval(alg, f, args) for xj, f in zip(x, eta_funcs) if xj != 0)
                ok &= (lhs == tf(args))
            res["solution re-verified at %d fresh samples" % verify] = ok
        out[tn] = res
    print("  %s: %s (%.0fs)" % (label, out, time.time() - t0), flush=True)
    return out


def closedness(alg, cochain, H, nsamp):
    vals = []
    for _ in range(nsamp):
        args = ce.random_args(alg, rng, 5, H)
        vals.append(ce.d_eval(alg, cochain, args))
    return all(v == 0 for v in vals)


def nonzero(alg, cochain, H, nsamp=6):
    return any(cochain(ce.random_args(alg, rng, 4, H)) != 0 for _ in range(nsamp))


print("\n=== 2. (b4) one-loop box at H = 3, weight (1,1,1) ===", flush=True)
mono111 = ce.monomials(2, (1, 1, 1))
for g in [le.sl(2), le.gl(2), le.sl(3), le.gl(3), le.sl(4)]:
    W1, W2, Wb = colour(g)
    wbox = ce.TensorCochain(PBOX, Wb)
    t0 = time.time()
    cl = closedness(g, wbox, 3, 6)
    nz = nonzero(g, wbox, 3)
    print(" %s: omega_box nonzero %s, d omega_box = 0 at 6 random samples: %s (%.0fs)" % (g.name, nz, cl, time.time() - t0), flush=True)
    etas = [ce.TraceCochain(g, mu, cyc, 3) for mu in mono111 for cyc in ce.trace_types(3)]
    exactness(g, etas, [("omega_box", wbox)], 3, 70, g.name + " exactness of omega_box in d(C^3_inv)")


# ------------------------------------------------------------------ part 3: (c4)
print("\n=== 3. (c4) trace-basis forms of W_133, W_223 (random integer matrices) ===", flush=True)


def T(cyc_word_list, X):
    """product of traces; cyc_word_list like [(0,1,2,3)] or [(0,1),(2,3)] (0-based)."""
    v = 1
    for w in cyc_word_list:
        m = np.eye(X[0].shape[0], dtype=object)
        for j in w:
            m = m.dot(X[j])
        v *= m.trace()
    return v


def W133_sl(X, N):
    return (2 * T([(0, 1, 2, 3)], X) + (2 + N * N) * T([(0, 1, 3, 2)], X) - 4 * T([(0, 2, 1, 3)], X) + (2 + N * N) * T([(0, 2, 3, 1)], X)
            - 4 * T([(0, 3, 1, 2)], X) + 2 * T([(0, 3, 2, 1)], X) + 6 * N * T([(0, 1), (2, 3)], X))


def W223_sl(X, N):
    return (-2 * T([(0, 1, 2, 3)], X) - 2 * T([(0, 1, 3, 2)], X) + 4 * T([(0, 2, 1, 3)], X) - 2 * T([(0, 2, 3, 1)], X)
            + 4 * T([(0, 3, 1, 2)], X) - 2 * T([(0, 3, 2, 1)], X) - 4 * N * T([(0, 1), (2, 3)], X)
            + 2 * N * T([(0, 2), (1, 3)], X) + 2 * N * T([(0, 3), (1, 2)], X))


def Wtot_gl(X, N):
    t = lambda *ws: T(list(ws), X)
    return (-10 * t((0, 1, 2, 3)) - (10 + 3 * N * N) * t((0, 1, 3, 2)) + 20 * t((0, 2, 1, 3)) - (10 + 3 * N * N) * t((0, 2, 3, 1))
            + 20 * t((0, 3, 1, 2)) - 10 * t((0, 3, 2, 1)) - 26 * N * t((0, 1), (2, 3)) + 4 * N * t((0, 2), (1, 3)) + 4 * N * t((0, 3), (1, 2))
            + 3 * N * (t((0,), (1, 2, 3)) + t((0,), (1, 3, 2)) + t((0, 1, 2), (3,)) + t((0, 1, 3), (2,)) + t((0, 2, 1), (3,))
                       + t((0, 2, 3), (1,)) + t((0, 3, 1), (2,)) + t((0, 3, 2), (1,)))
            + 20 * (t((0,), (1,), (2, 3)) + t((0, 1), (2,), (3,)))
            - 10 * (t((0,), (1, 2), (3,)) + t((0,), (1, 3), (2,)) + t((0, 2), (1,), (3,)) + t((0, 3), (1,), (2,))))


def Wtensor_eval(W, cs):
    r = W
    for c in reversed(cs):
        r = np.tensordot(r, np.array(c, dtype=object), axes=([r.ndim - 1], [0]))
    return r


for N in (2, 3, 4):
    for kind in ("sl", "gl"):
        g = le.sl(N) if kind == "sl" else le.gl(N)
        t0 = time.time()
        W1, W2, Wb = colour(g)
        sc = Fr(g.kinv_scale ** 7)
        ok1 = ok2 = okt = oks = True
        for _ in range(6):
            cs = [[rng.randint(-3, 3) for _ in range(g.dim)] for _ in range(4)]
            X = [g.mat(c).astype(object) for c in cs]
            w1 = Fr(int(Wtensor_eval(W1, cs))) / sc
            w2 = Fr(int(Wtensor_eval(W2, cs))) / sc
            if kind == "sl":
                ok1 &= (w1 == W133_sl(X, N))
                ok2 &= (w2 == W223_sl(X, N))
            okt &= (-3 * w1 + 2 * w2 == (Wtot_gl(X, N)))
            if N == 2 and kind == "sl":
                k = lambda i, j: (X[i].dot(X[j])).trace()
                oks &= (-3 * w1 + 2 * w2 == -104 * k(0, 1) * k(2, 3) + 16 * k(0, 2) * k(1, 3) + 40 * k(0, 3) * k(1, 2))
        msg = " %s: W_tot = -3W_133 + 2W_223 equals the gl_N trace formula: %s" % (g.name, okt)
        if kind == "sl":
            msg += "; W_133 = sl_N formula: %s; W_223 = sl_N formula: %s" % (ok1, ok2)
        if N == 2 and kind == "sl":
            msg += "; sl2 W_tot = -104 k12 k34 + 16 k13 k24 + 40 k14 k23: %s" % oks
        print(msg + " (%.0fs)" % (time.time() - t0), flush=True)


# ------------------------------------------------------------------ part 4: (c5), (c6), (c7)
print("\n=== 4. (c5) closedness, (c6) exactness, (c7) class relation at H = 2, weight (2,2) ===", flush=True)
mono22 = ce.monomials(2, (2, 2))
Bdims = {}
for g in [le.sl(2), le.gl(2), le.sl(3), le.gl(3), le.sl(4), le.gl(4), le.gl(5), le.sl(5)]:
    t0 = time.time()
    W1, W2, Wb = colour(g)
    w133 = ce.TensorCochain(P133, W1)
    w223 = ce.TensorCochain(P223, W2)
    F223 = ce.TensorCochain(PF, W2)
    F133 = ce.TensorCochain(PF, W1)
    wsum = lambda a: w133(a) + w223(a)
    comb = lambda a: 2 * w133(a) - 3 * w223(a)          # 2 (omega_133 - 3/2 omega_223)
    Fsum = lambda a: F133(a) + F223(a)
    ns = 4 if g.dim <= 9 else 2
    c133 = closedness(g, w133, 2, ns)
    c223 = closedness(g, w223, 2, ns)
    print(" %s: d omega_133 = 0: %s, d omega_223 = 0: %s (%d random 5-point samples each); omega_133, omega_223 nonzero: %s, %s (%.0fs)" % (
        g.name, c133, c223, ns, nonzero(g, w133, 2), nonzero(g, w223, 2), time.time() - t0), flush=True)
    etas = [ce.TraceCochain(g, mu, cyc, 3) for mu in mono22 for cyc in ce.trace_types(3)]
    out = exactness(g, etas, [("omega", wsum), ("omega_133", w133), ("omega_223", w223), ("omega_133-(3/2)omega_223", comb),
                              ("F_133+F_223", Fsum), ("F_223", F223)], 2, 90, g.name + " exactness")
    Bdims[g.name] = out["omega"]["rank D"]


# ------------------------------------------------------------------ part 5: (c8)
print("\n=== 5. (c8) invariant cohomology at weight (2,2), gl_5 and gl_6 (explicit matrices) ===", flush=True)
mono4 = ce.monomials(3, (2, 2))
for N in (5, 6):
    t0 = time.time()
    g = le.gl(N)
    funcs = [ce.TraceCochain(g, mu, cyc, 4) for mu in mono4 for cyc in ce.trace_types(4)]
    # value matrix on 4-point samples: dim C^4
    S4 = [ce.random_args(g, rng, 4, 2) for _ in range(70)]
    rows = [[f(a) for f in funcs] for a in S4]
    r1 = ce.rank_Q(rows)
    S4b = [ce.random_args(g, rng, 4, 2) for _ in range(20)]
    rows2 = rows + [[f(a) for f in funcs] for a in S4b]
    r1b = ce.rank_Q(rows2)
    # pick basis columns: greedy independent columns
    cols = ce.pivots_Q(rows)
    assert len(cols) == r1
    basis = [funcs[j] for j in cols]
    S5 = [ce.random_args(g, rng, 5, 2) for _ in range(60)]
    drows = [[ce.d_eval(g, f, a) for f in basis] for a in S5]
    rd = ce.rank_Q(drows)
    S5b = [ce.random_args(g, rng, 5, 2) for _ in range(15)]
    rdb = ce.rank_Q(drows + [[ce.d_eval(g, f, a) for f in basis] for a in S5b])
    print(" gl%d: dim C^4_inv(2,2) = %d (70 samples), %d (90 samples); rank d: C^4 -> C^5 = %d (60 samples), %d (75 samples); dim Z^4 = %d (%.0fs)" % (
        N, r1, r1b, rd, rdb, r1 - rdb, time.time() - t0), flush=True)
    if N == 5:
        etas = [ce.TraceCochain(g, mu, cyc, 3) for mu in mono22 for cyc in ce.trace_types(3)]
        rowsB = value_matrix(g, etas, [], 4, 70, 2)
        rB = ce.rank_Q(rowsB)
        print("   gl5: dim B^4 = rank d: C^3 -> C^4 = %d ; dim H^4 = %d" % (rB, r1 - rdb - rB), flush=True)
print("total %.0fs" % (time.time() - T00))
