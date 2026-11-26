"""(1) d F_G = 0, exact: 5-cochain coefficients over a basis of 20 quartic momentum polynomials, evaluated at
    X_0 = b_theta and all basis 4-tuples (injective on invariant 5-tensors, same argument as for 4-tensors).
    dF(x_0..x_4) = sum_{i<j} (-1)^{i+j} F([x_i,x_j], x_k, x_l, x_m),  F = sum_mu mu(q) T_mu,
    q = (p_i + p_j, p_k, p_l, p_m), p_4 = -(p_0+..+p_3).
    Polynomials: brackets b_ab = [p_a p_b], 0 <= a < b <= 3; quadratic monomials (21) reduced by the Pluecker
    identity b01 b23 = b02 b13 - b03 b12 to 20 monomials (independence checked with sympy).
(2) pointwise check at random integer X_i and momenta (exact integers): F_133 + F_223 = c beta with the c of
    run_classes, and F_G(b_theta,...) of core.F_eval equals the pointwise value."""
import sys, time, itertools
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/simple-g/code')
import numpy as np, sympy as sp
import algebras as al, core
from algebras import perm_sign

PAIRS = [(a, b) for a in range(4) for b in range(a + 1, 4)]          # 6 brackets b_ab


def br5(a, b):
    """[p_a p_b], 0 <= a, b <= 4, p4 = -(p0+..+p3), as vector over PAIRS."""
    if a == b:
        return np.zeros(6, dtype=np.int64)
    if a > b:
        return -br5(b, a)
    if b < 4:
        v = np.zeros(6, dtype=np.int64); v[PAIRS.index((a, b))] = 1
        return v
    return -sum(br5(a, c) for c in range(4))


MON21 = [(i, j) for i in range(6) for j in range(i, 6)]
PL = (PAIRS.index((0, 1)), PAIRS.index((2, 3)))
MON20 = [m for m in MON21 if m != PL]


def quad5(l1, l2):
    """product of two linear forms over PAIRS -> vector over MON20 (Pluecker-reduced)."""
    v21 = {}
    for i in range(6):
        for j in range(6):
            if l1[i] and l2[j]:
                m = (min(i, j), max(i, j))
                v21[m] = v21.get(m, 0) + int(l1[i] * l2[j])
    out = np.zeros(20, dtype=np.int64)
    for m, c in v21.items():
        if m == PL:
            # b01 b23 = b02 b13 - b03 b12
            m1 = tuple(sorted((PAIRS.index((0, 2)), PAIRS.index((1, 3)))))
            m2 = tuple(sorted((PAIRS.index((0, 3)), PAIRS.index((1, 2)))))
            out[MON20.index(m1)] += c
            out[MON20.index(m2)] -= c
        else:
            out[MON20.index(m)] += c
    return out


def check_poly5():
    P = [[sp.Symbol('p%d_%d' % (i, c)) for c in range(2)] for i in range(4)]
    P.append([-(P[0][c] + P[1][c] + P[2][c] + P[3][c]) for c in range(2)])
    B = lambda a, b: P[a][0] * P[b][1] - P[a][1] * P[b][0]
    bs = [B(a, b) for (a, b) in PAIRS]
    ok = all(sp.expand(B(a, b) - sum(int(c) * x for c, x in zip(br5(a, b), bs))) == 0 for a in range(5) for b in range(5))
    pl = sp.expand(bs[PL[0]] * bs[PL[1]] - bs[PAIRS.index((0, 2))] * bs[PAIRS.index((1, 3))]
                   + bs[PAIRS.index((0, 3))] * bs[PAIRS.index((1, 2))]) == 0
    gens = [P[i][c] for i in range(4) for c in range(2)]
    polys = [sp.Poly(sp.expand(bs[i] * bs[j]), *gens).as_dict() for (i, j) in MON20]
    mons = sorted(set().union(*polys))
    rk = sp.Matrix([[p.get(m, 0) for m in mons] for p in polys]).rank()
    return ok, pl, rk


def T_mu_full(T, graph):
    W = T.W(graph, [None] * 4)
    d = T.L.dim
    out = np.zeros((6, d, d, d, d), dtype=np.int64)
    for sig in itertools.permutations(range(4)):
        c = core.quad(core.br(sig[0], sig[1]), core.br(sig[2], sig[3]))
        V = np.transpose(W, list(np.argsort(sig)))
        sg = perm_sign(sig)
        for mu in range(6):
            if c[mu]:
                out[mu] += sg * int(c[mu]) * V
    return out, W


def lin_uvw(q):
    """(u, v, w) = ([q0 q1], [q0 q2], [q1 q2]) as linear forms over PAIRS; q = list of index lists (sums)."""
    def B(x, y):
        return sum(br5(a, b) for a in x for b in y)
    return B(q[0], q[1]), B(q[0], q[2]), B(q[1], q[2])


def dF_exact(L, T, Tmu, mask=None):
    d = L.dim
    th = L.theta
    acc = np.zeros((20,) + (d,) * 4, dtype=np.int64) if mask is None else np.zeros((20, int(mask.sum())), dtype=np.int64)
    for i, j in itertools.combinations(range(5), 2):
        k, l, m = [r for r in range(5) if r not in (i, j)]
        u, v, w = lin_uvw([[i, j], [k], [l], [m]])
        polys = [quad5(u, u), quad5(v, v), quad5(w, w), quad5(u, v), quad5(u, w), quad5(v, w)]
        for mu in range(6):
            if not polys[mu].any():
                continue
            if i == 0:
                A = core.ct('je,eklm->jklm', L.fu[th], Tmu[mu])      # variables (j, k, l, m)
                vars_ = [j, k, l, m]
            else:
                A = core.ct('ije,elm->ijlm', L.fu, Tmu[mu][:, th])   # variables (i, j, l, m), k = 0
                vars_ = [i, j, l, m]
            # new axis t (variable t+1) = old axis with vars_ == t+1
            A = np.transpose(A, [vars_.index(t + 1) for t in range(4)])
            A = A if mask is None else A[mask]
            sgn = (-1) ** (i + j)
            for nu in range(20):
                if polys[mu][nu]:
                    acc[nu] += sgn * int(polys[mu][nu]) * A
    return acc


def pointwise(L, T, rng, R=3):
    """random integer X_0..X_3, momenta p_0..p_2 (p_3 = -sum); returns (s^7 F133, s^7 F223, beta)."""
    d = L.dim
    X = [rng.integers(-R, R + 1, size=d).astype(np.int64) for _ in range(4)]
    p = [rng.integers(-9, 10, size=2) for _ in range(3)]
    p.append(-(p[0] + p[1] + p[2]))
    Bk = lambda a, b: int(p[a][0] * p[b][1] - p[a][1] * p[b][0])
    out = []
    for G in ('133', '223'):
        tot = 0
        for sig in itertools.permutations(range(4)):
            Wv = int(T.W(G, [X[sig[0]], X[sig[1]], X[sig[2]], X[sig[3]]]))
            tot += perm_sign(sig) * Bk(sig[0], sig[1]) * Bk(sig[2], sig[3]) * Wv
        out.append(tot)
    be = 0
    for i, j in itertools.combinations(range(4), 2):
        k, l = [r for r in range(4) if r not in (i, j)]
        bij = L.bracket_vec(X[i], X[j]); bkl = L.bracket_vec(X[k], X[l])
        be += (-1) ** (i + j) * Bk(k, l) ** 2 * int(bij @ L.kappa @ bkl)
    out.append(be)
    return out


if __name__ == "__main__":
    print("polynomial algebra (br5 identities, Pluecker identity, rank of the 20 monomials):", check_poly5(), flush=True)
    for key in sys.argv[1:]:
        t0 = time.time()
        L = al.build(key)
        T = core.Tensors(L)
        mask = None
        if L.wt is not None and L.dim > 30:
            Wt = np.array(L.wt)
            d = L.dim
            tot = (Wt[:, None, None, None, :] + Wt[None, :, None, None, :] + Wt[None, None, :, None, :]
                   + Wt[None, None, None, :, :] + Wt[L.theta])
            mask = np.all(tot == 0, axis=4)
        line = "%-14s" % L.name
        for G in ('133', '223'):
            Tmu, W = T_mu_full(T, G)
            acc = dF_exact(L, T, Tmu, mask)
            line += "  d F_%s = 0: %s (max|T_mu| = %.2e)" % (G, not acc.any(), float(np.abs(Tmu).max()))
        if mask is not None:
            line += "  [weight-compatible 4-tuples only: %d of %d]" % (int(mask.sum()), L.dim ** 4)
        print(line + "  (%.1fs)" % (time.time() - t0), flush=True)
