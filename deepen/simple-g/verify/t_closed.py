"""Closedness d F_133 = d F_223 = 0 (verifier's own d), exact.
(dX)(x0..x4) = sum_{i<j} (-1)^{i+j} X([x_i,x_j], x_rest),  momenta (p_i+p_j, p_rest), sum p = 0.
Mode 'full' (small algebras): all 5-tensor entries.  Mode 'slot': x0 = e_b fixed, all basis 4-tuples x1..x4
(injective on invariant 5-tensors by the ideal argument).  Momentum configurations Q_1..Q_24 (5 vectors, sum 0):
the 21 products of two brackets of p0..p3 evaluated at the Q_j have rank 20 (= dim of their span, Pluecker), so a
5-cochain in that span vanishes iff it vanishes at the Q_j.  Controls: d beta = 0 (d^2 = 0); d of X' = alt([p0p1]^2 W_223)
(another invariant cochain); d of alt([p0p1][p2p3] R) with R a random non-invariant integer tensor (must be nonzero)."""
import sys, time, itertools
import numpy as np
import sympy as sp
import valg, vcoch as vc


def alt_gen(W, P, poly):
    out = None
    for s in itertools.permutations(range(4)):
        k = vc.perm_sign(s) * poly([P[s[0]], P[s[1]], P[s[2]], P[s[3]]])
        if k:
            V = k * np.transpose(W, list(np.argsort(s)))
            out = V if out is None else out + V
    return np.zeros_like(W) if out is None else out


def poly_F(Q):
    return vc.brk(Q[0], Q[1]) * vc.brk(Q[2], Q[3])


def poly_bad(Q):
    return vc.brk(Q[0], Q[1]) ** 2


def d_of(c, tensor_at, Q, mode, b=0):
    """tensor_at(q) -> full n^4 coefficient tensor of the 4-cochain at the 4-momentum config q."""
    n = c.shape[0]
    shape = (n,) * 5 if mode == 'full' else (n,) * 4
    out = np.zeros(shape, dtype=np.int64)
    for i, j in itertools.combinations(range(5), 2):
        rest = [r for r in range(5) if r not in (i, j)]
        q = [tuple(int(Q[i][t]) + int(Q[j][t]) for t in range(2))] + [Q[r] for r in rest]
        T = tensor_at(q)
        sg = (-1) ** (i + j)
        if mode == 'full':
            # C[a_i, a_j, a_rest] = sum_e c[a_i,a_j,e] T[e, a_rest]
            C = vc.mm(c.reshape(n * n, n), T.reshape(n, n ** 3)).reshape((n,) * 5)
            cur = [i, j] + rest
            out += sg * np.transpose(C, list(np.argsort(cur)))
        else:
            if i == 0:
                # sum_e c[b, a_j, e] T[e, a_k, a_l, a_m] ; axes (j, k, l, m) -> slots 1..4
                C = vc.mm(np.ascontiguousarray(c[b]), T.reshape(n, n ** 3)).reshape((n,) * 4)
                cur = [j] + rest
            else:
                # rest = (0, l, m): T[e, b, a_l, a_m]
                V = np.ascontiguousarray(T[:, b, :, :]).reshape(n, n * n)
                C = vc.mm(c.reshape(n * n, n), V).reshape((n,) * 4)
                cur = [i, j] + rest[1:]
            cur = [x - 1 for x in cur]
            out += sg * np.transpose(C, list(np.argsort(cur)))
    return out


def configs5(rng, N=24):
    while True:
        Qs = []
        for _ in range(N):
            p = [rng.integers(-9, 10, 2) for _ in range(4)]
            p.append(-(p[0] + p[1] + p[2] + p[3]))
            Qs.append([tuple(int(t) for t in x) for x in p])
        pairs = list(itertools.combinations(range(4), 2))
        prods = list(itertools.combinations_with_replacement(range(6), 2))
        M = [[vc.brk(Q[pairs[a][0]], Q[pairs[a][1]]) * vc.brk(Q[pairs[b][0]], Q[pairs[b][1]]) for (a, b) in prods]
             for Q in Qs]
        r = valg.exact_rank(M)
        if r == 20:
            return Qs, r


if __name__ == '__main__':
    rng = np.random.default_rng(5)
    mode_full = {'sl2', 'sl3', 'so5', 'sp4', 'G2'}
    for name in sys.argv[1:]:
        t0 = time.time()
        A = valg.build(name)
        S = vc.Setup(A)
        W = {G: S.W_full(G) for G in ('133', '223')}
        mode = 'full' if name in mode_full else 'slot'
        Qs, r = configs5(rng)
        res = {'dF133': True, 'dF223': True, 'dbeta': True, 'd(bad)': True, 'd(random)': True}
        Rnd = rng.integers(-2, 3, (A.n,) * 4).astype(np.int64)
        for Q in Qs:
            D1 = d_of(A.c, lambda q: alt_gen(W['133'], q, poly_F), Q, mode)
            D2 = d_of(A.c, lambda q: alt_gen(W['223'], q, poly_F), Q, mode)
            res['dF133'] &= not D1.any()
            res['dF223'] &= not D2.any()
            if Q is Qs[0] or Q is Qs[1]:
                Db = d_of(A.c, lambda q: vc.beta_tensor(S.TT, q), Q, mode)
                res['dbeta'] &= not Db.any()
                Dx = d_of(A.c, lambda q: alt_gen(W['223'], q, poly_bad), Q, mode)
                res['d(bad)'] &= not Dx.any()
                Dr = d_of(A.c, lambda q: alt_gen(Rnd, q, poly_F), Q, mode)
                res['d(random)'] &= not Dr.any()
        print("%-4s mode %s, 24 configs (bracket-product rank %d): d F_133 = 0: %s ; d F_223 = 0: %s ; "
              "control d beta = 0: %s ; control d alt([p0p1]^2 W_223) = 0: %s ; "
              "control d alt([p0p1][p2p3] R), R random non-invariant, = 0: %s  (%.1fs)"
              % (name, mode, r, res['dF133'], res['dF223'], res['dbeta'], res['d(bad)'], res['d(random)'], time.time() - t0), flush=True)
    print("arith:", vc.STATS)
