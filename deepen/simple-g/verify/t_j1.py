"""J_1 = [q0 q1] T(Y0,Y1,Y2), T in (S^3 g*)^g.
Analytic: dJ(x0..x3) = sum_{i<j} (-1)^{i+j} [p_k p_l] T([X_i,X_j],X_k,X_l) (since [p_i+p_j, p_k] = [p_k p_l]);
the coefficient of u = [p0p1] is T([X0,X2],X1,X3) + T([X1,X2],X0,X3) - T([X2,X3],X0,X1), which is minus the
ad_{X2}-invariance identity of T applied to (X0,X1,X3), hence 0; same for v, w.
Numerical (exact): sl3, sl4, sl5 with T = d, d(X,Y,Z) = tr(X(YZ+ZY)) in the defining representation:
symmetric, invariant (checked on all basis entries), nonzero; dJ at 4 configurations (rank 3 on (u,v,w)): full
4-tensors.  Control: a random symmetric non-invariant T gives dJ != 0."""
import itertools, sys
import numpy as np
import valg, vcoch as vc

rng = np.random.default_rng(23)


def dJ(c, T, P):
    n = c.shape[0]
    TT4 = vc.ein('abe,ecd->abcd', c, T)          # T([e_a,e_b],e_c,e_d)
    out = np.zeros((n,) * 4, dtype=np.int64)
    for i, j in itertools.combinations(range(4), 2):
        k, l = [r for r in range(4) if r not in (i, j)]
        q0 = (P[i][0] + P[j][0], P[i][1] + P[j][1])
        cf = (-1) ** (i + j) * vc.brk(q0, P[k])
        if cf:
            out += cf * np.transpose(TT4, list(np.argsort([i, j, k, l])))
    return out


for N in (3, 4, 5):
    A = valg.sl(N)
    M = A.mats
    n = A.n
    d = np.zeros((n, n, n), dtype=np.int64)
    for a in range(n):
        for b in range(n):
            for c_ in range(n):
                d[a, b, c_] = int(np.trace(M[a] @ (M[b] @ M[c_] + M[c_] @ M[b])))
    sym = all(np.array_equal(d, np.transpose(d, p)) for p in itertools.permutations(range(3)))
    # invariance: d([e_a,X],Y,Z) + d(X,[e_a,Y],Z) + d(X,Y,[e_a,Z]) = 0 for all basis a
    inv = True
    for a in range(n):
        ad = A.c[a]                                   # [e_a, e_b] = sum_m ad[b, m] e_m
        t = np.einsum('bm,mcd->bcd', ad, d) + np.einsum('cm,bmd->bcd', ad, d) + np.einsum('dm,bcm->bcd', ad, d)
        inv &= not t.any()
    while True:
        Ps = [vc.config4(rng) for _ in range(4)]
        if valg.exact_rank([[vc.brk(P[0], P[1]), vc.brk(P[0], P[2]), vc.brk(P[1], P[2])] for P in Ps]) == 3:
            break
    zero = all(not dJ(A.c, d, P).any() for P in Ps)
    R = rng.integers(-2, 3, (n, n, n))
    Rs = sum(np.transpose(R, p) for p in itertools.permutations(range(3))).astype(np.int64)
    ctrl = all(not dJ(A.c, Rs, P).any() for P in Ps)
    print("sl%d: d symmetric %s, invariant %s, nonzero %s ; d J_1 = 0 (4 configs, all entries): %s ; "
          "control random symmetric T: d J = 0: %s" % (N, sym, inv, bool(d.any()), zero, ctrl))
