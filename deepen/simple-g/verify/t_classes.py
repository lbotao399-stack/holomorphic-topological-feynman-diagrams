"""Decisive computation (verifier's own code): classes of F_133, F_223, -3F_133 +- 2F_223 relative to beta.
Full coefficient 4-tensors (no fixed slot), at 6 momentum configurations P_1..P_6 whose 6 x 6 matrix
[mu(P_j)], mu in (u^2,v^2,w^2,uv,uw,vw), is invertible (exact det printed), so a cochain in
span(mu) (x) tensors vanishes iff its values at P_1..P_6 vanish.  A 7th configuration checks on 3000 random
entries that the values are interpolated by the 6 (i.e. the cochains lie in that span).
Exactness criterion (Lemma 1 + Lemma 2 of the report, re-derived in VERIFY.txt): X in Z^4_inv is exact iff X in C beta.
All arithmetic exact (vcoch.mm bound rule; final comparisons in Python integers)."""
import sys, time
from fractions import Fraction as Fr
import numpy as np
import sympy as sp
import valg, vcoch as vc


def as_obj(T):
    return T.astype(object) if T.dtype != object else T


def ratio_or_none(X, B, i0):
    """X = c B entrywise (exact)?  i0: an index with B[i0] != 0.  Returns (c, witness)."""
    c = Fr(int(X[i0]), int(B[i0]))
    bnd = float(np.abs(X).max()) * abs(c.denominator) + float(np.abs(B).max()) * abs(c.numerator)
    if X.dtype != object and B.dtype != object and bnd < 2.0 ** 62:
        D = X.astype(np.int64) * c.denominator - B.astype(np.int64) * c.numerator
    else:
        D = as_obj(X) * c.denominator - as_obj(B) * c.numerator
    nz = np.flatnonzero(D != 0)
    if nz.size == 0:
        return c, None
    j = int(nz[0])
    return None, (i0, j, int(B[i0]) * int(X[j]) - int(X[i0]) * int(B[j]))


rng = np.random.default_rng(2026)
for name in sys.argv[1:]:
    t0 = time.time()
    A = valg.build(name)
    S = vc.Setup(A)
    W = {G: S.W_full(G) for G in ('133', '223')}
    tW = time.time() - t0
    n = A.n
    # momentum configurations
    while True:
        Ps = [vc.config4(rng) for _ in range(7)]
        E = sp.Matrix([vc.quad_basis_eval(P) for P in Ps[:6]])
        if E.det() != 0:
            break
    detE = E.det()
    cols = {'beta': [], 'F133': [], 'F223': []}
    sub = rng.choice(n ** 4, size=min(3000, n ** 4), replace=False)
    subvals = {k: [] for k in cols}
    for j, P in enumerate(Ps):
        F1 = vc.alt_tensor(W['133'], P).reshape(-1)
        F2 = vc.alt_tensor(W['223'], P).reshape(-1)
        Bt = vc.beta_tensor(S.TT, P).reshape(-1)
        for k, T in (('beta', Bt), ('F133', F1), ('F223', F2)):
            subvals[k].append([int(T[i]) for i in sub])
            if j < 6:
                cols[k].append(T)
    # interpolation check on the subset: values at P_7 = sum_mu mu(P_7) C_mu, C = E^{-1} values
    Einv = E.inv()
    m7 = vc.quad_basis_eval(Ps[6])
    interp_ok = True
    for k in cols:
        V = sp.Matrix(subvals[k][:6])        # 6 x |sub|
        C = Einv * V
        pred = sp.Matrix([m7]) * C
        interp_ok &= all(pred[0, i] == subvals[k][6][i] for i in range(len(sub)))
    vec = {k: np.concatenate(cols[k]) for k in cols}
    s7 = S.s ** 7
    B = vec['beta']
    nzB = np.flatnonzero(B != 0)
    i0 = int(nzB[0])
    out = []
    out.append("%s: dim %d, kappa = Killing/%d, s = %d, det[mu(P_j)] = %s, interpolation at P_7 (3000 entries): %s"
               % (name, n, S.g, S.s, detE, interp_ok))
    out.append("   nonzero: beta %s, F133 %s, F223 %s" % (bool(B.any()), bool(vec['F133'].any()), bool(vec['F223'].any())))
    combos = [('F133', vec['F133']), ('F223', vec['F223']), ('F133+F223', vec['F133'] + vec['F223']),
              ('w133+w223 = -3F133+2F223', -3 * vec['F133'] + 2 * vec['F223']),
              ('w133-w223 = -3F133-2F223', -3 * vec['F133'] - 2 * vec['F223'])]
    for label, X in combos:
        c, wit = ratio_or_none(X, B, i0)
        if c is not None:
            cK = c / s7 / S.g ** 2          # X_K = X_kappa / g, beta_K = g beta_kappa
            out.append("   %-26s = c beta exactly on all %d entries: c(kappa) = %s, c(Killing) = %s"
                       % (label, X.size, c / s7, cK))
        else:
            out.append("   %-26s NOT in C beta: 2x2 minor beta[i0] X[j] - X[i0] beta[j] = %s (i0 = %d, j = %d)"
                       % (label, wit[2], wit[0], wit[1]))
    out.append("   arithmetic: %s ; time W %.1fs, total %.1fs" % (dict(vc.STATS), tW, time.time() - t0))
    print("\n".join(out), flush=True)
