"""Direct computation of C^4_inv and of d on it (sl2, sl3, so5, G2).
Basis of (g^{*(x)4})^g: rank-reduced pool of all slot permutations of kappa.kappa, kappa([.,.],[.,.]), W_133, W_223
and the trace tr(ad ad ad ad); its rank is compared with dim (g^{(x)4})^g from the Weyl formula (complete iff equal).
C^4_inv = span{alt(mu T)}: mu over the 21 products of two brackets of p0..p3, T over that basis; values at 6
configurations (injective, det printed in t_classes).  Then rank C^4_inv, rank together with beta, F_133, F_223,
and d of every spanning cochain (t_closed.d_of, full 5-tensors, 24 configurations)."""
import sys, itertools, time
import numpy as np
import sympy as sp
import valg, vcoch as vc
from t_closed import d_of, configs5

rng = np.random.default_rng(31)
WEYL_T4 = {'sl2': 3, 'sl3': 8, 'so5': 6, 'G2': 5}
for name in sys.argv[1:]:
    t0 = time.time()
    A = valg.build(name)
    S = vc.Setup(A)
    n = A.n
    kk = np.einsum('ab,cd->abcd', S.kap, S.kap)
    ad = A.ad
    box = np.einsum('aij,bjk,ckl,dli->abcd', ad, ad, ad, ad)
    pool0 = [kk, S.TT, S.W_full('133'), S.W_full('223'), box]
    if name.startswith('sl'):
        Mt = A.mats
        tr4 = np.array([[[[int(np.trace(Mt[a] @ Mt[b] @ Mt[c] @ Mt[d])) for d in range(n)] for c in range(n)]
                         for b in range(n)] for a in range(n)], dtype=np.int64)
        pool0.append(tr4)
    basis = []
    for T0 in pool0:
        for s in itertools.permutations(range(4)):
            T = np.transpose(T0, list(np.argsort(s)))
            cand = basis + [T]
            if valg.exact_rank([t.reshape(-1).tolist() for t in cand]) == len(cand):
                basis.append(T)
    rT = len(basis)
    while True:
        Ps = [vc.config4(rng) for _ in range(6)]
        if sp.Matrix([vc.quad_basis_eval(P) for P in Ps]).det() != 0:
            break
    brs = list(itertools.combinations(range(4), 2))
    prods = list(itertools.combinations_with_replacement(range(6), 2))

    def poly(m):
        m1, m2 = m
        return lambda Q: vc.brk(Q[brs[m1][0]], Q[brs[m1][1]]) * vc.brk(Q[brs[m2][0]], Q[brs[m2][1]])

    def alt_gen(T, P, pl):
        out = np.zeros_like(T)
        for s in itertools.permutations(range(4)):
            k = vc.perm_sign(s) * pl([P[s[0]], P[s[1]], P[s[2]], P[s[3]]])
            if k:
                out = out + k * np.transpose(T, list(np.argsort(s)))
        return out
    cols, gens = [], []
    for T in basis:
        for m in prods:
            v = np.concatenate([alt_gen(T, P, poly(m)).reshape(-1) for P in Ps])
            cand = cols + [v]
            if valg.exact_rank(np.stack(cand, axis=1).tolist()) == len(cand):
                cols.append(v)
                gens.append((T, m))
    dimC = len(cols)
    B = np.concatenate([vc.beta_tensor(S.TT, P).reshape(-1) for P in Ps])
    F1 = np.concatenate([vc.alt_tensor(S.W_full('133'), P).reshape(-1) for P in Ps])
    F2 = np.concatenate([vc.alt_tensor(S.W_full('223'), P).reshape(-1) for P in Ps])
    M = np.stack(cols, axis=1)
    rB = valg.exact_rank(np.concatenate([M, B[:, None]], axis=1).tolist())
    rF = valg.exact_rank(np.concatenate([M, B[:, None], F1[:, None], F2[:, None]], axis=1).tolist())
    rbF = valg.exact_rank(np.stack([B, F1, F2], axis=1).tolist())
    Qs, _ = configs5(rng)
    dvecs = []
    for (T, m) in gens:
        dvecs.append(np.concatenate([d_of(A.c, lambda q: alt_gen(T, q, poly(m)), Q, 'full').reshape(-1) for Q in Qs]))
    rd = valg.exact_rank(np.stack(dvecs, axis=1).tolist())
    dzero = (rd == 0)
    print("%-4s rank of invariant 4-tensor pool = %d (Weyl dim (g^(x)4)^g = %d) ; dim C^4_inv = %d ; rank[C^4_inv, beta] = %d ; "
          "rank[C^4_inv, beta, F133, F223] = %d ; rank[beta, F133, F223] = %d ; rank of d on C^4_inv (24 configs, full 5-tensors) = %d ; dim Z^4_inv = %d ; dim H^4_inv = %d  (%.1fs)"
          % (name, rT, WEYL_T4[name], dimC, rB, rF, rbF, rd, dimC - rd, dimC - rd - 1, time.time() - t0), flush=True)
