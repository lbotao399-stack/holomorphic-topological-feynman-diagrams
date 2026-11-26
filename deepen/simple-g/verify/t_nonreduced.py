"""Check of the reduction (Lemmas 1, 2) without the SL2 reduction: exactness of F_133, F_223, F_133+F_223 in the
full weight-(2,2) complex, eta over ALL G-invariant alternating 3-cochains of weight (2,2):
  eta = alt(mu(q0,q1) T(Y0,Y1,Y2)), mu over the 9 monomials a1^i b1^(2-i) a2^j b2^(2-j) (q0 = (a1,a2), q1 = (b1,b2)),
  T over a basis of (g^{*(x)3})^g:  so5, G2, sl2: {f} (dim 1, computed in t_inv_lin.py);  sl3: {f, d} with
  d(X,Y,Z) = tr(X(YZ+ZY)) (dim 2).
d eta(x0..x3) = sum_{i<j} (-1)^{i+j} Smu(p_i+p_j, p_k) T([X_i,X_j],X_k,X_l), Smu = sum_s mu(q_s0,q_s1) (T = f) or
sum_s sgn(s) mu(q_s0,q_s1) (T = d).  Values at 40 configurations whose evaluation matrix on the 36 weight-(2,2)
monomials of (p0,p1,p2) has rank 36 (so evaluation is injective on all weight-(2,2) 4-cochains).
Decision on a random subset of 2500 nonzero rows: if [d eta | X] has larger exact rank than [d eta] there, X is not
in span(d eta) (a subsystem is already inconsistent); if not, the solution is unique when rank = #columns and is
then verified on all rows."""
import sys, time, itertools
import numpy as np
import sympy as sp
import valg, vcoch as vc

rng = np.random.default_rng(99)


def S_poly(mu, q0, q1, sym):
    q2 = (-q0[0] - q1[0], -q0[1] - q1[1])
    qs = [q0, q1, q2]
    i, j = mu
    tot = 0
    for s in itertools.permutations(range(3)):
        a, b = qs[s[0]], qs[s[1]]
        v = a[0] ** i * b[0] ** (2 - i) * a[1] ** j * b[1] ** (2 - j)
        tot += v if sym else vc.perm_sign(s) * v
    return tot


def d_eta(TT4, P, mu, sym):
    out = None
    for i, j in itertools.combinations(range(4), 2):
        k, l = [r for r in range(4) if r not in (i, j)]
        q0 = (P[i][0] + P[j][0], P[i][1] + P[j][1])
        cf = (-1) ** (i + j) * S_poly(mu, q0, P[k], sym)
        if cf:
            V = cf * np.transpose(TT4, list(np.argsort([i, j, k, l])))
            out = V if out is None else out + V
    return np.zeros_like(TT4) if out is None else out


def mono36(P):
    # weight-(2,2) monomials in (p0,p1,p2): first components degree 2, second components degree 2
    e2 = [e for e in itertools.product(range(3), repeat=3) if sum(e) == 2]
    out = []
    for e in e2:
        for f in e2:
            v = 1
            for t in range(3):
                v *= P[t][0] ** e[t] * P[t][1] ** f[t]
            out.append(v)
    return out


for name in sys.argv[1:]:
    t0 = time.time()
    A = valg.build(name)
    S = vc.Setup(A)
    n = A.n
    W = {G: S.W_full(G) for G in ('133', '223')}
    tens = [('f', S.TT, True)]
    if name == 'sl3':
        M = A.mats
        d3 = np.zeros((n, n, n), dtype=np.int64)
        for a in range(n):
            for b in range(n):
                for c_ in range(n):
                    d3[a, b, c_] = int(np.trace(M[a] @ (M[b] @ M[c_] + M[c_] @ M[b])))
        TTd = vc.ein('abe,ecd->abcd', A.c, d3)
        tens.append(('d', TTd, False))
    while True:
        Ps = [vc.config4(rng) for _ in range(40)]
        if valg.exact_rank([mono36(P) for P in Ps]) == 36:
            break
    cols = []
    names = []
    for tname, TT4, sym in tens:
        for mu in itertools.product(range(3), repeat=2):
            cols.append(np.concatenate([d_eta(TT4, P, mu, sym).reshape(-1) for P in Ps]))
            names.append('d eta[%s,%s]' % (tname, mu))
    X = {}
    X['F133'] = np.concatenate([vc.alt_tensor(W['133'], P).reshape(-1) for P in Ps])
    X['F223'] = np.concatenate([vc.alt_tensor(W['223'], P).reshape(-1) for P in Ps])
    X['F133+F223'] = X['F133'] + X['F223']
    X['-3F133+2F223'] = -3 * X['F133'] + 2 * X['F223']
    X['-3F133-2F223'] = -3 * X['F133'] - 2 * X['F223']
    Dm = np.stack(cols, axis=1)
    allnz = np.flatnonzero(np.any(Dm != 0, axis=1) | (X['F133'] != 0) | (X['F223'] != 0))
    sub = np.sort(rng.choice(allnz, size=min(2500, allnz.size), replace=False))
    rD = valg.exact_rank(Dm[sub].tolist())
    line = ["%s: %d eta columns, rank d(eta) on subset = %d (config monomial rank 36)" % (name, len(cols), rD)]
    for k, v in X.items():
        r2 = valg.exact_rank(np.concatenate([Dm[sub], v[sub, None]], axis=1).tolist())
        if r2 > rD:
            line.append("   %-13s NOT exact in the full weight-(2,2) complex: rank [d eta | X] = %d > %d" % (k, r2, rD))
        else:
            # solve on subset, verify on all rows
            Ms = sp.Matrix(Dm[sub].tolist())
            sol, params = Ms.gauss_jordan_solve(sp.Matrix(v[sub].tolist()))
            sol = sol.subs({p_: 0 for p_ in params})
            den = sp.ilcm(*[sp.Rational(t).q for t in sol]) if len(sol) else 1
            coef = np.array([int(sp.Rational(t) * den) for t in sol], dtype=object)
            full_ok = bool(np.array_equal(Dm.astype(object).dot(coef), v.astype(object) * int(den)))
            line.append("   %-13s exact: X = sum a_mu d eta_mu on subset, unique: %s, verified on all %d rows: %s"
                        % (k, len(params) == 0, v.size, full_ok))
    line.append("   (%.1fs)" % (time.time() - t0))
    print("\n".join(line), flush=True)
