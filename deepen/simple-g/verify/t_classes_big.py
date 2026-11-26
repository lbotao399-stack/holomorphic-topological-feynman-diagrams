"""Same decisive computation as t_classes.py, processed configuration by configuration (memory-light; used for F4).
For each combination X: c fixed at the first configuration from an entry with beta != 0; then
X*den - beta*num is checked = 0 on every entry of every configuration, or a nonzero 2x2 minor is printed."""
import sys, time
from fractions import Fraction as Fr
import numpy as np
import sympy as sp
import valg, vcoch as vc

rng = np.random.default_rng(2026)
for name in sys.argv[1:]:
    t0 = time.time()
    A = valg.build(name)
    S = vc.Setup(A)
    tb = time.time() - t0
    W = {G: S.W_full(G) for G in ('133', '223')}
    tW = time.time() - t0
    n = A.n
    while True:
        Ps = [vc.config4(rng) for _ in range(7)]
        E = sp.Matrix([vc.quad_basis_eval(P) for P in Ps[:6]])
        if E.det() != 0:
            break
    sub = rng.choice(n ** 4, size=3000, replace=False)
    subvals = {k: [] for k in ('beta', 'F133', 'F223')}
    labels = ['F133', 'F223', 'F133+F223', 'w133+w223 = -3F133+2F223', 'w133-w223 = -3F133-2F223']
    state = {}
    i0 = None
    s7 = S.s ** 7
    for j, P in enumerate(Ps):
        F1 = vc.alt_tensor(W['133'], P).reshape(-1)
        F2 = vc.alt_tensor(W['223'], P).reshape(-1)
        B = vc.beta_tensor(S.TT, P).reshape(-1)
        for k, T in (('beta', B), ('F133', F1), ('F223', F2)):
            subvals[k].append([int(T[i]) for i in sub])
        if j == 6:
            break
        Xs = {'F133': F1, 'F223': F2, 'F133+F223': F1 + F2, labels[3]: -3 * F1 + 2 * F2, labels[4]: -3 * F1 - 2 * F2}
        if j == 0:
            i0 = int(np.flatnonzero(B != 0)[0])
            for l in labels:
                state[l] = ['prop', Fr(int(Xs[l][i0]), int(B[i0])), None]
            nonzero = (bool(B.any()), bool(F1.any()), bool(F2.any()))
            B0 = int(B[i0])
            X0 = {l: int(Xs[l][i0]) for l in labels}
        for l in labels:
            if state[l][0] != 'prop':
                continue
            c = state[l][1]
            X = Xs[l]
            bnd = float(np.abs(X).max()) * abs(c.denominator) + float(np.abs(B).max()) * abs(c.numerator)
            assert bnd < 2.0 ** 62
            D = X.astype(np.int64) * c.denominator - B.astype(np.int64) * c.numerator
            nz = np.flatnonzero(D)
            if nz.size:
                jj = int(nz[0])
                if j == 0:
                    minor = B0 * int(X[jj]) - X0[l] * int(B[jj])
                    state[l] = ['notprop', None, (j, jj, minor)]
                else:
                    state[l] = ['notprop', None, (j, jj, 'entry differs from c*beta with c fixed at config 0')]
        del F1, F2, B, Xs
    Einv = E.inv()
    m7 = vc.quad_basis_eval(Ps[6])
    interp_ok = True
    for k in subvals:
        V = sp.Matrix(subvals[k][:6])
        pred = sp.Matrix([m7]) * (Einv * V)
        interp_ok &= all(pred[0, i] == subvals[k][6][i] for i in range(len(sub)))
    print("%s: dim %d, kappa = Killing/%d, s = %d, det[mu(P_j)] = %s, interpolation at P_7 (3000 entries): %s"
          % (name, n, S.g, S.s, E.det(), interp_ok))
    print("   nonzero: beta %s, F133 %s, F223 %s" % nonzero)
    for l in labels:
        st = state[l]
        if st[0] == 'prop':
            print("   %-26s = c beta exactly on all %d entries: c(kappa) = %s, c(Killing) = %s"
                  % (l, 6 * n ** 4, st[1] / s7, st[1] / s7 / S.g ** 2))
        else:
            print("   %-26s NOT in C beta: config %d, entry %d, 2x2 minor with entry %d of config 0: %s"
                  % (l, st[2][0], st[2][1], i0, st[2][2]))
    print("   arithmetic: %s ; time build %.1fs, W %.1fs, total %.1fs" % (dict(vc.STATS), tb, tW - tb, time.time() - t0), flush=True)
