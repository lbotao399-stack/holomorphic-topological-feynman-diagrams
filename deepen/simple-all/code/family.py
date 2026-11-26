"""n4 = tr_{g(x)g}(P Om^4) for the classical families with kappa = tr_V (lie.so_tr, sp_tr, sl_tr), exact.

Degree bound.  V = defining representation.  so_N: g = Lambda^2 V, sp_N: g = S^2 V (as g-modules, X -> X J^{-1}),
sl_N: Om_sl = Om_gl on g (x) g and Om_gl = 0 on the centre components, gl_N = V (x) V^*.  ad(x) acts on V^{(x)2} (or
V (x) V^*) as x(x)1 + 1(x)x (resp. x(x)1 - 1(x)x^T), so Om = sum_{i in {1,2}, j in {3,4}} t_ij with
t = sum_a rho(e_a) (x) rho(e^a) on V (x) V (or the partial transpose); checked here for every N:
   so_N: t = (1/2)(P - E_J),  sp_N: t = (1/2)(P + E_J'),  gl_N: t = P   (E: rank-one contraction with the form J).
n4 = tr_{V^{(x)4}}(Pi S Om^4), Pi = projector onto g (x) g, S = exchange of the two factors: a linear combination,
with N-independent coefficients, of traces of products of permutations and at most 4 contraction elements E.
Each closed index loop of such a product either contains a cap of some E (at most 4 caps) or is a cycle of the
4 strands closed by the trace (at most 4): at most 8 loops, so n4(N) is a polynomial in N of degree <= 8
(each loop gives a factor +-N).  Interpolation through 9 values of N, checked at all further computed N.
Killing ratio r(N) = Killing / tr_V: computed exactly for each N.
"""
import sys, time, json
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/simple-all/code')
import numpy as np
import sympy as sp
from fractions import Fraction as Fr
import lie, t4


def defining_mats(L):
    """recover the N x N matrices of the basis from the construction (lie._mat_alg keeps no matrices): rebuild."""
    raise NotImplementedError


def casimir_tensor(kind, N):
    """t = sum_a rho(e_a) (x) rho(e^a) on V (x) V, with the basis/form of lie.<kind>_tr, compared with the Brauer form."""
    import itertools
    if kind == 'so':
        pr = lambda k: N - 1 - k
        mats = []
        for i in range(N):
            for j in range(N):
                if i + j < N - 1:
                    m = np.zeros((N, N), dtype=np.int64); m[i, j] += 1; m[pr(j), pr(i)] -= 1
                    mats.append(m)
        J = np.zeros((N, N), dtype=np.int64)
        for k in range(N):
            J[k, pr(k)] = 1
    elif kind == 'sp':
        n = N // 2
        pr = lambda k: N - 1 - k
        s = lambda k: 1 if k < n else -1
        mats = []
        for i in range(N):
            for j in range(N):
                if i + j < N - 1:
                    m = np.zeros((N, N), dtype=np.int64); m[i, j] += 1; m[pr(j), pr(i)] += s(pr(j)) * s(i)
                    mats.append(m)
                elif i + j == N - 1:
                    m = np.zeros((N, N), dtype=np.int64); m[i, j] = 1
                    mats.append(m)
        J = np.zeros((N, N), dtype=np.int64)
        for k in range(N):
            J[k, pr(k)] = s(k)
    else:
        mats = []
        for i in range(N):
            for j in range(N):
                if i != j:
                    m = np.zeros((N, N), dtype=np.int64); m[i, j] = 1; mats.append(m)
        for i in range(N - 1):
            m = np.zeros((N, N), dtype=np.int64); m[i, i] = 1; m[i + 1, i + 1] = -1; mats.append(m)
        J = None
    d = len(mats)
    K = sp.Matrix([[int(np.trace(a @ b)) for b in mats] for a in mats])
    Ki = K.inv()
    # t[(i,k),(j,l)] = sum_{a,b} (e_a)_{ij} Ki[a,b] (e_b)_{kl}  : operator on V(x)V, (v (x) w) -> sum e_a v (x) e^a w
    A = np.array([m.reshape(-1) for m in mats], dtype=object)            # d x N^2
    Kio = np.array(Ki.tolist(), dtype=object)
    T = A.T.dot(Kio).dot(A)                                              # (ij),(kl)
    T = T.reshape(N, N, N, N)                                            # [i,j,k,l] = sum (e_a)_ij (e^a)_kl
    # operator on V(x)V: out (i,k) <- in (j,l)
    P = np.zeros((N, N, N, N), dtype=object)                             # P[i,j,k,l]: (i,k) <- (j,l) : delta_il delta_kj
    E = np.zeros((N, N, N, N), dtype=object)
    for i in range(N):
        for k in range(N):
            P[i, k, k, i] = 1
    if J is not None:
        Jinv = np.array(sp.Matrix(J.tolist()).inv().tolist(), dtype=object)
        # E: (j,l) -> J_{jl}^{-1}-contraction ... E[i,j,k,l] = c_ik * b_jl with c = J^{-1}-type, b = J
        for i in range(N):
            for j in range(N):
                for k in range(N):
                    for l in range(N):
                        E[i, j, k, l] = Jinv[i, k] * J[j, l]
    # T as operator: out (i,k), in (j,l): T[i,j,k,l]; P as operator out (i,k), in (j,l) = delta_{i l} delta_{k j}
    Pop = np.zeros((N, N, N, N), dtype=object)
    for i in range(N):
        for k in range(N):
            Pop[i, k, k, i] = 1     # index order [i,j,k,l] with j = k_in, l = i_in
    if kind == 'so':
        ok = all(T[i, j, k, l] == Fr(1, 2) * (Pop[i, j, k, l] - E[i, j, k, l]) for i in range(N) for j in range(N) for k in range(N) for l in range(N))
        form = "(1/2)(P - E_J)"
    elif kind == 'sp':
        ok = all(T[i, j, k, l] == Fr(1, 2) * (Pop[i, j, k, l] + E[i, j, k, l]) for i in range(N) for j in range(N) for k in range(N) for l in range(N))
        if not ok:
            ok2 = all(T[i, j, k, l] == Fr(1, 2) * (Pop[i, j, k, l] - E[i, j, k, l]) for i in range(N) for j in range(N) for k in range(N) for l in range(N))
            form = "(1/2)(P - E_J)" if ok2 else "?"
            ok = ok2
        else:
            form = "(1/2)(P + E_J)"
    else:
        ok = all(T[i, j, k, l] == Pop[i, j, k, l] - Fr(1, N) * (1 if (i == j and k == l) else 0)
                 for i in range(N) for j in range(N) for k in range(N) for l in range(N))
        form = "P - (1/N) 1(x)1 (sl), i.e. P on gl"
    return ok, form


def n4_of(L):
    E = t4.Engine(L)
    Om2 = E.compose(E.Om, E.Om)
    Om4 = E.compose(Om2, Om2)
    kk_up = E.raise_all(E.kk)
    n4 = E.pair(Om4, t4.perm(kk_up, (0, 3, 1, 2)))
    m4 = E.pair(Om4, t4.perm(kk_up, (0, 2, 1, 3)))
    return n4, m4, E.Om.nnz(), Om4.nnz()


if __name__ == "__main__":
    kind = sys.argv[1]
    Ns = [int(x) for x in sys.argv[2:]]
    out = []
    for N in Ns:
        t0 = time.time()
        L = {'so': lie.so_tr, 'sp': lambda N: lie.sp_tr(N // 2), 'sl': lie.sl_tr}[kind](N)
        r, okr = L.killing_ratio()
        n4, m4, nO, nO4 = n4_of(L)
        ok, form = casimir_tensor(kind, N) if N <= 12 else (None, '-')
        print("%s N %2d dim %3d : Killing = %s tr_V (%s) ; n4(tr) = %s ; m4(tr) = %s ; t on V(x)V = %s: %s ; nnz Om %d Om4 %d (%.1fs)"
              % (kind, N, L.dim, r, okr, n4, m4, form, ok, nO, nO4, time.time() - t0), flush=True)
        out.append({'N': N, 'dim': L.dim, 'r': str(r), 'n4': str(n4), 'm4': str(m4)})
    json.dump(out, open((__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/simple-all/logs/family_%s.json' % kind, 'w'))
