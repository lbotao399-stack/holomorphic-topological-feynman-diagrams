"""Part (b): dense tensors modulo p < 2^19 (float64 BLAS, all partial sums < 2^52), Killing form.

W223[A1,A2,A3,A4] = - sum_{m,l} Y[A1,A2,m,l] Z[A3,A4,m,l]
   Y[A1,A2,m,l] = sum_{x,y,i,j} C[i,j,m] R[A1,i,x] R[A2,j,y] C[x,y,l],  R[A] = Ad_A K^{-1},  Z[A3,A4] = K Ad_A4 Ad_A3
W133[A1,A2,A3,A4] = - sum_{m,l} Y'[A1,A2,m,l] Z[A3,A4,m,l]
   Y'[A1,A2,m,l] = sum_{j,y} G[m,j,y,l] S[A1,A2,j,y],  G[m,j,y,l] = sum_{i,x} C[i,j,m] Kinv[i,x] C[x,y,l],  S[A1,A2] = Ad_A2 Ad_A1 K^{-1}
(both from W = -Phi(P1,P2,P3), Phi = sum_{x,y} kappa([P1 e^x, P2 e^y], P3 [e_x, e_y]); tested against vpoint.Engine).
Probes:  lambda_kk(X) = sum X^{uv}_{abce} k^{ab} k^{ce},  lambda_ff(X) = sum X^{uv}_{abce} kappa([e^a,e^b],[e^c,e^e]).
Pairings of a 4-tensor W:  w^kk_pi = <W, kk_pi>, w^ff_pi = <W, ff_pi>, pi in ((0123),(0213),(0312)),
   kk_(ijkl)(x) = kappa(x_i,x_j) kappa(x_k,x_l),  ff_(ijkl)(x) = kappa([x_i,x_j],[x_k,x_l]).
Om_mat[(c,g),(b,e)] = sum_x Ad_x[c,b] (Ad K^{-1})_x[g,e]  (Om = sum_x ad e_x (x) ad e^x),  m_k = tr Om^k,  n_k = tr P Om^k.
Mode 'random': C from a random totally antisymmetric f and a random symmetric K (no Jacobi): tests the pure
tensor-network identities (graph isomorphisms) W223.kk1 = W223.kk2 = -m3, W223.ff1 = W223.ff2 = n4,
W133.kk1 = -m3, W133.ff1 = -m4, W133.ff2 = -n4, and the transfer formulas for lambda_kk, lambda_ff.
"""
import sys, time, itertools
import numpy as np
from fractions import Fraction as Fr
import vlie, vpoint, vbk

P19 = [524287, 524269]


class Dense:
    def __init__(self, d, Cp, Kp, Kinvp, p):
        self.d, self.p = d, p
        self.C = Cp.astype(np.float64)        # C[a,b,c]
        self.K = Kp.astype(np.float64); self.Ki = Kinvp.astype(np.float64)
        self.Ad = np.transpose(self.C, (0, 2, 1)).copy()   # Ad[a][c,b] = C[a,b,c]

    def f(self, X):
        return np.fmod(X, self.p)

    def td(self, A, B, axes):
        with np.errstate(all='ignore'):
            R = np.tensordot(A, B, axes=axes)
        assert np.isfinite(R).all() and np.abs(R).max(initial=0) < 2 ** 52
        return self.f(R)

    def mm(self, A, B):
        p = self.p
        n = A.shape[1]
        out = np.zeros((A.shape[0], B.shape[1]))
        for s in range(0, n, 8192):
            with np.errstate(all='ignore'):
                R = A[:, s:s + 8192] @ B[s:s + 8192]
            assert np.isfinite(R).all() and R.max(initial=0) < 2 ** 52
            out = np.fmod(out + np.fmod(R, p), p)
        return out

    def build(self):
        d, p = self.d, self.p
        C, Ad, K, Ki = self.C, self.Ad, self.K, self.Ki
        R = self.td(Ad, Ki, ([2], [0]))                               # R[A,i,x]
        E1 = self.td(C, R, ([0], [1]))                                # (j,m,A1,x)
        E1m = np.transpose(E1, (2, 1, 0, 3)).reshape(d * d, d * d)    # (A1 m),(j x)
        E2 = self.td(R, C, ([2], [1]))                                # (A2,j,x,l)
        E2m = np.transpose(E2, (1, 2, 0, 3)).reshape(d * d, d * d)    # (j x),(A2 l)
        Y = self.mm(E1m, E2m).reshape(d, d, d, d)                     # A1,m,A2,l
        Ym = np.transpose(Y, (0, 2, 1, 3)).reshape(d * d, d * d)      # (A1 A2),(m l)
        KA = self.td(Ad, K, ([1], [1]))                               # Ad[A4][n,k] K[m,n] -> (A4,k,m)
        KA = np.transpose(KA, (0, 2, 1))                              # (A4,m,k) = (K Ad_A4)[m,k]
        Z = self.td(KA, Ad, ([2], [1]))                               # (A4,m,A3,l)
        Zm = np.transpose(Z, (2, 0, 1, 3)).reshape(d * d, d * d)      # (A3 A4),(m l)
        self.W223 = self.f(p - self.mm(Ym, Zm.T)).reshape(d, d, d, d)
        # W133
        G1 = self.td(C, Ki, ([0], [0]))                               # (j,m,x)
        G = self.td(G1, C, ([2], [0]))                                # (j,m,y,l)
        Gm = np.transpose(G, (0, 2, 1, 3)).reshape(d * d, d * d)      # (j y),(m l)
        AK = self.td(Ad, Ki, ([2], [0]))                              # (A1, k, y) = (Ad_A1 Ki)[k,y]
        S = self.td(Ad, AK, ([2], [1]))                               # (A2, j, A1, y)
        Sm = np.transpose(S, (2, 0, 1, 3)).reshape(d * d, d * d)      # (A1 A2),(j y)
        Yp = self.mm(Sm, Gm)                                          # (A1 A2),(m l)
        self.W133 = self.f(p - self.mm(Yp, Zm.T)).reshape(d, d, d, d)
        # ff lower and raised, kk raised
        Cm = self.C.reshape(d * d, d)
        self.FF = self.f(self.mm(self.mm(Cm, K), Cm.T)).reshape(d, d, d, d)          # kappa([e_a,e_b],[e_c,e_e])
        # C^{ab}_c = sum Ki[a,a'] Ki[b,b'] C[a',b',c]
        T1 = self.td(Ki, self.C, ([1], [0]))                          # (a, b', c)
        Cup = self.td(T1, Ki, ([1], [1])).transpose(0, 2, 1)          # (a, b, c)
        Cupm = np.ascontiguousarray(Cup).reshape(d * d, d)
        self.FFu = self.f(self.mm(self.mm(Cupm, K), Cupm.T)).reshape(d, d, d, d)    # kappa([e^a,e^b],[e^c,e^e])
        # Om
        Adu = self.td(Ki, Ad, ([1], [0]))                             # (x, g, e)
        Om = self.td(Ad, Adu, ([0], [0]))                             # (c,b,g,e)
        self.Om = np.transpose(Om, (0, 2, 1, 3)).reshape(d * d, d * d)   # (c g),(b e)

    def pair(self, X, Y):
        """full contraction sum X*Y mod p (no further raising)."""
        p = self.p
        a = X.reshape(-1).astype(np.int64); b = Y.reshape(-1).astype(np.int64)
        tot = 0
        for s in range(0, len(a), 1 << 22):
            tot = (tot + int((a[s:s + (1 << 22)] * b[s:s + (1 << 22)] % p).sum())) % p
        return tot

    def KK(self, perm):
        """raised kk_pi as 4-tensor: k^{x_i x_j} k^{x_k x_l}."""
        Ki = self.Ki
        T = np.einsum('ab,ce->abce', Ki, Ki)
        return self.f(np.transpose(T, np.argsort(perm)))

    def FFp(self, perm):
        return np.transpose(self.FFu, np.argsort(perm))

    def traces(self):
        p, d = self.p, self.d
        Om = self.Om
        O2 = self.mm(Om, Om)
        O3 = self.mm(O2, Om)
        Pidx = np.arange(d * d).reshape(d, d).T.reshape(-1)          # (c,g) -> (g,c)
        tr = lambda M: int(np.fmod(np.diag(M).sum(), p))
        PO2 = O2[Pidx, :]
        m2 = tr(O2); m3 = tr(O3); n2 = tr(PO2); n3 = tr(O3[Pidx, :])
        m4 = self.pair(O2, O2.T); n4 = self.pair(PO2, O2.T)
        return dict(m2=m2, m3=m3, m4=m4, n2=n2, n3=n3, n4=n4)


def permT(T, s):
    """(T o s)[i0..i3] = T[i_{s0}, .., i_{s3}]."""
    return np.transpose(T, np.argsort(s))


def alt_uv(Dn, W):
    p = Dn.p
    out = np.zeros_like(W)
    for s in vpoint.PERMS:
        sg, c = vpoint.CQ[s]
        if c[3]:
            out = np.fmod(out + ((sg * c[3]) % p) * permT(W, s), p)
    return out


def beta_uv(Dn):
    p = Dn.p
    out = np.zeros_like(Dn.FF)
    for i, j in itertools.combinations(range(4), 2):
        k, l = [r for r in range(4) if r not in (i, j)]
        c = vpoint.quad(vpoint.br(k, l), vpoint.br(k, l))
        if c[3]:
            out = np.fmod(out + ((-1) ** (i + j) * c[3] % p) * permT(Dn.FF, (i, j, k, l)), p)
    return out


def Gkk_uv(Dn):
    p = Dn.p
    KKl = np.fmod(np.einsum('ab,ce->abce', Dn.K, Dn.K), Dn.p)
    out = np.zeros_like(KKl)
    for s in vpoint.PERMS:
        sg, c = vpoint.CQ[s]
        if c[3]:
            out = np.fmod(out + ((sg * c[3]) % p) * permT(KKl, s), p)
    return out


PIS = [(0, 1, 2, 3), (0, 2, 1, 3), (0, 3, 1, 2)]


def r(x, p):
    x = Fr(x); return x.numerator % p * pow(x.denominator, p - 2, p) % p


def run_algebra(key, typ, n, p):
    t0 = time.time()
    L = vlie.build(key)
    d = L.d
    V, K, Kinv = L.mod(p)
    Cp = np.zeros((d, d, d), dtype=np.int64)
    np.add.at(Cp, (L.A, L.B, L.Cc), V)
    Cp %= p
    Dn = Dense(d, Cp, K, Kinv, p)
    Dn.build()
    # cross-check against the point engine at basis vectors
    E = vpoint.Engine(L, p)
    rng = np.random.default_rng(3)
    ok_eng = True
    for _ in range(3):
        idx = rng.integers(0, d, 4)
        X = [np.eye(d)[i] for i in idx]
        w1, w2 = E.W(X)
        ok_eng &= (w1 == int(Dn.W133[tuple(idx)])) and (w2 == int(Dn.W223[tuple(idx)]))
    tr = Dn.traces()
    bk = vbk.analyse(typ, n)
    pred = {'m2': d, 'm3': Fr(-d, 4), 'm4': bk['m'][4], 'n2': Fr(d, 2), 'n3': 0, 'n4': bk['n'][4]}
    ok_tr = all(tr[k] == r(v, p) for k, v in pred.items())
    w = {}
    for nm, W in (('W223', Dn.W223), ('W133', Dn.W133)):
        w[nm] = ([Dn.pair(W, Dn.KK(pi)) for pi in PIS], [Dn.pair(W, Dn.FFp(pi)) for pi in PIS])
    n4, m4, m3 = pred['n4'], pred['m4'], pred['m3']
    exp = {'W223': ([Fr(-d, 2), Fr(d, 4), Fr(d, 4)], [None, n4, n4]), 'W133': ([Fr(d), Fr(d, 4), Fr(0)], [Fr(-d, 8), -m4, -n4])}
    ok_w = all((e is None) or (r(e, p) == got) for nm in exp for kind in (0, 1) for e, got in zip(exp[nm][kind], w[nm][kind]))
    F223 = alt_uv(Dn, Dn.W223); F133 = alt_uv(Dn, Dn.W133); B = beta_uv(Dn); G = Gkk_uv(Dn)
    KK0 = Dn.KK(PIS[0]); FF0 = Dn.FFp(PIS[0])
    lam = {nm: (Dn.pair(X, KK0), Dn.pair(X, FF0)) for nm, X in (('F223', F223), ('F133', F133), ('beta', B), ('Gkk', G))}
    # transfer formulas
    tr_ok = all(lam[nm][0] == (8 * w[W][0][0] - 4 * w[W][0][1] - 4 * w[W][0][2]) % p and
                lam[nm][1] == (-4 * w[W][1][1] - 4 * w[W][1][2]) % p for nm, W in (('F223', 'W223'), ('F133', 'W133')))
    lp = {'F223': (-6 * d, -8 * n4), 'F133': (7 * d, 4 * (m4 + n4)), 'beta': (2 * d, d), 'Gkk': (8 * d * (d - 1), 8 * d)}
    ok_lam = all(lam[k][0] == r(v[0], p) and lam[k][1] == r(v[1], p) for k, v in lp.items())
    mu = lambda X: (lam[X][0] * lam['beta'][1] - lam[X][1] * lam['beta'][0]) % p
    ok_mu = mu('F223') == r(2 * d * (8 * n4 - 3 * d), p) and mu('Gkk') == r(8 * d * d * (d - 3), p) and (mu('F133') + mu('F223')) % p == 0
    print('%-16s p=%d d=%d  engine==dense %s ; m2,m3,m4,n2,n3,n4 = (d, -d/4, m4_BK, d/2, 0, n4_BK): %s ; '
          'W223 kk=(-d/2,d/4,d/4) ff=(*,n4,n4), W133 kk=(d,d/4,0) ff=(-d/8,-m4,-n4): %s ; transfer 8w0-4w1-4w2, -4w1\'-4w2\': %s ; '
          'lambda(F223)=(-6d,-8n4), lambda(F133)=(7d,4(m4+n4)), lambda(beta)=(2d,d), lambda(Gkk)=(8d(d-1),8d): %s ; '
          'mu(F223)=2d(8n4-3d), mu(Gkk)=8d^2(d-3), mu(F133)=-mu(F223): %s  (%.1fs)'
          % (L.name, p, d, ok_eng, ok_tr, ok_w, tr_ok, ok_lam, ok_mu, time.time() - t0), flush=True)


def run_random(d, p, seed):
    """generic antisymmetric f, generic symmetric K: graph-isomorphism identities."""
    rng = np.random.default_rng(seed)
    f = np.zeros((d, d, d), dtype=np.int64)
    for a, b, c in itertools.combinations(range(d), 3):
        v = int(rng.integers(-5, 6))
        for s in itertools.permutations(range(3)):
            idx = tuple((a, b, c)[k] for k in s)
            f[idx] = v * vlie.perm_sign(s)
    while True:
        M = rng.integers(-3, 4, (d, d))
        K = (M + M.T) % p
        try:
            Ki = vlie.inv_mod(K, p); break
        except ValueError:
            continue
    C = np.einsum('abe,ec->abc', f % p, Ki) % p                       # f = kappa([.,.],.)
    Dn = Dense(d, C, K, Ki, p)
    Dn.build()
    tr = Dn.traces()
    w = {}
    for nm, W in (('W223', Dn.W223), ('W133', Dn.W133)):
        w[nm] = ([Dn.pair(W, Dn.KK(pi)) for pi in PIS], [Dn.pair(W, Dn.FFp(pi)) for pi in PIS])
    m3, m4, n4, n3 = tr['m3'], tr['m4'], tr['n4'], tr['n3']
    checks = {'W223.kk1 = -m3': w['W223'][0][1] == (-m3) % p, 'W223.kk2 = -m3': w['W223'][0][2] == (-m3) % p,
              'W223.ff1 = n4': w['W223'][1][1] == n4, 'W223.ff2 = n4': w['W223'][1][2] == n4,
              'W133.kk1 = -m3': w['W133'][0][1] == (-m3) % p, 'W133.ff1 = -m4': w['W133'][1][1] == (-m4) % p,
              'W133.ff2 = -n4': w['W133'][1][2] == (-n4) % p}
    F223 = alt_uv(Dn, Dn.W223)
    KK0 = Dn.KK(PIS[0]); FF0 = Dn.FFp(PIS[0])
    lk, lf = Dn.pair(F223, KK0), Dn.pair(F223, FF0)
    checks['transfer kk'] = lk == (8 * w['W223'][0][0] - 4 * w['W223'][0][1] - 4 * w['W223'][0][2]) % p
    checks['transfer ff'] = lf == (-4 * w['W223'][1][1] - 4 * w['W223'][1][2]) % p
    # transfer for a completely random (non-invariant, no symmetry) 4-tensor
    Wr = rng.integers(0, p, (d, d, d, d)).astype(np.float64)
    Fr_ = alt_uv(Dn, Wr)
    wk = [Dn.pair(Wr, Dn.KK(pi)) for pi in PIS]; wf = [Dn.pair(Wr, Dn.FFp(pi)) for pi in PIS]
    checks['transfer kk (random W)'] = Dn.pair(Fr_, KK0) == (8 * wk[0] - 4 * wk[1] - 4 * wk[2]) % p
    checks['transfer ff (random W)'] = Dn.pair(Fr_, FF0) == (-4 * wf[1] - 4 * wf[2]) % p
    checks['n3 = 0'] = n3 == 0
    print('random antisymmetric f, d=%d, seed %d, p=%d: %s ; W133.kk2 = %d, n3 = %d, m3 = %d' % (d, seed, p, checks, w['W133'][0][2], n3, m3), flush=True)


if __name__ == '__main__':
    args = sys.argv[1:]
    if args[0] == 'random':
        for seed in (1, 2, 3):
            for p in P19:
                run_random(int(args[1]), p, seed)
    else:
        for p in P19[:int(args[0])]:
            for item in args[1:]:
                key, rt = item.split(':')
                run_algebra(key, rt[0], int(rt[1:]), p)
