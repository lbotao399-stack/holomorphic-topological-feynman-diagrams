"""Random-point evaluation (exact modulo primes p < 2^20) of W_133, W_223, F_133, F_223, beta, G_kk (Killing form).

W_133(A1..A4) = - f_xyz f_{A1 y'u} f_{A2 u'w} f_{A3 z's} f_{A4 s'r} f_{x'w'r'}
W_223(A1..A4) = - f_xyz f_{A1 x'u} f_{A2 y'w} f_{A3 z's} f_{A4 s'r} f_{u'w'r'}     (primed pairs: kappa^{-1})
Contracting the paths:  W = - Phi(P1, P2, P3),
   Phi(P1,P2,P3) = sum_{x,y,z} f(e_x,e_y,e_z) f(P1 e^x, P2 e^y, P3 e^z),   e^x = kappa^{xx'} e_x'
   W_133 = - Phi(1, ad A2 ad A1, ad A4 ad A3),   W_223 = - Phi(ad A1, ad A2, ad A4 ad A3).
Fast form:  Phi = sum_{a,x} B1[a,x] <Ad_a B2, (K P3) Ad_x>_F,  B_i = P_i K^{-1}  (checked against dense einsum).
F_G^mu(X) = sum_{s in S4} sgn(s) c_mu([p_s0 p_s1][p_s2 p_s3]) W_G(X_s0..X_s3),  mu in QB = (uu,vv,ww,uv,uw,vw),
u = [p0p1], v = [p0p2], w = [p1p2], p3 = -(p0+p1+p2).
beta^mu = sum_{i<j} (-1)^{i+j} c_mu([p_k p_l]^2) K([X_i,X_j],[X_k,X_l]);  G_kk^mu = alt([p0p1][p2p3] K_01 K_23)^mu.
Float64 arithmetic on residues: every product < 2^40, every partial sum < 2^52 (chunks of 2^12 terms), then fmod.
"""
import itertools
import numpy as np
import scipy.sparse as sps
from vlie import perm_sign

QB = ['uu', 'vv', 'ww', 'uv', 'uw', 'vw']
BR = {(0, 1): (1, 0, 0), (0, 2): (0, 1, 0), (1, 2): (0, 0, 1), (0, 3): (-1, -1, 0), (1, 3): (1, 0, -1), (2, 3): (0, 1, 1)}


def br(i, j):
    if i < j:
        return BR[(i, j)]
    return tuple(-x for x in BR[(j, i)])


def quad(l1, l2):
    return (l1[0] * l2[0], l1[1] * l2[1], l1[2] * l2[2], l1[0] * l2[1] + l1[1] * l2[0],
            l1[0] * l2[2] + l1[2] * l2[0], l1[1] * l2[2] + l1[2] * l2[1])


PERMS = list(itertools.permutations(range(4)))
CQ = {s: (perm_sign(s), quad(br(s[0], s[1]), br(s[2], s[3]))) for s in PERMS}


class Engine:
    def __init__(self, L, p):
        assert p < 2 ** 20
        self.L, self.p, d = L, p, L.d
        self.d = d
        V, K, Kinv = L.mod(p)
        self.Vp = V.astype(np.float64)
        self.K = K.astype(np.float64); self.Kinv = Kinv.astype(np.float64)
        # AdS: rows a*d + c, col b  (Ad_a[c,b] = C[a,b,c]);  AdH: rows c, cols x*d + b;  AdR: rows a, cols c*d+b
        self.AdS = sps.csr_matrix((self.Vp, (L.A * d + L.Cc, L.B)), shape=(d * d, d))
        self.AdH = sps.csr_matrix((self.Vp, (L.Cc, L.A * d + L.B)), shape=(d, d * d))
        self.AdR = sps.csr_matrix((self.Vp, (L.A, L.Cc * d + L.B)), shape=(d, d * d))

    def m(self, X):
        return np.fmod(X, self.p)

    def mm(self, A, B):
        """dense product mod p, chunked over the summation index."""
        p = self.p
        n = A.shape[1]
        if n <= 4096:
            with np.errstate(all='ignore'):
                R = A @ B
            assert np.isfinite(R).all() and R.max(initial=0) < 2 ** 53
            return np.fmod(R, p)
        out = np.zeros((A.shape[0], B.shape[1]))
        for s in range(0, n, 4096):
            with np.errstate(all='ignore'):
                R = A[:, s:s + 4096] @ B[s:s + 4096]
            assert np.isfinite(R).all() and R.max(initial=0) < 2 ** 53
            out = np.fmod(out + np.fmod(R, p), p)
        return out

    def ad(self, X):
        d = self.d
        return self.m(np.asarray(self.AdR.T @ X).reshape(d, d))

    def bracket(self, X, Y):
        R = self.ad(X).astype(np.int64) @ np.asarray(Y).astype(np.int64) % self.p
        return R.astype(np.float64)

    def kap(self, X, Y):
        p = self.p
        XK = np.asarray(X).astype(np.int64) @ self.K.astype(np.int64) % p
        return int(XK @ np.asarray(Y).astype(np.int64) % p)

    def Phi(self, P1, P2, P3):
        d, p = self.d, self.p
        B1 = self.mm(P1, self.Kinv) if P1 is not None else self.Kinv
        B2 = self.mm(P2, self.Kinv)
        D = self.mm(self.K, P3)
        U = self.m(np.asarray(self.AdS @ B2))                    # (a,c) x y
        Uf = U.reshape(d, d * d)                                  # a x (c,y)
        Vh = self.m(np.asarray(self.AdH.T @ D.T).T)                # c x (x,b)
        Vf = Vh.reshape(d, d, d).transpose(1, 0, 2).reshape(d, d * d)    # x x (c,b)
        T = self.mm(Uf, Vf.T)                                     # a x x
        return int((B1.astype(np.int64) * T.astype(np.int64)).sum() % p)

    def W(self, X4):
        """W_133, W_223 at (A1..A4) = X4."""
        p = self.p
        ads = [self.ad(X) for X in X4]
        Q1 = self.mm(ads[1], ads[0])      # ad A2 ad A1
        Q3 = self.mm(ads[3], ads[2])      # ad A4 ad A3
        w133 = (-self.Phi(None, Q1, Q3)) % p
        w223 = (-self.Phi(ads[0], ads[1], Q3)) % p
        return w133, w223

    def F(self, X):
        """F_133^mu, F_223^mu (mu in QB) at the 4 vectors X."""
        p = self.p
        F1 = [0] * 6; F2 = [0] * 6
        for s in PERMS:
            sg, c = CQ[s]
            if not any(c):
                continue
            w1, w2 = self.W([X[s[0]], X[s[1]], X[s[2]], X[s[3]]])
            for mu in range(6):
                if c[mu]:
                    F1[mu] = (F1[mu] + sg * c[mu] * w1) % p
                    F2[mu] = (F2[mu] + sg * c[mu] * w2) % p
        return F1, F2

    def beta(self, X):
        p = self.p
        out = [0] * 6
        for i, j in itertools.combinations(range(4), 2):
            k, l = [r for r in range(4) if r not in (i, j)]
            c = quad(br(k, l), br(k, l))
            val = self.kap(self.bracket(X[i], X[j]), self.bracket(X[k], X[l]))
            for mu in range(6):
                out[mu] = (out[mu] + (-1) ** (i + j) * c[mu] * val) % p
        return out

    def Gkk(self, X):
        p = self.p
        out = [0] * 6
        for s in PERMS:
            sg, c = CQ[s]
            val = self.kap(X[s[0]], X[s[1]]) * self.kap(X[s[2]], X[s[3]]) % p
            for mu in range(6):
                out[mu] = (out[mu] + sg * c[mu] * val) % p
        return out


def Phi_dense_exact(L, P1, P2, P3):
    """reference: exact Fractions, small d."""
    from fractions import Fraction as Fr
    import sympy as sp
    d = L.d
    F = L.f_dense()
    Kinv = sp.Matrix(L.K).inv()
    Ki = np.array(Kinv.tolist(), dtype=object)
    B1 = P1.dot(Ki); B2 = P2.dot(Ki); B3 = P3.dot(Ki)
    T = np.einsum('abc,ax->xbc', F, B1)
    T = np.einsum('xbc,by->xyc', T, B2)
    T = np.einsum('xyc,cz->xyz', T, B3)
    return sum(F[x, y, z] * T[x, y, z] for x in range(d) for y in range(d) for z in range(d))
