"""Root-data computations (independent of structure constants): Brauer-Klimyk / Racah-Speiser decompositions.

Weights in fundamental-weight coordinates; A_ij = <alpha_j, alpha_i^vee> = 2 B_ij / B_ii (B = Gram matrix of simple
roots, vlie.gram);  alpha_j = sum_i A_ij omega_i;  s_i(l) = l - l_i alpha_i;  rho = (1,..,1).
For a W-invariant weight multiset m:  sum_nu m_nu e^nu = sum_lambda c_lambda chi_lambda,
   c_lambda = sum_{nu, w : w(nu + rho) = lambda + rho} sgn(w) m_nu       (reflect nu + rho into the dominant chamber).
Brauer-Klimyk:  V_lam (x) V_mu = sum_{nu in wt(V_mu)} m_mu(nu) [reflect(lam + nu + rho) - rho].
Adams:  chi_{S^2 V} = (chi^2 + psi^2 chi)/2,  chi_{Lambda^2 V} = (chi^2 - psi^2 chi)/2,  psi^2: e^nu -> e^{2 nu}.
Casimir (Killing normalization, C(adjoint) = 1):  C(lam) = (lam, lam + 2 rho) / (theta, theta + 2 rho).
Om on V subset g (x) g:  w_V = (C(V) - 2)/2;  m_k = tr(Om^k) = sum_{V in g(x)g} dim V w_V^k;
n_k = tr(P Om^k) = sum_V c_V(psi^2 g) dim V w_V^k   (P = +1 on S^2, -1 on Lambda^2).
"""
import sys
from fractions import Fraction as Fr
import numpy as np
import vlie


class RootData:
    def __init__(self, typ, n):
        self.typ, self.n = typ, n
        B = vlie.gram(typ, n)
        self.B = [[Fr(B[i, j]) for j in range(n)] for i in range(n)]
        self.A = [[2 * self.B[i][j] / self.B[i][i] for j in range(n)] for i in range(n)]
        assert all(x.denominator == 1 for row in self.A for x in row)
        self.Ai = [[int(x) for x in row] for row in self.A]
        self.pos = vlie.positive_roots(B)                       # simple-root coordinates
        self.dvec = [self.B[i][i] / 2 for i in range(n)]
        self.posw = [self.to_w(r) for r in self.pos]
        ht = [sum(r) for r in self.pos]
        self.theta = self.posw[int(np.argmax(ht))]
        self.dim = 2 * len(self.pos) + n
        # inner product on omega coordinates: (omega_i, omega_j) = M B M^T, M = (A^T)^{-1}
        import sympy as sp
        At = sp.Matrix(self.Ai).T
        M = At.inv()
        G = M * sp.Matrix(self.B) * M.T
        self.G = [[Fr(int(sp.Rational(G[i, j]).p), int(sp.Rational(G[i, j]).q)) for j in range(n)] for i in range(n)]

    def to_w(self, c):
        n = self.n
        return tuple(sum(self.Ai[k][j] * c[j] for j in range(n)) for k in range(n))

    def ip(self, l, m):
        n = self.n
        return sum(l[i] * self.G[i][j] * m[j] for i in range(n) for j in range(n))

    def adj_weights(self):
        m = {}
        for r in self.posw:
            m[r] = m.get(r, 0) + 1
            nr = tuple(-x for x in r)
            m[nr] = m.get(nr, 0) + 1
        z = tuple([0] * self.n)
        m[z] = m.get(z, 0) + self.n
        return m

    def reflect(self, v):
        """reflect v into the dominant chamber; returns (sign, v) or (0, None) if singular."""
        v = list(v); sgn = 1
        n = self.n
        while True:
            neg = [i for i in range(n) if v[i] < 0]
            if not neg:
                break
            i = neg[0]
            li = v[i]
            for k in range(n):
                v[k] -= li * self.Ai[k][i]
            sgn = -sgn
        if any(x == 0 for x in v):
            return 0, None
        return sgn, tuple(v)

    def decompose(self, mult):
        out = {}
        for nu, m in mult.items():
            s, v = self.reflect(tuple(x + 1 for x in nu))
            if s:
                lam = tuple(x - 1 for x in v)
                out[lam] = out.get(lam, 0) + s * m
        return {k: v for k, v in out.items() if v}

    def tensor_with_adj(self, lam):
        """V_lam (x) g via Brauer-Klimyk."""
        out = {}
        for nu, m in self.adj_weights().items():
            s, v = self.reflect(tuple(l + x + 1 for l, x in zip(lam, nu)))
            if s:
                mu = tuple(x - 1 for x in v)
                out[mu] = out.get(mu, 0) + s * m
        return {k: v for k, v in out.items() if v}

    def weyl_dim(self, lam):
        num = Fr(1); den = Fr(1)
        for c in self.pos:
            num *= sum((lam[j] + 1) * c[j] * self.dvec[j] for j in range(self.n))
            den *= sum(c[j] * self.dvec[j] for j in range(self.n))
        return num / den

    def casimir(self, lam):
        rho = tuple([1] * self.n)
        t = self.theta
        return self.ip(lam, tuple(l + 2 * r for l, r in zip(lam, rho))) / self.ip(t, tuple(x + 2 for x in t))


def analyse(typ, n):
    R = RootData(typ, n)
    gg = R.tensor_with_adj(R.theta)
    psi2 = R.decompose({tuple(2 * x for x in nu): m for nu, m in R.adj_weights().items()})
    zero = tuple([0] * n)
    th = R.theta
    S2 = {k: (gg.get(k, 0) + psi2.get(k, 0)) // 2 for k in set(gg) | set(psi2)}
    L2 = {k: (gg.get(k, 0) - psi2.get(k, 0)) // 2 for k in set(gg) | set(psi2)}
    assert all((gg.get(k, 0) + psi2.get(k, 0)) % 2 == 0 for k in S2)
    S2 = {k: v for k, v in S2.items() if v}; L2 = {k: v for k, v in L2.items() if v}
    dims = {k: R.weyl_dim(k) for k in set(gg) | set(psi2)}
    w = {k: (R.casimir(k) - 2) / 2 for k in dims}
    d = R.dim
    chk_dim = (sum(v * dims[k] for k, v in gg.items()) == d * d, sum(v * dims[k] for k, v in S2.items()) == d * (d + 1) // 2)
    mk = {k: sum(v * dims[l] * w[l] ** k for l, v in gg.items()) for k in range(1, 5)}
    nk = {k: sum(v * dims[l] * w[l] ** k for l, v in psi2.items()) for k in range(1, 5)}
    res = dict(R=R, d=d, theta=th, gg=gg, S2=S2, L2=L2, dims=dims, w=w,
               triv_gg=gg.get(zero, 0), adj_gg=gg.get(th, 0), adj_S2=S2.get(th, 0), adj_L2=L2.get(th, 0),
               chk_dim=chk_dim, m=mk, n=nk,
               X2_w=[(k, w[k]) for k in L2 if k != th])
    return res


NAMES = {('A', 1): 'sl2', ('B', 2): 'so5', ('C', 2): 'sp4', ('G', 2): 'G2', ('F', 4): 'F4', ('E', 6): 'E6', ('E', 7): 'E7', ('E', 8): 'E8'}


def name(typ, n):
    if (typ, n) in NAMES:
        return NAMES[(typ, n)]
    return {'A': 'sl%d' % (n + 1), 'B': 'so%d' % (2 * n + 1), 'C': 'sp%d' % (2 * n), 'D': 'so%d' % (2 * n)}[typ]


if __name__ == '__main__':
    import time
    lst = []
    for a in sys.argv[1:]:
        lst.append((a[0], int(a[1:])))
    for typ, n in lst:
        t0 = time.time()
        r = analyse(typ, n)
        d = r['d']
        print('%-6s %s%d d=%d theta=%s  trivial in g(x)g %d  g in g(x)g %d [= dim (g*^(x)3)^g]  g in S^2 g %d [>= dim (S^3 g*)^g]  '
              'g in Lambda^2 g %d [>= dim (Lambda^3 g*)^g]  dims ok %s  m1..m4 %s  n1..n4 %s  w on Lambda^2 g minus g: %s  m4-n4-d/8 = %s  n4 = %s  8n4-3d = %s  (%.1fs)'
              % (name(typ, n), typ, n, d, r['theta'], r['triv_gg'], r['adj_gg'], r['adj_S2'], r['adj_L2'], r['chk_dim'],
                 [str(r['m'][k]) for k in range(1, 5)], [str(r['n'][k]) for k in range(1, 5)],
                 sorted(set(str(x[1]) for x in r['X2_w'])), r['m'][4] - r['n'][4] - Fr(d, 8), r['n'][4], 8 * r['n'][4] - 3 * d,
                 time.time() - t0), flush=True)
