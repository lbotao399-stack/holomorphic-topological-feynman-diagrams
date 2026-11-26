"""Direct numerical integration of K = int_{RP^{N-1}_{>0}} D iota_E vol for K4 subdivisions (no reduction formulas).

Coordinates: chain resistances R_E (6 K4 edges) and chain fractions x_{E,j}, t_{E,j} = R_E x_{E,j}
  d^N t = prod_E R_E^{a_E-1} dR_E d^{a_E-1}x_E  ->  K = int_{prod simplices} dx int_{RP^5} D(t(R,x)) prod_E R_E^{a_E-1}.
x-integral: Gauss product rules on the simplices (exact for polynomials; the x-degree is checked by m -> m+1).
R-integral: all 720 Hepp sectors R_{s0}=1 >= R_{s1} >= ... >= R_{s5}, R_{s,k} = u_0...u_{k-1},
  d^5R = u_0^4 u_1^3 u_2^2 u_3 du, tensor Gauss-Legendre with n points per u.
D(t): vk.NG (fresh Kirchhoff solution on the full graph, D = det[1;Psi]/(U^H sum t))."""
import itertools, time
import numpy as np
from vk import G, NG


def chains_of(nv, edges):
    deg = [0] * nv
    for a, b in edges:
        deg[a] += 1; deg[b] += 1
    adj = {v: [] for v in range(nv)}
    for e, (a, b) in enumerate(edges):
        adj[a].append((b, e)); adj[b].append((a, e))
    cubic = [v for v in range(nv) if deg[v] == 3]
    used = set(); out = []
    for p in cubic:
        for (w, e) in adj[p]:
            if e in used:
                continue
            es = [e]; used.add(e); cur = w; path = [p]
            while deg[cur] == 2:
                path.append(cur)
                (n1, e1), (n2, e2) = adj[cur]
                nxt, en = (n2, e2) if e1 == es[-1] else (n1, e1)
                used.add(en); es.append(en); cur = nxt
            path.append(cur)
            out.append((path, es))
    assert len(out) == 6 and sorted(sum([es for _, es in out], [])) == list(range(len(edges)))
    return cubic, out


def gl01(m):
    x, w = np.polynomial.legendre.leggauss(m)
    return (x + 1) / 2, w / 2


def simplex_rule(a, m):
    """nodes (Q, a) on the (a-1)-simplex sum x = 1, weights for Lebesgue measure in (x_2..x_a)."""
    if a == 1:
        return np.ones((1, 1)), np.ones(1)
    s, w = gl01(m)
    if a == 2:
        return np.stack([1 - s, s], 1), w
    if a == 3:
        U, V = np.meshgrid(s, s, indexing='ij'); WU, WV = np.meshgrid(w, w, indexing='ij')
        U, V = U.ravel(), V.ravel()
        return np.stack([U, (1 - U) * V, (1 - U) * (1 - V)], 1), (WU * WV).ravel() * (1 - U)
    raise NotImplementedError


def x_rule(chains, N, m):
    per = [simplex_rule(len(es), m) for _, es in chains]
    nodes, wts = [], []
    for combo in itertools.product(*[range(len(p[1])) for p in per]):
        xv = np.zeros(N); wt = 1.0
        for (path, es), p, k in zip(chains, per, combo):
            for j, e in enumerate(es):
                xv[e] = p[0][k, j]
            wt *= p[1][k]
        nodes.append(xv); wts.append(wt)
    return np.array(nodes), np.array(wts)


class K4Int:
    def __init__(self, nv, edges, lams, m=3):
        self.g = G(nv, edges)
        self.ng = NG(self.g, lams)
        self.cubic, self.chains = chains_of(nv, edges)
        self.N = len(edges)
        self.a = np.array([len(es) for _, es in self.chains])
        self.E_of_e = np.zeros(self.N, dtype=int)
        for k, (_, es) in enumerate(self.chains):
            for e in es:
                self.E_of_e[e] = k
        self.xn, self.xw = x_rule(self.chains, self.N, m)

    def f_R(self, R):
        """R: (P, 6) -> int dx D(t(R,x)) prod R^{a-1}: (P, S)"""
        P = R.shape[0]; Q = len(self.xw)
        t = R[:, None, self.E_of_e] * self.xn[None, :, :]          # (P, Q, N)
        Dv = self.ng.D(t.reshape(P * Q, self.N)).reshape(P, Q, -1)
        return np.einsum('pqs,q->ps', Dv, self.xw) * np.prod(R ** (self.a - 1), axis=1)[:, None]

    def sector_integral(self, n, chunk=20000, sectors=None):
        u, w = gl01(n)
        U = np.array(list(itertools.product(u, repeat=5)))
        W = np.prod(np.array(list(itertools.product(w, repeat=5))), axis=1)
        W = W * U[:, 0] ** 4 * U[:, 1] ** 3 * U[:, 2] ** 2 * U[:, 3]
        Rk = np.ones((len(U), 6))
        for k in range(1, 6):
            Rk[:, k] = Rk[:, k - 1] * U[:, k - 1]
        perms = list(itertools.permutations(range(6))) if sectors is None else sectors
        total = 0.0; per = []
        for sg in perms:
            R = np.empty_like(Rk)
            R[:, list(sg)] = Rk
            acc = 0.0
            for i in range(0, len(U), chunk):
                acc = acc + np.einsum('ps,p->s', self.f_R(R[i:i + chunk]), W[i:i + chunk])
            per.append(acc)
            total = total + acc
        return total, per
