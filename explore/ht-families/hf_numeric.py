"""Float brute-force density (vectorized) and QMC integration over ratio coordinates.

Density: D(t) of hf_brute (U^{-H} rho_n ^ prod psi = D iota_E vol), here in floats for batches of t,
T = 1, n = 4 (Pfaffian of a 4x4 matrix of 2-forms) or n = 2.
Ratio map: built from the plan only (series t_ij = x_ij R_i, parallel R_i = t_red / w_i, t_final = 1).
"""
import itertools
import numpy as np
from scipy.special import gammaincinv, gammaln
from scipy.stats import qmc


class RatioMap:
    """own construction of t(ratios) from a plan [(p, q, chains, signs)]"""

    def __init__(self, N0, plan):
        self.N0, self.plan = N0, plan
        self.blocks = []        # list of lists of ratio names; first one eliminated
        for mi, (p, q, chains, signs) in enumerate(plan):
            k = len(chains)
            if k >= 2:
                self.blocks.append(["w%d_%d" % (mi + 1, i + 1) for i in range(k)])
            for i, ch in enumerate(chains):
                if len(ch) >= 2:
                    self.blocks.append(["x%d_%d_%d" % (mi + 1, i + 1, j + 1) for j in range(len(ch))])
        self.names = [b for blk in self.blocks for b in blk]
        self.coords = [b for blk in self.blocks for b in blk[1:]]
        # exponents texp[e][name]
        new_final = N0 + len(plan) - 1
        texp = {new_final: {}}
        for mi in reversed(range(len(plan))):
            p, q, chains, signs = plan[mi]
            Rnew = texp[N0 + mi]
            k = len(chains)
            for i, ch in enumerate(chains):
                Ri = dict(Rnew)
                if k >= 2:
                    b = "w%d_%d" % (mi + 1, i + 1)
                    Ri[b] = Ri.get(b, 0) - 1
                for j, e in enumerate(ch):
                    te = dict(Ri)
                    if len(ch) >= 2:
                        b = "x%d_%d_%d" % (mi + 1, i + 1, j + 1)
                        te[b] = te.get(b, 0) + 1
                    texp[e] = te
        self.texp = [texp[e] for e in range(N0)]
        self.E = np.array([[self.texp[e].get(b, 0) for b in self.names] for e in range(N0)], dtype=float)
        # d b / d xi_c
        self.dbdxi = np.zeros((len(self.names), len(self.coords)))
        idx = {b: i for i, b in enumerate(self.names)}
        for c, name in enumerate(self.coords):
            blk = next(bl for bl in self.blocks if name in bl)
            self.dbdxi[idx[name], c] += 1
            self.dbdxi[idx[blk[0]], c] -= 1

    def t_and_jac(self, B):
        """B: (npts, nnames) barycentric values. Returns t (npts, N0) and |det[t; dt/dxi]|."""
        logt = np.log(B) @ self.E.T
        t = np.exp(logt)
        # dt_e/dxi_c = t_e sum_b E[e,b] (db/dxi_c)/b
        G = (self.E[None, :, :] / B[:, None, :]) @ self.dbdxi       # (npts, N0, ncoords)
        dT = t[:, :, None] * G
        Mat = np.concatenate([t[:, None, :], np.transpose(dT, (0, 2, 1))], axis=1)   # (npts, N0, N0)
        return t, np.abs(np.linalg.det(Mat))

    def sample(self, U, alpha=None):
        """U: (npts, nnames) uniforms -> Dirichlet(alpha_b) in each block (default alpha = 1/2);
        returns (B, density w.r.t. Lebesgue measure on the independent coordinates)."""
        B = np.empty_like(U)
        logdens = np.zeros(U.shape[0])
        col = 0
        for blk in self.blocks:
            m1 = len(blk)
            a = np.array([0.5 if alpha is None else alpha[b] for b in blk])
            g = gammaincinv(a[None, :], np.clip(U[:, col:col + m1], 1e-300, 1 - 1e-16))
            g = np.maximum(g, 1e-300)
            b = g / g.sum(axis=1, keepdims=True)
            B[:, col:col + m1] = b
            logdens += gammaln(a.sum()) - gammaln(a).sum() + ((a - 1)[None, :] * np.log(b)).sum(axis=1)
            col += m1
        return B, np.exp(logdens)


def density_batch(nv, edges, C, t, lam, T=1):
    """float D(t) for a batch t (npts, N); lam: (nv-1, H) array."""
    N = len(edges)
    g = nv - 1
    npts = t.shape[0]
    C = np.array(C, dtype=float)
    n = C.shape[0]
    H = lam.shape[1]
    Bm = np.zeros((N, g))
    for e, (a, b) in enumerate(edges):
        if a < g:
            Bm[e, a] += 1
        if b < g:
            Bm[e, b] -= 1
    L = np.einsum('ev,pe,ew->pvw', Bm, 1 / t, Bm)
    phi = np.linalg.solve(L, np.broadcast_to(lam, (npts, g, H)))
    y = np.einsum('ev,pvi->pei', Bm, phi) / t[:, :, None]                # (npts, N, H)
    M = np.einsum('re,pe,se->prs', C, t, C)
    Mi = np.linalg.inv(M)
    U = np.linalg.det(M)
    gens = list(range(N - 1))
    ng = N - 1
    Psi = np.einsum('re,pei->pire', C[:, :ng], y[:, :ng, :]).reshape(npts, H * n, ng)  # rows (i, r)
    if T == 1:
        Gr = np.einsum('re,ae->rae', C[:, :ng], C[:, :ng])               # (n, n, ng): dM_ra coeff on e
        # A_rs[e,f] = -1/2 sum_ab Mi_ab (G_r[a,e] G_s[b,f] - G_r[a,f] G_s[b,e])
        X = np.einsum('rae,pab,sbf->prsef', Gr, Mi, Gr)
        A = -0.5 * (X - np.transpose(X, (0, 1, 2, 4, 3)))
        if n == 2:
            subsets = list(itertools.combinations(gens, 2))
            rho = np.stack([A[:, 0, 1, a, b] for (a, b) in subsets], axis=1)
        elif n == 4:
            subsets = list(itertools.combinations(gens, 4))

            def w22(al, be):
                out = []
                for (a, b, c, d) in subsets:
                    out.append(al[:, a, b] * be[:, c, d] - al[:, a, c] * be[:, b, d] + al[:, a, d] * be[:, b, c]
                               + al[:, b, c] * be[:, a, d] - al[:, b, d] * be[:, a, c] + al[:, c, d] * be[:, a, b])
                return np.stack(out, axis=1)
            rho = (w22(A[:, 0, 1], A[:, 2, 3]) - w22(A[:, 0, 2], A[:, 1, 3]) + w22(A[:, 0, 3], A[:, 1, 2]))
        else:
            raise NotImplementedError
    else:
        subsets = [()]
        rho = np.ones((npts, 1))
    top = np.zeros(npts)
    for si, S in enumerate(subsets):
        rest = [e for e in gens if e not in S]
        inv = sum(1 for x in S for yv in rest if yv < x)
        sgn = -1.0 if inv % 2 else 1.0
        minor = np.linalg.det(Psi[:, :, rest]) if rest else np.ones(npts)
        top += sgn * rho[:, si] * minor
    D = (-1) ** (N - 1) * top / (t[:, N - 1] * U ** H)
    if T == 1:
        D = D / np.sqrt(U)
    return D


def qmc_integral(nv, edges, C, plan, lam, T=1, m=16, nrep=8, seed=0, batch=4096, alpha=None):
    rm = RatioMap(len(edges), plan)
    vals = []
    for rep in range(nrep):
        sob = qmc.Sobol(d=len(rm.names), scramble=True, seed=seed + rep)
        Uu = sob.random_base2(m)
        tot = 0.0
        for s in range(0, Uu.shape[0], batch):
            Bb, dens = rm.sample(Uu[s:s + batch], alpha)
            t, J = rm.t_and_jac(Bb)
            D = density_batch(nv, edges, C, t, lam, T)
            tot += np.sum(D * J / dens)
        vals.append(tot / Uu.shape[0])
    vals = np.array(vals)
    return vals.mean(), vals.std(ddof=1) / np.sqrt(nrep), rm


# ---------------------------------------------------------------------- general even n (vectorized)
_TAB = {}


def _wedge_table(m, p, q):
    """for (p+q)-subsets S of range(m): list of (index of P in p-subsets, index of Q in q-subsets, sign)"""
    key = (m, p, q)
    if key in _TAB:
        return _TAB[key]
    sp_ = list(itertools.combinations(range(m), p))
    sq = list(itertools.combinations(range(m), q))
    ip = {s: i for i, s in enumerate(sp_)}
    iq = {s: i for i, s in enumerate(sq)}
    spq = list(itertools.combinations(range(m), p + q))
    A, B, S, R = [], [], [], []
    for r, s in enumerate(spq):
        for P in itertools.combinations(s, p):
            Q = tuple(x for x in s if x not in P)
            inv = sum(1 for x in P for y in Q if y < x)
            A.append(ip[P]); B.append(iq[Q]); S.append(-1.0 if inv % 2 else 1.0); R.append(r)
    tab = (np.array(A), np.array(B), np.array(S), np.array(R), len(spq))
    _TAB[key] = tab
    return tab


def _wedge(f, g, m, p, q):
    A, B, S, R, nout = _wedge_table(m, p, q)
    prod = f[:, A] * g[:, B] * S[None, :]
    out = np.zeros((f.shape[0], nout))
    np.add.at(out, (slice(None), R), prod)
    return out


def density_batch_general(nv, edges, C, t, lam):
    """float D(t), T = 1, any even n (Pfaffian by expansion along the first row)."""
    N = len(edges)
    g = nv - 1
    npts = t.shape[0]
    C = np.array(C, dtype=float)
    n = C.shape[0]
    H = lam.shape[1]
    Bm = np.zeros((N, g))
    for e, (a, b) in enumerate(edges):
        if a < g:
            Bm[e, a] += 1
        if b < g:
            Bm[e, b] -= 1
    L = np.einsum('ev,pe,ew->pvw', Bm, 1 / t, Bm)
    phi = np.linalg.solve(L, np.broadcast_to(lam, (npts, g, H)))
    y = np.einsum('ev,pvi->pei', Bm, phi) / t[:, :, None]
    M = np.einsum('re,pe,se->prs', C, t, C)
    Mi = np.linalg.inv(M)
    U = np.linalg.det(M)
    ng = N - 1
    Psi = np.einsum('re,pei->pire', C[:, :ng], y[:, :ng, :]).reshape(npts, H * n, ng)
    Gr = np.einsum('re,ae->rae', C[:, :ng], C[:, :ng])
    X = np.einsum('rae,pab,sbf->prsef', Gr, Mi, Gr)
    A2 = -0.5 * (X - np.transpose(X, (0, 1, 2, 4, 3)))
    pairs = list(itertools.combinations(range(ng), 2))
    Af = np.stack([A2[:, :, :, a, b] for (a, b) in pairs], axis=3)     # (npts, n, n, npairs)

    def pf(idx):
        if len(idx) == 0:
            return np.ones((npts, 1)), 0
        out = None
        for j in range(1, len(idx)):
            rest = idx[1:j] + idx[j + 1:]
            sub, deg = pf(rest)
            term = _wedge(Af[:, idx[0], idx[j], :], sub, ng, 2, deg) if deg else Af[:, idx[0], idx[j], :] * 1.0
            term = term * (-1.0 if (j - 1) % 2 else 1.0)
            out = term if out is None else out + term
        return out, len(idx)

    rho, deg = pf(list(range(n)))
    subsets = list(itertools.combinations(range(ng), n))
    top = np.zeros(npts)
    for si, S in enumerate(subsets):
        rest = [e for e in range(ng) if e not in S]
        inv = sum(1 for x in S for yv in rest if yv < x)
        sgn = -1.0 if inv % 2 else 1.0
        top += sgn * rho[:, si] * np.linalg.det(Psi[:, :, rest])
    return (-1) ** (N - 1) * top / (t[:, N - 1] * U ** H * np.sqrt(U))


def qmc_integral_general(nv, edges, C, plan, lam, m=14, nrep=16, seed=0, batch=2048, alpha=None, absval=False):
    rm = RatioMap(len(edges), plan)
    vals, avals = [], []
    for rep in range(nrep):
        sob = qmc.Sobol(d=len(rm.names), scramble=True, seed=seed + rep)
        Uu = sob.random_base2(m)
        tot, atot = 0.0, 0.0
        for s in range(0, Uu.shape[0], batch):
            Bb, dens = rm.sample(Uu[s:s + batch], alpha)
            t, J = rm.t_and_jac(Bb)
            D = density_batch_general(nv, edges, C, t, lam)
            w = D * J / dens
            tot += np.sum(w)
            atot += np.sum(np.abs(w))
        vals.append(tot / Uu.shape[0])
        avals.append(atot / Uu.shape[0])
    vals = np.array(vals)
    return vals.mean(), vals.std(ddof=1) / np.sqrt(nrep), float(np.mean(avals))
