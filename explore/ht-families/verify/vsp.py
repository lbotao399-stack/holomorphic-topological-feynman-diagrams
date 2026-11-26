"""Own series-parallel decomposition and ratio coordinates (independent of htk.spred and of hf_numeric).

Reduction of the multigraph of super-edges: parallel merge of all super-edges with the same end points,
series merge at a vertex of degree 2 with two distinct neighbours; random order if rng is given.
Tree nodes: ('E', e) leaves, ('S', children) chains, ('P', children) parallel bundles; with
flatten=True nested nodes of the same kind are absorbed.
Ratios: S-node s with children c_1..c_a: R_{c_j} = x_{s,j} R_s; P-node p with children c_1..c_k:
R_{c_i} = R_p / w_{p,i}; root R = 1.  t_e = product along the path root -> e.
Blocks: the ratios of one node (sum = 1); the first ratio of each block is eliminated.
"""
from fractions import Fraction as Fr
import random
import numpy as np
from scipy.special import gammaincinv, gammaln
from scipy.stats import qmc


class Node:
    def __init__(self, kind, children=None, edge=None):
        self.kind, self.children, self.edge = kind, children or [], edge
        if kind == 'E':
            self.L = 0
        elif kind == 'S':
            self.L = sum(c.L for c in self.children)
        else:
            self.L = sum(c.L for c in self.children) + len(self.children) - 1


def sp_decompose(nv, edges, rng=None, flatten=True):
    alive = {i: (a, b, Node('E', edge=i)) for i, (a, b) in enumerate(edges)}
    nxt = len(edges)
    verts = set(range(nv))

    def mk(kind, kids):
        ch = []
        for c in kids:
            if flatten and c.kind == kind:
                ch.extend(c.children)
            else:
                ch.append(c)
        return Node(kind, ch)

    while len(alive) > 1:
        groups = {}
        for i, (a, b, nd) in alive.items():
            groups.setdefault(frozenset((a, b)), []).append(i)
        par = [g for g in groups.values() if len(g) >= 2]
        inc = {v: [] for v in verts}
        for i, (a, b, nd) in alive.items():
            inc[a].append(i)
            inc[b].append(i)
        ser = []
        if len(verts) > 2:
            for v in verts:
                if len(inc[v]) == 2:
                    i, j = inc[v]
                    u1 = alive[i][0] if alive[i][1] == v else alive[i][1]
                    u2 = alive[j][0] if alive[j][1] == v else alive[j][1]
                    if u1 != u2:
                        ser.append((v, i, j, u1, u2))
        # deterministic order: parallel merges first, then series merges with the lightest pair (sum of loop numbers)
        ser.sort(key=lambda s: (alive[s[1]][2].L + alive[s[2]][2].L, s[0]))
        cands = [('P', g) for g in par] + [('S', s) for s in ser]
        if not cands:
            return None
        if rng is None:
            kind, c = cands[0]
        else:
            kind, c = rng.choice(cands)
        if kind == 'P':
            a, b, _ = alive[c[0]]
            if rng is not None:
                # merge a random subset of size >= 2 (allows non-maximal bundles)
                k = rng.randint(2, len(c))
                c = rng.sample(c, k)
            nd = mk('P', [alive[i][2] for i in c])
            for i in c:
                del alive[i]
            alive[nxt] = (a, b, nd)
        else:
            v, i, j, u1, u2 = c
            nd = mk('S', [alive[i][2], alive[j][2]])
            del alive[i]
            del alive[j]
            verts.discard(v)
            alive[nxt] = (u1, u2, nd)
        nxt += 1
    (a, b, root), = alive.values()
    return root


class RatioCoords:
    def __init__(self, root, N):
        self.root, self.N = root, N
        self.blocks = []        # list of (node, [ratio names])
        self.texp = [dict() for _ in range(N)]
        self.Uexp = {}          # claimed exponents of U (c1)
        self.parity = {}        # for c2/c3
        cnt = [0]

        def walk(nd, path):
            if nd.kind == 'E':
                for b, ex in path:
                    self.texp[nd.edge][b] = self.texp[nd.edge].get(b, 0) + ex
                return
            m = cnt[0]
            cnt[0] += 1
            names = []
            for i, c in enumerate(nd.children):
                b = ('x' if nd.kind == 'S' else 'w') + '%d_%d' % (m, i)
                names.append(b)
                if nd.kind == 'S':
                    self.Uexp[b] = c.L
                    walk(c, path + [(b, 1)])
                else:
                    self.Uexp[b] = -1 - c.L
                    walk(c, path + [(b, -1)])
            self.blocks.append((nd, names))

        walk(root, [])
        self.names = [b for _, nm in self.blocks for b in nm]
        self.coords = [b for _, nm in self.blocks for b in nm[1:]]
        self.elim = {}
        for _, nm in self.blocks:
            for b in nm[1:]:
                self.elim[b] = nm[0]
        assert len(self.coords) == N - 1

    def ratios_of_t(self, t):
        """inverse map: node resistances bottom-up (series sum, parallel harmonic sum), then ratios"""
        R = {}

        def res(nd):
            if nd.kind == 'E':
                R[id(nd)] = t[nd.edge]
            elif nd.kind == 'S':
                R[id(nd)] = sum(res(c) for c in nd.children)
            else:
                R[id(nd)] = 1 / sum(1 / res(c) for c in nd.children)
            return R[id(nd)]

        res(self.root)
        out = {}
        for nd, nm in self.blocks:
            for c, b in zip(nd.children, nm):
                out[b] = R[id(c)] / R[id(nd)] if nd.kind == 'S' else R[id(nd)] / R[id(c)]
        return out, R[id(self.root)]

    def complete(self, xi):
        """xi: dict coord -> value; returns dict over all ratios"""
        out = dict(xi)
        for _, nm in self.blocks:
            out[nm[0]] = 1 - sum(xi[b] for b in nm[1:])
        return out

    def t_of(self, full):
        t = []
        for e in range(self.N):
            v = Fr(1) if isinstance(next(iter(full.values())), Fr) else 1.0
            for b, ex in self.texp[e].items():
                v = v * full[b] ** ex
            t.append(v)
        return t

    def jac_rows(self, full, t):
        """rows d t / d xi_c (exact or float)"""
        rows = []
        for c in self.coords:
            b0 = self.elim[c]
            rows.append([t[e] * (self.texp[e].get(c, 0) / full[c] - self.texp[e].get(b0, 0) / full[b0])
                         for e in range(self.N)])
        return rows

    def random_point(self, rng, den=97):
        xi = {}
        for _, nm in self.blocks:
            k = len(nm)
            cuts = sorted(rng.sample(range(1, den), k - 1))
            parts = [Fr(b - a, den) for a, b in zip([0] + cuts, cuts + [den])]
            for b, p in zip(nm, parts):
                xi[b] = p
        return {c: xi[c] for c in self.coords}

    # ---- numeric
    def sample(self, u, alpha):
        """u: (B, d) uniforms (d = len(names)); alpha: dict ratio -> Dirichlet parameter.
        Returns full ratio dict of arrays and log of the sampling density (w.r.t. d xi)."""
        full = {}
        logp = np.zeros(u.shape[0])
        col = 0
        for _, nm in self.blocks:
            al = np.array([alpha[b] for b in nm])
            g = np.stack([gammaincinv(al[j], np.clip(u[:, col + j], 1e-300, 1 - 1e-16)) for j in range(len(nm))], axis=1)
            col += len(nm)
            s = g.sum(axis=1, keepdims=True)
            bvals = g / s
            for j, b in enumerate(nm):
                full[b] = bvals[:, j]
            logp += gammaln(al.sum()) - gammaln(al).sum() + ((al - 1)[None, :] * np.log(bvals)).sum(axis=1)
        return full, logp

    def t_and_jac(self, full):
        Bn = next(iter(full.values())).shape[0]
        logt = np.zeros((Bn, self.N))
        for e in range(self.N):
            for b, ex in self.texp[e].items():
                logt[:, e] += ex * np.log(full[b])
        t = np.exp(logt)
        Jm = np.zeros((Bn, self.N, self.N))
        Jm[:, 0, :] = t
        for k, c in enumerate(self.coords):
            b0 = self.elim[c]
            for e in range(self.N):
                Jm[:, k + 1, e] = t[:, e] * (self.texp[e].get(c, 0) / full[c] - self.texp[e].get(b0, 0) / full[b0])
        return t, np.linalg.det(Jm)


def qmc_integrate(rc, func, alpha, m=14, nrep=16, seed=0, batch=4096, exact_cb=None, thresh=1e6):
    """int over the product of simplices of func(t, jacdet, full) d xi, scrambled Sobol, nrep replicas.
    Points with max(t)/min(t) > thresh (where the float density loses digits) are re-evaluated by
    exact_cb(xi) (xi: dict coord -> Fraction of the sampled float coordinates) when exact_cb is given.
    Returns mean, standard error, replica values, number of exactly re-evaluated points."""
    from fractions import Fraction as Fr
    d = len(rc.names)
    vals = []
    nex = 0
    for r in range(nrep):
        sob = qmc.Sobol(d, scramble=True, seed=seed + r)
        u = sob.random_base2(m)
        acc = 0.0
        for s in range(0, u.shape[0], batch):
            full, logp = rc.sample(u[s:s + batch], alpha)
            t, J = rc.t_and_jac(full)
            f = func(t, np.abs(J), full)
            if exact_cb is not None:
                bad = np.where(t.max(axis=1) / t.min(axis=1) > thresh)[0]
                for k in bad:
                    xi = {c: Fr(float(full[c][k])) for c in rc.coords}
                    f[k] = exact_cb(xi)
                nex += len(bad)
            w = f * np.exp(-logp)
            w[~np.isfinite(w)] = 0.0
            acc += w.sum()
        vals.append(acc / u.shape[0])
    vals = np.array(vals)
    return vals.mean(), vals.std(ddof=1) / np.sqrt(nrep), vals, nex


def probe_exponents(rc, gfunc, rng, eps=(1e-4, 1e-7), npts=3):
    """leading exponent of each ratio b in g (own float evaluation): b -> eps, the other ratios of its
    block rescaled to sum 1 - eps; slope of log|g|, rounded to a half-integer, minimum over npts points."""
    out = {}
    for _, nm in rc.blocks:
        for b in nm:
            vals = []
            for _ in range(npts):
                base = {}
                for _, nm2 in rc.blocks:
                    x = np.array([rng.random() + 0.2 for _ in nm2])
                    x /= x.sum()
                    for bb, v in zip(nm2, x):
                        base[bb] = v
                gs = []
                for e in eps:
                    full = dict(base)
                    rest = sum(base[c] for c in nm if c != b)
                    for c in nm:
                        full[c] = e if c == b else base[c] * (1 - e) / rest
                    full = {k: np.array([v]) for k, v in full.items()}
                    t, J = rc.t_and_jac(full)
                    gs.append(abs(gfunc(t, np.abs(J), full)[0]))
                if gs[0] > 0 and gs[1] > 0:
                    vals.append(round(2 * np.log(gs[1] / gs[0]) / np.log(eps[1] / eps[0])) / 2)
            out[b] = min(vals) if vals else 0.0
    return out
