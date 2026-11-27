"""Two-loop theta kernels at (H,1), any rooting, exact (Fractions).  Closed form (S3 thm-theta as
adapted in explore/msym-odd/code/theta_kernels.py, pointwise-verified there):

  theta: trivalent N, N'; arms P = 0,1,2 from N to N' with n_P edges, interior P_1..P_{n_P-1};
  y^P_j = w_P lam_N + sum_{Q,k} (w_P - delta_PQ) sigma^Q_k lam^Q_k + sum_{k<j} lam^P_k,
  sigma^Q_k = sum_{j>k} x^Q_j;   g = s (D/2) prod_P w_P^{H+1/2-n_P}  (Dirichlet measure),
  D = det[[Lam_0, Lam_1, 0], [0, Lam_1, Lam_2]],  Lam_P = (lam^P_1..lam^P_{n_P-1})  (H x (n_P-1)).

Rooting: root O = N (fan type) or an interior vertex r (tadpole type).  Survivor vertices = interior
vertices other than O (the bivalent interaction vertices); lam at N, N' (no survivors) set to 0:
  fan:     lam_N = -sum_interior lam_u;
  tadpole: lam_N = 0, lam_r = -sum_{u != r} lam_u.
Labelled graph (htk conventions, root = nv-1): BFS from O, labels nv-1, nv-2, ...; edges (v,w) with
v < w, ordered by (-v, w)  (= contraction order of hteng).  The sign s of each labelled graph is fixed
against htk.kernel_exact at a random rational point (magnitude asserted equal).

Moments: M_p(lam) = int g prod_{e at O} y_e^{p_e} / pi,  p_e an H-multi-index; returned as
{beta: Fraction}, beta = tuple over survivor vertices (increasing label) of H-tuples."""
import sys, random, itertools
from fractions import Fraction as Fr
from math import factorial

HTK = (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code"


# ------------------------------------------------------------------ sparse polynomials
def padd(P, m, c):
    if c == 0:
        return
    v = P.get(m, 0) + c
    if v == 0:
        P.pop(m, None)
    else:
        P[m] = v


def pmul(P, Q):
    out = {}
    for m1, c1 in P.items():
        for m2, c2 in Q.items():
            padd(out, tuple(a + b for a, b in zip(m1, m2)), c1 * c2)
    return out


def ppow(P, k, nv):
    out = {(0,) * nv: Fr(1)}
    for _ in range(k):
        out = pmul(out, P)
    return out


def psum(P, Q, c=1):
    out = dict(P)
    for m, v in Q.items():
        padd(out, m, c * v)
    return out


def gamma_half(a2):
    """Gamma(a2/2) for positive integer a2: (rational, power of sqrt(pi))"""
    if a2 % 2 == 0:
        return Fr(factorial(a2 // 2 - 1)), 0
    k = (a2 - 1) // 2  # Gamma(k + 1/2)
    return Fr(factorial(2 * k), 4 ** k * factorial(k)), 1


# ------------------------------------------------------------------ rooted theta
class RootedTheta:
    def __init__(self, H, arms, root, order_seed=None):
        """arms: (n0, n1, n2) edge counts; root: 'N' or (P, j) with 1 <= j <= n_P - 1.
        order_seed: None -> BFS labelling; int -> random connected search order (labelling test)."""
        self.H, self.arms, self.rootspec = H, tuple(arms), root
        # abstract vertices
        av = ['N', "N'"] + [(P, j) for P in range(3) for j in range(1, arms[P])]
        aedges = []  # (a, b, P, j) oriented along the path N -> N'
        for P in range(3):
            seq = ['N'] + [(P, j) for j in range(1, arms[P])] + ["N'"]
            for j in range(1, arms[P] + 1):
                aedges.append((seq[j - 1], seq[j], P, j))
        nv = len(av)
        self.nv = nv
        adj = {a: [] for a in av}
        for (a, b, P, j) in aedges:
            adj[a].append(b)
            adj[b].append(a)
        rootv = 'N' if root == 'N' else tuple(root)
        # BFS labelling
        order, seen = [rootv], {rootv}
        if order_seed is None:
            k = 0
            while k < len(order):
                for b in sorted(adj[order[k]], key=lambda t: str(t)):
                    if b not in seen:
                        seen.add(b)
                        order.append(b)
                k += 1
        else:
            rng = random.Random(order_seed)
            while len(order) < nv:
                front = sorted({b for a in order for b in adj[a] if b not in seen}, key=lambda t: str(t))
                b = rng.choice(front)
                seen.add(b)
                order.append(b)
        lab = {a: nv - 1 - i for i, a in enumerate(order)}
        self.lab, self.alab = lab, {v: a for a, v in lab.items()}
        E = []
        for (a, b, P, j) in aedges:
            va, vb = lab[a], lab[b]
            if va < vb:
                E.append(((va, vb), P, j, +1))
            else:
                E.append(((vb, va), P, j, -1))
        E.sort(key=lambda t: (-t[0][0], t[0][1]))
        self.edges = [t[0] for t in E]
        self.edge_path = [(t[1], t[2], t[3]) for t in E]
        self.root = nv - 1
        self.root_edges = [e for e, (v, w) in enumerate(self.edges) if w == self.root]
        self.deg = [0] * nv
        for (v, w) in self.edges:
            self.deg[v] += 1
            self.deg[w] += 1
        interior = [a for a in av if a not in ('N', "N'")]
        self.surv = sorted(lab[a] for a in interior if a != rootv)  # survivor vertices
        self.sidx = {v: i for i, v in enumerate(self.surv)}
        # polynomial variables: w0 w1 w2, x_{P,j} (j = 1..n_P), lam_{u,i}
        self.xidx = {}
        n = 3
        for P in range(3):
            for j in range(1, arms[P] + 1):
                self.xidx[(P, j)] = n
                n += 1
        self.nwx = n
        self.nvar = n + len(self.surv) * H
        self._build_currents()
        self._build_D()
        self.s = None
        self.cache = {}

    def lamvar(self, u, i):
        return self.nwx + self.sidx[u] * self.H + i

    def mono(self, **kw):
        m = [0] * self.nvar
        for k, v in kw.items():
            m[k] += v
        return tuple(m)

    def _var(self, k):
        m = [0] * self.nvar
        m[k] = 1
        return {tuple(m): Fr(1)}

    def _const(self, c):
        return {(0,) * self.nvar: Fr(c)} if c else {}

    def _build_currents(self):
        """coef[(P,j)][u] = polynomial in (w,x) multiplying lam_u (same for each component)."""
        arms = self.arms
        w = [self._var(P) for P in range(3)]
        x = {k: self._var(v) for k, v in self.xidx.items()}

        def sigma(Q, k):
            out = {}
            for j in range(k + 1, arms[Q] + 1):
                out = psum(out, x[(Q, j)])
            return out

        def raw(P, j, Q, k):
            """coefficient of lam^Q_k (formula variable) in y^P_j"""
            c = psum(pmul(w[P], sigma(Q, k)), sigma(Q, k), -1 if P == Q else 0)
            if Q == P and k < j:
                c = psum(c, self._const(1))
            return c

        coef = {}
        rs = self.rootspec
        for P in range(3):
            for j in range(1, arms[P] + 1):
                cj = {}
                for u in self.surv:
                    Q, k = self.alab[u]
                    if rs == 'N':
                        c = psum(raw(P, j, Q, k), w[P], -1)  # lam_N = -sum lam_u
                    else:
                        P0, k0 = rs
                        c = psum(raw(P, j, Q, k), raw(P, j, P0, k0), -1)  # lam_r = -sum lam_u
                    cj[u] = c
                coef[(P, j)] = cj
        self.coef = coef

    def current(self, e, i):
        """y_e^i as polynomial (linear in lam^i)"""
        P, j, sg = self.edge_path[e]
        out = {}
        for u, c in self.coef[(P, j)].items():
            out = psum(out, pmul(c, self._var(self.lamvar(u, i))), sg)
        return out

    def _lamform(self, a, i):
        """lambda^i at abstract interior vertex a as a linear polynomial in survivor lam's"""
        if a == self.rootspec or (self.rootspec != 'N' and a == tuple(self.rootspec)):
            out = {}
            for u in self.surv:
                out = psum(out, self._var(self.lamvar(u, i)), -1)
            return out
        return self._var(self.lamvar(self.lab[a], i))

    def _build_D(self):
        H, arms = self.H, self.arms
        cols = [[(P, j) for j in range(1, arms[P])] for P in range(3)]
        rows = []
        for i in range(H):
            rows.append([self._lamform(a, i) for a in cols[0]] + [self._lamform(a, i) for a in cols[1]]
                        + [{} for _ in cols[2]])
        for i in range(H):
            rows.append([{} for _ in cols[0]] + [self._lamform(a, i) for a in cols[1]]
                        + [self._lamform(a, i) for a in cols[2]])
        n = len(rows)
        assert n == 2 * H and all(len(r) == n for r in rows)
        D = {}
        for perm in itertools.permutations(range(n)):
            sg = 1
            p = list(perm)
            for a in range(n):
                for b in range(a + 1, n):
                    if p[a] > p[b]:
                        sg = -sg
            term = self._const(sg)
            for r in range(n):
                term = pmul(term, rows[r][perm[r]])
                if not term:
                    break
            if term:
                D = psum(D, term)
        self.D = D

    def integrate(self, P):
        """int prod_P w_P^{alpha_P} dw prod dx  P  -> polynomial in lam (exponent tuples over lam vars),
        divided by pi"""
        H, arms = self.H, self.arms
        out = {}
        for m, c in P.items():
            # w part: prod Gamma(H+3/2-n_P+m_P) / Gamma(H+3/2+sum m_P)
            q, kp = Fr(1), 0
            for Pp in range(3):
                g, k = gamma_half(2 * H + 3 - 2 * arms[Pp] + 2 * m[Pp])
                q *= g
                kp += k
            g, k = gamma_half(2 * H + 3 + 2 * sum(m[:3]))
            q /= g
            kp -= k
            assert kp == 2
            for Pp in range(3):
                ms = [m[self.xidx[(Pp, j)]] for j in range(1, arms[Pp] + 1)]
                num = 1
                for a in ms:
                    num *= factorial(a)
                q *= Fr(num, factorial(arms[Pp] - 1 + sum(ms)))
            padd(out, m[self.nwx:], c * q)
        return out

    def moment_raw(self, key):
        """key: tuple over root edges (edge order) of H-tuples p_e.  Returns s-less moment/pi as
        {lam exponent tuple: Fraction} (includes D/2)."""
        F = self._const(1)
        for e, p in zip(self.root_edges, key):
            for i in range(self.H):
                if p[i]:
                    F = pmul(F, ppow(self.current(e, i), p[i], self.nvar))
        lamD = {m[self.nwx:]: c for m, c in self.D.items()}
        I = self.integrate(F)
        out = {}
        for m1, c1 in I.items():
            for m2, c2 in lamD.items():
                padd(out, tuple(a + b for a, b in zip(m1, m2)), c1 * c2 / 2)
        return out

    def fix_sign(self, seed=1, check_moments=1):
        sys.path.insert(0, HTK)
        import htk
        from htk import graph
        rng = random.Random(seed)
        G = graph.Graph(self.nv, self.edges, "theta%s-%s" % (self.arms, self.rootspec))
        self.G = G
        lamv = {u: [Fr(rng.randint(-7, 7)) for _ in range(self.H)] for u in self.surv}
        lam = [lamv.get(v, [Fr(0)] * self.H) for v in range(self.nv - 1)]
        keys = [tuple((0,) * self.H for _ in self.root_edges)]
        for _ in range(check_moments):
            keys.append(tuple(tuple(rng.randint(0, 1) for _ in range(self.H)) for _ in self.root_edges))
        ratios = []
        for key in keys:
            raw = self.moment_raw(key)
            val = Fr(0)
            for m, c in raw.items():
                t = c
                for u in self.surv:
                    for i in range(self.H):
                        t *= lamv[u][i] ** m[self.sidx[u] * self.H + i]
                val += t
            powers = {}
            for e, p in zip(self.root_edges, key):
                for i in range(self.H):
                    if p[i]:
                        powers[(e, i)] = p[i]
            K = htk.kernel_exact(G, lam, T=1, powers=powers or None, check=1)
            nz = {k: v for k, v in K.c.items() if v}
            if not nz and val == 0:
                ratios.append(None)
                continue
            assert set(nz) == {1}, nz
            ratios.append(nz[1] / val)
        rs = [r for r in ratios if r is not None]
        assert rs and all(abs(r) == 1 for r in rs) and len(set(rs)) == 1, ratios
        self.s = rs[0]
        return self.s, ratios

    def moment(self, key):
        """M_key / pi as {beta: Fraction}; beta = tuple over survivors of H-tuples."""
        key = tuple(tuple(p) for p in key)
        if key in self.cache:
            return self.cache[key]
        assert self.s is not None
        raw = self.moment_raw(key)
        H = self.H
        out = {}
        for m, c in raw.items():
            beta = tuple(tuple(m[k * H:(k + 1) * H]) for k in range(len(self.surv)))
            out[beta] = c * self.s
        self.cache[key] = out
        return out


# ------------------------------------------------------------------ shapes
def theta_arms(H):
    D = H + 1
    out = []
    for a in range(1, 2 * D + 2):
        for b in range(a, 2 * D + 2):
            c = 2 * D + 1 - a - b
            if c < b:
                continue
            if a + b >= D + 1 and a + c >= D + 1 and b + c >= D + 1:
                out.append((a, b, c))
    return out


def automorphisms(nv, edges, fix):
    """number of vertex permutations preserving the (simple) edge set and fixing `fix`"""
    es = set(frozenset(e) for e in edges)
    deg = [0] * nv
    nb = [set() for _ in range(nv)]
    for (a, b) in edges:
        deg[a] += 1; deg[b] += 1
        nb[a].add(b); nb[b].add(a)
    count = 0
    verts = list(range(nv))

    def rec(i, img, used):
        nonlocal count
        if i == nv:
            count += 1
            return
        for t in verts:
            if t in used or deg[t] != deg[i]:
                continue
            if i == fix and t != fix:
                continue
            if i != fix and t == fix:
                continue
            ok = all((frozenset((img[j], t)) in es) == (frozenset((j, i)) in es) for j in range(i))
            if ok:
                img.append(t); used.add(t)
                rec(i + 1, img, used)
                img.pop(); used.discard(t)

    rec(0, [], set())
    return count


def rooted_shapes(H):
    """list of (arms, rootspec) orbit representatives with |Aut_O|"""
    out = []
    for arms in theta_arms(H):
        cands = ['N'] + [(P, j) for P in range(3) for j in range(1, arms[P])]
        seen = []
        for rs in cands:
            T = RootedTheta(H, arms, rs)
            # canonical invariant: sorted degree/arm data; dedupe by explicit isomorphism
            key = None
            dup = False
            for (T2, _) in seen:
                if _iso_rooted(T, T2):
                    dup = True
                    break
            if dup:
                continue
            seen.append((T, rs))
        for T, rs in seen:
            out.append((arms, rs, automorphisms(T.nv, T.edges, T.root)))
    return out


def _iso_rooted(T1, T2):
    if T1.nv != T2.nv or sorted(T1.deg) != sorted(T2.deg):
        return False
    nv = T1.nv
    e1 = set(frozenset(e) for e in T1.edges)
    e2 = set(frozenset(e) for e in T2.edges)
    for perm in itertools.permutations(range(nv - 1)):
        img = list(perm) + [nv - 1]
        if all(frozenset((img[a], img[b])) in e2 for (a, b) in [tuple(e) for e in e1]):
            return True
    return False
