"""Vertex-elimination reductions with exact current maps (T = 0 densities, T = 1 chain-move factors).

Plan: list of ('S', v) / ('Y', v); after every step all parallel edges are merged (P).
    S  v of degree 2, edges e1 = v-u, e2 = v-w:  t_e1 = a1 R, t_e2 = a2 R (block (a1, a2)), new edge u -> w,
       lambda_u += a2 lambda_v, lambda_w += a1 lambda_v,
       y(v->u) = a2 lambda_v - Y,  y(v->w) = a1 lambda_v + Y.
    P  k edges p-q:  t_i = T / w_i (block (w_1..w_k)), new edge p -> q, y_i(p->q) = w_i Y.
    Y  v of degree 3, neighbours u < w < x via a, b, c; new edges u->w, u->x, w->x with R_uw, R_ux, R_wx,
       P = R_uw + R_ux + R_wx, t_a = R_uw R_ux / P, t_b = R_uw R_wx / P, t_c = R_ux R_wx / P,
       g_u = R_wx/P, g_w = R_ux/P, g_x = R_uw/P,  lambda_a += g_a lambda_v  (rational unless lambda_v = 0),
       y(v->u) = g_u lambda_v - Y_uw - Y_ux,  y(v->w) = g_w lambda_v + Y_uw - Y_wx,
       y(v->x) = g_x lambda_v + Y_ux + Y_wx.
Coordinates: the independent barycentric ratios (first ratio of each block eliminated); final edge
time T_f (gauge 1) or the edge times T of a core (partial plan).
Density (T = 0, full reduction): g = orientation * sigma_T^H eps(n,H) det[d ell / d xi],
ell = chord currents of the Kruskal tree of G; K_Gamma(lambda; 0) = int g d xi.
"""
import os
import random
import sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
P = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(P, "code"))

from sympy.polys.domains import QQ
from sympy.polys.rings import PolyRing

from htk import brute
from htk.graph import Graph
from htk.linalg import det as fdet
from htk.conventions import sign_dy
from htk.spred import _poly_det


def poly_det_dp(rows, ring):
    """Laplace expansion over column subsets (dynamic programming), rows ordered by sparsity."""
    n = len(rows)
    if n == 0:
        return ring.one
    order = sorted(range(n), key=lambda i: sum(1 for x in rows[i] if x != 0))
    perm = list(order)
    sgn = 1
    for i in range(n):
        while perm[i] != i:
            j = perm[i]
            perm[i], perm[j] = perm[j], perm[i]
            sgn = -sgn
    dp = {0: ring.one}
    for i in order:
        row = rows[i]
        nz = [c for c in range(n) if row[c] != 0]
        new = {}
        for mask, val in dp.items():
            for c in nz:
                bit = 1 << c
                if mask & bit:
                    continue
                term = row[c] * val
                if bin(mask >> (c + 1)).count("1") & 1:
                    term = -term
                m2 = mask | bit
                new[m2] = new[m2] + term if m2 in new else term
        dp = {m: v for m, v in new.items() if v != 0}
        if not dp:
            return ring.zero
    res = dp.get((1 << n) - 1, ring.zero)
    return res if sgn > 0 else -res


class NotPolynomial(Exception):
    pass


# ---------------------------------------------------------------------------- exact forward AD

class Dual:
    __slots__ = ("v", "g")

    def __init__(self, v, g):
        self.v = Fr(v)
        self.g = g

    @staticmethod
    def const(v, n):
        return Dual(v, [Fr(0)] * n)

    @staticmethod
    def var(v, j, n):
        g = [Fr(0)] * n
        g[j] = Fr(1)
        return Dual(v, g)

    def _c(self, o):
        return o if isinstance(o, Dual) else Dual(o, [Fr(0)] * len(self.g))

    def __add__(self, o):
        o = self._c(o)
        return Dual(self.v + o.v, [a + b for a, b in zip(self.g, o.g)])

    __radd__ = __add__

    def __neg__(self):
        return Dual(-self.v, [-a for a in self.g])

    def __sub__(self, o):
        return self + (-self._c(o))

    def __rsub__(self, o):
        return self._c(o) - self

    def __mul__(self, o):
        if not isinstance(o, Dual):
            o = Fr(o)
            return Dual(self.v * o, [a * o for a in self.g])
        return Dual(self.v * o.v, [self.v * b + o.v * a for a, b in zip(self.g, o.g)])

    __rmul__ = __mul__

    def __truediv__(self, o):
        o = self._c(o)
        q = self.v / o.v
        return Dual(q, [(a - q * b) / o.v for a, b in zip(self.g, o.g)])

    def __rtruediv__(self, o):
        return self._c(o) / self

    def __eq__(self, o):
        return isinstance(o, Dual) and self.v == o.v and self.g == o.g


def _is_zero(x):
    if isinstance(x, Dual):
        return x.v == 0 and all(a == 0 for a in x.g)
    return x == 0


# ---------------------------------------------------------------------------- structure

class Step:
    __slots__ = ("kind", "v", "edges", "nbrs", "new", "block", "lam_v", "p", "q", "signs")


class Elimination:
    """Structural pass of a plan on the htk.Graph G."""

    def __init__(self, G, plan):
        self.G = G
        cur = {e: ab for e, ab in enumerate(G.edges)}
        self.steps, self.blocks = [], []
        nid = G.ne
        for kind, v in plan:
            inc = sorted(e for e, ab in cur.items() if v in ab)
            st = Step()
            st.kind, st.v = kind, v
            other = [cur[e][1] if cur[e][0] == v else cur[e][0] for e in inc]
            if kind == "S":
                if len(inc) != 2 or other[0] == other[1]:
                    raise ValueError("S at %d: degree %d" % (v, len(inc)))
                st.edges, st.nbrs = inc, other
                u, w = other
                st.new = [nid]
                cur[nid] = (u, w)
                nid += 1
                st.block = ["a%d_1" % (len(self.steps) + 1), "a%d_2" % (len(self.steps) + 1)]
                self.blocks.append(st.block)
            elif kind == "Y":
                if len(inc) != 3 or len(set(other)) != 3:
                    raise ValueError("Y at %d: degree %d" % (v, len(inc)))
                order = sorted(range(3), key=lambda i: other[i])
                st.edges = [inc[i] for i in order]
                st.nbrs = [other[i] for i in order]
                u, w, x = st.nbrs
                st.new = [nid, nid + 1, nid + 2]
                cur[nid], cur[nid + 1], cur[nid + 2] = (u, w), (u, x), (w, x)
                nid += 3
                st.block = None
            else:
                raise ValueError(kind)
            for e in inc:
                del cur[e]
            self.steps.append(st)
            # parallel merges
            groups = {}
            for e, (a, b) in cur.items():
                groups.setdefault((min(a, b), max(a, b)), []).append(e)
            for (p, q), es in sorted(groups.items()):
                if len(es) < 2:
                    continue
                ps = Step()
                ps.kind, ps.v, ps.p, ps.q = "P", None, p, q
                ps.edges = sorted(es)
                ps.signs = [1 if cur[e] == (p, q) else -1 for e in ps.edges]
                ps.new = [nid]
                m = len(self.steps) + 1
                ps.block = ["w%d_%d" % (m, i + 1) for i in range(len(es))]
                self.blocks.append(ps.block)
                for e in es:
                    del cur[e]
                cur[nid] = (p, q)
                nid += 1
                self.steps.append(ps)
        self.final = dict(cur)
        self.coords = [b for blk in self.blocks for b in blk[1:]]
        self.bary = [b for blk in self.blocks for b in blk]
        self.k = len(self.coords)

    # ------------------------------------------------------------------ generic passes
    def bary_values(self, xi, one):
        """xi: dict coord -> number (independent coordinates); returns all barycentric values."""
        val = {}
        for blk in self.blocks:
            s = one
            for b in blk[1:]:
                val[b] = xi[b]
                s = s - xi[b]
            val[blk[0]] = s
        return val

    def times(self, bv, Tcore):
        """Backward times: Tcore = dict edge id -> time of the final network."""
        t = dict(Tcore)
        for st in reversed(self.steps):
            if st.kind == "S":
                R = t[st.new[0]]
                t[st.edges[0]] = bv[st.block[0]] * R
                t[st.edges[1]] = bv[st.block[1]] * R
            elif st.kind == "P":
                T = t[st.new[0]]
                for e, b in zip(st.edges, st.block):
                    t[e] = T / bv[b]
            else:
                Ruw, Rux, Rwx = (t[e] for e in st.new)
                Ps = Ruw + Rux + Rwx
                a, b, c = st.edges
                t[a] = Ruw * Rux / Ps
                t[b] = Ruw * Rwx / Ps
                t[c] = Rux * Rwx / Ps
        return [t[e] for e in range(self.G.ne)], t

    def sources(self, lam0, bv, tall, zero, polynomial=False):
        """Forward source transfer. lam0: dict vertex -> H-vector (all vertices, sum zero).
        Records lam_v at each eliminated vertex in self._lamv[step index].  Returns final dict."""
        lam = {v: list(x) for v, x in lam0.items()}
        lamv = {}
        for i, st in enumerate(self.steps):
            if st.kind == "S":
                lv = lam.pop(st.v)
                lamv[i] = lv
                u, w = st.nbrs
                a1, a2 = bv[st.block[0]], bv[st.block[1]]
                lam[u] = [x + a2 * y for x, y in zip(lam[u], lv)]
                lam[w] = [x + a1 * y for x, y in zip(lam[w], lv)]
            elif st.kind == "Y":
                lv = lam.pop(st.v)
                lamv[i] = lv
                if all(_is_zero(y) for y in lv):
                    continue
                if polynomial:
                    raise NotPolynomial("sourced star-triangle at %d" % st.v)
                Ruw, Rux, Rwx = (tall[e] for e in st.new)
                Ps = Ruw + Rux + Rwx
                g = {st.nbrs[0]: Rwx / Ps, st.nbrs[1]: Rux / Ps, st.nbrs[2]: Ruw / Ps}
                for a in st.nbrs:
                    lam[a] = [x + g[a] * y for x, y in zip(lam[a], lv)]
        self._lamv = lamv
        return lam

    def currents(self, ycore, bv, tall, polynomial=False):
        """Backward currents from the currents ycore (dict edge id -> H-vector) of the final network."""
        y = dict(ycore)
        for i in range(len(self.steps) - 1, -1, -1):
            st = self.steps[i]
            if st.kind == "S":
                Y = y[st.new[0]]
                lv = self._lamv[i]
                a1, a2 = bv[st.block[0]], bv[st.block[1]]
                e1, e2 = st.edges
                vu = [a2 * l - Yc for l, Yc in zip(lv, Y)]
                vw = [a1 * l + Yc for l, Yc in zip(lv, Y)]
                y[e1] = vu if self.G_edge(e1)[0] == st.v else [-x for x in vu]
                y[e2] = vw if self.G_edge(e2)[0] == st.v else [-x for x in vw]
            elif st.kind == "P":
                Y = y[st.new[0]]
                for e, b, s in zip(st.edges, st.block, st.signs):
                    y[e] = [bv[b] * Yc if s > 0 else -(bv[b] * Yc) for Yc in Y]
            else:
                Yuw, Yux, Ywx = (y[e] for e in st.new)
                lv = self._lamv[i]
                H = len(Yuw)
                if all(_is_zero(l) for l in lv):
                    g = None
                else:
                    if polynomial:
                        raise NotPolynomial("sourced star-triangle")
                    Ruw, Rux, Rwx = (tall[e] for e in st.new)
                    Ps = Ruw + Rux + Rwx
                    g = (Rwx / Ps, Rux / Ps, Ruw / Ps)
                vu = [-Yuw[j] - Yux[j] for j in range(H)]
                vw = [Yuw[j] - Ywx[j] for j in range(H)]
                vx = [Yux[j] + Ywx[j] for j in range(H)]
                if g is not None:
                    vu = [vu[j] + g[0] * lv[j] for j in range(H)]
                    vw = [vw[j] + g[1] * lv[j] for j in range(H)]
                    vx = [vx[j] + g[2] * lv[j] for j in range(H)]
                for e, vec in zip(st.edges, (vu, vw, vx)):
                    y[e] = vec if self.G_edge(e)[0] == st.v else [-x for x in vec]
        return [y[e] for e in range(self.G.ne)]

    def G_edge(self, e):
        """orientation of edge e at the time it was removed (original or created)."""
        if e < self.G.ne:
            return self.G.edges[e]
        return self._created[e]

    def _record_created(self):
        cr = {}
        for st in self.steps:
            if st.kind == "S":
                cr[st.new[0]] = (st.nbrs[0], st.nbrs[1])
            elif st.kind == "P":
                cr[st.new[0]] = (st.p, st.q)
            else:
                u, w, x = st.nbrs
                cr[st.new[0]], cr[st.new[1]], cr[st.new[2]] = (u, w), (u, x), (w, x)
        self._created = cr


def lam_dict(G, lam):
    """sources on all vertices (root = -sum) from a list over non-root vertices."""
    lam = brute.sources(G, lam)
    H = len(lam[0]) if lam else 0
    d = {v: list(lam[v]) for v in range(G.nv - 1)}
    d[G.root] = [-sum((lam[v][i] for v in range(G.nv - 1)), Fr(0)) for i in range(H)]
    return d, H


# ---------------------------------------------------------------------------- full reduction (T = 0)

class FullReduction:
    """Plan eliminating all but two vertices.  Polynomial density at T = 0 if no sourced Y."""

    def __init__(self, G, plan):
        self.G = G
        self.el = Elimination(G, plan)
        self.el._record_created()
        if len(self.el.final) != 1:
            raise ValueError("plan leaves %d edges" % len(self.el.final))
        (self.fe, (self.pf, self.qf)), = self.el.final.items()
        if self.el.k != G.ne - 1:
            raise AssertionError("coordinates %d != |I| - 1" % self.el.k)

    def random_point(self, rng):
        pt = {}
        for blk in self.el.blocks:
            ns = [rng.randint(1, 40) for _ in blk]
            s = sum(ns)
            for b, k in zip(blk, ns):
                pt[b] = Fr(k, s)
        return {c: pt[c] for c in self.el.coords}

    def _eval(self, xi, lam0, duals):
        el = self.el
        k = el.k
        if duals:
            X = {c: Dual.var(xi[c], j, k) for j, c in enumerate(el.coords)}
            one = Dual.const(1, k)
        else:
            X, one = dict(xi), Fr(1)
        bv = el.bary_values(X, one)
        t, tall = el.times(bv, {self.fe: one})
        lamc = {v: [one * x for x in vec] for v, vec in lam0.items()}
        lam = el.sources(lamc, bv, tall, None)
        Yf = lam[self.pf]
        if any(not _is_zero(a + b) for a, b in zip(lam[self.pf], lam[self.qf])):
            raise AssertionError("final sources do not balance")
        y = el.currents({self.fe: Yf}, bv, tall)
        return t, y

    def check_point(self, lam, xi):
        """Exact pointwise check: currents = Kirchhoff; g(xi) = D(t) |det[t; dt/dxi]|.
        Returns (g, brute value)."""
        G = self.G
        lam0, H = lam_dict(G, lam)
        t, y = self._eval(xi, lam0, True)
        tv = [x.v for x in t]
        yb = brute.currents(G, tv, lam)
        yv = [[c.v for c in row] for row in y]
        if yb != yv:
            raise AssertionError("currents differ from Kirchhoff")
        Jt = fdet([[x.v for x in t]] + [[x.g[j] for x in t] for j in range(self.el.k)])
        Db = brute.density(G, tv, lam)
        cot = G.cotree()
        Jl = fdet([[y[c][i].g[j] for j in range(self.el.k)] for c in cot for i in range(H)])
        g = (1 if Jt > 0 else -1) * G.sigma() ** H * sign_dy(G.loops, H) * Jl
        return g, Db.coeff * abs(Jt), Jt

    # -------------------------------------------------------- polynomial density and integral
    def poly_density(self, lam):
        G, el = self.G, self.el
        lam0, H = lam_dict(G, lam)
        ring = PolyRing(el.coords, QQ)
        gens = dict(zip(el.coords, ring.gens))
        one = ring.one
        bv = el.bary_values(gens, one)
        lamc = {v: [one * QQ(int(x.numerator), int(x.denominator)) for x in vec] for v, vec in lam0.items()}
        lam = el.sources(lamc, bv, None, None, polynomial=True)
        Yf = lam[self.pf]
        y = el.currents({self.fe: Yf}, bv, None, polynomial=True)
        cot = G.cotree()
        J = [[y[c][i].diff(gens[x]) for x in el.coords] for c in cot for i in range(H)]
        dets = poly_det_dp(J, ring) if len(J) >= 8 else _poly_det(J, ring)
        rng = random.Random(5)
        signs = set()
        for _ in range(4):
            xi = self.random_point(rng)
            X = {c: Dual.var(xi[c], j, el.k) for j, c in enumerate(el.coords)}
            bvd = el.bary_values(X, Dual.const(1, el.k))
            t, _ = el.times(bvd, {self.fe: Dual.const(1, el.k)})
            Jt = fdet([[x.v for x in t]] + [[x.g[j] for x in t] for j in range(el.k)])
            signs.add(1 if Jt > 0 else -1)
        if len(signs) != 1:
            raise AssertionError("orientation not constant")
        self.orientation = signs.pop()
        s = self.orientation * G.sigma() ** H * sign_dy(G.loops, H)
        self.ring, self.y_poly = ring, y
        return dets * s

    def integrate(self, g, F=None):
        """int g F over the product of simplices (Dirichlet, integer exponents)."""
        from math import factorial
        el = self.el
        idx = {c: i for i, c in enumerate(el.coords)}
        blocks = [[idx[b] for b in blk[1:]] for blk in el.blocks]
        p = g if F is None else g * F
        tot = Fr(0)
        for mon, c in p.terms():
            val = Fr(int(c.numerator), int(c.denominator))
            for blk in blocks:
                ks = [mon[j] for j in blk]
                num = 1
                for kk in ks:
                    num *= factorial(kk)
                val *= Fr(num, factorial(len(blk) + sum(ks)))
            tot += val
        return tot


# ---------------------------------------------------------------------------- partial reduction to a core

class PartialReduction:
    """Plan leaving a core network; the core is relabelled as an htk.Graph (root = original root if it
    survives, else the largest remaining label)."""

    def __init__(self, G, plan):
        self.G = G
        self.el = Elimination(G, plan)
        self.el._record_created()
        fin = self.el.final
        verts = sorted({v for ab in fin.values() for v in ab})
        root = G.root if G.root in verts else verts[-1]
        order = [v for v in verts if v != root] + [root]
        self.vmap = {v: i for i, v in enumerate(order)}
        self.inv = order
        self.core_edges = sorted(fin)
        self.core = Graph(len(order), [(self.vmap[fin[e][0]], self.vmap[fin[e][1]]) for e in self.core_edges],
                          "core")

    def check(self, lam, xi, Tc, verbose=False):
        """T = 0 factorization identity at (xi, Tc):
           D_G(t) J_aff (-1)^k = sigma_G^H eps(n,H) sigma_c^H eps(n_c,H) J_poly D_core(Tc; lam_core(xi)).
        Returns (lhs, rhs)."""
        G, el, core = self.G, self.el, self.core
        lam0, H = lam_dict(G, lam)
        k, Nc = el.k, core.ne
        nvar = k + Nc
        # (xi, T) duals for t and J_aff
        X = {c: Dual.var(xi[c], j, nvar) for j, c in enumerate(el.coords)}
        one = Dual.const(1, nvar)
        bv = el.bary_values(X, one)
        Tmap = {e: Dual.var(Tc[i], k + i, nvar) for i, e in enumerate(self.core_edges)}
        t, tall = el.times(bv, Tmap)
        J_aff = fdet([[x.g[j] for j in range(nvar)] for x in t])
        tv = [x.v for x in t]
        D_G = brute.density(G, tv, lam).coeff
        # sources of the core (numbers)
        Xn = {c: Fr(xi[c]) for c in el.coords}
        bvn = el.bary_values(Xn, Fr(1))
        _, talln = el.times(bvn, {e: Fr(Tc[i]) for i, e in enumerate(self.core_edges)})
        lamc = el.sources(lam0, bvn, talln, None)
        lam_core = [lamc[self.inv[i]] for i in range(core.nv - 1)]
        D_c = brute.density(core, [Fr(x) for x in Tc], lam_core).coeff
        # J_poly: duals over (xi, ell_core)
        nc_ = core.loops
        nv2 = k + nc_ * H
        X2 = {c: Dual.var(xi[c], j, nv2) for j, c in enumerate(el.coords)}
        one2 = Dual.const(1, nv2)
        bv2 = el.bary_values(X2, one2)
        _, tall2 = el.times(bv2, {e: one2 * Fr(Tc[i]) for i, e in enumerate(self.core_edges)})
        lam2 = el.sources({v: [one2 * x for x in vec] for v, vec in lam0.items()}, bv2, tall2, None)
        lamcore2 = [lam2[self.inv[i]] for i in range(core.nv - 1)]
        tree = core.spanning_tree()
        cot = core.cotree(tree)
        B = core.incidence()
        BT = [[Fr(B[e][v]) for e in tree] for v in range(core.nv - 1)]      # (nv-1) x (nv-1): B_T^t
        from htk.linalg import inverse
        BTinv = inverse(BT)                                                 # y_T = BTinv lam
        C = core.fundamental_cycles(tree)
        ycore = {}
        for pos, e in enumerate(tree):
            vec = []
            for i in range(H):
                s = Dual.const(0, nv2)
                for v in range(core.nv - 1):
                    if BTinv[pos][v] != 0:
                        s = s + BTinv[pos][v] * lamcore2[v][i]
                vec.append(s)
            ycore[e] = vec
        for e in cot:
            ycore[e] = [Dual.const(0, nv2) for _ in range(H)]
        for r, c in enumerate(cot):
            for i in range(H):
                lv = Dual.var(0, k + r * H + i, nv2)
                for e in range(core.ne):
                    if C[r][e]:
                        ycore[e][i] = ycore[e][i] + C[r][e] * lv
        # set ell values to the actual core chord currents (point of evaluation)
        yc_num = brute.currents(core, [Fr(x) for x in Tc], lam_core)
        for r, c in enumerate(cot):
            for i in range(H):
                shift = yc_num[c][i]
                for e in range(core.ne):
                    if C[r][e]:
                        ycore[e][i] = ycore[e][i] + C[r][e] * shift
        y_in = {self.core_edges[e]: ycore[e] for e in range(core.ne)}
        yG = el.currents(y_in, bv2, tall2)
        # consistency of the reconstructed currents with Kirchhoff on G
        yb = brute.currents(G, tv, lam)
        if [[c.v for c in row] for row in yG] != yb:
            raise AssertionError("reconstructed currents differ from Kirchhoff")
        cotG = G.cotree()
        J_poly = fdet([[yG[c][i].g[j] for j in range(nv2)] for c in cotG for i in range(H)])
        s = (G.sigma() ** H * sign_dy(G.loops, H)) * (core.sigma() ** H * sign_dy(core.loops, H))
        lhs = D_G * J_aff * (-1) ** k
        rhs = s * J_poly * D_c
        return lhs, rhs, J_poly
