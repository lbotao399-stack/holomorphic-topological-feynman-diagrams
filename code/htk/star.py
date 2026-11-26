"""Star-triangle (Y-Delta) moves of the circuit reduction, T in {0, 1}, with exact pointwise
self-checks against htk.brute.

Sources: LLC = sources/LLC_main.tex, Theorem b-01159 (single source-free move), Theorem b-01173
(paired move), Theorem thm-sourced-star (sourced move).  Derivations and statements:
explore/circuit-general/ht-d2/REPORT.md.

Coordinates of a move: fractions f (2 for the single move, 4 for the paired move) and ALL reduced
edge times T (affine, no gauge).  The map (f, T) -> t is homogeneous of degree one in T, so

    iota_E prod_e dt_e = J_aff  df_1 ^ ... ^ df_k ^ iota_E' prod_e' dT_e'     (k even),
    J_aff = det d(t_0, ..., t_{N-1}) / d(f_1, ..., f_k, T_0, ..., T_{N'-1}),

and the density of Omega_Gamma in these coordinates is D'(f, T) = D_Gamma(t(f, T)) * J_aff, with
Omega_Gamma = D' df ^ iota_E' prod dT.  Every formula of this module is compared with that brute
force value (sign and square, exact Fractions).

Single move (LLC configuration): degree-3 vertex v with neighbours u, w, x; ux, wx present, uw absent.
  R_uw = T_uw, R_ux = T_ux/(1-r1), R_wx = T_wx/(1-r2), P = R_uw + R_ux + R_wx,
  t_vu = R_uw R_ux/P, t_vw = R_uw R_wx/P, t_vx = R_ux R_wx/P, t_ux = T_ux/r1, t_wx = T_wx/r2.
Paired move: degree-3 vertices v1, v2 joined to u, A, B; AB present, uA, uB absent.
  R1_uA = T_uA/r1, R2_uA = T_uA/(1-r1), R1_uB = T_uB/r2, R2_uB = T_uB/(1-r2),
  R1_AB = T_AB/r3, R2_AB = T_AB/r4, t_AB = T_AB/(1-r3-r4), star i from (Ri_uA, Ri_uB, Ri_AB).

T = 1 (any H): Kron decomposition of the topological Gaussian.  The vertex potentials X of the
topological direction are Gaussian with density exp(-X^t Lambda X); conditionally on the other
vertices X_v = Xbar_v + xi_v, xi_v ~ N(0, 1/(2 G_v)), G_v = sum_a 1/t_va, independent of the
reduced potentials, whose law is the reduced network's (Lambda_red = Schur complement).
Holomorphic potentials: Phi_v = Phibar_v + lambda_v/G_v, reduced sources lambda_a + lambda_v/(t_va G_v).
    Omega_Gamma = U^{-H-1/2} E_{X_red, xi}[prod_r psi^tau_r] ^ prod_i prod_r psi^i_r,
    psi_r = sum_e C_re (Phi_e/t_e) dt_e,   U_Gamma = U_red * prod_Gamma t / prod_red T * prod_v G_v.
"""
import itertools
import random
from fractions import Fraction as Fr

from . import brute
from .graph import Graph
from .linalg import det, inverse, solve, wedge, wedge_all, form_add, form_scale, one_form

__all__ = ["Dual", "SingleMove", "PairedMove", "find_single", "find_paired", "wick",
           "y_to_delta", "delta_to_y", "density_prime_brute"]


# ---------------------------------------------------------------------- exact forward-mode AD

class Dual:
    """value + sum_j grad[j] d xi_j, exact (Fractions)."""
    __slots__ = ("v", "g")

    def __init__(self, v, g=None, n=None):
        self.v = Fr(v)
        self.g = list(g) if g is not None else [Fr(0)] * n

    @staticmethod
    def var(v, j, n):
        g = [Fr(0)] * n
        g[j] = Fr(1)
        return Dual(v, g)

    def _c(self, o):
        return o if isinstance(o, Dual) else Dual(o, None, len(self.g))

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
        o = self._c(o)
        return Dual(self.v * o.v, [self.v * b + o.v * a for a, b in zip(self.g, o.g)])

    __rmul__ = __mul__

    def __truediv__(self, o):
        o = self._c(o)
        q = self.v / o.v
        return Dual(q, [(a - q * b) / o.v for a, b in zip(self.g, o.g)])

    def __rtruediv__(self, o):
        return self._c(o) / self

    def dform(self, offset=0):
        """the one-form d(self) = sum_j g_j d xi_j as an htk form (bit j + offset)."""
        return {1 << (j + offset): c for j, c in enumerate(self.g) if c != 0}


def delta_to_y(Ruw, Rux, Rwx):
    """Delta -> Y: star times (t_vu, t_vw, t_vx) at the centre v of the triangle with resistances
    R_uw, R_ux, R_wx (t_va = product of the two triangle resistances at a / sum)."""
    P = Ruw + Rux + Rwx
    return Ruw * Rux / P, Ruw * Rwx / P, Rux * Rwx / P


def y_to_delta(tvu, tvw, tvx):
    """Y -> Delta: (R_uw, R_ux, R_wx) = Sigma/(t_vx, t_vw, t_vu), Sigma = t_vu t_vw + t_vw t_vx + t_vx t_vu."""
    S = tvu * tvw + tvw * tvx + tvx * tvu
    return S / tvx, S / tvw, S / tvu


# ---------------------------------------------------------------------- Gaussian expectation (Wick)

def _pair(a, b, cov):
    """sum_{k,l} cov[k][l] a_k ^ b_l for items a, b: dicts {gaussian index: form}."""
    out = {}
    for k, fa in a.items():
        for l, fb in b.items():
            c = cov[k][l]
            if c != 0:
                out = form_add(out, form_scale(wedge(fa, fb), c))
    return out


def wick(items, degs, cov, const=None):
    """E[prod_j items_j] for items linear in a centred Gaussian vector g (covariance cov):
    item_j = sum_k g_k items_j[k] (forms of degree degs[j]).  Ordered product with Koszul signs."""
    if not items:
        return {0: Fr(1)}
    if len(items) % 2:
        return {}
    a, rest, drest = items[0], items[1:], degs[1:]
    out = {}
    acc_deg = 0
    for j in range(len(rest)):
        sign = -1 if (degs[1 + j] * acc_deg) % 2 else 1
        p = _pair(a, rest[j], cov)
        if p:
            sub = wick(rest[:j] + rest[j + 1:], drest[:j] + drest[j + 1:], cov)
            if sub:
                out = form_add(out, form_scale(wedge(p, sub), sign))
        acc_deg += degs[1 + j]
    return out


def _top_density(form, k, nT, T):
    """D' with form = D' df_1..df_k ^ iota_E' (dT_0..dT_{nT-1}) (f bits 0..k-1, T bits k..k+nT-1);
    checks that the form is basic (all e0 agree)."""
    full = (1 << (k + nT)) - 1
    vals = []
    for e0 in range(nT):
        c = form.get(full ^ (1 << (k + e0)), Fr(0))
        vals.append((-c if e0 % 2 else c) / Fr(T[e0]))
    if any(v != vals[0] for v in vals):
        raise ArithmeticError("form is not basic: %s" % vals)
    return vals[0]


def density_prime_brute(G, t, J_aff, lam, T):
    """D' = D_Gamma(t) * J_aff as a brute.Density."""
    D = brute.density(G, t, lam, T=T)
    return D.scaled(J_aff)


# ---------------------------------------------------------------------- moves

class _Move:
    """Common part: subclasses define self.G, self.centres, self.red (reduced Graph), self.k
    (number of fractions), self.vmap (old vertex -> new vertex, centres -> None), and
    tmap_dual(f, T) -> list of Dual (times of G), red_sources, current map."""

    # ---- coordinates and the time map
    def random_point(self, rng=None, den=40):
        rng = rng or random
        f = self._random_fractions(rng, den)
        T = [Fr(rng.randint(1, den), rng.randint(1, 9)) for _ in range(self.red.ne)]
        return f, T

    def tmap(self, f, T):
        return [x.v for x in self.tmap_dual(f, T, derivs=False)]

    def jacobian(self, f, T):
        """J_aff = det d t / d (f, T) (outputs in G's edge order, inputs f then T)."""
        tt = self.tmap_dual(f, T)
        return det([x.g for x in tt])

    def _duals(self, f, T, derivs=True):
        n = self.k + len(T)
        if not derivs:
            return [Dual(x, None, 0) for x in f], [Dual(x, None, 0) for x in T]
        fd = [Dual.var(x, j, n) for j, x in enumerate(f)]
        Td = [Dual.var(x, self.k + j, n) for j, x in enumerate(T)]
        return fd, Td

    # ---- U
    def U_ratio(self, f, T):
        """U_Gamma / U_red = prod_Gamma t / prod_red T * prod_v G_v (Kron reduction of the
        grounded Laplacian; matrix-tree)."""
        t = self.tmap(f, T)
        r = Fr(1)
        for x in t:
            r *= x
        for x in T:
            r /= x
        for v in self.centres:
            r *= sum((1 / t[e] for e in self._star_edges(v)), Fr(0))
        return r

    def _star_edges(self, v):
        return [e for e, (a, b) in enumerate(self.G.edges) if a == v or b == v]

    # ---- the reduced problem
    def red_sources(self, lam, t):
        """Reduced sources (H-vectors, list over the non-root vertices of red): transferred
        lambda_a + sum_v lambda_v / (t_va G_v)."""
        G = self.G
        lam = brute.sources(G, lam)
        H = len(lam[0]) if lam else 0
        full = lam + [[-sum((lam[v][i] for v in range(G.nv - 1)), Fr(0)) for i in range(H)]]
        new = [[Fr(0)] * H for _ in range(self.red.nv)]
        for v in range(G.nv):
            if self.vmap[v] is not None:
                for i in range(H):
                    new[self.vmap[v]][i] += full[v][i]
        for v in self.centres:
            Gv = sum((1 / t[e] for e in self._star_edges(v)), Fr(0))
            for e in self._star_edges(v):
                a = G.edges[e][1] if G.edges[e][0] == v else G.edges[e][0]
                for i in range(H):
                    new[self.vmap[a]][i] += full[v][i] / (t[e] * Gv)
        return new[:-1], H

    # ---- potentials of G from the reduced potentials (+ local fluctuations)
    def lifted_potentials(self, t, phi_red, extra):
        """Potentials of G's vertices (root = 0): reduced vertices keep phi_red; a centre v gets
        Phibar_v + extra[v] (extra: lambda_v/G_v, or the Gaussian fluctuation)."""
        G = self.G
        out = [None] * G.nv
        for v in range(G.nv):
            if self.vmap[v] is not None:
                out[v] = phi_red[self.vmap[v]]
        for v in self.centres:
            Gv = Fr(0)
            acc = None
            for e in self._star_edges(v):
                a = G.edges[e][1] if G.edges[e][0] == v else G.edges[e][0]
                c = 1 / t[e]
                Gv += c
                acc = _lin_add(acc, _lin_scale(out[a], c))
            out[v] = _lin_add(_lin_scale(acc, 1 / Gv), extra[v])
        return out

    # ---- T = 1 (any H) and T = 0 densities from the reduced problem (generic Kron path)
    def density_kron(self, f, T, lam, Tdir, xi_var=1):
        """D'(f, T) computed from the reduced network and the local data only:
        Kron decomposition of the topological Gaussian (T = 1) and source transfer."""
        G, red = self.G, self.red
        tt = self.tmap_dual(f, T)
        t = [x.v for x in tt]
        nvar = self.k + len(T)
        lamr, H = self.red_sources(lam, t)
        # holomorphic potentials (deterministic), as 'linear objects' over a dummy index
        phi_red = brute.potentials(red, T, lamr) if H else [[] for _ in range(red.nv)]
        lam_full = brute.sources(G, lam)
        lam_full = lam_full + [[-sum((lam_full[v][i] for v in range(G.nv - 1)), Fr(0)) for i in range(H)]]
        C = G.cycle_matrix()
        n = len(C)
        dt = [x.dform() for x in tt]
        hol = []
        for i in range(H):
            pr = [{"c": phi_red[w][i]} for w in range(red.nv)]
            extra = {}
            for v in self.centres:
                Gv = sum((1 / t[e] for e in self._star_edges(v)), Fr(0))
                extra[v] = {"c": lam_full[v][i] / Gv}
            pot = self.lifted_potentials(t, pr, extra)
            for r in range(n):
                form = {}
                for e, (a, b) in enumerate(G.edges):
                    if C[r][e]:
                        d = _lin_add(pot[a], _lin_scale(pot[b], -1)).get("c", Fr(0))
                        form = form_add(form, form_scale(dt[e], C[r][e] * d / t[e]))
                hol.append(form)
        U = brute.symanzik_U(red, T) * self.U_ratio(f, T)
        if Tdir == 0:
            top = wedge_all(hol)
            D = _top_density(top, self.k, len(T), T) / U ** H
            return brute.Density(D)
        # topological Gaussian: reduced potentials X_w (w < red root) and xi_v
        g = red.nv - 1
        Lr = brute.laplacian(red, T)
        covX = inverse(Lr)
        idx = {("X", w): w for w in range(g)}
        for j, v in enumerate(self.centres):
            idx[("xi", v)] = g + j
        m = g + len(self.centres)
        cov = [[Fr(0)] * m for _ in range(m)]
        for a in range(g):
            for b in range(g):
                cov[a][b] = covX[a][b] / 2
        for j, v in enumerate(self.centres):
            Gv = sum((1 / t[e] for e in self._star_edges(v)), Fr(0))
            cov[g + j][g + j] = Fr(xi_var) / (2 * Gv)
        pr = [{idx[("X", w)]: Fr(1)} if w < g else {} for w in range(red.nv)]
        extra = {v: {idx[("xi", v)]: Fr(1)} for v in self.centres}
        pot = self.lifted_potentials(t, pr, extra)
        items = []
        for r in range(n):
            item = {}
            for e, (a, b) in enumerate(G.edges):
                if C[r][e]:
                    d = _lin_add(pot[a], _lin_scale(pot[b], -1))
                    for kk, c in d.items():
                        item[kk] = form_add(item.get(kk, {}), form_scale(dt[e], C[r][e] * c / t[e]))
            items.append({kk: fm for kk, fm in item.items() if fm})
        pf = wick(items, [1] * n, cov)
        top = wedge(pf, wedge_all(hol))
        coeff = _top_density(top, self.k, len(T), T) / U ** H
        return brute.Density(coeff, U)

    def selfcheck(self, lam, H, Tdir, npoints=2, seed=0, methods=("kron",)):
        """At random rational points: (i) forward(tmap(f, T)) = (f, T); (ii) currents of G = current
        map of the reduced currents; (iii) U_Gamma = U_red * U_ratio; (iv) D' of each method equals
        the brute force D_Gamma(t) J_aff (sign and square)."""
        rng = random.Random(seed)
        out = []
        for _ in range(npoints):
            f, T = self.random_point(rng)
            t = self.tmap(f, T)
            f2, T2 = self.forward(t)
            assert f2 == list(f) and T2 == list(T), "forward map"
            assert brute.symanzik_U(self.G, t) == brute.symanzik_U(self.red, T) * self.U_ratio(f, T), "U ratio"
            if lam is not None:
                y = brute.currents(self.G, t, lam)
                assert y == self.currents(f, T, lam), "current map"
            ref = density_prime_brute(self.G, t, self.jacobian(f, T), lam, Tdir)
            for meth in methods:
                val = getattr(self, "density_" + meth)(f, T, lam, Tdir)
                if not val.same_as(ref):
                    raise AssertionError("%s: %r != brute %r at f=%s" % (meth, val, ref, f))
            out.append((f, T, ref))
        return out


def _lin_add(a, b):
    if a is None:
        return dict(b) if b is not None else None
    if b is None:
        return dict(a)
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, Fr(0)) + v
    return {k: v for k, v in out.items() if v != 0}


def _lin_scale(a, c):
    if a is None:
        return None
    return {k: v * c for k, v in a.items() if v * c != 0}


def _edge_between(G, p, q):
    es = [e for e, (a, b) in enumerate(G.edges) if {a, b} == {p, q}]
    return es


def _neighbours(G, v):
    out = []
    for a, b in G.edges:
        if a == v:
            out.append(b)
        elif b == v:
            out.append(a)
    return out


# ---------------------------------------------------------------------- single move

class SingleMove(_Move):
    """Source-free or sourced single star-triangle move (LLC Theorems b-01159, thm-sourced-star).
    Reduced graph: vertices of G without v (order kept, root last), edges
    (u->w, u->x, w->x) first, then the other edges of G in their order."""

    def __init__(self, G, v, u=None, w=None, x=None):
        if v == G.root:
            raise ValueError("the centre must not be the root")
        nb = _neighbours(G, v)
        if len(nb) != 3 or len(set(nb)) != 3:
            raise ValueError("centre must have degree three with distinct neighbours")
        if x is None:
            cand = [c for c in nb if all(len(_edge_between(G, c, o)) == 1 for o in nb if o != c)]
            others = lambda c: [o for o in nb if o != c]
            cand = [c for c in cand if len(_edge_between(G, *others(c))) == 0]
            if not cand:
                raise ValueError("no LLC configuration (ux, wx present, uw absent) at %d" % v)
            x = cand[0]
            u, w = sorted(o for o in nb if o != x)
        if sorted([u, w, x]) != sorted(nb):
            raise ValueError("u, w, x must be the neighbours of v")
        if len(_edge_between(G, u, x)) != 1 or len(_edge_between(G, w, x)) != 1 or _edge_between(G, u, w):
            raise ValueError("need ux, wx present (once) and uw absent")
        self.G, self.v, self.u, self.w, self.x = G, v, u, w, x
        self.centres = [v]
        self.k = 2
        E = lambda p, q: _edge_between(G, p, q)[0]
        self.e = {"vu": E(v, u), "vw": E(v, w), "vx": E(v, x), "ux": E(u, x), "wx": E(w, x)}
        # orientation signs relative to LLC (v->u, v->w, v->x, u->x, w->x)
        want = {"vu": (v, u), "vw": (v, w), "vx": (v, x), "ux": (u, x), "wx": (w, x)}
        self.s = {key: (1 if G.edges[e] == want[key] else -1) for key, e in self.e.items()}
        # reduced graph
        keep = [a for a in range(G.nv) if a != v]
        self.vmap = {a: (keep.index(a) if a in keep else None) for a in range(G.nv)}
        U_, W_, X_ = self.vmap[u], self.vmap[w], self.vmap[x]
        local = set(self.e.values())
        self.others = [e for e in range(G.ne) if e not in local]
        edges = [(U_, W_), (U_, X_), (W_, X_)] + [(self.vmap[G.edges[e][0]], self.vmap[G.edges[e][1]]) for e in self.others]
        self.red = Graph(G.nv - 1, edges, "%s / Y-Delta(%d)" % (G.name or "G", v))
        # red edge ids: 0 = uw, 1 = ux, 2 = wx, 3.. = others

    def _random_fractions(self, rng, den):
        return [Fr(rng.randint(1, den - 1), den), Fr(rng.randint(1, den - 1), den)]

    def tmap_dual(self, f, T, derivs=True):
        (r1, r2), Td = self._duals(f, T, derivs)
        Tuw, Tux, Twx = Td[0], Td[1], Td[2]
        Rux = Tux / (1 - r1)
        Rwx = Twx / (1 - r2)
        tvu, tvw, tvx = delta_to_y(Tuw, Rux, Rwx)
        t = [None] * self.G.ne
        t[self.e["vu"]], t[self.e["vw"]], t[self.e["vx"]] = tvu, tvw, tvx
        t[self.e["ux"]] = Tux / r1
        t[self.e["wx"]] = Twx / r2
        for j, e in enumerate(self.others):
            t[e] = Td[3 + j]
        return t

    def forward(self, t):
        t = [Fr(x) for x in t]
        tvu, tvw, tvx = t[self.e["vu"]], t[self.e["vw"]], t[self.e["vx"]]
        Ruw, Rux, Rwx = y_to_delta(tvu, tvw, tvx)
        Tux = 1 / (1 / t[self.e["ux"]] + 1 / Rux)
        Twx = 1 / (1 / t[self.e["wx"]] + 1 / Rwx)
        r1, r2 = Tux / t[self.e["ux"]], Twx / t[self.e["wx"]]
        return [r1, r2], [Ruw, Tux, Twx] + [t[e] for e in self.others]

    def J_closed(self, f, T):
        """Closed form of J_aff: J_loc = -T_uw T_ux^2 T_wx^2 / (r1^2 r2^2 (1-r1)^3 (1-r2)^3 P^3) for
        (r1, r2, T_uw, T_ux, T_wx) -> (t_ux, t_wx, t_vu, t_vw, t_vx), times the sign of the
        permutation placing these outputs at their edge positions in G."""
        r1, r2 = f
        Tuw, Tux, Twx = T[0], T[1], T[2]
        P = Tuw + Tux / (1 - r1) + Twx / (1 - r2)
        Jloc = -Tuw * Tux ** 2 * Twx ** 2 / (r1 ** 2 * r2 ** 2 * (1 - r1) ** 3 * (1 - r2) ** 3 * P ** 3)
        # columns: f (r1, r2), T (uw, ux, wx, others...); rows: G's edge order.
        # local block rows (ux, wx, vu, vw, vx); other rows are identity on T_others.
        order_rows = [self.e["ux"], self.e["wx"], self.e["vu"], self.e["vw"], self.e["vx"]] + self.others
        # permutation sign: position of row e in G's order vs canonical order
        perm = [order_rows.index(e) for e in range(self.G.ne)]
        from .linalg import perm_sign
        return Jloc * perm_sign(perm)

    def currents(self, f, T, lam):
        """Edge currents of G from the reduced currents (reverse map, LLC b-01166 and
        local-B-sourced-current-map)."""
        t = self.tmap(f, T)
        lamr, H = self.red_sources(lam, t)
        yr = brute.currents(self.red, T, lamr)
        lamf = brute.sources(self.G, lam)
        lv = lamf[self.v]
        Gv = 1 / t[self.e["vu"]] + 1 / t[self.e["vw"]] + 1 / t[self.e["vx"]]
        r1, r2 = f
        y = [None] * self.G.ne
        for i in range(H):
            pass
        def comb(*terms):
            return [sum((c * vec[i] for c, vec in terms), Fr(0)) for i in range(H)]
        q = lambda key: [lv[i] / (t[self.e[key]] * Gv) for i in range(H)]
        y_vu = comb((1, q("vu")), (-1, yr[0]), (-(1 - r1), yr[1]))
        y_vw = comb((1, q("vw")), (1, yr[0]), (-(1 - r2), yr[2]))
        y_vx = comb((1, q("vx")), (1 - r1, yr[1]), (1 - r2, yr[2]))
        y_ux = comb((r1, yr[1]))
        y_wx = comb((r2, yr[2]))
        for key, val in (("vu", y_vu), ("vw", y_vw), ("vx", y_vx), ("ux", y_ux), ("wx", y_wx)):
            y[self.e[key]] = [self.s[key] * c for c in val]
        for j, e in enumerate(self.others):
            y[e] = yr[3 + j]
        return y

    # ---- closed formulas (REPORT.md, Section 1)
    def cycle_basis(self):
        """B' = (gamma1, gamma2, c'_1, ...): gamma1 = ux - vx + vu, gamma2 = wx - vx + vw (as closed
        paths), c' = cycles of red minus uw (= G minus v) lifted.  Returns (Cprime, Cred) with
        Cred = (delta, c_1, ...) for red, delta = uw + wx - ux."""
        G, red = self.G, self.red
        s, e = self.s, self.e
        g1 = [0] * G.ne
        g1[e["ux"]] += s["ux"]; g1[e["vx"]] -= s["vx"]; g1[e["vu"]] += s["vu"]
        g2 = [0] * G.ne
        g2[e["wx"]] += s["wx"]; g2[e["vx"]] -= s["vx"]; g2[e["vw"]] += s["vw"]
        # red minus uw: edges 1.. of red
        sub = Graph(red.nv, red.edges[1:], "red-uw")
        Csub = sub.fundamental_cycles() if sub.loops else []
        Cred = [[1, -1, 1] + [0] * (red.ne - 3)]
        Cp = [g1, g2]
        for row in Csub:
            Cred.append([0] + row)
            lift = [0] * G.ne
            lift[e["ux"]] += s["ux"] * row[0]
            lift[e["wx"]] += s["wx"] * row[1]
            for j, ee in enumerate(self.others):
                lift[ee] += row[2 + j]
            Cp.append(lift)
        return Cp, Cred

    def density_closed(self, f, T, lam, Tdir, xi_var=1):
        """D'(f, T) from the closed transformation formulas (REPORT.md Section 1):
           psi_g1 = -a dr1/(r1(1-r1)) + p2 omega,  psi_g2 = -b dr2/(r2(1-r2)) - p3 omega,
           omega = -psi^red_delta + a dr1/(1-r1) - b dr2/(1-r2),
           psi_c' = psi^red_c - C_{c,ux} a dr1/r1 - C_{c,wx} b dr2/r2,
           dipole: psi_g1[d] = alpha = dlog(T_uw/T_wx) - dr2/(1-r2), psi_g2[d] = beta = dlog(T_uw/T_ux) - dr1/(1-r1),
           a = Phi_u - Phi_x, b = Phi_w - Phi_x (reduced potentials), p2 = R_ux/P, p3 = R_wx/P,
           U = U_red t_vx/(r1 r2), Omega = U^{-H-T/2} E[prod psi^tau] ^ prod_i prod psi^i."""
        G, red = self.G, self.red
        r1, r2 = [Fr(x) for x in f]
        T = [Fr(x) for x in T]
        k, nT = 2, len(T)
        t = self.tmap(f, T)
        Tuw, Tux, Twx = T[0], T[1], T[2]
        Rux, Rwx = Tux / (1 - r1), Twx / (1 - r2)
        P = Tuw + Rux + Rwx
        p2, p3 = Rux / P, Rwx / P
        tvx = Rux * Rwx / P
        Gv = 1 / t[self.e["vu"]] + 1 / t[self.e["vw"]] + 1 / t[self.e["vx"]]
        Cp, Cred = self.cycle_basis()
        n = len(Cp)
        dr1, dr2 = {1: Fr(1)}, {2: Fr(1)}
        dlogT = [{1 << (k + j): 1 / T[j]} for j in range(nT)]
        U_red = brute.symanzik_U(red, T)
        U = U_red * tvx / (r1 * r2)
        lamr, H = self.red_sources(lam, t)
        lamf = brute.sources(G, lam)
        lv = lamf[self.v] if H else []
        alpha = form_add(form_add(dlogT[0], form_scale(dlogT[2], -1)), form_scale(dr2, -1 / (1 - r2)))
        beta = form_add(form_add(dlogT[0], form_scale(dlogT[1], -1)), form_scale(dr1, -1 / (1 - r1)))

        def forms_from_potentials(pot):
            """pot: list over red vertices of linear objects (dict index -> coefficient).
            Returns the list over B' of linear objects index -> 1-form."""
            def diff(p, q):
                return _lin_add(pot[p], _lin_scale(pot[q], -1)) or {}
            red_psi = []
            for row in Cred:
                item = {}
                for j, (p, q) in enumerate(red.edges):
                    if row[j]:
                        for kk, c in diff(p, q).items():
                            item[kk] = form_add(item.get(kk, {}), form_scale(dlogT[j], row[j] * c))
                red_psi.append(item)
            U_, W_, X_ = self.vmap[self.u], self.vmap[self.w], self.vmap[self.x]
            a, b = diff(U_, X_), diff(W_, X_)
            keys = set(a) | set(b) | set().union(*[set(it) for it in red_psi])
            out = []
            omega = {}
            for kk in keys:
                om = form_scale(red_psi[0].get(kk, {}), -1)
                om = form_add(om, form_scale(dr1, a.get(kk, Fr(0)) / (1 - r1)))
                om = form_add(om, form_scale(dr2, -b.get(kk, Fr(0)) / (1 - r2)))
                omega[kk] = om
            g1, g2 = {}, {}
            for kk in keys:
                g1[kk] = form_add(form_scale(dr1, -a.get(kk, Fr(0)) / (r1 * (1 - r1))), form_scale(omega[kk], p2))
                g2[kk] = form_add(form_scale(dr2, -b.get(kk, Fr(0)) / (r2 * (1 - r2))), form_scale(omega[kk], -p3))
            out = [g1, g2]
            for ci, row in enumerate(Cred[1:]):
                item = {}
                for kk in keys:
                    fm = red_psi[1 + ci].get(kk, {})
                    fm = form_add(fm, form_scale(dr1, -row[1] * a.get(kk, Fr(0)) / r1))
                    fm = form_add(fm, form_scale(dr2, -row[2] * b.get(kk, Fr(0)) / r2))
                    item[kk] = fm
                out.append(item)
            return [{kk: fm for kk, fm in it.items() if fm} for it in out]

        # holomorphic factors (transferred sources + dipole lambda_v/G)
        phi_red = brute.potentials(red, T, lamr) if H else None
        hol = []
        for i in range(H):
            pot = [{"c": phi_red[wv][i]} for wv in range(red.nv)]
            fl = forms_from_potentials(pot)
            psi = [it.get("c", {}) for it in fl]
            q = lv[i] / Gv
            psi[0] = form_add(psi[0], form_scale(alpha, q))
            psi[1] = form_add(psi[1], form_scale(beta, q))
            hol.extend(psi)
        # sign between B' and G's positive basis: Omega[C'] = det(A)^D Omega[C]
        Cpos = G.cycle_matrix()
        A = _basis_change(Cp, Cpos)
        dA = det(A)
        D = H + Tdir
        if Tdir == 0:
            top = wedge_all(hol)
            val = _top_density(top, k, nT, T) / U ** H
            return brute.Density(val * dA ** D)
        g = red.nv - 1
        covX = inverse(brute.laplacian(red, T))
        m = g + 1
        cov = [[Fr(0)] * m for _ in range(m)]
        for i_ in range(g):
            for j_ in range(g):
                cov[i_][j_] = covX[i_][j_] / 2
        cov[g][g] = Fr(xi_var) / (2 * Gv)
        pot = [{wv: Fr(1)} if wv < g else {} for wv in range(red.nv)]
        items = forms_from_potentials(pot)
        items[0] = dict(items[0]); items[0][g] = alpha
        items[1] = dict(items[1]); items[1][g] = beta
        pf = wick(items, [1] * n, cov)
        top = wedge(pf, wedge_all(hol))
        coeff = _top_density(top, k, nT, T) / U ** H * dA ** D
        return brute.Density(coeff, U)


    def check_components(self, f, T, lam):
        """H-independent transformation laws, compared with the pull-back of the brute-force forms
        in the basis B': (i) psi^i_r[lambda] (each component, sourced centre) from the closed
        formulas; (ii) Pf(-1/2 dM M^{-1} dM) = E_{X_red, xi}[prod psi^tau] (closed formulas);
        (iii) M_Gamma: det M_gg = M^red_dd t_vx/(r1 r2), Schur(M_Gamma; gamma) = Schur(M_red; delta).
        Returns True or raises AssertionError."""
        G, red = self.G, self.red
        T = [Fr(x) for x in T]
        tt = self.tmap_dual(f, T)
        t = [x.v for x in tt]
        k, nT = 2, len(T)
        Cp, Cred = self.cycle_basis()
        n = len(Cp)
        dt = [x.dform() for x in tt]
        # brute psi in B' (pulled back)
        lamf = brute.sources(G, lam)
        H = len(lamf[0]) if lamf else 0
        y = brute.currents(G, t, lam)
        psi_b = [[{}] * n for _ in range(H)]
        for i in range(H):
            for r in range(n):
                fm = {}
                for e in range(G.ne):
                    if Cp[r][e]:
                        fm = form_add(fm, form_scale(dt[e], Cp[r][e] * y[e][i]))
                psi_b[i][r] = fm
        # closed psi: reuse density_closed internals through a light re-implementation
        r1, r2 = [Fr(x) for x in f]
        Tuw, Tux, Twx = T[0], T[1], T[2]
        Rux, Rwx = Tux / (1 - r1), Twx / (1 - r2)
        P = Tuw + Rux + Rwx
        p2, p3 = Rux / P, Rwx / P
        Gv = 1 / t[self.e["vu"]] + 1 / t[self.e["vw"]] + 1 / t[self.e["vx"]]
        dr1, dr2 = {1: Fr(1)}, {2: Fr(1)}
        dlogT = [{1 << (k + j): 1 / T[j]} for j in range(nT)]
        alpha = form_add(form_add(dlogT[0], form_scale(dlogT[2], -1)), form_scale(dr2, -1 / (1 - r2)))
        beta = form_add(form_add(dlogT[0], form_scale(dlogT[1], -1)), form_scale(dr1, -1 / (1 - r1)))
        lamr, _ = self.red_sources(lam, t)
        phi = brute.potentials(red, T, lamr)
        U_, W_, X_ = self.vmap[self.u], self.vmap[self.w], self.vmap[self.x]
        for i in range(H):
            a = phi[U_][i] - phi[X_][i]
            b = phi[W_][i] - phi[X_][i]
            red_psi = []
            for row in Cred:
                fm = {}
                for j, (p, q) in enumerate(red.edges):
                    if row[j]:
                        fm = form_add(fm, form_scale(dlogT[j], row[j] * (phi[p][i] - phi[q][i])))
                red_psi.append(fm)
            omega = form_add(form_add(form_scale(red_psi[0], -1), form_scale(dr1, a / (1 - r1))), form_scale(dr2, -b / (1 - r2)))
            g1 = form_add(form_scale(dr1, -a / (r1 * (1 - r1))), form_scale(omega, p2))
            g2 = form_add(form_scale(dr2, -b / (r2 * (1 - r2))), form_scale(omega, -p3))
            q = lamf[self.v][i] / Gv
            g1 = form_add(g1, form_scale(alpha, q))
            g2 = form_add(g2, form_scale(beta, q))
            closed = [g1, g2]
            for ci, row in enumerate(Cred[1:]):
                closed.append(form_add(form_add(red_psi[1 + ci], form_scale(dr1, -row[1] * a / r1)), form_scale(dr2, -row[2] * b / r2)))
            for r in range(n):
                if closed[r] != psi_b[i][r]:
                    raise AssertionError("psi transformation, component %d, cycle %d" % (i, r))
        # topological Pfaffian in B'
        top_b, _ = brute.topological_form(G, t, Cp)
        # pull back the n-form top_b (in dt) to (f, T): expand each monomial
        pulled = {}
        for mask, c in top_b.items():
            es = [e for e in range(G.ne) if mask >> e & 1]
            fm = {0: Fr(c)}
            for e in es:
                fm = wedge(fm, dt[e])
            pulled = form_add(pulled, fm)
        g = red.nv - 1
        covX = inverse(brute.laplacian(red, T))
        cov = [[Fr(0)] * (g + 1) for _ in range(g + 1)]
        for i_ in range(g):
            for j_ in range(g):
                cov[i_][j_] = covX[i_][j_] / 2
        cov[g][g] = 1 / (2 * Gv)
        items = []
        # topological psi in B' via the closed formulas with potentials X (linear objects)
        def lin_diff(p, q):
            out = {}
            if p < g:
                out[p] = out.get(p, Fr(0)) + 1
            if q < g:
                out[q] = out.get(q, Fr(0)) - 1
            return out
        a_l, b_l = lin_diff(U_, X_), lin_diff(W_, X_)
        red_items = []
        for row in Cred:
            it = {}
            for j, (p, q_) in enumerate(red.edges):
                if row[j]:
                    for kk, c in lin_diff(p, q_).items():
                        it[kk] = form_add(it.get(kk, {}), form_scale(dlogT[j], row[j] * c))
            red_items.append(it)
        keys = set(a_l) | set(b_l) | set().union(*[set(it) for it in red_items])
        g1, g2 = {}, {}
        for kk in keys:
            om = form_add(form_add(form_scale(red_items[0].get(kk, {}), -1), form_scale(dr1, a_l.get(kk, 0) / (1 - r1))),
                          form_scale(dr2, -b_l.get(kk, 0) / (1 - r2)))
            g1[kk] = form_add(form_scale(dr1, -Fr(a_l.get(kk, 0)) / (r1 * (1 - r1))), form_scale(om, p2))
            g2[kk] = form_add(form_scale(dr2, -Fr(b_l.get(kk, 0)) / (r2 * (1 - r2))), form_scale(om, -p3))
        g1[g], g2[g] = alpha, beta
        items = [{kk: v for kk, v in g1.items() if v}, {kk: v for kk, v in g2.items() if v}]
        for ci, row in enumerate(Cred[1:]):
            it = {}
            for kk in keys:
                fm = form_add(form_add(red_items[1 + ci].get(kk, {}), form_scale(dr1, -row[1] * Fr(a_l.get(kk, 0)) / r1)),
                              form_scale(dr2, -row[2] * Fr(b_l.get(kk, 0)) / r2))
                if fm:
                    it[kk] = fm
            items.append(it)
        pf = wick(items, [1] * n, cov)
        if pf != pulled:
            raise AssertionError("Pfaffian (Kron decomposition) mismatch")
        # period matrix statements
        M = brute.period_matrix(G, t, Cp)
        Mr = brute.period_matrix(red, T, Cred)
        Mgg = [row[:2] for row in M[:2]]
        tvx = Rux * Rwx / P
        if det(Mgg) != Mr[0][0] * tvx / (r1 * r2):
            raise AssertionError("det M_gg")
        Mi = inverse(Mgg)
        S = [[M[2 + i][2 + j] - sum(M[2 + i][a_] * Mi[a_][b_] * M[b_][2 + j] for a_ in range(2) for b_ in range(2))
              for j in range(n - 2)] for i in range(n - 2)]
        Sr = [[Mr[1 + i][1 + j] - Mr[1 + i][0] * Mr[0][1 + j] / Mr[0][0] for j in range(n - 2)] for i in range(n - 2)]
        if S != Sr:
            raise AssertionError("Schur complements")
        return True

    def density_llc_T0(self, f, T, lam, Tdir=0):
        """T = 0, H = 2, source-free: D' = (-1)^(n-1) r1 r2 det(y^red_ux, y^red_wx) D_red(T)
        (LLC eq-kernel-star-triangle in the htk orientation)."""
        if Tdir != 0:
            raise ValueError("T = 0 only")
        t = self.tmap(f, T)
        lamr, H = self.red_sources(lam, t)
        assert H == 2
        yr = brute.currents(self.red, T, lamr)
        Dred = brute.density(self.red, T, lamr, T=0).rational()
        r1, r2 = f
        n = self.G.loops
        dety = yr[1][0] * yr[2][1] - yr[1][1] * yr[2][0]
        return brute.Density((-1) ** (n - 1) * r1 * r2 * dety * Dred)

    def xi0_T1_identity(self, f, T, lam):
        """H = 1, T = 1, source-free: the xi^0 part (Gaussian over the reduced potentials only) of
        D' equals (-1)^(n-1) (r1 r2 t_vx)^(1/2) D'_red[w_ux y_wx - y_ux w_wx], the reduced T = 1
        density with the mixed current insertion.  Returns (lhs, rhs) as Densities."""
        G, red = self.G, self.red
        r1, r2 = [Fr(x) for x in f]
        T = [Fr(x) for x in T]
        t = self.tmap(f, T)
        lamr, H = self.red_sources(lam, t)
        assert H == 1
        # lhs: closed path with the xi-variance set to zero
        lhs = self.density_closed(f, T, lam, 1, xi_var=0)
        # rhs: reduced density with insertion, computed by Wick on the reduced graph
        g = red.nv - 1
        covX = inverse(brute.laplacian(red, T))
        cov = [[covX[i][j] / 2 for j in range(g)] for i in range(g)]
        Cr = red.cycle_matrix()
        nr = len(Cr)
        nT = len(T)
        dT = [{1 << j: Fr(1)} for j in range(nT)]
        yr = brute.currents(red, T, lamr)
        # topological currents as linear objects: w_e = (X_p - X_q)/T_e
        def wlin(e):
            p, q = red.edges[e]
            out = {}
            if p < g:
                out[p] = out.get(p, Fr(0)) + 1 / T[e]
            if q < g:
                out[q] = out.get(q, Fr(0)) - 1 / T[e]
            return out
        items = []
        for row in Cr:
            item = {}
            for e in range(red.ne):
                if row[e]:
                    for kk, c in wlin(e).items():
                        item[kk] = form_add(item.get(kk, {}), form_scale(dT[e], row[e] * c))
            items.append(item)
        hol = [form_add({}, {1 << e: row[e] * yr[e][0] for e in range(red.ne) if row[e]}) for row in Cr]
        # insertion F = w_ux y_wx - y_ux w_wx  (red edges 1 = ux, 2 = wx), a degree-0 item
        F = {}
        for kk, c in wlin(1).items():
            F[kk] = F.get(kk, Fr(0)) + c * yr[2][0]
        for kk, c in wlin(2).items():
            F[kk] = F.get(kk, Fr(0)) - c * yr[1][0]
        Fitem = {kk: {0: c} for kk, c in F.items() if c != 0}
        ex = wick([Fitem] + items, [0] + [1] * nr, cov)
        top = wedge(ex, wedge_all(hol))
        Ured = brute.symanzik_U(red, T)
        full = (1 << nT) - 1
        c0 = top.get(full ^ (1 << (nT - 1)), Fr(0))
        Dred = (-c0 if (nT - 1) % 2 else c0) / T[nT - 1] / Ured
        # D'_red = Dred * U_red^{-1/2};  rhs = (-1)^(n-1) (r1 r2 tvx)^(1/2) D'_red
        Rux, Rwx = T[1] / (1 - r1), T[2] / (1 - r2)
        P = T[0] + Rux + Rwx
        tvx = Rux * Rwx / P
        n = G.loops
        # (r1 r2 tvx)^(1/2) U_red^(-1/2) = coefficient * radicand^(-1/2) with radicand U_red/(r1 r2 tvx)
        rhs = brute.Density((-1) ** (n - 1) * Dred, Ured / (r1 * r2 * tvx))
        return lhs, rhs


def reduced_insertion_density(red, T, lamr, insertion):
    """D'_red[F] for the reduced graph at T = 1, H = 1: Omega_red[F] = U_red^{-3/2}
    E_X[F prod_r psi^tau_r] ^ prod_r psi_r[lambda] = D'_red[F] iota_E' prod dT, with F a product of
    topological currents: insertion = list of (edge, holomorphic factor) meaning
    prod_j (sum over its terms) ... given as a list of linear items {X index: coefficient}."""
    g = red.nv - 1
    covX = inverse(brute.laplacian(red, T))
    cov = [[covX[i][j] / 2 for j in range(g)] for i in range(g)]
    Cr = red.cycle_matrix()
    nr, nT = len(Cr), len(T)
    yr = brute.currents(red, T, lamr)
    items = []
    for row in Cr:
        item = {}
        for e in range(red.ne):
            if row[e]:
                for kk, c in w_linear(red, T, e).items():
                    item[kk] = form_add(item.get(kk, {}), form_scale({1 << e: Fr(1)}, row[e] * c))
        items.append(item)
    hol = [{1 << e: row[e] * yr[e][0] for e in range(red.ne) if row[e]} for row in Cr]
    ins = [{kk: {0: c} for kk, c in lin.items() if c != 0} for lin in insertion]
    ex = wick(ins + items, [0] * len(ins) + [1] * nr, cov)
    top = wedge(ex, wedge_all(hol))
    Ured = brute.symanzik_U(red, T)
    full = (1 << nT) - 1
    c0 = top.get(full ^ (1 << (nT - 1)), Fr(0))
    return (-c0 if (nT - 1) % 2 else c0) / T[nT - 1] / Ured, Ured


def w_linear(red, T, e):
    """Topological current w_e = (X_p - X_q)/T_e as a linear object over the non-root vertices."""
    g = red.nv - 1
    p, q = red.edges[e]
    out = {}
    if p < g:
        out[p] = out.get(p, Fr(0)) + 1 / Fr(T[e])
    if q < g:
        out[q] = out.get(q, Fr(0)) - 1 / Fr(T[e])
    return out


def _lin_comb(*terms):
    out = {}
    for c, lin in terms:
        for k, v in lin.items():
            out[k] = out.get(k, Fr(0)) + c * v
    return {k: v for k, v in out.items() if v != 0}


def _basis_change(Cnew, Cold):
    """Integer matrix A with Cnew = A Cold (rows = cycles)."""
    n = len(Cold)
    # solve A Cold = Cnew: Cold has full row rank; pick n independent columns
    cols = []
    for e in range(len(Cold[0])):
        trial = cols + [e]
        sub = [[Cold[r][c] for c in trial] for r in range(n)]
        # rank check via det of a square selection when complete
        if len(trial) <= n:
            M = [[Fr(x) for x in row] for row in sub]
            if _rank(M) == len(trial):
                cols = trial
        if len(cols) == n:
            break
    S = [[Fr(Cold[r][c]) for c in cols] for r in range(n)]       # n x n
    Bm = [[Fr(Cnew[r][c]) for c in cols] for r in range(n)]
    # A S = Bm  ->  S^t A^t = Bm^t
    St = [list(x) for x in zip(*S)]
    Bt = [list(x) for x in zip(*Bm)]
    At = solve(St, Bt)
    A = [list(x) for x in zip(*At)]
    # verify on all columns
    for r in range(n):
        for e in range(len(Cold[0])):
            if sum(A[r][s] * Cold[s][e] for s in range(n)) != Cnew[r][e]:
                raise ArithmeticError("not a change of cycle basis")
    return A


def _rank(M):
    M = [list(r) for r in M]
    rk, rows, cols = 0, len(M), len(M[0]) if M else 0
    for c in range(cols):
        p = next((r for r in range(rk, rows) if M[r][c] != 0), None)
        if p is None:
            continue
        M[rk], M[p] = M[p], M[rk]
        for r in range(rows):
            if r != rk and M[r][c] != 0:
                fct = M[r][c] / M[rk][c]
                M[r] = [a - fct * b for a, b in zip(M[r], M[rk])]
        rk += 1
    return rk


# ---------------------------------------------------------------------- paired move

class PairedMove(_Move):
    """Paired source-free (or sourced) star-triangle move (LLC Theorem b-01173).  Reduced graph:
    vertices of G without v1, v2 (order kept, root last); edges (u->A, u->B, A->B) first, then the
    other edges of G in their order.  Fractions (r1, r2, r3, r4) as in LLC b-01180."""

    def __init__(self, G, v1, v2, u=None, A=None, B=None):
        for v in (v1, v2):
            if v == G.root:
                raise ValueError("centres must not be the root")
        n1, n2 = _neighbours(G, v1), _neighbours(G, v2)
        if len(n1) != 3 or sorted(n1) != sorted(n2) or len(set(n1)) != 3:
            raise ValueError("v1, v2 must have the same three distinct neighbours")
        if u is None:
            for c in n1:
                o = [x for x in n1 if x != c]
                if not _edge_between(G, c, o[0]) and not _edge_between(G, c, o[1]) and len(_edge_between(G, *o)) == 1:
                    u = c
                    A, B = sorted(o)
            if u is None:
                raise ValueError("no paired configuration")
        self.G, self.v1, self.v2, self.u, self.A, self.B = G, v1, v2, u, A, B
        self.centres = [v1, v2]
        self.k = 4
        E = lambda p, q: _edge_between(G, p, q)[0]
        self.e = {"1u": E(v1, u), "1A": E(v1, A), "1B": E(v1, B), "2u": E(v2, u), "2A": E(v2, A),
                  "2B": E(v2, B), "AB": E(A, B)}
        want = {"1u": (v1, u), "1A": (v1, A), "1B": (v1, B), "2u": (v2, u), "2A": (v2, A),
                "2B": (v2, B), "AB": (A, B)}
        self.s = {key: (1 if G.edges[e] == want[key] else -1) for key, e in self.e.items()}
        keep = [a for a in range(G.nv) if a not in (v1, v2)]
        self.vmap = {a: (keep.index(a) if a in keep else None) for a in range(G.nv)}
        U_, A_, B_ = self.vmap[u], self.vmap[A], self.vmap[B]
        local = set(self.e.values())
        self.others = [e for e in range(G.ne) if e not in local]
        edges = [(U_, A_), (U_, B_), (A_, B_)] + [(self.vmap[G.edges[e][0]], self.vmap[G.edges[e][1]]) for e in self.others]
        self.red = Graph(G.nv - 2, edges, "%s / paired(%d,%d)" % (G.name or "G", v1, v2))

    def _random_fractions(self, rng, den):
        r1 = Fr(rng.randint(1, den - 1), den)
        r2 = Fr(rng.randint(1, den - 1), den)
        a, b, c = rng.randint(1, 20), rng.randint(1, 20), rng.randint(1, 20)
        return [r1, r2, Fr(a, a + b + c), Fr(b, a + b + c)]

    def tmap_dual(self, f, T, derivs=True):
        (r1, r2, r3, r4), Td = self._duals(f, T, derivs)
        TuA, TuB, TAB = Td[0], Td[1], Td[2]
        R1 = (TuA / r1, TuB / r2, TAB / r3)
        R2 = (TuA / (1 - r1), TuB / (1 - r2), TAB / r4)
        t = [None] * self.G.ne
        for i, R in ((1, R1), (2, R2)):
            tu, tA, tB = R[0] * R[1], R[0] * R[2], R[1] * R[2]
            P = R[0] + R[1] + R[2]
            t[self.e["%du" % i]] = tu / P       # star edge to u: adjacent triangle edges uA, uB
            t[self.e["%dA" % i]] = tA / P       # to A: uA, AB
            t[self.e["%dB" % i]] = tB / P       # to B: uB, AB
        t[self.e["AB"]] = TAB / (1 - r3 - r4)
        for j, e in enumerate(self.others):
            t[e] = Td[3 + j]
        return t

    def forward(self, t):
        t = [Fr(x) for x in t]
        tri = []
        for i in (1, 2):
            tu, tA, tB = t[self.e["%du" % i]], t[self.e["%dA" % i]], t[self.e["%dB" % i]]
            S = tu * tA + tA * tB + tB * tu
            tri.append((S / tB, S / tA, S / tu))      # (R_uA, R_uB, R_AB)
        TuA = 1 / (1 / tri[0][0] + 1 / tri[1][0])
        TuB = 1 / (1 / tri[0][1] + 1 / tri[1][1])
        TAB = 1 / (1 / tri[0][2] + 1 / tri[1][2] + 1 / t[self.e["AB"]])
        f = [TuA / tri[0][0], TuB / tri[0][1], TAB / tri[0][2], TAB / tri[1][2]]
        return f, [TuA, TuB, TAB] + [t[e] for e in self.others]

    def currents(self, f, T, lam):
        """LLC b-01180 (+ source-transfer terms lambda_vi/(t G_i) for sourced centres)."""
        t = self.tmap(f, T)
        lamr, H = self.red_sources(lam, t)
        yr = brute.currents(self.red, T, lamr)
        lamf = brute.sources(self.G, lam)
        r1, r2, r3, r4 = f
        G1 = sum((1 / t[self.e[k_]] for k_ in ("1u", "1A", "1B")), Fr(0))
        G2 = sum((1 / t[self.e[k_]] for k_ in ("2u", "2A", "2B")), Fr(0))
        l1, l2 = lamf[self.v1], lamf[self.v2]

        def comb(*terms):
            return [sum((c * vec[i] for c, vec in terms), Fr(0)) for i in range(H)]
        q = lambda key, lv, Gi: [lv[i] / (t[self.e[key]] * Gi) for i in range(H)]
        vals = {
            "1u": comb((1, q("1u", l1, G1)), (-r1, yr[0]), (-r2, yr[1])),
            "2u": comb((1, q("2u", l2, G2)), (-(1 - r1), yr[0]), (-(1 - r2), yr[1])),
            "1A": comb((1, q("1A", l1, G1)), (r1, yr[0]), (-r3, yr[2])),
            "2A": comb((1, q("2A", l2, G2)), (1 - r1, yr[0]), (-r4, yr[2])),
            "1B": comb((1, q("1B", l1, G1)), (r2, yr[1]), (r3, yr[2])),
            "2B": comb((1, q("2B", l2, G2)), (1 - r2, yr[1]), (r4, yr[2])),
            "AB": comb((1 - r3 - r4, yr[2])),
        }
        y = [None] * self.G.ne
        for key, val in vals.items():
            y[self.e[key]] = [self.s[key] * c for c in val]
        for j, e in enumerate(self.others):
            y[e] = yr[3 + j]
        return y


    def J_abs_closed(self, f, T):
        """|J_aff| = T_uA T_uB T_AB^2 / (r1^2 (1-r1)^2 r2^2 (1-r2)^2 r3^2 r4^2 (1-r3-r4)^2)
        * prod_i R^i_uA R^i_uB R^i_AB / P_i^3."""
        r1, r2, r3, r4 = [Fr(x) for x in f]
        TuA, TuB, TAB = [Fr(x) for x in T[:3]]
        s = 1 - r3 - r4
        J = TuA * TuB * TAB ** 2 / (r1 ** 2 * (1 - r1) ** 2 * r2 ** 2 * (1 - r2) ** 2 * r3 ** 2 * r4 ** 2 * s ** 2)
        for R in ((TuA / r1, TuB / r2, TAB / r3), (TuA / (1 - r1), TuB / (1 - r2), TAB / r4)):
            J *= R[0] * R[1] * R[2] / (R[0] + R[1] + R[2]) ** 3
        return J

    def U_ratio_closed(self, f, T):
        """U_Gamma/U_red = T_uA T_uB T_AB^2 / (r1 (1-r1) r2 (1-r2) r3 r4 (1-r3-r4) P1 P2)."""
        r1, r2, r3, r4 = [Fr(x) for x in f]
        TuA, TuB, TAB = [Fr(x) for x in T[:3]]
        P1 = TuA / r1 + TuB / r2 + TAB / r3
        P2 = TuA / (1 - r1) + TuB / (1 - r2) + TAB / r4
        return TuA * TuB * TAB ** 2 / (r1 * (1 - r1) * r2 * (1 - r2) * r3 * r4 * (1 - r3 - r4) * P1 * P2)

    def density_llc_T0(self, f, T, lam, Tdir=0):
        """T = 0, H = 2, source-free: D' = -(r1 - r2)(1 - r3 - r4) det(y_uA, y_AB) det(y_uB, y_AB) D_red."""
        t = self.tmap(f, T)
        lamr, H = self.red_sources(lam, t)
        assert H == 2
        yr = brute.currents(self.red, T, lamr)
        Dred = brute.density(self.red, T, lamr, T=0).rational()
        r1, r2, r3, r4 = f
        d = lambda a, b: a[0] * b[1] - a[1] * b[0]
        return brute.Density(-(r1 - r2) * (1 - r3 - r4) * d(yr[0], yr[2]) * d(yr[1], yr[2]) * Dred)

    def xi0_T1_identity(self, f, T, lam):
        """H = 1, T = 1, source-free centres: the xi^0 part of D' (both fluctuations switched off)
        equals -(U/U_red)^{1/2} (r1 - r2)(1 - r3 - r4) D'_red[F],
        F = (w_uA y_AB - y_uA w_AB)(w_uB y_AB - y_uB w_AB).  Returns (lhs, rhs)."""
        T = [Fr(x) for x in T]
        r1, r2, r3, r4 = [Fr(x) for x in f]
        lhs = self.density_kron(f, T, lam, 1, xi_var=0)
        t = self.tmap(f, T)
        lamr, H = self.red_sources(lam, t)
        yr = brute.currents(self.red, T, lamr)
        F1 = _lin_comb((yr[2][0], w_linear(self.red, T, 0)), (-yr[0][0], w_linear(self.red, T, 2)))
        F2 = _lin_comb((yr[2][0], w_linear(self.red, T, 1)), (-yr[1][0], w_linear(self.red, T, 2)))
        Dred, Ured = reduced_insertion_density(self.red, T, lamr, [F1, F2])
        ratio = self.U_ratio(f, T)
        coeff = -(r1 - r2) * (1 - r3 - r4) * Dred
        # (U/U_red)^{1/2} U_red^{-1/2} = (U_red/ratio)^{-1/2}
        rhs = brute.Density(coeff, Ured / ratio)
        return lhs, rhs


# ---------------------------------------------------------------------- detection

def find_single(G):
    """All (v, u, w, x) with v a non-root degree-3 vertex in the LLC configuration."""
    out = []
    for v in range(G.nv - 1):
        try:
            mv = SingleMove(G, v)
        except ValueError:
            continue
        out.append((v, mv.u, mv.w, mv.x))
    return out


def find_paired(G):
    out = []
    for v1, v2 in itertools.combinations(range(G.nv - 1), 2):
        try:
            mv = PairedMove(G, v1, v2)
        except ValueError:
            continue
        out.append((v1, v2, mv.u, mv.A, mv.B))
    return out
