"""Layer d): series-parallel reduction of the network (CONVENTIONS.md, Section 5).

A reduction is a list of chain moves.  A chain move replaces k >= 1 chains between two vertices
p, q (interior vertices of degree two, sources allowed on all of them) by one edge p -> q:

    series    R_i = sum_j t_{i,j},  x_{i,j} = t_{i,j} / R_i               (simplex, sum_j x_{i,j} = 1)
    parallel  1/t_red = sum_i 1/R_i,  w_i = t_red / R_i                   (simplex, sum_i w_i = 1)
    sources   sigma_{i,k} = sum_{j>k} x_{i,j},  mhat_i = sum_k sigma_{i,k} lam_{i,k},
              lam_p += sum_i mhat_i,  lam_q += sum_{i,k} (1 - sigma_{i,k}) lam_{i,k}
    currents  y_{i,j} = w_i y_red - mhat_i + Q_{i,j},  Q_{i,j} = sum_{k<j} lam_{i,k}   (direction p -> q)

After the last move one edge is left; its resistance is the gauge t_final = 1.
All Schwinger times are Laurent monomials in the ratios and all currents are polynomials.
"""
import itertools
import random
from fractions import Fraction as Fr

import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.rings import PolyRing
try:
    from flint import fmpq, fmpq_mpoly_ctx
except ImportError:
    fmpq_mpoly_ctx = None

from . import brute
from .linalg import det as fdet, shuffle_sign
from .simplex import dirichlet, PiNumber

__all__ = ["SPReduction", "RatioDensity", "sp_reduce", "find_chains", "plan_chain_moves",
           "plan_generic", "NotSeriesParallel"]


class NotSeriesParallel(ValueError):
    pass


# ---------------------------------------------------------------------- planning

def _adjacency(cur):
    adj = {}
    for eid in sorted(cur):
        a, b = cur[eid]
        adj.setdefault(a, []).append((eid, b, +1))
        adj.setdefault(b, []).append((eid, a, -1))
    return adj


def find_chains(cur, p, q):
    """Chains from p to q in the multigraph cur = {edge id: (tail, head)} whose interior vertices
    have degree two.  Each chain is the list of its edge ids in the order from p to q."""
    adj = _adjacency(cur)
    out = []
    for (eid, other, _s) in adj.get(p, []):
        chain, prev, v, ok = [eid], eid, other, True
        while v != q:
            if v == p or len(adj[v]) != 2:
                ok = False
                break
            e2, o2, _ = next(x for x in adj[v] if x[0] != prev)
            chain.append(e2)
            prev, v = e2, o2
        if ok:
            out.append(chain)
    return out


def _apply(cur, chains, p, q, new_id):
    new = dict(cur)
    for ch in chains:
        for e in ch:
            del new[e]
    new[new_id] = (p, q)
    return new


def plan_chain_moves(G, D, prefer="small"):
    """Sequence of tight chain moves (sum_i a_i = D (k - 1) + 1, k >= 2 chains) reducing G to one
    edge, found by depth-first search; None if there is none.  prefer = 'small': try small k
    first; 'large': large k first.  Returns [(p, q, [chain, ...]), ...]; reduced edges get the ids
    G.ne, G.ne + 1, ... in the order of the moves."""
    def rec(cur, new_id):
        if len(cur) == 1:
            return []
        verts = sorted({v for ab in cur.values() for v in ab})
        for p, q in itertools.combinations(verts, 2):
            chains = find_chains(cur, p, q)
            if len(chains) < 2:
                continue
            ks = range(2, len(chains) + 1)
            for k in (ks if prefer == "small" else reversed(ks)):
                for sub in itertools.combinations(chains, k):
                    if sum(len(c) for c in sub) != D * (k - 1) + 1:
                        continue
                    r = rec(_apply(cur, sub, p, q, new_id), new_id + 1)
                    if r is not None:
                        return [(p, q, [list(c) for c in sub])] + r
        return None

    return rec({e: ab for e, ab in enumerate(G.edges)}, G.ne)


def plan_generic(G):
    """Series-parallel reduction without the tightness condition: all chains between a pair of
    vertices are merged whenever there are at least two; otherwise a chain with interior vertices
    is replaced by one edge.  Raises NotSeriesParallel if the graph does not reduce to one edge."""
    cur = {e: ab for e, ab in enumerate(G.edges)}
    new_id = G.ne
    plan = []
    while len(cur) > 1:
        verts = sorted({v for ab in cur.values() for v in ab})
        adj = _adjacency(cur)
        move = None
        for p, q in itertools.combinations(verts, 2):
            chains = find_chains(cur, p, q)
            if len(chains) >= 2:
                move = (p, q, chains)
                break
        if move is None:
            branch = [v for v in verts if len(adj[v]) != 2]
            for p, q in itertools.combinations(branch, 2):
                chains = [c for c in find_chains(cur, p, q) if len(c) >= 2]
                if chains:
                    move = (p, q, chains[:1])
                    break
        if move is None:
            raise NotSeriesParallel("no series or parallel move available")
        plan.append(move)
        cur = _apply(cur, move[2], move[0], move[1], new_id)
        new_id += 1
    return plan


# ---------------------------------------------------------------------- helpers

def _q(x):
    if isinstance(x, Fr):
        return QQ(x.numerator, x.denominator)
    if isinstance(x, int):
        return QQ(x)
    x = sp.Rational(x)
    return QQ(int(x.p), int(x.q))


def _fr(x):
    return Fr(int(x.numerator), int(x.denominator))


def _evalp(p, vals):
    """Evaluate a PolyElement at the tuple vals (Fractions) of all generators."""
    tot = Fr(0)
    for mon, c in p.terms():
        term = _fr(c)
        for v, k in zip(vals, mon):
            if k:
                term *= v ** k
        tot += term
    return tot


def _poly_det(rows, ring):
    """Polynomial determinant by column subsets for small minors and Bareiss elimination for larger ones."""
    n = len(rows)
    if n == 0:
        return ring.one
    if n == 1:
        return rows[0][0]
    if n == 2:
        return rows[0][0] * rows[1][1] - rows[0][1] * rows[1][0]
    if n >= 8:
        matrix = [list(row) for row in rows]
        previous = ring.one
        native = fmpq_mpoly_ctx is not None and ring.domain == QQ
        if native:
            context = fmpq_mpoly_ctx.get([str(symbol) for symbol in ring.symbols])
            matrix = [[context.from_dict({m: fmpq(int(c.numerator), int(c.denominator))
                                          for m, c in poly.items()}) for poly in row]
                      for row in matrix]
            previous = context.constant(1)
        def exact_divide(value, divisor):
            return value / divisor if native else value.exquo(divisor)
        sign = 1
        while len(matrix) > 1:
            size = len(matrix)
            # A short polynomial pivot limits growth of the remaining minors.
            pivot = min(((len(value), i, j) for i, row in enumerate(matrix)
                         for j, value in enumerate(row) if value), default=None)
            if pivot is None:
                return ring.zero
            _, i, j = pivot
            if i:
                matrix[0], matrix[i] = matrix[i], matrix[0]
                sign = -sign
            if j:
                for row in matrix:
                    row[0], row[j] = row[j], row[0]
                sign = -sign
            top = matrix[0]
            matrix = [[exact_divide(top[0]*row[k] - row[0]*top[k], previous)
                       for k in range(1, size)] for row in matrix[1:]]
            previous = top[0]
        answer = sign * matrix[0][0]
        if native:
            return ring.from_dict({tuple(map(int, m)): QQ(int(c.numerator), int(c.denominator))
                                   for m, c in answer.terms()})
        return answer
    order = sorted(range(n), key=lambda i: sum(1 for x in rows[i] if x != 0))
    # sign of the row permutation
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
                if m2 in new:
                    new[m2] = new[m2] + term
                else:
                    new[m2] = term
        dp = {m: v for m, v in new.items() if v != 0}
        if not dp:
            return ring.zero
    res = dp.get((1 << n) - 1, ring.zero)
    return res if sgn > 0 else -res


class Block:
    """A simplex of barycentric ratios bary[0], ..., bary[m] (sum = 1); bary[0] is eliminated."""
    __slots__ = ("kind", "move", "chain", "bary")

    def __init__(self, kind, move, chain, bary):
        self.kind, self.move, self.chain, self.bary = kind, move, chain, list(bary)

    @property
    def indep(self):
        return self.bary[1:]

    def __repr__(self):
        return "Block(%s, move %d, %s)" % (self.kind, self.move, self.bary)


class Move:
    __slots__ = ("index", "p", "q", "chains", "signs", "interiors", "new_edge", "wblock", "xblocks",
                 "mhat", "Q", "yred")

    def __repr__(self):
        return "Move(%d: %d->%d, chains %s -> edge %d)" % (self.index, self.p, self.q, self.chains,
                                                            self.new_edge)


class LF:
    """num * prod_b b^{exp_b}: polynomial in the independent ratios times a Laurent monomial
    (rational exponents) in all barycentric ratios."""
    __slots__ = ("num", "exp")

    def __init__(self, num, exp):
        self.num, self.exp = num, tuple(exp)


# ---------------------------------------------------------------------- the reduction

class SPReduction:
    """Series-parallel reduction of G.

    Attributes
      G, D, n         graph, dimension D = H + T used for the plan, loop number
      moves           list of Move
      blocks          list of Block (simplices); the integration domain is their product
      bary            all barycentric ratio names;  coords: the independent ones, in order
      ring            QQ[coords]
      texp[e]         {ratio name: integer exponent}:  t_e = prod ratio^exponent  (gauge t_final = 1)
      resp[e][v]      PolyElement, y_e = sum_{v != root} resp[e][v] lambda_v
      orientation     sign of det[t ; dt/dcoords] (ratio coordinates versus iota_E prod dt_e)
    """

    def __init__(self, G, D=None, plan=None, prefer="small"):
        self.G = G
        self.n = G.loops
        if D is None:
            if self.n and (G.ne - 1) % self.n == 0:
                D = (G.ne - 1) // self.n
        self.D = D
        self.tight_plan = True
        if plan is None:
            plan = plan_chain_moves(G, D, prefer) if D is not None else None
            if plan is None:
                plan = plan_generic(G)
                self.tight_plan = False
        self._build(plan)

    # -------------------------------------------------------------- construction
    def _build(self, plan):
        G = self.G
        cur = {e: ab for e, ab in enumerate(G.edges)}
        new_id = G.ne
        self.moves, self.blocks = [], []
        for mi, (p, q, chains) in enumerate(plan):
            adj = _adjacency(cur)
            mv = Move()
            mv.index, mv.p, mv.q, mv.new_edge = mi, p, q, new_id
            mv.chains, mv.signs, mv.interiors = [], [], []
            used = set()
            for ch in chains:
                v, signs, ints = p, [], []
                for pos, e in enumerate(ch):
                    if e not in cur or e in used:
                        raise ValueError("move %d: edge %s not available" % (mi, e))
                    used.add(e)
                    a, b = cur[e]
                    if a == v:
                        signs.append(+1)
                        v = b
                    elif b == v:
                        signs.append(-1)
                        v = a
                    else:
                        raise ValueError("move %d: edge %s does not continue the chain" % (mi, e))
                    if pos < len(ch) - 1:
                        if len(adj[v]) != 2 or v in (p, q):
                            raise ValueError("move %d: interior vertex %d must have degree two" % (mi, v))
                        ints.append(v)
                if v != q:
                    raise ValueError("move %d: chain %s does not end at %d" % (mi, ch, q))
                mv.chains.append(list(ch))
                mv.signs.append(signs)
                mv.interiors.append(ints)
            k = len(mv.chains)
            mv.wblock = None
            if k >= 2:
                mv.wblock = Block("w", mi, None, ["w%d_%d" % (mi + 1, i + 1) for i in range(k)])
                self.blocks.append(mv.wblock)
            mv.xblocks = []
            for i, ch in enumerate(mv.chains):
                if len(ch) >= 2:
                    blk = Block("x", mi, i, ["x%d_%d_%d" % (mi + 1, i + 1, j + 1) for j in range(len(ch))])
                    self.blocks.append(blk)
                    mv.xblocks.append(blk)
                else:
                    mv.xblocks.append(None)
            self.moves.append(mv)
            cur = _apply(cur, mv.chains, p, q, new_id)
            new_id += 1
        if len(cur) != 1:
            raise ValueError("the plan leaves %d edges" % len(cur))
        (self.final_edge, (pf, qf)), = cur.items()
        self.final = (pf, qf)

        self.bary = [b for blk in self.blocks for b in blk.bary]
        self.coords = [b for blk in self.blocks for b in blk.indep]
        self._bidx = {b: i for i, b in enumerate(self.bary)}
        self._cidx = {c: i for i, c in enumerate(self.coords)}
        if len(self.coords) != G.ne - 1:
            raise AssertionError("number of ratios %d != |I| - 1" % len(self.coords))
        self.ring = PolyRing(self.coords, QQ) if self.coords else PolyRing(["_dummy"], QQ)
        self._bpoly, self._dep = self._bary_polys(self.ring)
        # d b / d coord_j : list of (bary index, +-1) for each coordinate j
        self._dcoord = []
        for c in self.coords:
            blk = next(b for b in self.blocks if c in b.bary)
            self._dcoord.append([(self._bidx[c], 1), (self._bidx[blk.bary[0]], -1)])

        # resistances, backwards
        self.texp = {self.final_edge: {}}
        for mv in reversed(self.moves):
            Rnew = self.texp[mv.new_edge]
            for i, ch in enumerate(mv.chains):
                Ri = dict(Rnew)
                if mv.wblock is not None:
                    b = mv.wblock.bary[i]
                    Ri[b] = Ri.get(b, 0) - 1
                for j, e in enumerate(ch):
                    te = dict(Ri)
                    if mv.xblocks[i] is not None:
                        b = mv.xblocks[i].bary[j]
                        te[b] = te.get(b, 0) + 1
                    self.texp[e] = te

        self.resp = self._reconstruct(self.ring, self._bpoly)
        pt = self.random_point(random.Random(12345))
        J = self.jacobian_t(pt)
        if J == 0:
            raise AssertionError("degenerate ratio coordinates")
        self.orientation = 1 if J > 0 else -1

    def _bary_polys(self, ring):
        gens = dict(zip([str(g) for g in ring.gens], ring.gens))
        bpoly, dep = {}, {}
        for blk in self.blocks:
            s = ring.one
            for b in blk.indep:
                bpoly[b] = gens[b]
                s = s - gens[b]
                dep[b] = blk.bary[0]
            bpoly[blk.bary[0]] = s
        return bpoly, dep

    def _reconstruct(self, ring, bpoly):
        """Forward source transfer and backward reconstruction of the currents.
        Vectors are lists over the non-root vertices of the original graph."""
        G = self.G
        g = G.nv - 1
        zero = [ring.zero] * g

        def add(a, b):
            return [x + y for x, y in zip(a, b)]

        def scal(c, a):
            return [c * x for x in a]

        S = {}
        for v in range(g):
            S[v] = [ring.one if u == v else ring.zero for u in range(g)]
        S[G.root] = [-ring.one] * g
        for mv in self.moves:
            mv.mhat, mv.Q = [], []
            addp, addq = zero, zero
            for i, ch in enumerate(mv.chains):
                a = len(ch)
                xs = [bpoly[b] for b in mv.xblocks[i].bary] if a >= 2 else [ring.one]
                mh, Qs, Qacc = zero, [zero], zero
                for kk, u in enumerate(mv.interiors[i]):       # vertex between edges kk and kk+1 (0-based)
                    lam_u = S.pop(u)
                    sigma = ring.zero
                    for j in range(kk + 1, a):
                        sigma = sigma + xs[j]
                    mh = add(mh, scal(sigma, lam_u))
                    addq = add(addq, scal(ring.one - sigma, lam_u))
                    Qacc = add(Qacc, lam_u)
                    Qs.append(Qacc)
                mv.mhat.append(mh)
                mv.Q.append(Qs)
                addp = add(addp, mh)
            S[mv.p] = add(S[mv.p], addp)
            S[mv.q] = add(S[mv.q], addq)
        pf, qf = self.final
        Y = {self.final_edge: S[pf]}
        for mv in reversed(self.moves):
            ynew = Y[mv.new_edge]
            mv.yred = ynew
            for i, ch in enumerate(mv.chains):
                Yi = scal(bpoly[mv.wblock.bary[i]], ynew) if mv.wblock is not None else ynew
                for j, e in enumerate(ch):
                    ydir = add(add(Yi, scal(-ring.one, mv.mhat[i])), mv.Q[i][j])
                    Y[e] = ydir if mv.signs[i][j] > 0 else scal(-ring.one, ydir)
        return [Y[e] for e in range(G.ne)]

    # -------------------------------------------------------------- description
    def describe(self):
        lines = ["%s, D = %s, %d moves, coordinates %s" % (self.G.name or "graph", self.D,
                                                             len(self.moves), self.coords)]
        for mv in self.moves:
            lines.append("  move %d: %d -> %d, chains %s (signs %s, interior %s) -> edge %d" % (
                mv.index + 1, mv.p, mv.q, mv.chains, mv.signs, mv.interiors, mv.new_edge))
        lines.append("  final edge %d = (%d -> %d), t = 1; orientation %+d" % (
            self.final_edge, self.final[0], self.final[1], self.orientation))
        return "\n".join(lines)

    # -------------------------------------------------------------- points and the map t(ratios)
    def random_point(self, rng=None, den=None):
        """Random rational point of the domain: {ratio name: Fraction} for all barycentric ratios."""
        rng = rng or random
        pt = {}
        for blk in self.blocks:
            ns = [rng.randint(1, 30) for _ in blk.bary]
            tot = sum(ns)
            for b, k in zip(blk.bary, ns):
                pt[b] = Fr(k, tot)
        return pt

    def complete_point(self, pt):
        """Add the eliminated ratios b_0 = 1 - sum of the others to a dict of independent ratios."""
        pt = {k: Fr(v) for k, v in pt.items()}
        for blk in self.blocks:
            pt[blk.bary[0]] = 1 - sum(pt[b] for b in blk.indep)
            if any(pt[b] <= 0 for b in blk.bary):
                raise ValueError("point outside the domain")
        return pt

    def _vals(self, pt):
        return tuple(Fr(pt[c]) for c in self.coords)

    def tmap(self, pt):
        """Schwinger times t_e (Fractions, gauge t_final = 1) at the point pt."""
        out = []
        for e in range(self.G.ne):
            v = Fr(1)
            for b, k in self.texp[e].items():
                v *= Fr(pt[b]) ** k
            out.append(v)
        return out

    def t_expr(self):
        """t_e as sympy expressions in the barycentric ratios (symbols named as in self.bary)."""
        return [sp.Mul(*[sp.Symbol(b) ** k for b, k in self.texp[e].items()]) for e in range(self.G.ne)]

    def ratios_of_t(self, t):
        """Inverse map: all barycentric ratios at given Schwinger times (any scale).
        t: Fractions/ints, or sympy expressions (then the ratios are sympy expressions)."""
        R = {e: (x if isinstance(x, sp.Basic) else Fr(x)) for e, x in enumerate(t)}
        pt = {}
        for mv in self.moves:
            Rs = []
            for i, ch in enumerate(mv.chains):
                Ri = sum(R[e] for e in ch)
                if mv.xblocks[i] is not None:
                    for b, e in zip(mv.xblocks[i].bary, ch):
                        pt[b] = R[e] / Ri
                Rs.append(Ri)
            tred = 1 / sum(1 / r for r in Rs)
            if mv.wblock is not None:
                for b, r in zip(mv.wblock.bary, Rs):
                    pt[b] = tred / r
            R[mv.new_edge] = tred
        return pt

    def jacobian_t(self, pt):
        """det[t ; dt/dxi_1 ; ... ; dt/dxi_{N-1}] (rows) at pt: the value of iota_E prod dt_e on
        the coordinate vector fields of the ratios."""
        t = self.tmap(pt)
        rows = [t]
        for j, c in enumerate(self.coords):
            row = []
            for e in range(self.G.ne):
                s = Fr(0)
                for (bi, sg) in self._dcoord[j]:
                    b = self.bary[bi]
                    k = self.texp[e].get(b, 0)
                    if k:
                        s += sg * k / Fr(pt[b])
                row.append(t[e] * s)
            rows.append(row)
        return fdet(rows)

    # -------------------------------------------------------------- rings with parameters
    def _setup(self, objs):
        """Ring containing the coordinates and all free symbols of objs; returns
        (ring, bary polynomials, response matrix, converter)."""
        syms = set()
        for o in objs:
            if isinstance(o, sp.Basic):
                syms |= o.free_symbols
        if not syms:
            ring, bpoly, resp = self.ring, self._bpoly, self.resp
        else:
            names = self.coords + sorted(str(s) for s in syms)
            if len(set(names)) != len(names):
                raise ValueError("parameter names clash with the ratio names")
            ring = PolyRing(names, QQ)
            bpoly, _ = self._bary_polys(ring)
            resp = [[p.set_ring(ring) for p in row] for row in self.resp]

        def conv(x):
            if isinstance(x, sp.Basic):
                if x.free_symbols:
                    return ring.from_expr(x)
                return ring.ground_new(_q(x))
            return ring.ground_new(_q(x))

        return ring, bpoly, resp, conv

    @staticmethod
    def _flat(lam):
        if lam is None:
            return []
        if isinstance(lam, dict):
            return [x for v in lam.values() for x in v]
        return [x for v in lam for x in v]

    def _lam(self, lam, conv, H=None):
        g = self.G.nv - 1
        if lam is None:
            return [[] for _ in range(g)], 0
        if isinstance(lam, dict):
            if H is None:
                H = len(next(iter(lam.values())))
            rows = [list(lam.get(v, [0] * H)) for v in range(g)]
        else:
            rows = [list(v) for v in lam]
            if len(rows) == self.G.nv:
                Hh = len(rows[0])
                for i in range(Hh):
                    tot = sum(sp.sympify(r[i]) for r in rows)
                    if sp.simplify(tot) != 0:
                        raise ValueError("sources do not sum to zero")
                rows = rows[:g]          # the root source is minus the sum of the others
            if len(rows) != g:
                raise ValueError("need sources at the %d non-root vertices" % g)
        H = len(rows[0]) if rows else 0
        return [[conv(x) for x in r] for r in rows], H

    # -------------------------------------------------------------- currents, displacement map
    def currents(self, lam):
        """Polynomial currents y[e][i] (PolyElements in the ratios) for sources lam.
        Entries of lam may be numbers or sympy expressions in parameters."""
        ring, _bp, resp, conv = self._setup(self._flat(lam))
        lamr, H = self._lam(lam, conv)
        return self._currents(ring, resp, lamr, H)

    def _currents(self, ring, resp, lamr, H):
        g = self.G.nv - 1
        out = []
        for e in range(self.G.ne):
            row = []
            for i in range(H):
                s = ring.zero
                for v in range(g):
                    if lamr[v][i] != 0 and resp[e][v] != 0:
                        s = s + resp[e][v] * lamr[v][i]
                row.append(s)
            out.append(row)
        return out

    def currents_at(self, lam, pt):
        """Currents at the point pt as Fractions (lam numeric)."""
        vals = self._vals(pt)
        lam = brute.sources(self.G, lam)
        H = len(lam[0]) if lam else 0
        c = [[_evalp(p, vals) for p in row] for row in self.resp]
        return [[sum((c[e][v] * lam[v][i] for v in range(self.G.nv - 1)), Fr(0)) for i in range(H)]
                for e in range(self.G.ne)]

    def displacement(self):
        """X[v][e] (PolyElements) with  - sum_e y_e . z_e = sum_{v != root} lambda_v . X_v,
        X_v = sum_e X[v][e] z_e,  i.e. X[v][e] = - d y_e / d lambda_v."""
        g = self.G.nv - 1
        return [[-self.resp[e][v] for e in range(self.G.ne)] for v in range(g)]

    def displacement_expr(self):
        """The same as sympy expressions in the independent ratios."""
        return [[p.as_expr() for p in row] for row in self.displacement()]

    # -------------------------------------------------------------- densities
    def _sign_dy(self, H, T):
        """K integrand = sign * rho^T ^ wedge_r d^H y_{c_r} for the fundamental cycles of the
        Kruskal tree:  sigma_T^{H+T} (-1)^{nH + C(n,2) C(H,2)}."""
        n = self.n
        s = (-1) ** (n * H + (n * (n - 1) // 2) * (H * (H - 1) // 2))
        return s * self.G.sigma() ** (H + T)

    def _jac(self, ring, resp, lamr, H):
        """Rows (r, i): d y^i_{c_r} / d xi_j as PolyElements, c_r = cotree edges of the Kruskal tree."""
        y = self._currents(ring, resp, lamr, H)
        gens = ring.gens[:len(self.coords)]
        return [[y[c][i].diff(x) for x in gens] for c in self.G.cotree() for i in range(H)]

    def density(self, lam, T=0, method="auto", check=2, seed=0):
        """Density g of the reduced kernel in the ratio coordinates at z = 0:

            K(lam; 0) = int_{prod simplices} g(xi) d xi,     g = prod_b b^{exp_b} * polynomial.

        T = 0: exp = 0.  T = 1: half-integer exponents.  method: 'auto', 'local' (product of the
        local determinants of tight chain moves, T = 0) or 'generic' (Jacobian of the polynomial
        currents; T = 1: Laplace expansion with rho_n).  check: number of random rational points at
        which g is compared with the brute-force layer (0: no check)."""
        ring, bpoly, resp, conv = self._setup(self._flat(lam))
        lamr, H = self._lam(lam, conv)
        n, N = self.n, self.G.ne
        if n * (H + T) != N - 1:
            raise ValueError("degree mismatch: n (H + T) = %d, |I| - 1 = %d" % (n * (H + T), N - 1))
        zero_exp = tuple(Fr(0) for _ in self.bary)
        if T == 0:
            if method == "auto":
                method = "local" if self._local_ok(H) else "generic"
            if method == "local":
                dens = self._density_local(ring, bpoly, resp, lamr, H)
            else:
                J = self._jac(ring, resp, lamr, H)
                num = _poly_det(J, ring) * (self.orientation * self._sign_dy(H, 0))
                dens = RatioDensity(self, ring, zero_exp, [num], H, 0)
        elif T == 1:
            dens = self._density_T1(ring, bpoly, resp, lamr, H)
        else:
            raise ValueError("T must be 0 or 1")
        dens.method = method
        if check:
            self.selfcheck(lam, T, npoints=check, seed=seed, dens=dens)
        return dens

    def _local_ok(self, H):
        return all(mv.wblock is not None and sum(len(c) for c in mv.chains) == H * (len(mv.chains) - 1) + 1
                   for mv in self.moves)

    def _local_matrix(self, mv, ring, bpoly, lamr, H):
        """Local matrix of a tight chain move at T = 0.  Rows (i, a), i = 2..k, a = 1..H; columns
        w_2..w_k, then x_{j,l} (j = 1..k, l = 2..a_j):
            d(l_i - w_i sum_j l_j) = y_red dw_i + sum_j (w_i - delta_ij) Q_{j,l} dx_{j,l}."""
        g = self.G.nv - 1
        k = len(mv.chains)

        def vec(resp_row):
            out = []
            for a in range(H):
                s = ring.zero
                for v in range(g):
                    if lamr[v][a] != 0 and resp_row[v] != 0:
                        s = s + resp_row[v] * lamr[v][a]
                out.append(s)
            return out

        conv = (lambda p: p) if ring is self.ring else (lambda p: p.set_ring(ring))
        yred = vec([conv(p) for p in mv.yred])
        w = [bpoly[b] for b in mv.wblock.bary]
        Q = [[vec([conv(p) for p in Qjl]) for Qjl in Qj] for Qj in mv.Q]
        rows = []
        for i in range(1, k):
            for a in range(H):
                row = [yred[a] if ii == i else ring.zero for ii in range(1, k)]
                for j in range(k):
                    for l in range(1, len(mv.chains[j])):
                        row.append((w[i] - (1 if i == j else 0)) * Q[j][l][a])
                rows.append(row)
        return rows

    def _density_local(self, ring, bpoly, resp, lamr, H):
        if not self._local_ok(H):
            raise ValueError("the plan is not a sequence of tight chain moves for H = %d" % H)
        factors = [_poly_det(self._local_matrix(mv, ring, bpoly, lamr, H), ring) for mv in self.moves]
        # overall sign: compare with the Jacobian of the polynomial currents at rational points
        rng = random.Random(777)
        sign = None
        npar = len(ring.gens) - len(self.coords)
        J = self._jac(ring, resp, lamr, H)
        for _ in range(20):
            pt = self.random_point(rng)
            vals = self._vals(pt) + tuple(Fr(rng.randint(-9, 9)) for _ in range(npar))
            prod = Fr(1)
            for f in factors:
                prod *= _evalp(f, vals)
            ref = fdet([[_evalp(p, vals) for p in row] for row in J]) * self.orientation * self._sign_dy(H, 0)
            if prod == 0 and ref == 0:
                continue
            if prod == ref:
                sign = 1
            elif prod == -ref:
                sign = -1
            else:
                raise AssertionError("local determinants do not reproduce the Jacobian of the currents")
            break
        if sign is None:
            sign = 1          # the density vanishes identically at the points tried
        factors = [ring.ground_new(QQ(sign))] + factors
        return RatioDensity(self, ring, tuple(Fr(0) for _ in self.bary), factors, H, 0)

    # ---- Laurent-fraction algebra (polynomial times Laurent monomial in the barycentric ratios)
    def _lf_ops(self, ring, bpoly):
        red = self
        nb = len(self.bary)
        P = [bpoly[b] for b in self.bary]

        def mul(a, b):
            return LF(a.num * b.num, [x + y for x, y in zip(a.exp, b.exp)])

        def scale(a, c):
            return LF(a.num * _q(c), a.exp)

        def add(a, b):
            if a.num == 0:
                return b
            if b.num == 0:
                return a
            e = [min(x, y) for x, y in zip(a.exp, b.exp)]
            na, nb_ = a.num, b.num
            for i in range(nb):
                da, db = a.exp[i] - e[i], b.exp[i] - e[i]
                if da.denominator != 1 or db.denominator != 1:
                    raise ArithmeticError("incompatible exponents")
                if da:
                    na = na * P[i] ** int(da)
                if db:
                    nb_ = nb_ * P[i] ** int(db)
            return LF(na + nb_, e)

        def diff(a, j):
            out = LF(a.num.diff(ring.gens[j]), a.exp)
            for (bi, sg) in red._dcoord[j]:
                if a.exp[bi] != 0:
                    e = list(a.exp)
                    e[bi] -= 1
                    out = add(out, LF(a.num * _q(a.exp[bi] * sg), e))
            return out

        def normalize(a):
            if a.num == 0:
                return LF(ring.zero, [Fr(0)] * nb)
            num, e = a.num, list(a.exp)
            for i in range(nb):
                while True:
                    qq, rr = num.div(P[i])
                    if rr != 0:
                        break
                    num = qq
                    e[i] += 1
            return LF(num, e)

        return mul, scale, add, diff, normalize

    def _rho(self, ring, bpoly):
        """rho_n(M) in the ratio coordinates for the fundamental cycles of the Kruskal tree:
        {sorted tuple of coordinate indices: LF}."""
        mul, scale, add, diff, normalize = self._lf_ops(ring, bpoly)
        G, n = self.G, self.n
        nb = len(self.bary)
        zero = LF(ring.zero, [Fr(0)] * nb)
        one = LF(ring.one, [Fr(0)] * nb)
        C = G.fundamental_cycles()
        tlf = [LF(ring.one, [Fr(self.texp[e].get(b, 0)) for b in self.bary]) for e in range(G.ne)]
        M = [[zero] * n for _ in range(n)]
        for r in range(n):
            for s in range(n):
                acc = zero
                for e in range(G.ne):
                    c = C[r][e] * C[s][e]
                    if c:
                        acc = add(acc, scale(tlf[e], c))
                M[r][s] = acc

        def minor_det(rows, cols):
            if not rows:
                return one
            r0 = rows[0]
            acc = zero
            for k, c in enumerate(cols):
                if M[r0][c].num == 0:
                    continue
                term = mul(M[r0][c], minor_det(rows[1:], cols[:k] + cols[k + 1:]))
                acc = add(acc, scale(term, -1 if k % 2 else 1))
            return acc

        detM = normalize(minor_det(list(range(n)), list(range(n))))
        if not detM.num.is_ground:
            raise ArithmeticError("det M is not a monomial in the ratios")
        c = _fr(detM.num.LC) if detM.num != 0 else Fr(0)
        if c != 1:
            raise ArithmeticError("det M has coefficient %s" % c)
        self._U_exp = dict(zip(self.bary, detM.exp))
        inv_det = LF(ring.one, [-x for x in detM.exp])
        Minv = [[zero] * n for _ in range(n)]
        idx = list(range(n))
        for a in range(n):
            for b in range(n):
                cof = minor_det(idx[:b] + idx[b + 1:], idx[:a] + idx[a + 1:])     # cofactor (b, a)
                Minv[a][b] = scale(mul(cof, inv_det), -1 if (a + b) % 2 else 1)
        nc = len(self.coords)
        dM = [[[diff(M[r][a], j) for j in range(nc)] for a in range(n)] for r in range(n)]
        mcoords = [j for j in range(nc) if any(dM[r][a][j].num != 0 for r in range(n) for a in range(n))]
        A = [[{} for _ in range(n)] for _ in range(n)]
        for r in range(n):
            for s in range(r + 1, n):
                form = {}
                for x, j in enumerate(mcoords):
                    for l in mcoords[x + 1:]:
                        acc = zero
                        for a in range(n):
                            for b in range(n):
                                if Minv[a][b].num == 0:
                                    continue
                                t1 = mul(dM[r][a][j], dM[s][b][l])
                                t2 = mul(dM[r][a][l], dM[s][b][j])
                                acc = add(acc, mul(Minv[a][b], add(t1, scale(t2, -1))))
                        if acc.num != 0:
                            form[(j, l)] = scale(acc, Fr(-1, 2))
                A[r][s] = form
                A[s][r] = {key: scale(v, -1) for key, v in form.items()}

        def fwedge(f, h):
            out = {}
            for ka, va in f.items():
                for kb, vb in h.items():
                    if set(ka) & set(kb):
                        continue
                    key = tuple(sorted(ka + kb))
                    inv = sum(1 for x in ka for y in kb if x > y)
                    term = scale(mul(va, vb), -1 if inv % 2 else 1)
                    out[key] = add(out[key], term) if key in out else term
            return {k_: v for k_, v in out.items() if v.num != 0}

        def pf(ix):
            if not ix:
                return {(): one}
            if len(ix) % 2:
                return {}
            out = {}
            for j in range(1, len(ix)):
                term = fwedge(A[ix[0]][ix[j]], pf(ix[1:j] + ix[j + 1:]))
                for key, v in term.items():
                    v = scale(v, -1 if (j - 1) % 2 else 1)
                    out[key] = add(out[key], v) if key in out else v
            return {k_: v for k_, v in out.items() if v.num != 0}

        Uhalf = LF(ring.one, [-x / 2 for x in detM.exp])
        return {key: mul(Uhalf, v) for key, v in pf(list(range(n))).items()}

    def _density_T1(self, ring, bpoly, resp, lamr, H):
        mul, scale, add, diff, normalize = self._lf_ops(ring, bpoly)
        nb, nc, n = len(self.bary), len(self.coords), self.n
        rho = self._rho(ring, bpoly)
        J = self._jac(ring, resp, lamr, H)
        total = LF(ring.zero, [Fr(0)] * nb)
        for S, coef in rho.items():
            rest = [j for j in range(nc) if j not in S]
            minor = _poly_det([[row[j] for j in rest] for row in J], ring) if rest else ring.one
            if minor == 0:
                continue
            term = LF(coef.num * minor, coef.exp)
            total = add(total, scale(term, shuffle_sign(S, nc)))
        total = normalize(scale(total, self.orientation * self._sign_dy(H, 1)))
        return RatioDensity(self, ring, total.exp, [total.num], H, 1)

    # -------------------------------------------------------------- self-check against layer b)
    def selfcheck(self, lam, T=0, npoints=3, seed=0, dens=None, kernel_check=True):
        """Pointwise comparison with the brute-force layer at random rational points:
        (i) t(xi) reproduces xi under ratios_of_t; (ii) polynomial currents = Kirchhoff currents;
        (iii) g(xi) = D(t(xi)) |det[t ; dt/dxi]|, sign included.  Symbolic parameters in lam are
        replaced by random integers.  Returns True or raises AssertionError."""
        rng = random.Random(seed)
        flat = self._flat(lam)
        syms = sorted({s for o in flat if isinstance(o, sp.Basic) for s in o.free_symbols}, key=str)
        for _ in range(npoints):
            sub = {s: rng.randint(-7, 7) for s in syms}

            def num(x):
                if isinstance(x, sp.Basic):
                    x = sp.Rational(x.subs(sub))
                    return Fr(int(x.p), int(x.q))
                return Fr(x)

            if lam is None:
                lamn = None
            elif isinstance(lam, dict):
                lamn = {v: [num(x) for x in vec] for v, vec in lam.items()}
            else:
                lamn = [[num(x) for x in vec] for vec in lam]
                if len(lamn) == self.G.nv:
                    lamn = lamn[:-1]
            pt = self.random_point(rng)
            t = self.tmap(pt)
            back = self.ratios_of_t(t)
            if any(back[b] != pt[b] for b in self.bary):
                raise AssertionError("ratios_of_t(tmap(xi)) != xi")
            yb = brute.currents(self.G, t, lamn)
            if lamn is not None and yb != self.currents_at(lamn, pt):
                raise AssertionError("polynomial currents differ from the Kirchhoff solution")
            if dens is not None:
                Dt = brute.density(self.G, t, lamn, T=T)
                J = self.jacobian_t(pt)
                if (J > 0) != (self.orientation > 0):
                    raise AssertionError("orientation of the ratio coordinates is not constant")
                lhs = Dt.scaled(abs(J))
                vals = self._vals(pt) + tuple(Fr(sub[s]) for s in syms)
                rat, rad = dens.evaluate(vals)
                rhs = brute.Density(rat, 1 / rad if rad != 0 else 1)
                if not lhs.same_as(rhs):
                    raise AssertionError("density mismatch at %s: brute %r, reduced %r" % (pt, lhs, rhs))
        return True

    # -------------------------------------------------------------- integration
    def integrate(self, lam, T=0, powers=None, insertion=None, method="auto", check=2, dens=None):
        """int g(xi) F(y(xi)) d xi over the product of simplices, exactly.
        powers: {(e, i): p} inserts prod (y^i_e)^p (i 0-based); insertion: callable taking the
        polynomial currents y[e][i] and returning a PolyElement.  Returns a Fraction (T = 0, numeric
        lam), a PiNumber (T = 1), or a sympy expression if lam contains parameters."""
        if dens is None:
            dens = self.density(lam, T, method=method, check=check)
        F = None
        if powers or insertion is not None:
            ring, _bp, resp, conv = self._setup(self._flat(lam))
            lamr, H = self._lam(lam, conv)
            y = self._currents(ring, resp, lamr, H)
            F = ring.one
            if powers:
                for (e, i), p in powers.items():
                    F = F * y[e][i] ** p
            if insertion is not None:
                F = F * insertion(y)
        return dens.integrate(F)

    def kernel(self, lam, T=0, **kw):
        """K(lam; 0)."""
        return self.integrate(lam, T, **kw)

    def taylor(self, lam, z, order, T=0, dens=None, check=2):
        """[K_0, ..., K_order] with K(lam; z) = sum_m K_m, K_m = int g (-sum_e y_e.z_e)^m / m!.
        z: list over edges of H-vectors, or dict {edge: vector}; entries numbers or sympy expressions."""
        zflat = [x for v in (z.values() if isinstance(z, dict) else z) for x in v]
        ring, _bp, resp, conv = self._setup(self._flat(lam) + zflat)
        lamr, H = self._lam(lam, conv)
        if dens is None:
            dens = self.density(lam, T, check=check)
        if dens.ring is not ring:
            dens = dens.in_ring(ring)
        y = self._currents(ring, resp, lamr, H)
        E = ring.zero
        items = z.items() if isinstance(z, dict) else enumerate(z)
        for e, ze in items:
            for i, zz in enumerate(ze):
                E = E - y[e][i] * conv(zz)
        out, F, fact = [], ring.one, 1
        for m in range(order + 1):
            if m:
                F = F * E
                fact *= m
            val = dens.integrate(F)
            out.append(val / fact)
        return out


class RatioDensity:
    """g(xi) = prod_b b^{exps[b]} * prod(factors): density of K in the ratio coordinates."""

    def __init__(self, red, ring, exps, factors, H, T):
        self.red, self.ring, self.H, self.T = red, ring, H, T
        self.exps = tuple(Fr(x) for x in exps)
        self.factors = list(factors)
        self._num = None
        self.method = None

    @property
    def num(self):
        """Expanded polynomial part."""
        if self._num is None:
            p = self.ring.one
            for f in self.factors:
                p = p * f
            self._num = p
        return self._num

    def exponents(self):
        """{ratio name: exponent} of the monomial prefactor (nonzero entries)."""
        return {b: e for b, e in zip(self.red.bary, self.exps) if e != 0}

    def in_ring(self, ring):
        d = RatioDensity(self.red, ring, self.exps, [f.set_ring(ring) for f in self.factors], self.H, self.T)
        d.method = self.method
        return d

    def as_expr(self, factored=True):
        """sympy expression in the ratio symbols (eliminated ratios appear as symbols b_0 in the
        monomial prefactor only)."""
        mono = sp.Mul(*[sp.Symbol(b) ** sp.Rational(e.numerator, e.denominator)
                        for b, e in zip(self.red.bary, self.exps) if e != 0])
        if factored:
            return mono * sp.Mul(*[sp.factor(f.as_expr()) for f in self.factors])
        return mono * self.num.as_expr()

    def evaluate(self, vals):
        """Value at vals (tuple over the ring generators) as (rational, radicand):
        g = rational * sqrt(radicand)."""
        pt = dict(zip(self.red.coords, vals))
        pt = self.red.complete_point(pt)
        rat = Fr(1)
        for f in self.factors:
            rat *= _evalp(f, vals)
        rad = Fr(1)
        for b, e in zip(self.red.bary, self.exps):
            if e == 0:
                continue
            if e.denominator == 1:
                rat *= pt[b] ** int(e)
            elif e.denominator == 2:
                k = (e.numerator - 1) // 2
                rat *= pt[b] ** k
                rad *= pt[b]
            else:
                raise ValueError("unsupported exponent %s" % e)
        return rat, rad

    def _integrate_polynomial(self, factors):
        """Integrate each ratio simplex as soon as its dependent factors are multiplied."""
        ring = self.ring
        blocks = [tuple(self.red._cidx[b] for b in block.indep) for block in self.red.blocks]
        if any(not factor for factor in factors):
            return ring.zero
        # The supports identify exactly which factors can be pulled out of an integral.
        def support(poly):
            return {j for j, degree in enumerate(poly.degrees()) if degree > 0}
        active = [(factor, support(factor)) for factor in factors]
        while blocks:
            def cost(block):
                terms = 1
                for factor, variables in active:
                    if variables.intersection(block):
                        terms *= len(factor)
                return terms, -len(block)
            block = min(blocks, key=cost)
            blocks.remove(block)
            dependent = [(p, v) for p, v in active if v.intersection(block)]
            active = [(p, v) for p, v in active if not v.intersection(block)]
            combined = ring.one
            variables = set()
            for factor, indices in sorted(dependent, key=lambda item: len(item[0])):
                combined *= factor
                variables.update(indices)
            moments, coefficients = {}, {}
            for mon, coefficient in combined.items():
                powers = tuple(mon[j] for j in block)
                if powers not in moments:
                    value, _ = dirichlet((0,) + powers)
                    moments[powers] = QQ(value.numerator, value.denominator)
                rest = list(mon)
                for j in block:
                    rest[j] = 0
                key = tuple(rest)
                coefficients[key] = coefficients.get(key, QQ.zero) + coefficient * moments[powers]
            integrated = ring.from_dict(coefficients)
            if not integrated:
                return ring.zero
            active.append((integrated, variables.difference(block)))
        result = ring.one
        for factor, _ in active:
            result *= factor
        return result

    def integrate(self, F=None):
        """int g F d xi over the product of simplices (Dirichlet integrals, exact)."""
        red, ring = self.red, self.ring
        if self.T == 0 and not any(self.exps) and red.blocks:
            factors = self.factors if F is None else [*self.factors, ring(F)]
            result = self._integrate_polynomial(factors)
            if len(ring.gens) == len(red.coords):
                return _fr(result.LC)
            return result.as_expr()
        p = self.num if F is None else self.num * F
        nc = len(red.coords)
        npar = len(ring.gens) - nc
        layout = []
        for blk in red.blocks:
            e0 = self.exps[red._bidx[blk.bary[0]]]
            layout.append((e0, [(red._cidx[b], self.exps[red._bidx[b]]) for b in blk.indep]))
        acc = {}
        for mon, c in p.terms():
            q, kpi = _fr(c), 0
            for e0, items in layout:
                dq, dk = dirichlet((e0,) + tuple(e + mon[j] for j, e in items))
                q *= dq
                kpi += dk
            key = (mon[nc:], kpi)
            acc[key] = acc.get(key, Fr(0)) + q
        if npar == 0:
            res = PiNumber({k: v for (_m, k), v in acc.items()})
            return res.rational() if self.T == 0 else res
        pars = [sp.Symbol(str(g)) for g in ring.gens[nc:]]
        expr = sp.Integer(0)
        for (mon, k), v in acc.items():
            if v:
                expr += sp.Rational(v.numerator, v.denominator) * sp.pi ** k * sp.Mul(*[s ** m for s, m in zip(pars, mon)])
        return sp.expand(expr)


def sp_reduce(G, D=None, plan=None, prefer="small"):
    """Series-parallel reduction of G for dimension D = H + T (default: (|I| - 1) / n).
    plan: explicit list of chain moves [(p, q, [chain, ...]), ...], chain = list of edge ids from p
    to q; reduced edges are numbered G.ne, G.ne + 1, ... in the order of the moves."""
    return SPReduction(G, D=D, plan=plan, prefer=prefer)
