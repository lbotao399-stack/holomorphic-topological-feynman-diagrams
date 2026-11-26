"""Linear system of exact, regulator-free relations among convergent vacuum periods of one
core at fixed q (normalized p = P / prod Gamma(a_e)):

 (T) triangle rule, used only when all five periods converge (zero powers = contraction)
 (I) parametric IBP:  sum_{C not containing i} I(a+1_C) = ((q-1-a_i)/a_i) sum_{C containing i} I(a+1_C)
     for a base a at level q-1 (all a_e>=1, a_i<=q-2), written for P (not p);
 (B) boundary:  I_{G/i}(a_{-i}; q-1) = (q-1) sum_{C containing i} I(a+1_C; q),  a_i = 0.
Unknowns: orbits of convergent power vectors of the core at level q.  Contracted graphs are
evaluated by vac.full_period (series/parallel + wheels).  Solves modulo a prime for the
rank, then exactly over Q for the determined unknowns.
"""
import sys
import time
from fractions import Fraction as F
from math import factorial, prod

import gutil
import vac


def contract_value(edges, pw, q):
    """P for a power vector with zeros (contracted graph), or None if divergent"""
    ne, npw = list(edges), list(pw)
    while 0 in npw:
        i = npw.index(0)
        u, v = ne[i]
        if u == v:
            return None
        ne, npw, _ = vac.contract(ne, npw, [0] * len(npw), i)
    if any(a == b for a, b in ne) or not gutil.convergent(tuple(ne), npw, q):
        return None
    try:
        return vac.full_period(ne, npw, q, mark=None)
    except (vac.Stuck, RecursionError):
        return 'stuck'


class System:
    def __init__(self, edges, q):
        self.edges = tuple(edges)
        self.q = q
        self.autos = gutil.edge_automorphisms(self.edges)
        self.n = gutil.cycle_rank(self.edges)
        self.cot = gutil.cotrees(self.edges)
        self.index = {}
        self.reps = []
        self.rows = []      # list of (dict unknown->coef, const dict label->coef)  meaning sum coef*P_unknown = const
        self.stuck_consts = 0
        self.use_star = True
        self.star_known = 0

    def var(self, pw):
        rep = gutil.orbit_rep(pw, self.autos)
        if rep not in self.index:
            self.index[rep] = len(self.reps)
            self.reps.append(rep)
        return self.index[rep]

    def term(self, pw, coef, row, const):
        """add coef * P(pw) to the relation (unknown if all >= 1, else known boundary)"""
        q = self.q
        if any(p < 0 for p in pw):
            return False
        if 0 in pw:
            val = contract_value(self.edges, pw, q)
            if val is None or val == 'stuck':
                return False
            for k, v in val.items():
                const[k] = const.get(k, 0) - coef * v
            return True
        if not gutil.convergent(self.edges, list(pw), q):
            return False
        j = self.var(tuple(pw))
        row[j] = row.get(j, 0) + coef
        return True

    def add_row(self, items):
        row, const = {}, {}
        for pw, coef in items:
            if coef == 0:
                continue
            if not self.term(pw, F(coef), row, const):
                return False
        row = {k: v for k, v in row.items() if v}
        const = {k: v for k, v in const.items() if v}
        if row or const:
            self.rows.append((row, const))
        return True

    # ---- relation generators ----
    def triangle_rules(self, pw):
        """triangle rule written for P = p * prod Gamma: p(a') = P(a')/prod Gamma(a'_e)"""
        edges, q = self.edges, self.q
        verts = sorted({x for e in edges for x in e})
        deg = {v: sum(1 for e in edges if v in e) for v in verts}
        out = []
        for i, (A0, B0) in enumerate(edges):
            for A, B in ((A0, B0), (B0, A0)):
                if deg[A] != 3 or deg[B] != 3:
                    continue
                incA = [t for t, e in enumerate(edges) if A in e and t != i]
                incB = [t for t, e in enumerate(edges) if B in e and t != i]
                for t3 in incA:
                    C = edges[t3][0] if edges[t3][1] == A else edges[t3][1]
                    for t2 in incB:
                        C2 = edges[t2][0] if edges[t2][1] == B else edges[t2][1]
                        if C2 != C or C in (A, B):
                            continue
                        t5 = [t for t in incA if t != t3][0]
                        t4 = [t for t in incB if t != t2][0]
                        a1, a2, a3 = pw[i], pw[t2], pw[t3]
                        c0 = 2 * q - 2 * a1 - a2 - a3
                        items = [(list(pw), c0)]
                        for (x, dx), (y, dy), coef in (((i, -1), (t2, 1), -a2), ((t2, 1), (t4, -1), a2),
                                                       ((i, -1), (t3, 1), -a3), ((t3, 1), (t5, -1), a3)):
                            np_ = list(pw)
                            np_[x] += dx
                            np_[y] += dy
                            items.append((np_, coef))
                        # convert p-relation to P-relation: p(b) = P(b)/G(b), G(b) = prod (b_e-1)!
                        conv = []
                        base = prod(factorial(x - 1) for x in pw)
                        ok = True
                        for np_, coef in items:
                            if min(np_) < 0:
                                ok = False
                                break
                            g = prod(factorial(x - 1) for x in np_ if x > 0)
                            conv.append((np_, F(coef) * base / g))
                        if ok:
                            out.append(conv)
        return out

    def ibp_rules(self, base):
        """base at level q-1 (sum = n(q-1))"""
        q = self.q
        out = []
        for i in range(len(self.edges)):
            ai = base[i]
            if ai == 0:
                # boundary relation: I_{G/i}(base_{-i}; q-1) = (q-1) sum_{C ∋ i} I(base+1_C; q)
                items = []
                for C in self.cot:
                    if i in C:
                        np_ = list(base)
                        for e in C:
                            np_[e] += 1
                        items.append((np_, q - 1))
                # constant from the contracted graph at level q-1
                out.append(('B', i, items))
            elif ai <= q - 2:
                items = []
                for C in self.cot:
                    np_ = list(base)
                    for e in C:
                        np_[e] += 1
                    coef = F(1) if i not in C else -F(q - 1 - ai, ai)
                    items.append((np_, coef))
                out.append(('I', i, items))
        return out

    def build(self, use_ibp=True, use_tri=True, verbose=True):
        q, N = self.q, len(self.edges)
        t = time.time()
        level = [c for c in gutil.compositions(self.n * q, N, 1, q - 1) if gutil.convergent(self.edges, list(c), q)]
        seen = set()
        reps = []
        for c in level:
            r = gutil.orbit_rep(c, self.autos)
            if r not in seen:
                seen.add(r)
                reps.append(r)
        for r in reps:
            self.var(r)
        if use_tri:
            for r in reps:
                for items in self.triangle_rules(r):
                    self.add_row(items)
        if self.use_star:
            import startri
            for r in reps:
                us = startri.unique_vertices(self.edges, list(r), q)
                if not us:
                    continue
                ne, npw = startri.star_to_triangle(self.edges, list(r), us[0], q)
                try:
                    val = vac.full_period(ne, npw, q, mark=None)
                except (vac.Stuck, RecursionError):
                    continue
                const = {k: v for k, v in val.items()}
                self.rows.append(({self.var(r): F(1)}, const))
                self.star_known += 1
        if use_ibp:
            bases = set()
            for c in gutil.compositions(self.n * (q - 1), N, 0, q - 1):
                bases.add(gutil.orbit_rep(c, self.autos))
            for b in bases:
                for kind, i, items in self.ibp_rules(b):
                    if kind == 'I':
                        if not all(x >= 1 for x in b):
                            continue
                        self.add_row(items)
                    else:
                        if any(x < 1 for j, x in enumerate(b) if j != i) or sum(1 for x in b if x == 0) != 1:
                            continue
                        val = contract_value(self.edges, list(b), q - 1)
                        if val is None or val == 'stuck':
                            continue
                        row, const = {}, {}
                        ok = True
                        for np_, coef in items:
                            if not self.term(np_, F(coef), row, const):
                                ok = False
                                break
                        if not ok:
                            continue
                        for k, v in val.items():
                            const[k] = const.get(k, 0) + v
                        row = {k: v for k, v in row.items() if v}
                        const = {k: v for k, v in const.items() if v}
                        if row or const:
                            self.rows.append((row, const))
        if verbose:
            print('q=%d unknowns %d (convergent orbits %d) rows %d (star-triangle values %d)  %.1fs' % (q, len(self.reps), len(reps), len(self.rows), self.star_known, time.time() - t), flush=True)
        return reps

    def solve(self, p=2147483629, verbose=True):
        """Gaussian elimination over Q (sparse, Fraction). Returns dict unknown -> value dict or None."""
        t = time.time()
        labels = sorted({k for r, c in self.rows for k in c} | {'1'})
        # augmented rows: unknown coefs + constants per label
        pivots = {}   # col -> (row dict, const dict)
        order = sorted(range(len(self.rows)), key=lambda i: len(self.rows[i][0]))
        for idx in order:
            row, const = dict(self.rows[idx][0]), dict(self.rows[idx][1])
            # reduce by existing pivots
            changed = True
            while True:
                cols = [c for c in row if c in pivots]
                if not cols:
                    break
                c = cols[0]
                prow, pconst = pivots[c]
                f = row[c]
                for k, v in prow.items():
                    row[k] = row.get(k, 0) - f * v
                    if row[k] == 0:
                        del row[k]
                for k, v in pconst.items():
                    const[k] = const.get(k, 0) - f * v
                    if const[k] == 0:
                        del const[k]
            if not row:
                if const:
                    raise ArithmeticError('inconsistent relation: %s' % const)
                continue
            c = min(row, key=lambda k: (k))
            f = row[c]
            row = {k: v / f for k, v in row.items()}
            const = {k: v / f for k, v in const.items()}
            # eliminate c from other pivots (keep reduced)
            for c2, (r2, k2) in list(pivots.items()):
                if c in r2:
                    g = r2[c]
                    for k, v in row.items():
                        r2[k] = r2.get(k, 0) - g * v
                        if r2[k] == 0:
                            del r2[k]
                    for k, v in const.items():
                        k2[k] = k2.get(k, 0) - g * v
                        if k2[k] == 0:
                            del k2[k]
            pivots[c] = (row, const)
        solved = {}
        for c, (row, const) in pivots.items():
            if len(row) == 1:
                solved[c] = const
        if verbose:
            print('rank %d of %d unknowns, fully determined %d  %.1fs' % (len(pivots), len(self.reps), len(solved), time.time() - t), flush=True)
        return solved, pivots
