"""[copy of explore/physics/6d-nilpotency/blocktri.py, H general]\nDensity of a series-parallel graph with a non-tight plan as a block lower-triangular determinant.

Move m of the plan merges chains i = 0..k-1 between p and q into the edge e_m.  For every chain i except a
reference chain r (the chain without original edges, if any; else r = 0) let the loop variable be the current
of the first original edge j of chain i,
    l_i = w_i y_red - mhat_i + Q_{i,j}            (htk CONVENTIONS section 5, direction p -> q).
In G_m the cotree formed by the loop variables of the later moves leaves the spanning tree T_m, and
    y_red^{(m)} = (constant combination of later loop variables) + T^{(m)},
T^{(m)} the flow through e_m in T_m from the sources of G_m (sum of the sources on the p side).  In the top
form prod_{m,i,a} d l^{(m)}_{i,a} the later differentials drop out, leaving the 1-forms
    R^{(m)}_{i,a} = y_red,a dw_i + w_i dT^{(m)}_a - d mhat_{i,a} + d Q_{i,j,a},   dw_0 = - sum_{i>=1} dw_i,
which contain only differentials of the coordinates of moves <= m.  The density is +- det[R] (sign fixed
against the brute-force Kirchhoff density at random rational points); det[R] is evaluated by the
generalized Laplace expansion along the row blocks.
"""
import itertools
import random
from fractions import Fraction as Fr

import numpy as np
import sympy as sp

from htk import brute
from htk.spred import _apply

import kern_h as kern
from kern_h import SymPoly, xp_add, xp_mul, xp_scale, xp_zero, xp_const, P


def _perm_sign(seq):
    inv = sum(1 for a in range(len(seq)) for b in range(a + 1, len(seq)) if seq[a] > seq[b])
    return -1 if inv & 1 else 1


class BlockTriangular:
    def __init__(self, gk, npoints=6):
        self.gk = gk
        red = gk.red
        G = red.G
        g = G.nv - 1
        ring0 = red.ring
        bp = red._bpoly
        one, zero = ring0.one, ring0.zero
        # forward source transfer with snapshots
        S = {v: [one if u == v else zero for u in range(g)] for v in range(g)}
        S[G.root] = [-one] * g
        snaps = []
        for mv in red.moves:
            addp, addq = [zero] * g, [zero] * g
            for i, ch in enumerate(mv.chains):
                a = len(ch)
                xs = [bp[b] for b in mv.xblocks[i].bary] if a >= 2 else [one]
                for kk, u in enumerate(mv.interiors[i]):
                    lam_u = S.pop(u)
                    sigma = zero
                    for j in range(kk + 1, a):
                        sigma = sigma + xs[j]
                    addp = [x + sigma * y for x, y in zip(addp, lam_u)]
                    addq = [x + (one - sigma) * y for x, y in zip(addq, lam_u)]
            S[mv.p] = [x + y for x, y in zip(S[mv.p], addp)]
            S[mv.q] = [x + y for x, y in zip(S[mv.q], addq)]
            snaps.append({v: list(vec) for v, vec in S.items()})
        # edge sets of G_m and the loop-variable edges
        cur = {e: ab for e, ab in enumerate(G.edges)}
        graphs_m = []
        loopedges = []
        self.choice = []
        for mv in red.moves:
            allmerged = [i for i, ch in enumerate(mv.chains) if all(e >= G.ne for e in ch)]
            if len(allmerged) > 1:
                raise ValueError('two chains without original edges')
            ref = allmerged[0] if allmerged else 0
            le, ch_sel = [], []
            for i, ch in enumerate(mv.chains):
                if i == ref:
                    continue
                j = next(jj for jj, e in enumerate(ch) if e < G.ne)
                le.append(ch[j])
                ch_sel.append((i, j))
            loopedges.append(le)
            self.choice.append((ref, ch_sel))
            cur = _apply(cur, mv.chains, mv.p, mv.q, mv.new_edge)
            graphs_m.append(dict(cur))
        cidx = {c: i for i, c in enumerate(red.coords)}
        self.colblocks = []
        for mv in red.moves:
            cols = []
            if mv.wblock is not None:
                cols += [cidx[b] for b in mv.wblock.indep]
            for blk in mv.xblocks:
                if blk is not None:
                    cols += [cidx[b] for b in blk.indep]
            self.colblocks.append(cols)
        assert sorted(c for cb in self.colblocks for c in cb) == list(range(len(red.coords)))
        # rows
        gens = ring0.gens
        nc = len(red.coords)
        rows = []          # list over (m, i, a) of (list over columns of vectors over vertices)
        self.rowblocks = []
        for m, mv in enumerate(red.moves):
            later = set(e for le in loopedges[m + 1:] for e in le)
            Gm = graphs_m[m]
            tree = {e: ab for e, ab in Gm.items() if e not in later}
            verts = sorted({x for ab in Gm.values() for x in ab})
            if len(tree) != len(verts) - 1:
                raise AssertionError('T_m is not a spanning tree (move %d)' % m)
            adj = {v: [] for v in verts}
            for e, (a, b) in tree.items():
                if e == mv.new_edge:
                    continue
                adj[a].append(b)
                adj[b].append(a)
            side, stack = {mv.p}, [mv.p]
            while stack:
                v = stack.pop()
                for w in adj[v]:
                    if w not in side:
                        side.add(w)
                        stack.append(w)
            if mv.q in side:
                raise AssertionError('T_m is not a tree')
            snap = snaps[m]
            Tm = [zero] * g
            for v in side:
                Tm = [x + y for x, y in zip(Tm, snap[v])]
            yred = mv.yred
            w = [bp[b] for b in mv.wblock.bary] if mv.wblock is not None else [one]
            rb = []
            ref, ch_sel = self.choice[m]
            windep = [cidx[b] for b in mv.wblock.indep]
            for i, jsel in ch_sel:
                coeffs = []
                for j in range(nc):
                    vec = []
                    for v in range(g):
                        e = (w[i] * Tm[v].diff(gens[j]) - mv.mhat[i][v].diff(gens[j])
                             + mv.Q[i][jsel][v].diff(gens[j]))
                        if i >= 1 and j == cidx[mv.wblock.bary[i]]:
                            e = e + yred[v]
                        elif i == 0 and j in windep:
                            e = e - yred[v]
                        vec.append(e)
                    coeffs.append(vec)
                for a in range(self.gk.H):
                    rb.append(len(rows))
                    rows.append((coeffs, a))
            self.rowblocks.append(rb)
        # block triangularity
        for m, rb in enumerate(self.rowblocks):
            later_cols = [c for cb in self.colblocks[m + 1:] for c in cb]
            for r in rb:
                for c in later_cols:
                    assert all(x == 0 for x in rows[r][0][c]), 'not block triangular'
        self.rows = rows
        self.g = g
        self.nc = nc
        # Laplace terms
        self.terms = []

        def rec(m, avail, chosen):
            if m == len(red.moves):
                if not avail:
                    order = [c for U in chosen for c in U]
                    self.terms.append((_perm_sign(order), [tuple(U) for U in chosen]))
                return
            pool = sorted(avail + self.colblocks[m])
            r = len(self.rowblocks[m])
            for U in itertools.combinations(pool, r):
                rest = [c for c in pool if c not in U]
                rec(m + 1, rest, chosen + [list(U)])
        rec(0, [], [])
        # symbolic entries in the parameter ring
        outs = gk.outputs
        lam_sym = gk.lam_sym
        ring, _bp, _resp, _conv = red._setup(red._flat(lam_sym))
        self.ring = ring
        names = [str(x) for x in ring.gens]
        pg = dict(zip(names, ring.gens))
        self.entries = {}
        for r, (coeffs, a) in enumerate(rows):
            for c in range(nc):
                e = ring.zero
                for v in outs:
                    x = coeffs[c][v]
                    if x != 0:
                        e = e + x.set_ring(ring) * pg['l%d_%d' % (v, a)]
                self.entries[(r, c)] = e
        self.sign = self._fix_sign(npoints)
        self.sym = {k: SymPoly(v, nc, gk.lam_index) for k, v in self.entries.items() if v != 0}

    # ------------------------------------------------------------------ exact check
    def _numeric_det(self, vals, lamn):
        nc = self.nc
        M = [[Fr(0)] * nc for _ in range(len(self.rows))]
        for r, (coeffs, a) in enumerate(self.rows):
            for c in range(nc):
                s = Fr(0)
                for v in range(self.g):
                    x = coeffs[c][v]
                    if x != 0 and lamn[v][a] != 0:
                        s += kern_eval(x, vals) * lamn[v][a]
                M[r][c] = s
        from htk.linalg import det as fdet
        return fdet(M)

    def _fix_sign(self, npoints):
        red = self.gk.red
        rng = random.Random(4242)
        sign = None
        for _ in range(npoints):
            lamn = [[Fr(rng.randint(-5, 5)) if v in self.gk.outputs else Fr(0) for _a in range(self.gk.H)]
                    for v in range(self.g)]
            pt = red.random_point(rng)
            vals = red._vals(pt)
            t = red.tmap(pt)
            Dt = brute.density(red.G, t, lamn, T=0)
            ref = Dt.coeff * abs(red.jacobian_t(pt))
            mine = self._numeric_det(vals, lamn)
            if ref == 0 and mine == 0:
                continue
            s = Fr(ref) / mine if mine else None
            if s not in (1, -1):
                raise AssertionError('block determinant differs from the Kirchhoff density: ratio %s' % s)
            if sign is None:
                sign = int(s)
            elif sign != s:
                raise AssertionError('sign of the block determinant is not constant')
        if sign is None:
            raise AssertionError('no nonzero check point')
        self.checked = npoints
        return sign

    # ------------------------------------------------------------------ evaluation
    def evaluate(self, arr, L):
        """list of XP factors: [common block determinants..., sum over Laplace terms of the rest]."""
        ent = {k: v.evaluate(arr, L) for k, v in self.sym.items()}

        def det(rows, cols):
            rows, cols = list(rows), list(cols)
            dp = {(): xp_const(1, L)}
            for i, r in enumerate(rows):
                new = {}
                for used, val in dp.items():
                    for c in cols:
                        if c in used or (r, c) not in ent:
                            continue
                        pos = sum(1 for u in used if cols.index(u) > cols.index(c))
                        term = xp_mul(val, ent[(r, c)])
                        if pos & 1:
                            term = xp_scale(term, P - 1)
                        key = tuple(sorted(used + (c,), key=cols.index))
                        new[key] = xp_add(new[key], term) if key in new else term
                dp = {k: v for k, v in new.items() if len(v)}
                if not dp:
                    return xp_zero(L)
            return dp.get(tuple(cols), xp_zero(L))
        # factor out blocks common to all terms
        nb = len(self.rowblocks)
        common = [m for m in range(nb) if len({t[1][m] for t in self.terms}) == 1]
        factors = []
        for m in common:
            factors.append(det(self.rowblocks[m], self.terms[0][1][m]))
        cache = {}
        parts = []
        for sgn, Us in self.terms:
            prod = None
            for m in range(nb):
                if m in common:
                    continue
                key = (m, Us[m])
                if key not in cache:
                    cache[key] = det(self.rowblocks[m], Us[m])
                prod = cache[key] if prod is None else xp_mul(prod, cache[key])
            if prod is None:
                prod = xp_const(1, L)
            s = sgn * self.sign
            parts.append(prod if s > 0 else xp_scale(prod, P - 1))
        parts = [p for p in parts if len(p)]
        factors.append(xp_add(*parts) if parts else xp_zero(L))
        return factors


def kern_eval(poly, vals):
    tot = Fr(0)
    for mon, c in poly.terms():
        term = Fr(int(c.numerator), int(c.denominator))
        for v, k in zip(vals, mon):
            if k:
                term *= v ** k
        tot += term
    return tot
