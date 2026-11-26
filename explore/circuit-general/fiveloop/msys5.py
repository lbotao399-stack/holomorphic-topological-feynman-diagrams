"""[fiveloop copy of cores/msys.py; change: planar-duality hook in _classify (DUAL_HOOK): a reduced
graph G that is not a family member, planar with rank(G*) = |V|-1 < rank(G), is replaced by its dual
with powers q - a (I_q(G; a) = I_q(G*; q - a o pi), symmetry/d0_planar_duality.py, check_rel5.log (P)),
so that duals of family members become unknowns.]
Multi-family linear system of exact, regulator-free relations among convergent
projective periods  P(G; a) = int prod t^(a-1)/U^q  at fixed q.

Unknowns: orbits of convergent power vectors on a list of 'family' graphs (e.g. prism,
K33, K5-e).  Relations (each used only if every period in it converges):
 (T) triangle rule for p = P/prod Gamma(a)  (momentum IBP of the cut integral, eps -> 0);
 (I) parametric IBP at base b (level q-1, all b_e>=1, b_i<=q-2):
       sum_{C not∋ i} P(b+1_C) = ((q-1-b_i)/b_i) sum_{C∋ i} P(b+1_C);
 (B) boundary: P(G/i; b_{-i}, q-1) = (q-1) sum_{C∋ i} P(b+1_C),  b_i = 0;
 (S) star-triangle at a unique vertex (sum of the three powers = q): P(G;a) = P(G';a').
Every graph met is reduced by contraction of zero powers, series and parallel push-forward
(normalized factor G_2q) and wheel recognition; what remains must be a family member.
"""
import time
from fractions import Fraction as F
from itertools import permutations
from math import factorial, prod

import gutil
import halfq
import startri
import vac
import wheel


def struct_key(edges):
    verts = sorted({x for e in edges for x in e})
    best = None
    for perm in permutations(range(len(verts))):
        mp = dict(zip(verts, perm))
        k = tuple(sorted((min(mp[u], mp[v]), max(mp[u], mp[v])) for u, v in edges))
        if best is None or k < best:
            best = k
    return best


class Unknown(Exception):
    pass


DUAL_HOOK = [True]


class MultiSystem:
    def __init__(self, families, q):
        self.q = q
        self.fam = []
        self.fam_by_key = {}
        for edges in families:
            edges = tuple(edges)
            d = dict(edges=edges, autos=gutil.edge_automorphisms(edges), n=gutil.cycle_rank(edges),
                     cot=gutil.cotrees(edges), skey=struct_key(edges))
            self.fam_by_key[d['skey']] = len(self.fam)
            self.fam.append(d)
        self.fam_edges = {d['edges']: i for i, d in enumerate(self.fam)}
        self.fam_sigs = set()
        for d in self.fam:
            vs = sorted({x for e in d['edges'] for x in e})
            dg = {v: sum(1 for e in d['edges'] if v in e) for v in vs}
            self.fam_sigs.add((len(vs), len(d['edges']), tuple(sorted(dg.values()))))
        self.class_cache = {}
        self.skey_cache = {}
        self.index = {}
        self.reps = []
        self.rows = []
        self.iso_cache = {}
        self.counts = {'T': 0, 'D': 0, 'I': 0, 'B': 0, 'S': 0}
        self.tuple_labels = F(q).denominator != 1
        self.lower = None          # callable(edges, powers, q) -> value dict or None
        self.use_diamond = True

    def var(self, f, pw):
        rep = gutil.orbit_rep(tuple(pw), self.fam[f]['autos'])
        key = (f, rep)
        if key not in self.index:
            self.index[key] = len(self.reps)
            self.reps.append(key)
        return self.index[key]

    # ---------- reduction of an arbitrary graph to (known value) or (factor * unknown) ----------
    def classify(self, edges, pw, q):
        """returns (const_dict, unknown_list) with P(edges;pw) = const + sum c*P_unknown,
        or None if divergent/unrecognized.  Memoized; fast path for family members."""
        te = tuple(edges)
        if te in self.fam_edges and min(pw) >= 1:
            f = self.fam_edges[te]
            if not gutil.convergent(te, list(pw), q):
                return None
            return ({}, [(self.var(f, pw), F(1))])
        key = (te, tuple(pw), q)
        if key in self.class_cache:
            return self.class_cache[key]
        res = self._classify(edges, pw, q)
        self.class_cache[key] = res
        return res

    def _classify(self, edges, pw, q):
        edges = list(edges)
        pw = list(pw)
        factor = F(1)          # P = factor * (prod Gamma(pw)/prod Gamma(final))-style bookkeeping
        sqk = 0                # extra power of sqrt(pi) (half-integer q)
        # work with normalized p: P = p * prod Gamma(a)
        norm0 = prod(factorial(x - 1) for x in pw if x > 0)
        if any(x < 0 for x in pw):
            return None
        while True:
            if any(x == 0 for x in pw):
                i = pw.index(0)
                u, v = edges[i]
                if u == v:
                    return None
                edges, pw, _ = vac.contract(edges, pw, [0] * len(pw), i)
                continue
            if any(u == v for u, v in edges):
                return None
            verts = sorted({x for e in edges for x in e})
            if len(verts) == 2:
                # banana of k edges: P = prod Gamma(q - s_i) / Gamma(q), sum s = (k-1) q
                k = len(edges)
                if sum(pw) != (k - 1) * q or any(not (0 < x < q) for x in pw):
                    return None
                num, nk = F(1), 0
                for x in pw:
                    r, kk = halfq.gamma(q - x)
                    num *= r
                    nk += kk
                    r, kk = halfq.gamma(x)
                    num /= r
                    nk -= kk
                r, kk = halfq.gamma(q)
                num /= r
                nk -= kk
                val = factor * norm0 * num          # P_orig = factor * norm0 * p_final
                nk += sqk
                if not self.tuple_labels and nk == 0:
                    return ({'1': val}, [])
                return ({((('sqrtpi', nk),) if nk else ()): val}, [])
            deg = {v: 0 for v in verts}
            for u, v in edges:
                deg[u] += 1
                deg[v] += 1
            done = False
            for v in verts:
                if deg[v] == 2:
                    i, j = [t for t, e in enumerate(edges) if v in e]
                    a = edges[i][0] if edges[i][1] == v else edges[i][1]
                    b = edges[j][0] if edges[j][1] == v else edges[j][1]
                    if a == b:
                        return None
                    keep = [t for t in range(len(edges)) if t not in (i, j)]
                    edges = [edges[t] for t in keep] + [(a, b)]
                    pw = [pw[t] for t in keep] + [pw[i] + pw[j]]
                    done = True
                    break
            if done:
                continue
            seen = {}
            for i, (u, v) in enumerate(edges):
                pr = (min(u, v), max(u, v))
                if pr in seen:
                    j = seen[pr]
                    s1, s2 = pw[j], pw[i]
                    if not (0 < s1 < q and 0 < s2 < q and s1 + s2 > q):
                        return None
                    keep = [t for t in range(len(edges)) if t not in (i, j)]
                    edges = [edges[t] for t in keep] + [pr]
                    pw = [pw[t] for t in keep] + [s1 + s2 - q]
                    if F(q).denominator == 1:
                        factor *= vac.G2q(s1, s2, q)
                    else:
                        r, k = halfq.G2q(s1, s2, q)
                        factor *= r
                        sqk += k
                    done = True
                    break
                seen[pr] = i
            if done:
                if len(edges) == 1:
                    return None
                continue
            break
        if not gutil.convergent(tuple(edges), pw, q):
            return None
        verts = sorted({x for e in edges for x in e})
        deg = {v: 0 for v in verts}
        for u, v in edges:
            deg[u] += 1
            deg[v] += 1
        if all(F(x).denominator == 1 for x in pw):
            normf = F(prod(factorial(int(x) - 1) for x in pw))
        else:
            normf = F(1)
            for x in pw:
                r, k = halfq.gamma(x)
                normf *= r
                sqk -= k
        # P_original = p_original * norm0 ; p_original = factor * p_final ; p_final = P_final/normf
        scale = factor * F(norm0) / normf
        w = vac.wheel_form(edges, pw, verts, deg)
        if w is not None:
            if F(q).denominator == 1 and not self.tuple_labels:
                val = wheel.integrate(w[0], w[1], q)
                return ({('1' if k == 0 else 'z%d' % k): v * scale for k, v in val.items() if v}, [])
            val = halfq.wheel_value(w[0], w[1], q)
            return (halfq.scale_value(val, scale, sqk), [])
        te2 = tuple(sorted((min(u, v), max(u, v)) for u, v in edges))
        sk = self.skey_cache.get(te2)
        if sk is None:
            sig = (len(verts), len(edges), tuple(sorted(deg.values())))
            if sig in self.fam_sigs:
                sk = struct_key(edges)
            else:
                sk = ('nofamily', sig)          # no family has this degree signature
            self.skey_cache[te2] = sk
        if sk not in self.fam_by_key:
            if DUAL_HOOK[0] and F(q).denominator == 1 and not sqk and 2 * len(verts) - 2 < len(edges):
                import rel5
                dm = rel5.dual_map([tuple(e) for e in edges])
                if dm is not None and all(a != b for a, b in dm):
                    res = self.classify(dm, [q - x for x in pw], q)
                    if res is not None:
                        c, unk = res
                        return ({k: v * scale for k, v in c.items() if v}, [(j, s_ * scale) for j, s_ in unk])
            # lower cores (prism, K33) with known tables
            if self.lower is not None:
                val = self.lower(edges, pw, q)
                if val is None:
                    return None
                if self.tuple_labels:
                    return (halfq.scale_value(val, scale, sqk), [])
                if sqk:
                    return None
                return ({k: v * scale for k, v in val.items() if v}, [])
            return None
        if sqk:
            raise ValueError('family member reached through a sqrt(pi) bubble factor')
        f = self.fam_by_key[sk]
        mapped = self.map_to_family(edges, pw, f)
        if mapped is None:
            return None
        return ({}, [(self.var(f, mapped), scale)])

    def map_to_family(self, edges, pw, f):
        fe = self.fam[f]['edges']
        verts = sorted({x for e in edges for x in e})
        fv = sorted({x for e in fe for x in e})
        idx = {frozenset(e): i for i, e in enumerate(fe)}
        key = (tuple(edges), f)
        perm = self.iso_cache.get(key)
        if perm is None:
            for p in permutations(fv):
                mp = dict(zip(verts, p))
                img = []
                ok = True
                for u, v in edges:
                    j = idx.get(frozenset((mp[u], mp[v])))
                    if j is None:
                        ok = False
                        break
                    img.append(j)
                if ok and len(set(img)) == len(img):
                    perm = img
                    break
            if perm is None:
                return None
            self.iso_cache[key] = perm
        out = [0] * len(fe)
        for t, j in enumerate(perm):
            out[j] = pw[t]
        return out

    # ---------- relations ----------
    def add_relation(self, items, kind):
        """items: list of (edges, powers, coef) meaning sum coef * P = 0"""
        row, const = {}, {}
        for edges, pw, coef in items:
            if coef == 0:
                continue
            res = self.classify(edges, pw, self.q)
            if res is None:
                return False
            c, unk = res
            for k, v in c.items():
                const[k] = const.get(k, 0) - coef * v
            for j, s in unk:
                row[j] = row.get(j, 0) + coef * s
        row = {k: v for k, v in row.items() if v}
        const = {k: v for k, v in const.items() if v}
        if row or const:
            self.rows.append((row, const))
            self.counts[kind] += 1
        return True

    def triangle_rows(self, f, pw):
        edges, q = self.fam[f]['edges'], self.q
        verts = sorted({x for e in edges for x in e})
        deg = {v: sum(1 for e in edges if v in e) for v in verts}
        base = prod(factorial(x - 1) for x in pw)
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
                        if any(min(np_) < 0 for np_, c in items):
                            continue
                        conv = [(edges, np_, F(coef) * base / prod(factorial(x - 1) for x in np_ if x > 0))
                                for np_, coef in items]
                        self.add_relation(conv, 'T')

    def diamond_rows(self, f, pw):
        """diamond rule (Ruijl-Ueda-Vermaseren) at the period level, d = 2q:
        B: all neighbours X_i cubic, each X_i adjacent to a common T (U_i = X_i T),
        third line E_i at X_i; spectators = direct B-T edges (power s_j).
        E p = sum_i a_i [p(U_i+, D_i-) - p(U_i+, E_i-)],
        E = (L+S) 2q - sum_i (a_i + 2 b_i) - 2 sum_j s_j,  L+1 = number of X_i."""
        edges, q = self.fam[f]['edges'], self.q
        verts = sorted({x for e in edges for x in e})
        nbr = {v: [] for v in verts}
        for t, (u, v) in enumerate(edges):
            nbr[u].append((v, t))
            nbr[v].append((u, t))
        base = prod(factorial(x - 1) for x in pw)
        for B in verts:
            for T in verts:
                if T == B:
                    continue
                spect = [t for (w, t) in nbr[B] if w == T]
                Xs = [(w, t) for (w, t) in nbr[B] if w != T]
                if len(Xs) < 2:
                    continue
                ok = True
                paths = []
                for X, tD in Xs:
                    if len(nbr[X]) != 3:
                        ok = False
                        break
                    ups = [t for (w, t) in nbr[X] if w == T]
                    others = [t for (w, t) in nbr[X] if w not in (T, B)]
                    if len(ups) != 1 or len(others) != 1:
                        ok = False
                        break
                    paths.append((ups[0], tD, others[0]))
                if not ok:
                    continue
                L = len(paths) - 1
                S = len(spect)
                E = (L + S) * 2 * q - sum(pw[u] + 2 * pw[dd] for u, dd, e in paths) - 2 * sum(pw[t] for t in spect)
                items = [(list(pw), E)]
                for u, dd, e in paths:
                    for low, sign in ((dd, -1), (e, +1)):
                        np_ = list(pw)
                        np_[u] += 1
                        np_[low] -= 1
                        items.append((np_, sign * pw[u]))
                if any(min(np_) < 0 for np_, c in items):
                    continue
                conv = [(edges, np_, F(coef) * base / prod(factorial(x - 1) for x in np_ if x > 0))
                        for np_, coef in items]
                self.add_relation(conv, 'D')

    def build(self, verbose=True):
        q = self.q
        t0 = time.time()
        level_reps = []
        for f, d in enumerate(self.fam):
            N = len(d['edges'])
            reps = set()
            amax = -(-F(q).numerator // F(q).denominator) - 1          # largest integer < q
            for c in gutil.compositions(int(d['n'] * q), N, 1, amax):
                if gutil.convergent(d['edges'], list(c), q):
                    reps.add(gutil.orbit_rep(c, d['autos']))
            reps = sorted(reps)
            for r in reps:
                self.var(f, r)
            level_reps.append(reps)
        # star-triangle links
        for f, reps in enumerate(level_reps):
            edges = self.fam[f]['edges']
            for r in reps:
                for u in startri.unique_vertices(edges, list(r), q):
                    ne, npw = startri.star_to_triangle(edges, list(r), u, q)
                    self.add_relation([(edges, list(r), 1), (ne, npw, -1)], 'S')
        # triangle and diamond rules
        for f, reps in enumerate(level_reps):
            for r in reps:
                self.triangle_rows(f, list(r))
                if self.use_diamond:
                    self.diamond_rows(f, list(r))
        # IBP and boundary relations
        for f, d in enumerate(self.fam):
            edges, N, n, cot = d['edges'], len(d['edges']), d['n'], d['cot']
            bases = set()
            amax1 = -(-F(q - 1).numerator // F(q - 1).denominator) - 1
            for c in gutil.compositions(int(n * (q - 1)), N, 0, amax1 + 1):
                bases.add(gutil.orbit_rep(c, d['autos']))
            for b in bases:
                zeros = [i for i in range(N) if b[i] == 0]
                if len(zeros) > 1:
                    continue
                if len(zeros) == 1:
                    i = zeros[0]
                    items = [(edges, list(b), F(-1, 1))]   # P(G/i; b, q-1) handled below
                    lower = self.classify(edges, list(b), q - 1)
                    if lower is None or lower[1]:
                        continue
                    row, const = {}, dict(lower[0])
                    ok = True
                    for C in cot:
                        if i in C:
                            np_ = list(b)
                            for e in C:
                                np_[e] += 1
                            res = self.classify(edges, np_, q)
                            if res is None:
                                ok = False
                                break
                            c, unk = res
                            for k, v in c.items():
                                const[k] = const.get(k, 0) - (q - 1) * v
                            for j, s in unk:
                                row[j] = row.get(j, 0) + (q - 1) * s
                    if not ok:
                        continue
                    row = {k: v for k, v in row.items() if v}
                    const = {k: v for k, v in const.items() if v}
                    if row or const:
                        self.rows.append((row, const))
                        self.counts['B'] += 1
                    continue
                for i in range(N):
                    if b[i] >= q - 1:
                        continue
                    items = []
                    for C in cot:
                        np_ = list(b)
                        for e in C:
                            np_[e] += 1
                        coef = F(1) if i not in C else -F(q - 1 - b[i], b[i])
                        items.append((edges, np_, coef))
                    self.add_relation(items, 'I')
        if verbose:
            print('q=%d unknowns %d rows %d %s  %.1fs' % (q, len(self.reps), len(self.rows), self.counts,
                                                          time.time() - t0), flush=True)

    def solve(self, verbose=True):
        t = time.time()
        pivots = {}
        order = sorted(range(len(self.rows)), key=lambda i: len(self.rows[i][0]))
        for idx in order:
            row, const = dict(self.rows[idx][0]), dict(self.rows[idx][1])
            while True:
                cols = [c for c in row if c in pivots]
                if not cols:
                    break
                c = cols[0]
                prow, pconst = pivots[c]
                f = row[c]
                for k, v in prow.items():
                    nv = row.get(k, 0) - f * v
                    if nv:
                        row[k] = nv
                    else:
                        row.pop(k, None)
                for k, v in pconst.items():
                    nv = const.get(k, 0) - f * v
                    if nv:
                        const[k] = nv
                    else:
                        const.pop(k, None)
            if not row:
                if const:
                    raise ArithmeticError('inconsistent relation: %s' % const)
                continue
            c = min(row)
            f = row[c]
            row = {k: v / f for k, v in row.items()}
            const = {k: v / f for k, v in const.items()}
            for c2, (r2, k2) in pivots.items():
                if c in r2:
                    g = r2[c]
                    for k, v in row.items():
                        nv = r2.get(k, 0) - g * v
                        if nv:
                            r2[k] = nv
                        else:
                            r2.pop(k, None)
                    for k, v in const.items():
                        nv = k2.get(k, 0) - g * v
                        if nv:
                            k2[k] = nv
                        else:
                            k2.pop(k, None)
            pivots[c] = (row, const)
        solved = {c: const for c, (row, const) in pivots.items() if len(row) == 1}
        if verbose:
            print('rank %d of %d unknowns, determined %d  %.1fs' % (len(pivots), len(self.reps), len(solved),
                                                                   time.time() - t), flush=True)
        return solved, pivots
