"""Recursive exact evaluation of single periods of a five-loop core with triangles (Z5) at any q,
without the full table: memoized DFS over
  direct maps: (Y) unique triangle -> star, (S) unique vertex -> triangle, (C) completion decompletions,
               (P) planar self-duality (only to vectors of smaller potential);
  recursion:   (T) triangle rule; every T rule lowers  Phi = sum of the edges {AB, third at A, third at B}
               over all triangle rules (checked: raised edges never in that set), terms with a zero
               power are contractions (per5.value).
Targets outside the family are evaluated by per5.value (K5-e = dual prism, 2-sums, W_n, 4-loop tables,
tabulated five-loop cores)."""
import sys
from fractions import Fraction as F
from math import factorial, prod

import per5
import gutil
import rel5
import startri


class Rec:
    def __init__(self, name, q, use_table=True):
        self.name, self.q = name, q
        self.E = list(per5.SKEL[name])
        self.autos = gutil.edge_automorphisms(tuple(self.E))
        self.memo = {}
        self.table = per5.table(name, q) if use_table else None
        self.tris = rel5.triangles(self.E)
        self.sdp = rel5.self_dual_perms(self.E)
        self.trules = self._trules()
        low = set()
        for i, t2, t3, t4, t5 in self.trules:
            low |= {i, t4, t5}
        for i, t2, t3, t4, t5 in self.trules:
            assert t2 not in low and t3 not in low, 'no common potential'
        self.low = sorted(low)
        self.stats = {'Y': 0, 'S': 0, 'C': 0, 'P': 0, 'T': 0, 'fail': 0, 'table': 0}

    def _trules(self):
        E = self.E
        verts = sorted({x for e in E for x in e})
        deg = {v: sum(1 for e in E if v in e) for v in verts}
        out = []
        for i, (A0, B0) in enumerate(E):
            for A, B in ((A0, B0), (B0, A0)):
                if deg[A] != 3 or deg[B] != 3:
                    continue
                incA = [t for t, e in enumerate(E) if A in e and t != i]
                incB = [t for t, e in enumerate(E) if B in e and t != i]
                for t3 in incA:
                    C = E[t3][0] if E[t3][1] == A else E[t3][1]
                    for t2 in incB:
                        C2 = E[t2][0] if E[t2][1] == B else E[t2][1]
                        if C2 != C or C in (A, B):
                            continue
                        t5 = [t for t in incA if t != t3][0]
                        t4 = [t for t in incB if t != t2][0]
                        out.append((i, t2, t3, t4, t5))
        return out

    def phi(self, a):
        return sum(a[e] for e in self.low)

    def other(self, edges, pw):
        return per5.value(edges, pw, self.q)

    def val(self, a, depth=0):
        q, E = self.q, self.E
        a = list(a)
        if min(a) < 1 or not gutil.convergent(tuple(E), a, q):
            return None
        rep = gutil.orbit_rep(tuple(a), self.autos)
        if rep in self.memo:
            return self.memo[rep]
        if self.table is not None and rep in self.table:
            self.stats['table'] += 1
            self.memo[rep] = self.table[rep]
            return self.table[rep]
        self.memo[rep] = None            # cycle guard
        res = self._direct(a)
        if res is None:
            res = self._trec(a, depth)
        if res is None:
            self.stats['fail'] += 1
        self.memo[rep] = res
        return res

    def _direct(self, a):
        q, E = self.q, self.E
        for tri in self.tris:
            if a[tri[3]] + a[tri[4]] + a[tri[5]] == 2 * q:
                ne, npw = rel5.tri_to_star(E, a, tri, q)
                v = self.other(ne, npw)
                if v is not None:
                    self.stats['Y'] += 1
                    return v
        for u in startri.unique_vertices(E, a, q):
            try:
                ne, npw = startri.star_to_triangle(E, a, u, q)
            except ValueError:
                continue
            v = self.other(ne, npw)
            if v is not None:
                self.stats['S'] += 1
                return v
        for items in rel5.completion_items(E, a, q):
            ne, npw, c = items[1]
            v = self.other(ne, npw)
            if v is not None:
                self.stats['C'] += 1
                return per5.vscale(v, -c)
        for perm in self.sdp:
            b = [q - a[perm[k]] for k in range(len(a))]
            if min(b) >= 1 and self.phi(b) < self.phi(a):
                rep = gutil.orbit_rep(tuple(b), self.autos)
                if rep in self.memo and self.memo[rep] is not None:
                    self.stats['P'] += 1
                    return self.memo[rep]
        return None

    def _trec(self, a, depth):
        q, E = self.q, self.E
        base = prod(factorial(x - 1) for x in a)
        for i, t2, t3, t4, t5 in self.trules:
            a1, a2, a3 = a[i], a[t2], a[t3]
            c0 = 2 * q - 2 * a1 - a2 - a3
            if c0 == 0:
                continue
            terms = []
            for (x, dx), (y, dy), coef in (((i, -1), (t2, 1), a2), ((t2, 1), (t4, -1), -a2),
                                           ((i, -1), (t3, 1), a3), ((t3, 1), (t5, -1), -a3)):
                b = list(a)
                b[x] += dx
                b[y] += dy
                terms.append((b, coef))
            # p-normalized: c0 p(a) = sum coef p(b);  I = p prod Gamma
            out = {}
            ok = True
            for b, coef in terms:
                if min(b) < 0:
                    ok = False
                    break
                if min(b) == 0:
                    v = self.other(E, b)
                else:
                    v = self.val(b, depth + 1)
                if v is None:
                    ok = False
                    break
                nb = prod(factorial(x - 1) for x in b if x > 0)
                out = per5.vadd(out, v, F(coef) * base / (nb * c0))
            if ok:
                self.stats['T'] += 1
                return out
        return None


if __name__ == '__main__':
    import random
    import time
    name, q = sys.argv[1], int(sys.argv[2])
    tab = per5.table(name, q)
    per5._tables[(name, q)] = None          # hide the family table from per5 (test of the recursion alone)
    R = Rec(name, q, use_table=False)
    rng = random.Random(1)
    reps = sorted(tab) if tab else []
    rng.shuffle(reps)
    t = time.time()
    ok = bad = fail = 0
    for r in reps[:int(sys.argv[3]) if len(sys.argv) > 3 else 300]:
        v = R.val(r)
        if v is None:
            fail += 1
        elif v == tab[r]:
            ok += 1
        else:
            bad += 1
    print('%s q=%d: %d agree with table, %d disagree, %d not reached; stats %s  %.1fs' % (
        name, q, ok, bad, fail, R.stats, time.time() - t))
