"""Lean builder for prism / K33 period tables (subclass of cores/msys.MultiSystem, unchanged relations):
optionally without the parametric IBP rows (I), which carry up to |cotrees| terms each and dominate
memory; S (star-triangle), T (triangle) / D (diamond), B (boundary, level q-1) are kept.
  python3 lean_tables.py NAME Q [--noI] [--check]   (--check: compare with an existing table)"""
import json
import os
import resource
import sys
import time
from fractions import Fraction as F

import periods
import gutil
import msys
import startri


class LeanSystem(msys.MultiSystem):
    def __init__(self, families, q, use_I=True):
        super().__init__(families, q)
        self.use_I = use_I

    def build(self, verbose=True):
        q = self.q
        t0 = time.time()
        level_reps = []
        for f, d in enumerate(self.fam):
            N = len(d['edges'])
            reps = set()
            amax = q - 1
            for c in gutil.compositions(int(d['n'] * q), N, 1, amax):
                if gutil.convergent(d['edges'], list(c), q):
                    reps.add(gutil.orbit_rep(c, d['autos']))
            reps = sorted(reps)
            for r in reps:
                self.var(f, r)
            level_reps.append(reps)
        for f, reps in enumerate(level_reps):
            edges = self.fam[f]['edges']
            for r in reps:
                for u in startri.unique_vertices(edges, list(r), q):
                    ne, npw = startri.star_to_triangle(edges, list(r), u, q)
                    self.add_relation([(edges, list(r), 1), (ne, npw, -1)], 'S')
        for f, reps in enumerate(level_reps):
            for r in reps:
                self.triangle_rows(f, list(r))
                if self.use_diamond:
                    self.diamond_rows(f, list(r))
        self.class_cache.clear()
        for f, d in enumerate(self.fam):
            edges, N, n, cot = d['edges'], len(d['edges']), d['n'], d['cot']
            bases = set()
            for c in gutil.compositions(int(n * (q - 1)), N, 0, q - 1):
                bases.add(gutil.orbit_rep(c, d['autos']))
            for b in bases:
                zeros = [i for i in range(N) if b[i] == 0]
                if len(zeros) > 1:
                    continue
                if len(zeros) == 1:
                    i = zeros[0]
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
                if not self.use_I:
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
        self.class_cache.clear()
        if verbose:
            print('q=%d unknowns %d rows %d %s  %.1fs  maxrss %.2f GB' % (q, len(self.reps), len(self.rows), self.counts,
                  time.time() - t0, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9), flush=True)


def main(argv):
    name, q = argv[0], int(argv[1])
    use_I = '--noI' not in argv
    t = time.time()
    S = LeanSystem([periods.FAMILIES[name]], q, use_I=use_I)
    S.build()
    solved, piv = S.solve()
    print('%s q=%d: determined %d of %d (%.1fs, maxrss %.2f GB)' % (name, q, len(solved), len(S.reps), time.time() - t,
          resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9), flush=True)
    tab = {S.reps[j][1]: v for j, v in solved.items()}
    if '--check' in argv:
        ref = periods.table(name, q)
        bad = sum(1 for r, v in tab.items() if ref.get(r) != v)
        print('  compared with existing table: %d differences' % bad, flush=True)
    elif len(solved) == len(S.reps):
        path = os.path.join(periods.TABLES_OUT, '%s_q%d.json' % (name, q))
        if not os.path.exists(path):
            json.dump({'core': name, 'edges': periods.FAMILIES[name], 'q': str(q), 'seconds': round(time.time() - t, 1),
                       'relations': 'lean (no I)' if not use_I else 'all',
                       'values': [[list(r), {k: str(x) for k, x in v.items()}] for r, v in sorted(tab.items())]},
                      open(path, 'w'))
            print('  written', path, flush=True)


if __name__ == '__main__':
    main(sys.argv[1:])
