"""Joint relation system with different family sets per level, dimension links between levels for
families present at both, and the q=2 link I_2(G; 1) = sum_C I_3(G; 1 + 1_C) for families whose
all-ones vector appears as an unknown (B relations of larger cores at level 2).
usage: python3 joint5.py 'Q:NAME,NAME;Q:NAME' [--k=K] [--z5q2]   (--z5q2: add P(Z5) = 441/8 zeta7 as input)"""
import resource
import sys
import time
from fractions import Fraction as F

import per5
import gutil
import modsolve
import msys5 as msys
import multi5
import rel5


def run(spec, kmax=2, z5q2=False, verbose=True, save=True, extra_known=None):
    t = time.time()
    levels = []
    for part in spec.split(';'):
        q, names = part.split(':')
        levels.append((int(q), names.split(',')))
    systems, offset, rows, total = {}, {}, [], 0
    for q, names in levels:
        S = msys.MultiSystem([per5.SKEL[nm] for nm in names], q)
        S.lower = per5.lower
        S.build(verbose=False)
        rel5.add_extra(S)
        if kmax:
            rel5.add_syzygy_rows(S, tuple(range(1, kmax + 1)))
        systems[q] = (S, names)
        offset[q] = total
        for row, const in S.rows:
            rows.append(({offset[q] + j: v for j, v in row.items()}, const))
        if verbose:
            print('  level q=%d %s: unknowns %d rows %d %s  %.1fs' % (q, names, len(S.reps), len(S.rows), S.counts,
                                                                    time.time() - t), flush=True)
        total += len(S.reps)
    nlink = 0
    qs = sorted(systems)
    for q in qs:
        S, names = systems[q]
        # q-1 -> q links inside the same system (unknowns with lower total, e.g. all-ones at level 2)
        for j, (f, rep) in enumerate(list(S.reps)):
            E = S.fam[f]['edges']
            n = S.fam[f]['n']
            lev = sum(rep) // n
            if lev == q:
                continue
            row = {offset[q] + j: F(1)}
            const = {}
            ok = True
            for C in gutil.cotrees(tuple(E)):
                b = list(rep)
                for e in C:
                    b[e] += 1
                res = S.classify(list(E), b, lev + 1) if lev + 1 == q else None
                if res is None:
                    ok = False
                    break
                c, unk = res
                for k, v in c.items():
                    const[k] = const.get(k, 0) + v
                for jj, s in unk:
                    row[offset[q] + jj] = row.get(offset[q] + jj, 0) - s
            if ok:
                rows.append(({k: v for k, v in row.items() if v}, {k: v for k, v in const.items() if v}))
                nlink += 1
            if z5q2 and names[f] == 'Z5' and lev == 2 and all(x == 1 for x in rep):
                rows.append(({offset[q] + j: F(1)}, {'z7': F(441, 8)}))
        if q - 1 in systems:
            lo, lnames = systems[q - 1]
            for j, (f, rep) in enumerate(lo.reps):
                nm = lnames[f]
                if nm not in names:
                    continue
                if sum(rep) != lo.fam[f]['n'] * (q - 1):
                    continue
                E = lo.fam[f]['edges']
                row = {offset[q - 1] + j: F(1)}
                const = {}
                ok = True
                for C in gutil.cotrees(tuple(E)):
                    b = list(rep)
                    for e in C:
                        b[e] += 1
                    res = S.classify(list(E), b, q)
                    if res is None:
                        ok = False
                        break
                    c, unk = res
                    for k, v in c.items():
                        const[k] = const.get(k, 0) + v
                    for jj, s in unk:
                        row[offset[q] + jj] = row.get(offset[q] + jj, 0) - s
                if ok:
                    rows.append(({k: v for k, v in row.items() if v}, {k: v for k, v in const.items() if v}))
                    nlink += 1
    level = {}
    for q, (S, names) in systems.items():
        for j, (f, rep) in enumerate(S.reps):
            level[offset[q] + j] = (q, f, max(rep), sum(x * x for x in rep), rep)
    if verbose:
        print('  joint: %d unknowns, %d rows (%d links)  %.1fs' % (total, len(rows), nlink, time.time() - t), flush=True)
    sol, rank = multi5.frac_solve(rows, total, key=lambda c: level[c])
    out = {}
    for q, (S, names) in systems.items():
        for f, nm in enumerate(names):
            js = [j for j in range(len(S.reps)) if S.reps[j][0] == f and sum(S.reps[j][1]) == S.fam[f]['n'] * q]
            tab = {S.reps[j][1]: sol[offset[q] + j] for j in js if offset[q] + j in sol}
            out[(q, nm)] = (len(tab), len(js))
            if save and len(tab) == len(js) and per5.table(nm, q) is None:
                per5.save_table(nm, q, tab, time.time() - t, extra={'joint': spec, 'z5q2_input': z5q2})
        for j, (f, rep) in enumerate(S.reps):
            if sum(rep) != S.fam[f]['n'] * q and offset[q] + j in sol:
                print('  lower-level unknown %s %s = %s' % (names[f], rep, per5.fmt(sol[offset[q] + j])), flush=True)
    if verbose:
        print('%s: rank %d of %d; determined %s  %.1fs maxrss %.2f GB' % (
            spec, rank, total, {'%s q=%d' % (k[1], k[0]): '%d/%d' % v for k, v in out.items()}, time.time() - t,
            resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9), flush=True)
    return systems, sol, out, offset


if __name__ == '__main__':
    kmax = 2
    for a in sys.argv:
        if a.startswith('--k='):
            kmax = int(a[4:])
    run(sys.argv[1], kmax=kmax, z5q2='--z5q2' in sys.argv, save='--nosave' not in sys.argv)
