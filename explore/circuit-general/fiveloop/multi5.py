"""Joint relation system of one five-loop core over levels q0..q1, linked by the dimension
relation  I_{q-1}(a) = sum_C I_q(a + 1_C)  (C over cotrees).  Saves the table of every level
that closes completely.  usage: python3 multi5.py Q0 Q1 NAME [NAME ...] [--mod]"""
import resource
import sys
import time
from fractions import Fraction as F

import per5
import gutil
import modsolve
import msys5 as msys
import rel5


KMAX = [0]


def frac_solve(rows, nunk, key):
    """Gauss-Jordan over Q, pivot = most complex column (key), rows sorted by length"""
    pivots = {}
    order = sorted(range(len(rows)), key=lambda i: len(rows[i][0]))
    for idx in order:
        row, const = dict(rows[idx][0]), dict(rows[idx][1])
        while True:
            cols = [c for c in row if c in pivots]
            if not cols:
                break
            for c in cols:
                if c not in row:
                    continue
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
        c = max(row, key=key)
        f = row[c]
        row = {k: v / f for k, v in row.items()}
        const = {k: v / f for k, v in const.items()}
        pivots[c] = (row, const)
    # back substitution
    sol = {}
    changed = True
    while changed:
        changed = False
        for c, (row, const) in pivots.items():
            if c in sol:
                continue
            others = [k for k in row if k != c]
            if all(k in sol for k in others):
                val = dict(const)
                for k in others:
                    for l, x in sol[k].items():
                        val[l] = val.get(l, 0) - row[k] * x
                sol[c] = {l: x for l, x in val.items() if x}
                changed = True
    return sol, len(pivots)


def build_multi(names, q0, q1, use_mod=False, verbose=True, save=True):
    t = time.time()
    fams = [per5.SKEL[nm] for nm in names]
    systems, offset, rows = {}, {}, []
    total = 0
    for q in range(q0, q1 + 1):
        S = msys.MultiSystem(fams, q)
        S.lower = per5.lower
        S.build(verbose=False)
        rel5.add_extra(S)
        if KMAX[0]:
            rel5.add_syzygy_rows(S, tuple(range(1, KMAX[0] + 1)))
        systems[q] = S
        offset[q] = total
        for row, const in S.rows:
            rows.append(({offset[q] + j: v for j, v in row.items()}, const))
        total += len(S.reps)
        if verbose:
            print('  level q=%d: unknowns %d rows %d %s  %.1fs' % (q, len(S.reps), len(S.rows), S.counts,
                                                                time.time() - t), flush=True)
    nlink = 0
    for q in range(q0 + 1, q1 + 1):
        lo, hi = systems[q - 1], systems[q]
        for j, (f, rep) in enumerate(lo.reps):
            E = fams[f]
            row = {offset[q - 1] + j: F(1)}
            const = {}
            ok = True
            for C in gutil.cotrees(tuple(E)):
                b = list(rep)
                for e in C:
                    b[e] += 1
                res = hi.classify(list(E), b, q)
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
    for q, S in systems.items():
        for j, (f, rep) in enumerate(S.reps):
            level[offset[q] + j] = (q, f, max(rep), sum(x * x for x in rep), rep)
    if verbose:
        print('  joint: %d unknowns, %d rows (%d links)  %.1fs' % (total, len(rows), nlink, time.time() - t), flush=True)
    if use_mod:
        sol, rank = modsolve.solve(rows, total, rank_key=lambda c: level[c], verbose=verbose)
    else:
        sol, rank = frac_solve(rows, total, key=lambda c: level[c])
    out = {}
    for q, S in systems.items():
        for f, nm in enumerate(names):
            js = [j for j in range(len(S.reps)) if S.reps[j][0] == f]
            tab = {S.reps[j][1]: sol[offset[q] + j] for j in js if offset[q] + j in sol}
            out[(q, nm)] = (len(tab), len(js))
            if save and len(tab) == len(js):
                per5.save_table(nm, q, tab, time.time() - t, extra={'joint_levels': [q0, q1], 'joint': names})
    if verbose:
        print('%s q=%d..%d: rank %d of %d; determined %s  %.1fs maxrss %.2f GB' % (
            '+'.join(names), q0, q1, rank, total, {'%s q=%d' % (k[1], k[0]): '%d/%d' % v for k, v in out.items()}, time.time() - t,
            resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9), flush=True)
    return systems, sol, out


if __name__ == '__main__':
    q0, q1 = int(sys.argv[1]), int(sys.argv[2])
    names = [a for a in sys.argv[3:] if not a.startswith('--')]
    for a in sys.argv:
        if a.startswith('--k='):
            KMAX[0] = int(a[4:])
    build_multi(names, q0, q1, use_mod='--mod' in sys.argv)
