"""Level-q tables from the single-level system (all relation types, (L) up to degree K) plus the
dimension relation with known lower table as rows:  sum_C I_q(a + 1_C) = I_{q-1}(a)  (a over the level q-1
table).  usage: python3 buildup.py Q NAME [NAME ...] [--k=K]"""
import resource
import sys
import time
from fractions import Fraction as F

import per5
import gutil
import msys5 as msys
import multi5
import rel5


def build(names, q, kmax=2, verbose=True, save=True):
    t = time.time()
    fams = [per5.SKEL[nm] for nm in names]
    S = msys.MultiSystem(fams, q)
    S.lower = per5.lower
    S.build(verbose=False)
    rel5.add_extra(S)
    if kmax:
        rel5.add_syzygy_rows(S, tuple(range(1, kmax + 1)))
    S.counts['R'] = 0
    for f, nm in enumerate(names):
        low = per5.table(nm, q - 1)
        if low is None:
            continue
        E = list(fams[f])
        cot = gutil.cotrees(tuple(E))
        for a, val in low.items():
            # row . x = I_{q-1}(a) - (known parts of sum_C I_q(a + 1_C))
            row, const = {}, dict(val)
            ok = True
            for C in cot:
                b = list(a)
                for e in C:
                    b[e] += 1
                res = S.classify(E, b, q)
                if res is None:
                    ok = False
                    break
                c, unk = res
                for k, v in c.items():
                    const[k] = const.get(k, 0) - v
                for j, s in unk:
                    row[j] = row.get(j, 0) + s
            if ok:
                S.rows.append(({k: v for k, v in row.items() if v}, {k: v for k, v in const.items() if v}))
                S.counts['R'] += 1
    nunk = len(S.reps)
    if verbose:
        print('%s q=%d unknowns %d rows %d %s  %.1fs' % ('+'.join(names), q, nunk, len(S.rows), S.counts,
                                                         time.time() - t), flush=True)
    level = {j: (S.reps[j][0], max(S.reps[j][1]), sum(x * x for x in S.reps[j][1]), S.reps[j][1]) for j in range(nunk)}
    solved, rank = multi5.frac_solve(S.rows, nunk, key=lambda c: level[c])
    res = {}
    for f, nm in enumerate(names):
        js = [j for j in range(nunk) if S.reps[j][0] == f and sum(S.reps[j][1]) == S.fam[f]['n'] * q]
        tab = {S.reps[j][1]: solved[j] for j in js if j in solved}
        res[nm] = (len(tab), len(js))
        if save and len(tab) == len(js) and per5.table(nm, q) is None:
            per5.save_table(nm, q, tab, time.time() - t, extra={'rows': len(S.rows), 'counts': S.counts, 'joint': names,
                                                                 'method': 'buildup'})
    if verbose:
        print('q=%d: rank %d of %d; determined %s  %.1fs  maxrss %.2f GB' % (
            q, rank, nunk, {k: '%d/%d' % v for k, v in res.items()}, time.time() - t,
            resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9), flush=True)
    return S, solved, res


if __name__ == '__main__':
    q = int(sys.argv[1])
    names = [a for a in sys.argv[2:] if not a.startswith('--')]
    kmax = 2
    for a in sys.argv:
        if a.startswith('--k='):
            kmax = int(a[4:])
    build(names, q, kmax)
