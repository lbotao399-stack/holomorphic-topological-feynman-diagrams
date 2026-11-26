"""build5 + (L) syzygy relations of degrees 1, 2.  usage: python3 build5L.py Q NAME [NAME ...]"""
import resource
import sys
import time

import per5
import msys5 as msys
import modsolve
import multi5
import rel5


def build(names, q, degrees=(1, 2), verbose=True, save=True, use_mod=None):
    t = time.time()
    fams = [per5.SKEL[nm] for nm in names]
    S = msys.MultiSystem(fams, q)
    S.lower = per5.lower
    S.build(verbose=False)
    rel5.add_extra(S)
    rel5.add_syzygy_rows(S, degrees)
    nunk = len(S.reps)
    if verbose:
        print('%s q=%d unknowns %d rows %d %s  %.1fs' % ('+'.join(names), q, nunk, len(S.rows), S.counts,
                                                         time.time() - t), flush=True)
    level = {j: (S.reps[j][0], max(S.reps[j][1]), sum(x * x for x in S.reps[j][1]), S.reps[j][1]) for j in range(nunk)}
    if use_mod is None:
        use_mod = nunk > 20000
    if use_mod:
        solved, rank = modsolve.solve(S.rows, nunk, rank_key=lambda c: level[c], verbose=verbose)
    else:
        solved, rank = multi5.frac_solve(S.rows, nunk, key=lambda c: level[c])
    res = {}
    for f, nm in enumerate(names):
        js = [j for j in range(nunk) if S.reps[j][0] == f]
        tab = {S.reps[j][1]: solved[j] for j in js if j in solved}
        res[nm] = (len(tab), len(js))
        if save and len(tab) == len(js):
            per5.save_table(nm, q, tab, time.time() - t, extra={'rows': len(S.rows), 'counts': S.counts, 'joint': names})
    if verbose:
        print('q=%d: rank %d of %d; determined %s  %.1fs  maxrss %.2f GB' % (
            q, rank, nunk, {k: '%d/%d' % v for k, v in res.items()}, time.time() - t,
            resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9), flush=True)
    return S, solved, res


if __name__ == '__main__':
    q = int(sys.argv[1])
    names = [a for a in sys.argv[2:] if not a.startswith('--')]
    kmax = 3
    for a in sys.argv:
        if a.startswith('--k='):
            kmax = int(a[4:])
    build(names, q, degrees=tuple(range(1, kmax + 1)))
