"""Build period tables of five-loop cores at level q from one joint regulator-free relation system
over the listed cores (msys.MultiSystem: (T) triangle, (D) diamond, (I) parametric IBP, (B)
contraction boundary, (S) star-triangle; rel5: (Y) triangle-star, (X) position-space vertex IBP,
(P) planar self-duality, (C) conformal completion), boundary data from per5.lower (2-sums,
K5-e = dual prism, prism/K33 tables, W_n, lower five-loop tables).
usage: python3 build5.py Q NAME [NAME ...] [--mod]"""
import resource
import sys
import time

import per5
import msys5 as msys
import modsolve
import rel5


def build(names, q, use_mod=None, verbose=True, save=True):
    t = time.time()
    fams = [per5.SKEL[nm] for nm in names]
    S = msys.MultiSystem(fams, q)
    S.lower = per5.lower
    S.build(verbose=False)
    rel5.add_extra(S)
    if verbose:
        print('%s q=%d unknowns %d rows %d %s  %.1fs' % ('+'.join(names), q, len(S.reps), len(S.rows), S.counts,
                                                         time.time() - t), flush=True)
    nunk = len(S.reps)
    if use_mod is None:
        use_mod = nunk > 6000
    level = {j: (S.reps[j][0], max(S.reps[j][1]), sum(x * x for x in S.reps[j][1]), S.reps[j][1]) for j in range(nunk)}
    if use_mod:
        solved, rank = modsolve.solve(S.rows, nunk, rank_key=lambda c: level[c], verbose=verbose)
    else:
        import multi5
        solved, rank = multi5.frac_solve(S.rows, nunk, key=lambda c: level[c])
    res = {}
    for f, nm in enumerate(names):
        js = [j for j in range(nunk) if S.reps[j][0] == f]
        tab = {S.reps[j][1]: solved[j] for j in js if j in solved}
        res[nm] = (len(tab), len(js))
        if save and len(tab) == len(js):
            per5.save_table(nm, q, tab, time.time() - t, extra={'rows': len(S.rows), 'counts': S.counts,
                                                                 'joint': names})
    if verbose:
        print('q=%d: rank %d of %d; determined %s  %.1fs  maxrss %.2f GB' % (
            q, rank, nunk, {k: '%d/%d' % v for k, v in res.items()}, time.time() - t,
            resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9), flush=True)
    return S, solved, res


if __name__ == '__main__':
    q = int(sys.argv[1])
    names = [a for a in sys.argv[2:] if not a.startswith('--')]
    build(names, q, True if '--mod' in sys.argv else None)
