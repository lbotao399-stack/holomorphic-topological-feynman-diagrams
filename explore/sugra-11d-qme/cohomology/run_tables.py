"""Driver: dim H(C_quot) and dim H(V) on all torus weights t (one representative per S_n-orbit, t non-increasing)
for w fixed and G in [Gmin, Gmax].  Exact ranks (jets.Echelon).
usage: python3 run_tables.py n w Gmin Gmax [HC|HV|both] [maxslot] [all|diag|supp] [tag]
  diag: only t = (k,..,k);  supp: only t with t + e_K = (k,..,k) for some K
appends lines to tables_n{n}_w{w}{tag}.log:
  HC G t N dimC dimH
  HV G t N dimC dimV dimH koszulV
  SKIP kind G t total_dim
"""
import sys, time, itertools, json
from cohom import Space, Cohom, all_weights, koszul_dimV

n, w, Gmin, Gmax = (int(x) for x in sys.argv[1:5])
what = sys.argv[5] if len(sys.argv) > 5 else 'both'
maxslot = int(sys.argv[6]) if len(sys.argv) > 6 else 10 ** 9
tfilter = sys.argv[7] if len(sys.argv) > 7 else 'all'
tag = sys.argv[8] if len(sys.argv) > 8 else ''
logf = open('tables_n%d_w%d%s.log' % (n, w, tag), 'a')


def in_supp(t):
    for r in range(n + 1):
        for K in itertools.combinations(range(n), r):
            u = [t[j] + (j in K) for j in range(n)]
            if len(set(u)) == 1:
                return True
    return False


def log(s):
    print(s, flush=True)
    logf.write(s + '\n')
    logf.flush()


S = Space(n, w, Gmax + n + 1)
C = Cohom(S)
log('# run n=%d w=%d G=%d..%d %s maxslot=%d filter=%s  %s' % (n, w, Gmin, Gmax, what, maxslot, tfilter, time.ctime()))
for G in range(Gmin, Gmax + 1):
    t0 = time.time()
    ts = all_weights(S, G)
    reps = sorted({tuple(sorted(t, reverse=True)) for t in ts}, reverse=True)
    if tfilter == 'diag':
        reps = [t for t in reps if len(set(t)) == 1]
    elif tfilter == 'supp':
        reps = [t for t in reps if in_supp(t)]
    for t in reps:
        tot = sum(len(v) for v in S.slot(t).values())
        if what in ('HC', 'both'):
            if tot > maxslot:
                log('SKIP HC %d %s %d' % (G, t, tot))
            else:
                ts0 = time.time()
                h = C.HC(t)
                for N in sorted(h):
                    d, hh = h[N]
                    log('HC %d %s %d %d %d' % (G, t, N, d, hh))
                log('# HC %d %s done %.1fs' % (G, t, time.time() - ts0))
        if what in ('HV', 'both'):
            src = sum(sum(len(v) for v in S.slot(C.shiftw(t, k)).values()) for k in range(n))
            if tot + src > maxslot:
                log('SKIP HV %d %s %d' % (G, t, tot + src))
            else:
                ts0 = time.time()
                h = C.HV(t)
                for N in sorted(h):
                    d, dV, hh = h[N]
                    kz = koszul_dimV(S, t, N)
                    log('HV %d %s %d %d %d %d %d' % (G, t, N, d, dV, hh, kz))
                log('# HV %d %s done %.1fs' % (G, t, time.time() - ts0))
    # free memory of slots two levels below
    for t in list(S.cache):
        if S.G_of(t) < G - n - 1:
            S.drop(t)
    log('# G=%d done %.1fs' % (G, time.time() - t0))
