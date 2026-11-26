"""Reads tables_n{n}_w{w}*.log and produces: lemma comparison (T4), nonzero tables and SL_n-invariant parts (T6).
usage: python3 analyze.py n w  -> analysis_n{n}_w{w}.log"""
import sys, glob, itertools, re
from collections import defaultdict

n, w = int(sys.argv[1]), int(sys.argv[2])
out = open('analysis_n%d_w%d.log' % (n, w), 'w')


def log(s):
    print(s)
    out.write(s + '\n')


HC, HV, SKIP = {}, {}, []
for fn in sorted(glob.glob('tables_n%d_w%d*.log' % (n, w))):
    for line in open(fn):
        if line.startswith('HC '):
            m = re.match(r'HC (-?\d+) \((.*?)\) (-?\d+) (\d+) (-?\d+)', line)
            G = int(m.group(1)); t = tuple(int(x) for x in m.group(2).split(',') if x.strip()); N = int(m.group(3))
            HC[(t, N)] = (int(m.group(4)), int(m.group(5)))
        elif line.startswith('HV '):
            m = re.match(r'HV (-?\d+) \((.*?)\) (-?\d+) (\d+) (-?\d+) (-?\d+) (-?\d+)', line)
            G = int(m.group(1)); t = tuple(int(x) for x in m.group(2).split(',') if x.strip()); N = int(m.group(3))
            HV[(t, N)] = (int(m.group(4)), int(m.group(5)), int(m.group(6)), int(m.group(7)))
        elif line.startswith('SKIP'):
            SKIP.append(line.strip())


def canon(t):
    return tuple(sorted(t, reverse=True))


def G_of(t):
    return -w - sum(t)


tsC = {canon(t) for (t, N) in HC}
tsV = {canon(t) for (t, N) in HV}


GsC = {G_of(t) for t in tsC}
GminC = min(GsC) if GsC else 0


def hc(t, N):
    """returns (value, status): status 'chk' computed, 'cartan' zero by the Cartan identity (t not prop. to 1), None unknown."""
    tc = canon(t)
    if (tc, N) in HC:
        return HC[(tc, N)][1], 'chk'
    if tc in tsC or G_of(t) in GsC or G_of(t) < GminC:
        return 0, 'chk'   # slot computed (N outside range) or weight absent at a fully enumerated G => dim C = 0
    if len(set(t)) > 1:
        return 0, 'cartan'
    return None, None


log('n=%d w=%d: %d HC entries on %d weights, %d HV entries on %d weights' % (n, w, len(HC), len(tsC), len(HV), len(tsV)))
if SKIP:
    log('skipped: ' + '; '.join(SKIP))

# Cartan support check
off = [(t, N, v) for (t, N), (d, v) in HC.items() if v and len(set(t)) > 1]
log('H(C_quot) nonzero off t=(k,..,k): %s' % (off if off else 'none'))

# Koszul check
if HV:
    bad = [(t, N) for (t, N), (d, dV, h, kz) in HV.items() if dV != kz]
    log('dim V == Koszul alternating sum on all HV slots: %s%s' % (not bad, '' if not bad else ' ; failures: %s' % sorted(bad)[:20]))

# lemma comparison
log('translation lemma: H(V)_{t,N} =?= sum_K H(C)_{t+e_K, N+|K|}')
nchk = ncart = nunk = 0
fails = []
for (t, N), (d, dV, h, kz) in sorted(HV.items(), key=lambda kv: (G_of(kv[0][0]), kv[0])):
    pred = 0
    status = 'chk'
    for r in range(n + 1):
        for K in itertools.combinations(range(n), r):
            tt = tuple(t[j] + (j in K) for j in range(n))
            v, st = hc(tt, N + r)
            if v is None:
                status = None
                break
            pred += v
            if st == 'cartan':
                status = 'cartan' if status == 'chk' else status
        if status is None:
            break
    if status is None:
        nunk += 1
        continue
    if status == 'chk':
        nchk += 1
    else:
        ncart += 1
    if pred != h:
        fails.append((G_of(t), t, N, h, pred, status))
log('  slots compared: %d with all H(C) inputs computed, %d using H(C)=0 off-diagonal (Cartan), %d not comparable' % (nchk, ncart, nunk))
log('  mismatches: %s' % (fails if fails else 'none'))

# tables
log('nonzero H(C_quot) (G, t, N: dim H, dim C):')
for (t, N), (d, v) in sorted(HC.items(), key=lambda kv: (G_of(kv[0][0]), kv[0])):
    if v:
        log('  G=%3d t=%s N=%2d (%s): H=%d  (dim C=%d)' % (G_of(t), t, N, 'odd' if N % 2 else 'even', v, d))
log('nonzero H(V) (G, t up to S_n, N: dim H, dim V):')
for (t, N), (d, dV, h, kz) in sorted(HV.items(), key=lambda kv: (G_of(kv[0][0]), kv[0])):
    if h:
        log('  G=%3d t=%s N=%2d (%s): H=%d  (dim V=%d, dim C=%d)' % (G_of(t), t, N, 'odd' if N % 2 else 'even', h, dV, d))

# SL_n invariants of H(V): multiplicity of det^k = sum_s sgn(s) m(lam + rho - s rho)
rho = tuple(n - 1 - j for j in range(n))


def sgn(p):
    s = 1
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]:
                s = -s
    return s


def hv(t, N):
    tc = canon(t)
    if (tc, N) in HV:
        return HV[(tc, N)][2]
    if tc in tsV:
        return 0
    return None


log('SL_n-invariant part of H(V) (weight (k,..,k) highest-weight multiplicity, Weyl formula):')
for k in range(-8, 3):
    lam = (k,) * n
    if canon(lam) not in tsV:
        continue
    Ns = sorted(N for (t, N) in HV if t == lam)
    for N in Ns:
        tot = 0
        okk = True
        for p in itertools.permutations(range(n)):
            mu = tuple(lam[j] + rho[j] - rho[p[j]] for j in range(n))
            v = hv(mu, N)
            if v is None:
                if len(set(mu)) > 0 and sum(mu) == sum(lam):
                    v = None
                okk = False
                break
            tot += sgn(p) * v
        if not okk:
            log('  t=%s N=%d: incomplete' % (lam, N))
        elif tot:
            log('  G=%d t=%s N=%d: %d   [H(V) at t: %d]' % (G_of(lam), lam, N, tot, hv(lam, N)))
