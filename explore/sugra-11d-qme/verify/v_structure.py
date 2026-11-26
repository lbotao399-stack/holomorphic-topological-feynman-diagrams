"""Claims a1, a2, a3: independent enumeration by transfer matrices along the three paths of theta_pqr.

Edge label = (field at e(0), field at e(1)) in {(mu,ga),(ga,mu),(nu,be),(be,nu)}   (S3 eq-Delta).
Vertex types from W (S3 R-action, f(nu) = sum_k nu^k):
   a   : beta d_i mu^i                  letters [be, mu]
   b_k : nu^k mu^i mu^j d_i gamma_j     letters [ga, mu, mu] + [nu]*k
   c   : (g/6) eps gamma d gamma d gamma letters [ga, ga, ga]
Weight of a vertex with given port fields = number of injective letter maps (letters distinguishable),
monomial in (m = #c, n_a = #a, external mu, ga, nu, be).
Rule 'Dtheta': every interior (degree-2) vertex keeps >= 1 external letter (S3 eq-insertion, eq-Dtheta).
"""
import sys, re, os
from collections import Counter
from itertools import product
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vcommon import theta_graph, PQR

LABELS = [('mu', 'ga'), ('ga', 'mu'), ('nu', 'be'), ('be', 'nu')]
FIELDS = ['mu', 'ga', 'nu', 'be']


def vtypes(kmax):
    T = [('a', ['be', 'mu'])]
    for k in range(kmax + 1):
        T.append(('b%d' % k, ['ga', 'mu', 'mu'] + ['nu'] * k))
    T.append(('c', ['ga', 'ga', 'ga']))
    return T


def falling(n, k):
    r = 1
    for i in range(k):
        r *= (n - i)
    return r


def vweight(ports, interior, kmax, rule):
    out = Counter()
    pc = Counter(ports)
    for name, L in vtypes(kmax):
        lc = Counter(L)
        mult = 1
        for f in FIELDS:
            mult *= falling(lc[f], pc[f])
        if mult == 0:
            continue
        ext = tuple(lc[f] - pc[f] for f in FIELDS)
        if interior and rule and sum(ext) == 0:
            continue
        key = (1 if name == 'c' else 0, 1 if name == 'a' else 0) + ext
        out[key] += mult
    return out


def pmul(A, B):
    C = Counter()
    for ka, va in A.items():
        for kb, vb in B.items():
            C[tuple(x + y for x, y in zip(ka, kb))] += va * vb
    return C


ONE = Counter({(0, 0, 0, 0, 0, 0): 1})


def path_transfer(n_edges, kmax, rule):
    """T[(l_first, l_last)] = Counter over the interior vertices of a path with n_edges edges."""
    T = {}
    if n_edges == 1:
        for l in range(4):
            T[(l, l)] = ONE
        return T
    # state: (l_first, l_current) -> Counter
    S = {(l, l): ONE for l in range(4)}
    for j in range(n_edges - 1):
        S2 = {}
        for (lf, lc), P in S.items():
            for ln in range(4):
                w = vweight([LABELS[lc][1], LABELS[ln][0]], True, kmax, rule)
                if not w:
                    continue
                Q = pmul(P, w)
                key = (lf, ln)
                S2[key] = S2.get(key, Counter()) + Q
        S = S2
    return S


def enumerate_theta(pqr, kmax, rule):
    Ts = [path_transfer(n, kmax, rule) for n in pqr]
    tot = Counter()
    keys = [list(T.keys()) for T in Ts]
    for kA, kB, kC in product(*keys):
        portsN = [LABELS[kA[0]][0], LABELS[kB[0]][0], LABELS[kC[0]][0]]
        portsNp = [LABELS[kA[1]][1], LABELS[kB[1]][1], LABELS[kC[1]][1]]
        wN = vweight(portsN, False, kmax, rule)
        if not wN:
            continue
        wNp = vweight(portsNp, False, kmax, rule)
        if not wNp:
            continue
        P = pmul(pmul(wN, wNp), pmul(pmul(Ts[0][kA], Ts[1][kB]), Ts[2][kC]))
        tot += P
    return tot


def read_explorer_log(path):
    """Parse the explorer's brute-force lines of structure.log (data only)."""
    ref = {}
    pat = re.compile(r"^\((\d+), (\d+), (\d+)\) m = (\d+) n_a = (\d+) external (\{.*\}) maps: (\d+)")
    for line in open(path):
        mm = pat.match(line.strip())
        if mm:
            pqr = tuple(int(mm.group(i)) for i in (1, 2, 3))
            ext = eval(mm.group(6), {})
            key = (int(mm.group(4)), int(mm.group(5)), ext.get('mu', 0), ext.get('gamma', 0), ext.get('nu', 0), 0)
            ref[(pqr, key)] = int(mm.group(7))
    return ref


if __name__ == '__main__':
    here = os.path.dirname(os.path.abspath(__file__))
    ref = read_explorer_log(os.path.join(here, '..', 'structure.log'))
    print('key = (m, n_a, ext mu, ext gamma, ext nu, ext beta)')
    for pqr in PQR:
        nv, edges, paths = theta_graph(*pqr)
        # (i) no Dtheta rule, kmax = 1: all admissible letter maps
        R0 = enumerate_theta(pqr, 1, False)
        m0 = sum(v for k, v in R0.items() if k[0] == 0)
        na_pos = sum(v for k, v in R0.items() if k[1] > 0)
        print(pqr, 'no Dtheta rule, kmax=1: total maps', sum(R0.values()), '; with m=0:', m0,
              '; with n_a>0:', na_pos, '; min m =', min(k[0] for k in R0))
        for kmax in (1, 2):
            R = enumerate_theta(pqr, kmax, True)
            m0 = sum(v for k, v in R.items() if k[0] == 0)
            na_pos = sum(v for k, v in R.items() if k[1] > 0)
            beta_ext = sum(v for k, v in R.items() if k[5] > 0)
            bad = [k for k in R if not (k[2] == 11 - 2 * k[0] and k[3] == 2 * k[0] - 1 and k[5] == 0)]
            ms = sorted(set(k[0] for k in R))
            print(pqr, 'Dtheta rule, kmax=%d:' % kmax, 'maps', sum(R.values()), '; m=0:', m0, '; n_a>0:', na_pos,
                  '; ext beta>0:', beta_ext, '; m values', ms,
                  '; keys violating (mu,ga)=(11-2m,2m-1):', len(bad),
                  '; max K by m:', {m: max(k[4] for k in R if k[0] == m) for m in ms})
            if kmax == 1:
                K0 = {k[0]: v for k, v in R.items() if k[4] == 0}
                print('   K=0 maps by m:', dict(sorted(K0.items())))
                agree = [R.get(key, 0) == val for (pq, key), val in ref.items() if pq == pqr]
                if agree:
                    extra = [k for k in R if (pqr, k) not in ref]
                    print('   comparison with structure.log lines for this graph:', len(agree), 'lines,',
                          'all equal' if all(agree) else 'MISMATCH', '; keys of mine not in log:', len(extra))
