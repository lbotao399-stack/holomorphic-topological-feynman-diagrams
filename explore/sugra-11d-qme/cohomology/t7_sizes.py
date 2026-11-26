"""T7: exact dimensions of C_quot and V slots by characters (no enumeration).

chi_X = (sum_i z_i - 1) prod_j (1 - z_j^{-1})^{-1}          (X jets mod divergence)
chi_g = (sum_i z_i^{-1}) prod_j (1 - z_j^{-1})^{-1},  chi_b = prod_j (1 - z_j^{-1})^{-1}
Lambda^a(chi_X) = e_a[chi_X] = sum_{lambda |- a} eps_lambda / z_lambda prod_r psi_{lambda_r}(chi_X),
psi_k(chi_X) = (sum_i z_i^k - 1) prod_j (1 - z_j^{-k})^{-1}.
[z^t] e_a[chi_X] * chi_w = sum_lambda eps/z sum_{bins} sign * prod_j p_{lambda (+1)}(c_j - s_j - t_j),
p_mu(b) = #{beta in N^len(mu): sum_r mu_r beta_r = b}.
dim C[#X=a, gamma] at t, dim C[#X=a, beta] at t;  N = a + #beta.
dim V_t = sum_K (-1)^|K| dim C_{t+e_K}  (Koszul; exact because C_quot^{w=1} is free over C[D_1..D_n]).
Validated against explicit enumeration (n = 2, 3).  Writes t7.log."""
import sys, itertools, time
from fractions import Fraction
from math import factorial, comb
from functools import lru_cache
from collections import defaultdict

out = None


def log(s):
    print(s, flush=True)
    if out is not None:
        out.write(s + '\n')
        out.flush()


def partitions(a, maxpart=None):
    if maxpart is None:
        maxpart = a
    if a == 0:
        yield ()
        return
    for k in range(min(a, maxpart), 0, -1):
        for rest in partitions(a - k, k):
            yield (k,) + rest


@lru_cache(maxsize=None)
def pcount(mu, b):
    """# nonnegative solutions of sum_r mu_r beta_r = b, mu a sorted tuple."""
    if b < 0:
        return 0
    if not mu:
        return 1 if b == 0 else 0
    ways = [0] * (b + 1)
    ways[0] = 1
    for k in mu:
        for s in range(k, b + 1):
            ways[s] += ways[s - k]
    return ways[b]


@lru_cache(maxsize=None)
def bins(lam, n):
    """dict (c tuple in N^n, sign) -> multinomial weight, distributing the parts of lam over bins 0 ('-1') and 1..n."""
    mult = defaultdict(int)
    for k in lam:
        mult[k] += 1
    cur = {((0,) * n, 1): 1}
    for k, m in mult.items():
        new = defaultdict(int)
        # compositions of m into n+1 bins
        for comp in itertools.product(range(m + 1), repeat=n):
            s = sum(comp)
            if s > m:
                continue
            n0 = m - s
            w = factorial(m)
            for x in comp:
                w //= factorial(x)
            w //= factorial(n0)
            sg = -1 if n0 % 2 else 1
            for (c, sgn), v in cur.items():
                nc = tuple(c[j] + k * comp[j] for j in range(n))
                new[(nc, sgn * sg)] += v * w
        cur = dict(new)
    return cur


def zlam(lam):
    mult = defaultdict(int)
    for k in lam:
        mult[k] += 1
    z = 1
    for k, m in mult.items():
        z *= k ** m * factorial(m)
    return z


@lru_cache(maxsize=None)
def dimC(n, a, wtype, t):
    """dim of C_quot with #X = a, w-letter wtype in {None, 'g', 'b'}, torus weight t."""
    tot = Fraction(0)
    for lam in partitions(a):
        eps = -1 if (a - len(lam)) % 2 else 1
        mu = tuple(sorted(lam + ((1,) if wtype else ())))
        s = 0
        for (c, sg), wgt in bins(lam, n).items():
            if wtype == 'g':
                for ip in range(n):
                    pr = 1
                    for j in range(n):
                        pr *= pcount(mu, c[j] - (j == ip) - t[j])
                        if not pr:
                            break
                    s += sg * wgt * pr
            else:
                pr = 1
                for j in range(n):
                    pr *= pcount(mu, c[j] - t[j])
                    if not pr:
                        break
                s += sg * wgt * pr
        tot += Fraction(eps * s, zlam(lam))
    assert tot.denominator == 1
    return int(tot)


def dimC_N(n, w, t, N):
    """dim C^N_t split by content; w = 0 or 1."""
    if w == 0:
        return {(N, 0, 0): dimC(n, N, None, t)}
    return {(N, 1, 0): dimC(n, N, 'g', t), (N - 1, 0, 1): dimC(n, N - 1, 'b', t) if N >= 1 else 0}


def dimV_N(n, w, t, N):
    res = defaultdict(int)
    for r in range(n + 1):
        for K in itertools.combinations(range(n), r):
            tt = tuple(t[j] + (j in K) for j in range(n))
            for key, v in dimC_N(n, w, tt, N).items():
                res[key] += (-1) ** r * v
    return dict(res)


def weyl_mult(n, w, lam_t, N, which='C'):
    """multiplicity of the irreducible gl_n module of highest weight lam_t in C^N (or V^N) at fixed (w,G):
    sum_{s in S_n} sgn(s) m(lam + rho - s rho)."""
    rho = tuple(n - 1 - j for j in range(n))
    tot = 0
    for perm in itertools.permutations(range(n)):
        sg = 1
        for i in range(n):
            for j in range(i + 1, n):
                if perm[i] > perm[j]:
                    sg = -sg
        mu = tuple(lam_t[j] + rho[j] - rho[perm[j]] for j in range(n))
        d = dimC_N(n, w, mu, N) if which == 'C' else dimV_N(n, w, mu, N)
        tot += sg * sum(d.values())
    return tot


def validate():
    sys.path.insert(0, '.')
    from cohom import Space, koszul_dimV
    ok = True
    cnt = 0
    for n, w, Gs in ((2, 1, range(-3, 6)), (2, 0, range(-2, 6)), (3, 1, range(-4, 0)), (3, 0, range(-3, 2))):
        S = Space(n, w, max(Gs) + n + 1)
        from cohom import all_weights
        for G in Gs:
            for t in all_weights(S, G):
                sl = S.slot(t)
                for N in set(sl) | {max(sl) + 1 if sl else 0}:
                    d = len(sl.get(N, []))
                    if w == 1:
                        dd = dimC(n, N, 'g', t) + (dimC(n, N - 1, 'b', t) if N >= 1 else 0)
                    else:
                        dd = dimC(n, N, None, t)
                    cnt += 1
                    if d != dd:
                        ok = False
                        log('  MISMATCH n=%d w=%d t=%s N=%d enum=%d char=%d' % (n, w, t, N, d, dd))
            S.cache.clear()
    return ok, cnt


if __name__ == '__main__':
    out = open('t7.log', 'w')
    t0 = time.time()
    log('T7 (t7_sizes.py)')
    ok, cnt = validate()
    log('validation: character formula == explicit enumeration on %d slots (n=2: w=1 G<=5, w=0 G<=5; n=3: w=1 G<=-1, w=0 G<=1): %s' % (cnt, ok))
    n = 5
    t = (-3,) * 5
    log('n=5, w=1, t=(-3)^5, G=14.  dim C_quot^N_t by content (#X,#gamma,#beta):')
    for N in range(6, 13):
        d = dimC_N(n, 1, t, N)
        log('  N=%2d (%s)  %s  total %d' % (N, 'odd' if N % 2 else 'even', d, sum(d.values())))
    log('neighbouring C slots t+e_K (G = 14-|K|), N = 8, 9, 10 (by S_5 symmetry only |K| matters):')
    for r in range(6):
        tt = tuple(-3 + (j < r) for j in range(5))
        for N in (8, 9, 10):
            d = dimC_N(n, 1, tt, N)
            log('  |K|=%d t=%s N=%2d  %s  total %d' % (r, tt, N, d, sum(d.values())))
    log('V = C_quot / sum_k D_k C_quot at t=(-3)^5 (Koszul alternating sum), by content:')
    for N in range(6, 13):
        d = dimV_N(n, 1, t, N)
        log('  N=%2d  %s  total %d' % (N, d, sum(d.values())))
    log('SL_5-invariant parts (multiplicity of det^{-3}, Weyl character formula on weight multiplicities):')
    for N in (8, 9, 10):
        log('  N=%2d  C_quot: %d   V: %d' % (N, weyl_mult(n, 1, t, N, 'C'), weyl_mult(n, 1, t, N, 'V')))
    log('slot t=(-2)^5 (G=9), feeding H(V)_{(-3)^5,N} via K={1..5}, N+5:')
    for N in (13, 14, 15):
        d = dimC_N(n, 1, (-2,) * 5, N)
        log('  N=%2d  %s  total %d   SL_5-inv %d' % (N, d, sum(d.values()), weyl_mult(n, 1, (-2,) * 5, N, 'C')))
    log('time %.1fs' % (time.time() - t0))
