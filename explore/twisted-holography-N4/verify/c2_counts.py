"""C2: unmarked and marked face counts, reduced core algorithm, large-m quasi-polynomials.

Conventions (derived from the union-find gluing of c1_ribbon.py):
second-vertex corners i = 0..m-1 (corner i after slot i), alpha(i) = i+1, beta = m-cycle,
psi(i) = beta(alpha(i)); E in corner m-1, D in corner m-1-k (k=0: D,E in corner m-1, order DE);
coefficient of tr(DEF):  P_m(N) = sum_beta N^{c(psi)-1} [ l(m-1) + sum_{0<k<m} (-1)^k C(m,k) dist_psi(m-1 -> m-1-k) ].
"""
import sys
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import permutations, combinations
from math import comb, factorial, prod


def m_cycles(m):
    if m == 1:
        yield (0,)
        return
    for p in permutations(range(1, m)):
        beta = [0] * m
        seq = (0,) + p
        for a, b in zip(seq, seq[1:] + (0,)):
            beta[a] = b
        yield tuple(beta)


def cycles_of(perm):
    m = len(perm)
    seen = [False] * m
    out = []
    for i in range(m):
        if not seen[i]:
            cyc = []
            j = i
            while not seen[j]:
                seen[j] = True
                cyc.append(j)
                j = perm[j]
            out.append(cyc)
    return out


def unmarked(m):
    """{c: number of m-cycles beta with c(beta alpha) = c}."""
    res = Counter()
    for beta in m_cycles(m):
        psi = tuple(beta[(i + 1) % m] for i in range(m))
        res[len(cycles_of(psi))] += 1
    return res


def unmarked_formula(m, N):
    return Fraction(prod(N + j for j in range(m + 1)) - prod(N - j for j in range(m + 1)), m * (m + 1))


def elem_sym(k, xs):
    e = [1] + [0] * k
    for x in xs:
        for j in range(k, 0, -1):
            e[j] += e[j - 1] * x
    return e[k]


def marked_face(m):
    """{h: P_{m,h}} via the face formula, enumerating all m-cycles beta."""
    res = Counter()
    binom = [(-1)**k * comb(m, k) for k in range(m + 1)]
    for beta in m_cycles(m):
        psi = [beta[(i + 1) % m] for i in range(m)]
        c = len(cycles_of(psi))
        w = 0
        x, dist = psi[m - 1], 1
        while x != m - 1:
            k = m - 1 - x
            w += binom[k] * dist
            x, dist = psi[x], dist + 1
        w += dist                                  # cycle length l(m-1), k = 0 term
        h2 = m - c
        assert h2 % 2 == 0
        res[h2 // 2] += w
    return dict(sorted(res.items()))


# ---------------- reduced core algorithm ----------------
def fpf_perms_with_cycles(s, c):
    """fixed-point-free permutations of range(s) with exactly c cycles."""
    def rec(remaining, cycles_left, perm):
        if not remaining:
            if cycles_left == 0:
                yield tuple(perm)
            return
        if cycles_left == 0:
            return
        first = remaining[0]
        rest = remaining[1:]
        # cycle containing `first`, of length L >= 2; leave at least 2*(cycles_left-1) points
        for L in range(2, len(remaining) - 2 * (cycles_left - 1) + 1):
            for others in combinations(rest, L - 1):
                left = [x for x in rest if x not in others]
                for order in permutations(others):
                    cyc = (first,) + order
                    for a, b in zip(cyc, cyc[1:] + (first,)):
                        perm[a] = b
                    yield from rec(left, cycles_left - 1, perm)
                    for a in cyc:
                        perm[a] = -1
    yield from rec(list(range(s)), c, [-1] * s)


def is_single_cycle(perm):
    s = len(perm)
    x, n = perm[0], 1
    while x != 0:
        x, n = perm[x], n + 1
    return n == s


def core_statistics(h, smax=None):
    """For each support size s: (#cores, sum of l(s-1), distance vector D[j] = sum dist(s-1 -> j))."""
    stats = {}
    if h == 0:
        return {0: (1, 0, [])}
    for s in range(2 * h + 1, (4 * h if smax is None else min(4 * h, smax)) + 1):
        count, L, D = 0, 0, [0] * max(s - 1, 0)
        for psi in fpf_perms_with_cycles(s, s - 2 * h):
            # beta = psi o alpha^{-1}: beta(y) = psi(y-1)
            beta = tuple(psi[(y - 1) % s] for y in range(s))
            if not is_single_cycle(beta):
                continue
            count += 1
            x, dist = psi[s - 1], 1
            while x != s - 1:
                D[x] += dist
                x, dist = psi[x], dist + 1
            L += dist
        stats[s] = (count, L, D)
    return stats


def P_reduced(m, stats, closed=False):
    """P_{m,h} from core statistics; direct alternating sum unless closed=True."""
    total = 0
    for s, (count, L, D) in stats.items():
        if s == 0:
            total += count * comb(m - 1, 0)
            continue
        if s > m:
            continue
        total += count * comb(m - 1, s) + L * comb(m - 1, s - 1)
        for j, dj in enumerate(D):
            if not dj:
                continue
            if closed:
                # finite-difference evaluation (valid for s-2 < m)
                total -= dj * ((-1)**(s - 2 - j) * comb(m - 1, j) + (-1)**(m + j) * comb(m - 1, s - 2 - j))
            else:
                acc = 0
                for k in range(1, m):
                    a, b = m - 1 - k, k - 1
                    if a >= j and b >= s - 2 - j:
                        acc += (-1)**k * comb(m, k) * comb(a, j) * comb(b, s - 2 - j)
                total += dj * acc
    return total


def B_formula(m, h):
    return Fraction(elem_sym(2 * h + 1, range(1, m + 1)), comb(m + 1, 2))


if __name__ == '__main__':
    what = sys.argv[1] if len(sys.argv) > 1 else 'all'
    if what in ('unmarked', 'all'):
        mmax = int(sys.argv[2]) if len(sys.argv) > 2 else 9
        for m in range(1, mmax + 1):
            u = unmarked(m)
            ok = all(sum(v * N**c for c, v in u.items()) == unmarked_formula(m, N) for N in range(-3, 12))
            okB = all(u.get(m - 2 * h, 0) == B_formula(m, h) for h in range(0, m // 2 + 1))
            print(f'unmarked m={m}: {dict(sorted(u.items(), reverse=True))}  formula match={ok}  B=e_(2h+1)/C(m+1,2) match={okB}', flush=True)
    if what in ('marked', 'all'):
        mmax = int(sys.argv[3]) if len(sys.argv) > 3 else 9
        for m in range(2, mmax + 1):
            print(f'marked face formula m={m} (n={m-1}): {marked_face(m)}', flush=True)
