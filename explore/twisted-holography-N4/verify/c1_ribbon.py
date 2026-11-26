"""Symbolic-in-N fundamental gluing by union-find on matrix index variables (independent of code/).

Trace k-letter word tr(L_0...L_{k-1}) = sum prod (L_t)_{v(t), v(t+1)}.  Completeness
sum_c (e^c)_{pq} (e_c)_{rs} = d_{ps} d_{qr}  (tau = 1).  A gluing returns {(N-power, letter cycles): coeff}.
"""
import sys
from collections import Counter
from itertools import permutations, product
from math import comb


def canon(cycle):
    k = len(cycle)
    return min(tuple(cycle[i:] + cycle[:i]) for i in range(k))


def glue(trace1, trace2, pairs):
    """trace1, trace2: lists of letters; letter = ('g', label) or ('m', name).
    pairs: dict label -> None; each generator label occurs once in each trace (e^c in trace1, e_c in trace2)."""
    parent = {}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    letters = []            # (trace, pos, letter, rowvar, colvar)
    for t, tr in enumerate((trace1, trace2)):
        k = len(tr)
        for pos, L in enumerate(tr):
            row, col = (t, pos), (t, (pos + 1) % k)
            parent.setdefault(row, row)
            parent.setdefault(col, col)
            letters.append((t, pos, L, row, col))
    gen = {}
    for t, pos, L, row, col in letters:
        if L[0] == 'g':
            gen.setdefault(L[1], {})[t] = (row, col)
    for lab, d in gen.items():
        (p, q), (r, s) = d[0], d[1]
        union(p, s)
        union(q, r)
    mats = [(L[1], find(row), find(col)) for t, pos, L, row, col in letters if L[0] == 'm']
    touched = set()
    rowclass = {}
    for name, rc, cc in mats:
        touched.add(rc); touched.add(cc)
        assert rc not in rowclass
        rowclass[rc] = name
    classes = {find(x) for x in parent}
    free = len(classes - touched)
    succ = {name: rowclass[cc] for name, rc, cc in mats}
    seen, cycles = set(), []
    for name, _, _ in mats:
        if name in seen:
            continue
        cyc, x = [], name
        while x not in seen:
            seen.add(x); cyc.append(x); x = succ[x]
        cycles.append(canon(cyc))
    return free, tuple(sorted(cycles))


def nested_form(m, spectator=True):
    """sum_sigma sum_signs: tr(Y_1..Y_m F) with tr(ad_{X_sigma(m)}...ad_{X_sigma(1)}(D) E), X_j glued to Y_sigma(j)."""
    res = Counter()
    t1 = [('g', i) for i in range(m)] + ([('m', 'F')] if spectator else [])
    for sig in permutations(range(m)):
        for signs in product((1, -1), repeat=m):
            word = [('m', 'D')]
            for j in range(m):
                L = ('g', sig[j])
                word = [L] + word if signs[j] == 1 else word + [L]
            s = 1
            for x in signs:
                s *= x
            res[glue(t1, word + [('m', 'E')], None)] += s
    return {k: v for k, v in res.items() if v}


def k_form(m, spectator=True):
    """sum_k (-1)^k C(m,k) sum_tau gluing of tr(Y_1..Y_m F) with tr(X^{m-k} D X^k E)."""
    res = Counter()
    t1 = [('g', i) for i in range(m)] + ([('m', 'F')] if spectator else [])
    for k in range(m + 1):
        w = (-1)**k * comb(m, k)
        for tau in permutations(range(m)):
            g = [('g', tau[i]) for i in range(m)]
            word = g[:m - k] + [('m', 'D')] + g[m - k:] + [('m', 'E')]
            res[glue(t1, word, None)] += w
    return {k: v for k, v in res.items() if v}


def poly(res, structure):
    """coefficient polynomial {N-power: coeff} of a given letter-cycle structure."""
    out = {}
    for (free, cyc), v in res.items():
        if cyc == structure:
            out[free] = out.get(free, 0) + v
    return dict(sorted(out.items(), reverse=True))


def one_handle(n):
    return comb(n + 2, 4) + n * (n - 1) + comb(n, 3) + 2 + (-1)**n * (n - 1)


if __name__ == '__main__':
    mmax_nested = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    mmax_k = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    DEF, DFE = ('D', 'E', 'F'), ('D', 'F', 'E')
    # cubic trace without spectator (n=2): colour tensor must vanish identically
    r = nested_form(3, spectator=False)
    print('n=2, tr(Y1Y2Y3) without spectator, nested form: nonzero terms =', r)
    for m in range(3, mmax_nested + 1):
        a, b = nested_form(m), k_form(m)
        assert a == b, (m, a, b)
        others = {k: v for k, v in a.items() if k[1] not in ((DEF,), (DFE,))}
        print(f'm={m} (n={m-1}) nested==k_form; DEF:{poly(a,(DEF,))} DFE:{poly(a,(DFE,))} other structures:{others}', flush=True)
    for m in range(mmax_nested + 1, mmax_k + 1):
        b = k_form(m)
        others = {k: v for k, v in b.items() if k[1] not in ((DEF,), (DFE,))}
        print(f'm={m} (n={m-1}) k_form; DEF:{poly(b,(DEF,))} DFE:{poly(b,(DFE,))} other structures:{others}', flush=True)
        n = m - 1
        pd = poly(b, (DEF,))
        print('   one-handle coefficient', pd.get(n - 2), ' formula', one_handle(n), flush=True)
