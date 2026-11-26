"""Coverage of the canonicalizer on five-loop skeletons and K5.

For each core G and level q: uniform samples of integer weight vectors a in [1, q-1]^E with
sum a = n q (exact sampling of bounded compositions), kept if convergent (gutil.convergent);
engine.period(G, a, q) with search budget B: reached iff a chain to tabulated periods is found.
Prints fraction, 95% Wilson interval, route statistics, one exact value with its chain.
usage: coverage.py q N kmax [names...]   (kmax: two-point terminal with numerator k <= kmax; 0 = moves only)"""
import json
import math
import os
import random
import sys
import time
from collections import Counter
from functools import lru_cache

import engine as en
import gutil

HERE = os.path.dirname(os.path.abspath(__file__))
SK = {r['name']: [tuple(e) for e in r['edges']] for r in json.load(open(os.path.join(HERE, '..', 'fiveloop', 'skel5.json')))}
SK['K5'] = [(i, j) for i in range(5) for j in range(i + 1, 5)]
ORDER = ['K5-e', 'W5', 'K33+e', 'Z5', 'K4+K4', 'V7_0', 'V7_1', 'V7_2', 'V7_3', 'V8_0', 'V8_1', 'V8_2',
         'V8_3', 'V8_4', 'K5']
ALIAS = {'V8_3': 'cube', 'V8_4': 'Moebius-Wagner V8'}


@lru_cache(maxsize=None)
def count(t, p, hi):
    if p == 0:
        return 1 if t == 0 else 0
    return sum(count(t - x, p - 1, hi) for x in range(1, hi + 1) if t - x >= 0)


def sample(t, p, hi, rng):
    out = []
    for i in range(p):
        r = rng.randrange(count(t, p - i, hi))
        for x in range(1, hi + 1):
            c = count(t - x, p - i - 1, hi) if t - x >= 0 else 0
            if r < c:
                out.append(x)
                t -= x
                break
            r -= c
    return out


def wilson(k, n, z=1.96):
    if n == 0:
        return (0, 0)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0, c - h), min(1, c + h))


def route_kind(chain):
    kinds = []
    s = ' '.join(chain)
    for tag in ('PR', 'FT', 'TW', 'M2', 'M1'):
        if tag in s:
            kinds.append(tag)
    term = [x for x in chain[-1:]]
    return '+'.join(kinds) or 'direct'


if __name__ == '__main__':
    q = int(sys.argv[1])
    N = int(sys.argv[2])
    en.NUMER[0] = int(sys.argv[3])
    names = sys.argv[4:] or ORDER
    budget = 60
    out = {}
    for name in names:
        G = SK[name]
        n = gutil.cycle_rank(G)
        if name == 'K5' and q < 8:
            pass
        rng = random.Random(12345 + q)
        t = time.time()
        tried = acc = hit = 0
        kinds = Counter()
        terms = Counter()
        example = None
        while acc < N:
            a = sample(n * q, len(G), q - 1, rng)
            tried += 1
            if not gutil.convergent(tuple(G), a, q):
                continue
            acc += 1
            r = en.period(G, a, q, budget=budget)
            if r is None:
                continue
            hit += 1
            kinds[route_kind(r[1])] += 1
            tm = r[1][-1].split(':')[-1].rstrip(']') if r[1][-1].startswith('T:') else 'PR'
            terms[tm] += 1
            if (example is None or len(r[1]) > len(example[2])) and len(r[1]) >= 2 and len(' '.join(r[1])) < 300:
                example = (a, r[0], r[1])
        lo, hi = wilson(hit, acc)
        dt = time.time() - t
        print('%-6s %-18s q=%d  convergent/sampled %d/%d  reached %d/%d = %.3f  [%.3f, %.3f]  %.1fs'
              % (name, ALIAS.get(name, ''), q, acc, tried, hit, acc, hit / acc, lo, hi, dt), flush=True)
        print('       routes %s' % dict(kinds.most_common()), flush=True)
        if example:
            print('       example I_%d%s = %s\n         chain %s' % (q, tuple(example[0]), en.vstr(example[1]),
                                                                    ' '.join(example[2])), flush=True)
        out[name] = dict(q=q, sampled=tried, convergent=acc, reached=hit, frac=hit / acc, ci=[lo, hi],
                         routes=dict(kinds), example=[list(example[0]), en.vstr(example[1]), list(example[2])]
                         if example else None, seconds=round(dt, 1))
    json.dump(out, open(os.path.join(HERE, 'coverage_q%d_k%d_%s.json' % (q, en.NUMER[0], '_'.join(names) if len(names) < 4 else 'all')),
                        'w'), indent=1)
