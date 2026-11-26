"""Compare tower values (X + DS from the level below) with independent direct evaluations
(symmetry engine E, per5 P = K33+e tables / dual prism, graphical functions S) at level q.
usage: python3 check_direct.py Q [methods] [max_per_kind] [seed]"""
import random
import sys
import time

import bval
import k6
from tower import load_values, vstr

q = int(sys.argv[1])
methods = tuple(sys.argv[2]) if len(sys.argv) > 2 else ('E', 'P', 'S')
mx = int(sys.argv[3]) if len(sys.argv) > 3 else 10 ** 9
seed = int(sys.argv[4]) if len(sys.argv) > 4 else 1
vals = load_values(q, 'tower')
keys = sorted(vals)
random.Random(seed).shuffle(keys)
stats = {}
cnt = {}
t0 = time.time()
for b in keys:
    kd = k6.kind(b, q)
    dec = bool(k6.k5_decompletions(b, q))
    tag = kd + ('d' if dec else 'a')
    if cnt.get(tag, 0) >= mx:
        continue
    cnt[tag] = cnt.get(tag, 0) + 1
    for m in methods:
        t = time.time()
        try:
            r = bval.value(b, q, methods=(m,))
        except Exception as ex:
            r = None
        dt = time.time() - t
        key = (tag, m)
        s = stats.setdefault(key, [0, 0, 0])
        if r is None:
            s[2] += 1
            continue
        if r[0] == vals[b]:
            s[0] += 1
        else:
            s[1] += 1
            print('MISMATCH', q, b, m, r, 'tower:', vstr(vals[b]), flush=True)
print('q=%d checked classes per kind: %s' % (q, cnt))
for (tag, m), (ok, bad, na) in sorted(stats.items()):
    print('  kind %-3s method %s: agree %d, differ %d, not evaluable %d' % (tag, m, ok, bad, na))
print('time %.0fs' % (time.time() - t0))
