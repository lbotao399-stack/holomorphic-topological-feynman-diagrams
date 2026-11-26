"""Assembled first and second moments vs the existing exact route cores.kernel(..., insertions)
(interpolation of U^q D prod(U y) on the whole graph + core periods), C^3 three-loop K4 subdivisions,
sources MOM of cores/check_c3_three_loop.py; bare values vs the k4_periods fixtures."""
import random
import sys
import time
sys.path.insert(0, '../cores')
import assembly
import cores
from check_c3_three_loop import K4SUB, MOM

rng = random.Random(5)
tot = bad = 0
for name, edges, (r, z) in K4SUB:
    t = time.time()
    res = assembly.kernel_moments(8, edges, MOM, order=2, verbose=False)
    ta = time.time() - t
    ref = {k: x for k, x in (('1', r), ('z3', z)) if x}
    okb = res['bare'] == ref
    ins = [(rng.randrange(10), rng.randrange(3)) for _ in range(5)]
    prs = [tuple(sorted([(rng.randrange(10), rng.randrange(3)), (rng.randrange(10), rng.randrange(3))])) for _ in range(3)]
    nb = 0
    t = time.time()
    for a in ins:
        v = cores.kernel(8, edges, MOM, [a])
        nb += v != res.get(('y',) + a, {})
    for a, b in prs:
        v = cores.kernel(8, edges, MOM, [a, b])
        nb += v != res.get(('yy',) + a + b, {})
    tot += 1 + len(ins) + len(prs)
    bad += nb + (not okb)
    print('%-18s assembly (bare + 30 first + 465 second) %.1fs; bare fixture %s; %d first + %d second vs cores.kernel: %d mismatches (%.0fs)'
          % (name, ta, 'equal' if okb else 'DIFFERENT', len(ins), len(prs), nb, time.time() - t), flush=True)
print('total compared %d, mismatches %d' % (tot, bad))
