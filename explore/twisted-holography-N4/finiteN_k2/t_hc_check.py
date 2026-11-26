"""Check: Q-cohomology of graded-cyclic words of the cobar algebra at weight (a;0,0) = binom(s(a)-1, n) at gh -n."""
import sys, itertools
from math import comb
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
from holography import bg_words as bw, bg_fast as bf, rank_mod
lam = (3, 3, 2, 0, 0)
alg, _ = bw.n4_letters(2 * sum(lam))
tb = bf.sector_tables(alg, lam)
cyc, CQ = tb["cyc"], tb["CQ"]
bywt = {}
for i, cw in enumerate(cyc):
    key = (alg.wword(cw), alg.ghword(cw))
    bywt.setdefault(key, []).append(i)
bad = 0; tested = 0
for a in itertools.product(*[range(x + 1) for x in lam[:3]]):
    if not sum(a): continue
    om = tuple(a) + (0, 0)
    s = sum(1 for x in a if x)
    ghs = sorted(set(g for (w, g) in bywt if w == om))
    rk = {}
    for g in ghs:
        sp = rank_mod.SpanP()
        for i in bywt[(om, g)]:
            v = {j: c for j, c in CQ[i]}
            if v: sp.add(v)
        rk[g] = sp.dim()
    for g in ghs:
        h = len(bywt[(om, g)]) - rk.get(g, 0) - rk.get(g - 1, 0)
        exp = comb(s - 1, -g) if -g >= 0 and -g <= s - 1 else 0
        tested += 1
        if h != exp:
            bad += 1; print('MISMATCH', a, g, h, exp)
print('weights tested', tested, 'mismatches', bad)
