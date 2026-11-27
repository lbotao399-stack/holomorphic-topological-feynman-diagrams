"""2d (H=1) nilpotency: Q = Q0 + hbar Q1 (hatQ_2 = 0: the n=2 banana has no output source [drv]), so
(Q^2)_hbar = {Q0, hatQ1} = 0 and (Q^2)_hbar^2 = hatQ1^2 = 0, checked exactly (mod P) on random jet monomials
with hatQ1 O != 0 (resp. hatQ1 hatQ1 O computed with hatQ1 O != 0 having >= 2 letters).
usage: python3 c_nilp2d.py THEORY GAUGE MATTER NOPS SEED JETMAX
"""
import random
import sys
import time

import polyops as po
import theories
from wick import Algebra
from explicit import GraphOp
import rgraphs


def main(argv):
    th, gname, matter = argv[0], argv[1], (None if argv[2] == '-' else argv[2])
    nops, seed, jetmax = int(argv[3]), int(argv[4]), int(argv[5])
    H, sl = theories.build(th, gname, matter)
    alg = Algebra(sl, H)
    assert H == 1
    g = dict(rgraphs.rooted_trees(1, 1)[0])
    g['name'] = 'B1'
    Q1 = GraphOp(alg, g)
    rng = random.Random(seed)
    t0 = time.time()
    print('theory %s gauge %s matter %s H=%d N=%d nops=%d seed=%d jetmax=%d P=%d'
          % (th, gname, matter, H, alg.N, nops, seed, jetmax, po.P), flush=True)
    # n = 2 banana: no output vertex
    g2 = rgraphs.rooted_trees(1, 2)[0]
    print('n=2 banana %s: outputs %s -> hatQ_2 = 0' % (g2['edges'], rgraphs.outputs(g2['nb'], g2['edges'])), flush=True)
    stats = {'anti': [0, 0, 0], 'sq': [0, 0, 0]}
    tries = 0
    while (stats['anti'][0] < nops or stats['sq'][0] < nops) and tries < 20000:
        tries += 1
        k = rng.randint(2, 4)
        lets = [(rng.randrange(alg.N), (rng.randint(0, jetmax),)) for _ in range(k)]
        s, m = po.canon(lets, alg.pi)
        if not s:
            continue
        O = {m: 1}
        img = Q1.apply(O)
        if not img:
            continue
        if stats['anti'][0] < nops:
            a = po.psum(alg.q0(img), Q1.apply(alg.q0(O)))
            stats['anti'][0] += 1
            stats['anti'][1] += (len(a) == 0)
            stats['anti'][2] = max(stats['anti'][2], len(a))
        if stats['sq'][0] < nops and k >= 3:
            sq = Q1.apply(img)
            stats['sq'][0] += 1
            stats['sq'][1] += (len(sq) == 0)
            stats['sq'][2] = max(stats['sq'][2], len(sq))
            if 'nontriv' not in stats:
                stats['nontriv'] = 0
            # nontrivial: some monomial of hatQ1 O with >= 2 letters on which hatQ1 acts nonzero
            if any(Q1.apply({mm: 1}) for mm in img if len(mm) >= 2):
                stats['nontriv'] += 1
    print('{Q0, hatQ1} O = 0 on %d/%d operators (max terms of residual %d)' % (
        stats['anti'][1], stats['anti'][0], stats['anti'][2]), flush=True)
    print('hatQ1^2 O = 0 on %d/%d operators with >= 3 letters (max terms %d); hatQ1 acts nonzero on a term of '
          'hatQ1 O in %d of them' % (stats['sq'][1], stats['sq'][0], stats['sq'][2], stats.get('nontriv', 0)),
          flush=True)
    print('time %.1fs' % (time.time() - t0), flush=True)


if __name__ == '__main__':
    main(sys.argv[1:])
