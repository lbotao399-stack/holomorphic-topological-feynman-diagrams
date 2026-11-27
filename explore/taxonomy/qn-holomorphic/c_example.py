"""Representative cocycles of lifted (or chain-level nonzero) 2d sectors: z in Z_S (gauge invariant, c_0-free), Q_1 z,
and the check Q_1 z notin im Q0 (target).  usage: python3 c_example.py THEORY GAUGE MATTER q1.q2.. d k
Q_1 = hatQ_1 (kappa_1 = 1, normalization not fixed by nilpotency at H = 1).  Log: c_example_*.log"""
import sys
from fractions import Fraction as Fr

import numpy as np

import polyops as po
import theories
from wick import Algebra
from explicit import GraphOp
import rgraphs
import lifts as lf


def main(argv):
    th, gname, matter = argv[0], argv[1], (None if argv[2] == '-' else argv[2])
    q = tuple(int(x) for x in argv[3].split('.'))
    d = (int(argv[4]),)
    k = int(argv[5])
    H, sl = theories.build(th, gname, matter)
    alg = Algebra(sl, H)
    g = theories.gauge(gname)
    cartan = [a for a, m in enumerate(g.mats) if not (m - np.diag(np.diag(m))).any()]
    wts = [[int(g.f[h, a, a]) for h in cartan] for a in range(g.dim)]
    pos = [a for a in range(g.dim) if a not in cartan and wts[a][0] > 0]
    gidx = sl.blocks['g']
    lth = lf.LiftTheory(alg, cartan=[gidx[a] for a in cartan], posroots=[gidx[a] for a in pos])
    gr = dict(rgraphs.rooted_trees(1, 1)[0])
    gr['name'] = 'B1'
    Q1 = GraphOp(alg, gr)
    names = {i: '%s[%d]' % (b, j) for b, ids in lth.part.items() for j, i in enumerate(ids)}
    print('theory %s %s %s, sector q=%s d=%s k=%d; charges %s' % (th, gname, matter, q, d, k, lth.charges))

    def inv_sector(qq, dd, kk):
        monos = lf.sector_monomials(lth, qq, dd, kk)
        inv = [lf.vec_to_poly(v, monos) for v in lf.invariants(lth, monos)] if monos else []
        return inv

    inv = inv_sector(q, d, k)
    imgs = [lf.drop_c0(lth, alg.q0(p)) for p in inv]
    idx = {}
    cols = [{idx.setdefault(m, len(idx)): v for m, v in im.items()} for im in imgs]
    Z = []
    for kv in lf.kernel(cols, len(inv)):
        z = {}
        for j, c in kv.items():
            for m, v in inv[j].items():
                po.padd(z, m, c * v)
        Z.append(z)
    dq = (-1, -1) if len(q) == 2 else None
    for zi, z in enumerate(Z):
        print('cocycle z_%d (%d terms):' % (zi, len(z)))
        print(po.show(z, names, H, maxterms=30))
        im = lf.drop_c0(lth, Q1.apply(z))
        print('Q_1 z (%d terms):' % len(im))
        print(po.show(im, names, H, maxterms=30))
    print('letter blocks:', {b: '%s, letter parity %d' % (len(ids), alg.pi[ids[0]]) for b, ids in lth.part.items()})


if __name__ == '__main__':
    main(sys.argv[1:])
