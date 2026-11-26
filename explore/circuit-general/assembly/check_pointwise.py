"""Pointwise checks of the assembled integrand against htk.brute (exact rational layer):
 (a) chain decomposition of all currents y^i_e (Network.check_point, exact rationals, 3 points);
 (b) interpolated numerator P(t) = U^(n+H) D(t) mod p  vs  brute.density(G, t, lam) * U^(n+H) at
     random rational t (integer-scaled; P homogeneous of degree n^2-1);
 (c) assembled insertion integrand at the same t: P(t) * y^i_e(t) [chain formula] vs brute, mod p.
usage: python3 check_pointwise.py SET IDX [IDX ...]"""
import random
import sys
from fractions import Fraction as F

import numpy as np

import sets
import assembly
from numerator import BareNumerator
from htk import brute
from run_batch import HSET

P = 2147483629


def main(argv):
    name = argv[0]
    H = HSET[name]
    its = dict(sets.items(name))
    for idx in [int(x) for x in argv[1:]]:
        G, O, S0 = its[idx]
        Gh, lab, lam, full = sets.instance(name, idx, G, O, S0, H)
        net = assembly.Network(Gh.nv, Gh.edges, full)
        for s in (1, 2, 3):
            net.check_point(seed=s)
        num = BareNumerator(net.G, net.lam)
        coef = num.coefficients_mod(P, check=False)
        ex = num.exponents().astype(np.int64)
        rng = random.Random(idx)
        nbad = 0
        for trial in range(3):
            t = [rng.randint(1, 50) for _ in range(net.N)]
            Db = brute.density(net.G, t, net.lam[:-1])
            U = brute.symanzik_U(net.G, t)
            ref = Db.rational() * U ** (net.n + net.H)
            val = np.ones(len(coef), dtype=np.int64)
            for e in range(net.N):
                pw = [1]
                for _ in range(int(ex[:, e].max())):
                    pw.append(pw[-1] * t[e] % P)
                val = val * np.array(pw, dtype=np.int64)[ex[:, e]] % P
            Pt = int((val * coef % P).sum() % P)
            r = F(ref)
            if Pt != r.numerator % P * pow(r.denominator % P, P - 2, P) % P:
                nbad += 1
            y = brute.currents(net.G, t, net.lam[:-1])
            for e in range(net.N):
                for i in range(net.H):
                    s, A1, A0 = net.insertion_coefficients(e, i)
                    R = [sum(t[f] for f in c['edges']) for c in net.chains]
                    xv = [F(1)] + [F(t[net.chains[ci]['edges'][j]], R[ci]) for ci, j in net.xw]
                    cotm = [np.prod([R[f] for f in range(net.Ns) if net.cot[ti, f]]) for ti in range(len(net.trees))]
                    Us = sum(cotm)
                    ya = s * (sum(A1[ti, b] * xv[b] * cotm[ti] for ti in range(len(net.trees)) for b in range(net.nw)) / Us
                              + sum(A0[b] * xv[b] for b in range(net.nw)))
                    lhs = F(Pt) * 0 + ya * r          # exact: assembled current times exact P
                    rhs = y[e][i] * r
                    if lhs != rhs:
                        nbad += 1
        print('%s #%d: N=%d, currents exact at 3 points, P(t) and P(t) y^i_e(t) at 3 rational points x %d (e,i): %d failures'
              % (name, idx, net.N, net.N * net.H, nbad), flush=True)


if __name__ == '__main__':
    main(sys.argv[1:])
