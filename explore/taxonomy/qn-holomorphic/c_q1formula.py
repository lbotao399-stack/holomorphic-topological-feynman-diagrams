"""2d one-loop bracket.  (i) banana kernel with input jets:  I(lambda; a, b) = -(a! b!/(a+b+1)!) lambda^{1+a+b}
(spoke jets a, b; kernel_symbolic) ;  (ii) colour:  hatQ_1(phi^A phi^B) = eps(A,B) * (-1) * d omega(Phi, [T^A, T^B])
with T^A = sum_A' winv[A,A'] T_A' (omega(T^A, T_B) = delta), i.e. the coefficient of d phi^C equals
eps(A,B) * (-1) * omega(T_C, [T^A, T^B]) / norm ; eps tabulated against the parities; (iii) block table.
Log: c_q1formula.log"""
import itertools
import math
from fractions import Fraction as Fr

import numpy as np

import polyops as po
import theories
import superalg as sa
from wick import Algebra
from explicit import GraphOp
import rgraphs
import kern_h


def banana_kernel():
    g = rgraphs.rooted_trees(1, 1)[0]
    gk = kern_h.GraphKernel(1, g['nb'], g['edges'])
    print('(i) banana n=1 edges %s, spokes %s, output vertex %s' % (g['edges'], gk.spokes, gk.outputs))
    ok = True
    for a in range(6):
        for b in range(6):
            lexp, coef = gk.kernel_symbolic(((0,) * a, (0,) * b))
            terms = {int(lexp[t, 0]): int(coef[t]) % po.P for t in range(lexp.shape[0]) if coef[t] % po.P}
            ref = {1 + a + b: po.modq(Fr(-math.factorial(a) * math.factorial(b), math.factorial(a + b + 1)))}
            if terms != ref:
                ok = False
                print('  mismatch a=%d b=%d: %s vs %s' % (a, b, terms, ref))
    print('  I(lambda; a, b) = -(a! b!/(a+b+1)!) lambda^(1+a+b) for 0 <= a, b <= 5: %s' % ok)


def colour_rule(name, th, gname, matter):
    H, sl = theories.build(th, gname, matter)
    alg = Algebra(sl, H)
    g = dict(rgraphs.rooted_trees(1, 1)[0])
    g['name'] = 'B1'
    Q1 = GraphOp(alg, g)
    N = alg.N
    winv = np.array(sa._inverse_fr(sl.w), dtype=object)
    f = sl.f.astype(object)
    w = sl.w.astype(object)
    # [T^A, T^B] components, omega(T_C, .)
    epsrule = {}
    bad = 0
    from blocksutil import blocks_for
    part = blocks_for(sl)
    blk = {i: b for b, ids in part.items() for i in ids}
    table = {}
    for A in range(N):
        for B in range(N):
            s, m = po.canon([(A, (0,)), (B, (0,))], alg.pi)
            if not s:
                continue
            img = Q1.apply({m: s})          # hatQ1(phi^A phi^B) in the word order (A, B)
            X = np.zeros(N, dtype=object)
            for A2 in range(N):
                if winv[A, A2] == 0:
                    continue
                for B2 in range(N):
                    if winv[B, B2] == 0:
                        continue
                    X = X + winv[A, A2] * winv[B, B2] * f[A2, B2]
            ref = {C: sum(w[C, D] * X[D] for D in range(N)) for C in range(N)}
            ref = {C: v for C, v in ref.items() if v != 0}
            got = {}
            for mm, v in img.items():
                assert len(mm) == 1 and mm[0][1] == (1,)
                got[mm[0][0]] = po.ratrec(v)
            if not ref and not got:
                continue
            if set(ref) != set(got):
                bad += 1
                continue
            ratios = {Fr(got[C]) / Fr(ref[C].numerator, ref[C].denominator) if isinstance(ref[C], Fr) else
                      Fr(got[C]) / ref[C] for C in ref}
            if len(ratios) != 1:
                bad += 1
                continue
            r = ratios.pop()
            key = (int(sl.par[A]), int(sl.par[B]))
            epsrule.setdefault(key, set()).add(r)
            outb = sorted({blk[C] for C in ref})
            table.setdefault((blk[A], blk[B]), set()).update(outb)
    print('(ii) %s: hatQ1(phi^A phi^B) = r(p_A, p_B) * d omega(Phi,[T^A,T^B]); r by parities (p_A, p_B): %s; '
          'mismatches %d' % (name, {k: sorted(v) for k, v in sorted(epsrule.items())}, bad))
    print('     block table (A-block, B-block) -> output blocks: %s' % (
        {'%s.%s' % k: v for k, v in sorted(table.items())}))


if __name__ == '__main__':
    banana_kernel()
    for name, spec in [('(0,2) sl2+fund*4', ('2d02', 'sl2', 'fund*4')), ('(0,2) sl2+adj', ('2d02', 'sl2', 'adj')),
                       ('(2,2) sl2+fund', ('2d22', 'sl2', 'fund')), ('(4,4) sl2+fund', ('2d44', 'sl2', 'fund')),
                       ('(8,8) sl2', ('2d88', 'sl2', None))]:
        colour_rule(name, *spec)
