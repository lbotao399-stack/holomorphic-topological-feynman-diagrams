"""Fork B: 10d hCS on C^5 (g = su2, f = eps, kappa = delta): explicit {Q0, hatQ_1} on physical invariant operators
O = tr(d_a c d_b c) and the C^{5|2} control (su2 (x) Lambda[th1,th2], same H = 5 hexagon) on tr(d_a b d_b b) type.
Log: b_10d_anomaly.log"""
import itertools
import polyops as po
import superalg as sa
import theories
from wick import Algebra
from explicit import graph_ops


def run(label, sl, H, ops_letters):
    alg = Algebra(sl, H)
    G1 = graph_ops(alg, 1)[0]
    for name, O in ops_letters(alg):
        q1 = G1.apply(O)
        R = po.psum(alg.q0(q1), G1.apply(alg.q0(O)))
        lets = {}
        for m in R:
            key = tuple(sorted(sum(al) for _, al in m))
            lets[key] = lets.get(key, 0) + 1
        print('%s: O = %s: |hatQ1 O| = %d monomials, |{Q0, hatQ1} O| = %d monomials; derivative multiset %s'
              % (label, name, len(q1), len(R), sorted(lets.items())[:6]), flush=True)


def tenD(alg):
    n = 3
    for a, b in ((0, 1), (0, 0)):
        O = {}
        for A in range(n):
            ea = tuple(1 if i == a else 0 for i in range(5))
            eb = tuple(1 if i == b else 0 for i in range(5))
            s, m = po.canon([(A, ea), (A, eb)], alg.pi)
            if s:
                po.padd(O, m, s)
        if O:
            yield 'sum_A d%d c^A d%d c^A' % (a + 1, b + 1), O
    # cubic invariant eps_ABC d1 c^A d2 c^B d3 c^C
    O = {}
    for A, B, C in itertools.permutations(range(3)):
        s0 = 1 if (A, B, C) in ((0, 1, 2), (1, 2, 0), (2, 0, 1)) else -1
        s, m = po.canon([(A, (1, 0, 0, 0, 0)), (B, (0, 1, 0, 0, 0)), (C, (0, 0, 1, 0, 0))], alg.pi)
        if s:
            po.padd(O, m, s * s0)
    yield 'eps_ABC d1 c^A d2 c^B d3 c^C', O


def control(alg):
    n = 3
    top = 3 * n                     # mask th1th2 block offset
    O = {}
    for A in range(n):
        s, m = po.canon([(top + A, (0,) * 5), (top + A, (1, 0, 0, 0, 0))], alg.pi)
        if s:
            po.padd(O, m, s)
    yield 'sum_A b^A d1 b^A', O
    O = {}
    for A in range(n):
        s, m = po.canon([(A, (1, 0, 0, 0, 0)), (top + A, (0, 1, 0, 0, 0))], alg.pi)
        if s:
            po.padd(O, m, s)
    yield 'sum_A d1 c^A d2 b^A', O


run('10d su2', sa.lie_even(sa.su2_eps()), 5, tenD)
run('C^{5|2} su2 (x) Lambda2 (control)', sa.berezin_tensor(sa.su2_eps(), 2), 5, control)
