"""Closed forms of hat Q_1 on two letters without derivatives (checked against ops_c41.Theory.Q1):
  hat Q_1(gamma^A gamma^B) = (1/12) eps^{abcd} F_{A D1 D2 D3 D4 B} sum_q (-1)^{q-1} d_a C_1 d_b C_2 d_c C_3 d_d C_4 |_{C_q = gamma, others c}
  hat Q_1(gamma^A c^B)     = -(1/12) eps^{abcd} F_{A D1 D2 D3 D4 B} d_a c^{D1} d_b c^{D2} d_c c^{D3} d_d c^{D4},
  F_{A D1 D2 D3 D4 B} = f_{A D1 R1} f_{R1 D2 R2} f_{R2 D3 R3} f_{R3 D4 B}."""
import sys
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/msym-D8/code')
from fractions import Fraction as Fr
import ops_c41 as O

for alg in ('su2', 'so5'):
    n, f = O.structure_constants(alg)
    T = O.Theory(n, f)
    e = lambda a: tuple(1 if i == a else 0 for i in range(4))
    for (A, B) in ((0, 1), (0, 0), (1, 2)):
        gg, gc = {}, {}
        for (D1, R1, v1) in T.fby.get(A, []):
            for (D2, R2, v2) in T.fby.get(R1, []):
                for (D3, R3, v3) in T.fby.get(R2, []):
                    for (D4, BB, v4) in T.fby.get(R3, []):
                        if BB != B:
                            continue
                        F = v1 * v2 * v3 * v4
                        for (a, b, c, d), s in O.EPS4.items():
                            Ds = [(D1, a), (D2, b), (D3, c), (D4, d)]
                            cterm = O.const(1)
                            for (D, x) in Ds:
                                cterm = O.mul(cterm, O.L('c', D, e(x)))
                            gc = O.add(gc, O.scale(cterm, Fr(-1, 12) * F * s))
                            for q in range(4):
                                term = O.const(1)
                                for i, (D, x) in enumerate(Ds):
                                    term = O.mul(term, O.L('g' if i == q else 'c', D, e(x)))
                                gg = O.add(gg, O.scale(term, Fr(1, 12) * F * s * (-1) ** q))
        r1 = T.Q1(O.mul(O.L('g', A), O.L('g', B)))
        r2 = T.Q1(O.mul(O.L('g', A), O.L('c', B)))
        print('%s A=%d B=%d: hatQ1(gamma gamma) closed form: %s (%d terms); hatQ1(gamma c) closed form: %s (%d terms)'
              % (alg, A, B, r1 == gg, len(r1), r2 == gc, len(r2)))
