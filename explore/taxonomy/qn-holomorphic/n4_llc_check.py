"""4d N=4 (g (x) Lambda[th1,th2,th3], su(2), f = eps = C of LLC) against LLC eq-sum-N4 / eq-2-26 / eq-2-27 / eq-2-28.
Letters phi_m (Phi = sum_m m phi_m = C):  phi_0 = i g c, phi_{thI} = gamma^I, phi_{th2th3} = -beta_1, phi_{th1th3} = +beta_2,
phi_{th1th2} = -beta_3, phi_{th1th2th3} = b/(ig).
LLC in phi letters [drv from src LLC eq-2-26..2-28, Q_0 C = (i/2) C C C]:
  Q0 phi_0^D = a0 C_DBC phi_0^B phi_0^C,  a0 = i/2
  Q1(phi_th1^A phi_th2th3^B)          contains a1 (CC_DE + CC_ED) d1 phi_0^D d2 phi_0^E, CC_DE = C_ARD C_BRE, a1 = 1/(64 pi^2)
  Q1(phi_th2th3^A phi_th1th3^B)       contains b1 (C_ARD C_BRE + C_BRD C_ARE) d1 phi_0^D d2 phi_th3^E,  b1 = -1/(64 pi^2)
  Q2(phi_th2th3^A1 phi_th1th3^A2 phi_th1th2^A3) contains a2 (T_DE - T_ED) d1d1 phi_0^D d2d2 phi_0^E,
       T_DE = sum_sigma C_{A_s1 D R} C_{A_s2 R S} C_{A_s3 S E},  a2 = i/(48 (32 pi^2)^2)
Invariants under Q -> l Q, Q_n -> m^n Q_n, Phi -> v Phi:  a2 a0 / a1^2 = -1/24,  b1 / a1 = -1.
Here Q_n = kappa_n hatQ_n with kappa_1 = kappa_2 = 1/2 (nilp_explicit), hatQ_0 = Q0.  Log: n4_llc_check.log"""
import itertools
from fractions import Fraction as Fr

import numpy as np

import polyops as po
import superalg as sa
from wick import Algebra
from explicit import graph_ops

g = sa.su2_eps()
sl = sa.berezin_tensor(g, 3)
alg = Algebra(sl, 2)
eps = g.f
n = g.dim
masks = {'0': 0, 't1': 1, 't2': 2, 't3': 4, 't2t3': 6, 't1t3': 5, 't1t2': 3, 't1t2t3': 7}


def L(m, a, jet=(0, 0)):
    return (masks[m] * n + a, jet)


def coeff(poly, word):
    s, mono = po.canon(word, alg.pi)
    return s * poly.get(mono, 0) % po.P


def rat(x):
    return po.ratrec(x % po.P)


k1, k2 = po.modq(Fr(1, 2)), po.modq(Fr(1, 2))
G1 = graph_ops(alg, 1)[0]
G2 = graph_ops(alg, 2)[0]
# a0
q0 = alg.q0({(L('0', 0),): 1})
c12 = coeff(q0, [L('0', 1), L('0', 2)])           # coefficient of phi_0^1 phi_0^2 in Q0 phi_0^0
# Q0 phi^0 = a0 (C_012 phi^1 phi^2 + C_021 phi^2 phi^1) = 2 a0 eps_012 phi^1 phi^2
a0 = Fr(rat(c12)) / (2 * int(eps[0, 1, 2]))
print('Q0 phi_0^0 = %s * phi_0^1 phi_0^2  ->  a0 = %s' % (rat(c12), a0))
# a1
ratios = set()
for A in range(n):
    for B in range(n):
        img = po.pscale(G1.apply({po.canon([L('t1', A), L('t2t3', B)], alg.pi)[1]:
                                  po.canon([L('t1', A), L('t2t3', B)], alg.pi)[0]}), k1)
        for D in range(n):
            for E in range(n):
                cc = sum(int(eps[A, R, D]) * int(eps[B, R, E]) + int(eps[A, R, E]) * int(eps[B, R, D]) for R in range(n))
                v = coeff(img, [L('0', D, (1, 0)), L('0', E, (0, 1))])
                if cc:
                    ratios.add(Fr(rat(v)) / cc)
                else:
                    assert v == 0
print('Q1(phi_t1^A phi_t2t3^B): coefficient / (CC_DE + CC_ED) of d1 phi_0^D d2 phi_0^E:', ratios)
a1 = ratios.pop() if len(ratios) == 1 else None
ratios = set()
for A in range(n):
    for B in range(n):
        s, m = po.canon([L('t2t3', A), L('t1t3', B)], alg.pi)
        img = po.pscale(G1.apply({m: s}), k1)
        for D in range(n):
            for E in range(n):
                cc = sum(int(eps[A, R, D]) * int(eps[B, R, E]) + int(eps[B, R, D]) * int(eps[A, R, E]) for R in range(n))
                v = coeff(img, [L('0', D, (1, 0)), L('t3', E, (0, 1))])
                if cc:
                    ratios.add(Fr(rat(v)) / cc)
                else:
                    assert v == 0
print('Q1(phi_t2t3^A phi_t1t3^B): coefficient / (CC + CC) of d1 phi_0^D d2 phi_t3^E:', ratios)
b1 = ratios.pop() if len(ratios) == 1 else None
ratios = set()
for A1, A2, A3 in itertools.product(range(n), repeat=3):
    s, m = po.canon([L('t2t3', A1), L('t1t3', A2), L('t1t2', A3)], alg.pi)
    if not s:
        continue
    img = po.pscale(G2.apply({m: s}), k2)
    As = (A1, A2, A3)
    for D in range(n):
        for E in range(n):
            T = 0
            for sg in itertools.permutations(range(3)):
                for R in range(n):
                    for S in range(n):
                        T += int(eps[As[sg[0]], D, R]) * int(eps[As[sg[1]], R, S]) * int(eps[As[sg[2]], S, E])
                        T -= int(eps[As[sg[0]], E, R]) * int(eps[As[sg[1]], R, S]) * int(eps[As[sg[2]], S, D])
            v = coeff(img, [L('0', D, (2, 0)), L('0', E, (0, 2))])
            if T:
                ratios.add(Fr(rat(v)) / T)
            else:
                assert v == 0, (A1, A2, A3, D, E, rat(v))
print('Q2(phi_t2t3^A1 phi_t1t3^A2 phi_t1t2^A3): coefficient / (T_DE - T_ED) of d1d1 phi_0^D d2d2 phi_0^E:', ratios)
a2 = ratios.pop() if len(ratios) == 1 else None
print('a0 = %s, a1 = %s, b1 = %s, a2 = %s' % (a0, a1, b1, a2))
print('a2 a0 / a1^2 = %s   (LLC: -1/24)' % (a2 * a0 / a1 ** 2))
print('b1 / a1 = %s   (LLC: -1)' % (b1 / a1))
