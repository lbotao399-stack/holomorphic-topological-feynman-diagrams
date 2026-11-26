"""Tests of vlib against identities that hold independently of any convention choice:
 d^2 = 0;  {s_0, w} = -d w;  Moyal CME: d j_M = 0, (1/2){j_M,j_M} = d i_M (gl_N);  {s_0,s_0} = 0;
 bracket symmetry {a,b} = {b,a} for arity-3 ghost-0 cochains; derivation d{a,b} = {da,b} - {a,db} (a arity 3);
 formal trace-basis bracket == explicit dual-basis contraction at N = 3 (gl_3, sl_3) on random integer points."""
import random, itertools, time
from fractions import Fraction as Fr
import numpy as np
from vlib import *
import vexp

random.seed(1)
t0 = time.time()


def rnd_cochain(k, deg, sl, nterms=6):
    perms = all_perms(k, sl)
    monos = momentum_monos(k, deg)
    w = {}
    for _ in range(nterms):
        w[(random.choice(perms), random.choice(monos))] = LP({random.choice([0, 1, -1]): Fr(random.randint(-5, 5))})
    return alt(clean(w), k, sl)


for sl in (False, True):
    tag = 'sl' if sl else 'gl'
    s0 = cubic('f', 0, sl)
    for k, deg in [(2, (1, 0)), (3, (1, 1)), (3, (2, 2)), (4, (1, 1))]:
        w = rnd_cochain(k, deg, sl)
        if not w:
            continue
        dd = d(d(w, k, sl), k + 1, sl)
        br = bracket(s0, 3, w, k, sl)
        print("%s k=%d deg=%s: d^2 w = 0: %s;  {s0,w} = -dw: %s" % (tag, k, deg, is_zero(dd), is_zero(lin((1, br), (1, d(w, k, sl))))), flush=True)
    print("%s {s0,s0} = 0: %s" % (tag, is_zero(bracket(s0, 3, s0, 3, sl))))
    a = rnd_cochain(3, (1, 1), sl)
    b = rnd_cochain(3, (1, 0), sl)
    print("%s {a,b} = {b,a} (arity 3,3): %s" % (tag, is_zero(lin((1, bracket(a, 3, b, 3, sl)), (-1, bracket(b, 3, a, 3, sl))))))
    c = rnd_cochain(4, (1, 0), sl)
    lhs = d(bracket(a, 3, c, 4, sl), 5, sl)
    rhs = lin((1, bracket(d(a, 3, sl), 4, c, 4, sl)), (-1, bracket(a, 3, d(c, 4, sl), 5, sl)))
    print("%s d{a,c} = {da,c} - {a,dc} (a arity 3, c arity 4): %s" % (tag, is_zero(lin((1, lhs), (-1, rhs)))), flush=True)

# Moyal (gl_N): s_e = tr(123) e^{(e/2)[p1p2]} - tr(132) e^{-(e/2)[p1p2]};  j_M = (1/2)[p1p2] d,  i_M = (1/8)[p1p2]^2 f
jM = scal(Fr(1, 2), cubic('d', 1, False))
iM = scal(Fr(1, 8), cubic('f', 2, False))
print("Moyal: d j_M = 0: %s" % is_zero(d(jM, 3, False)))
print("Moyal: (1/2){j_M,j_M} = d i_M: %s;  {j_M,j_M} != 0: %s" %
      (is_zero(lin((Fr(1, 2), bracket(jM, 3, jM, 3, False)), (-1, d(iM, 3, False)))), not is_zero(bracket(jM, 3, jM, 3, False))))
# order e^3 of the Moyal CME: 2{s0, k_M} + 2{j_M, i_M} = 0, k_M = (1/48)[p1p2]^3 d
kM = scal(Fr(1, 48), cubic('d', 3, False))
s0 = cubic('f', 0, False)
print("Moyal order e^3: {s0,k_M} + {j_M,i_M} = 0: %s" % is_zero(lin((1, bracket(s0, 3, kM, 3, False)), (1, bracket(jM, 3, iM, 3, False)))))

# formal vs explicit at N = 3
for sl in (False, True):
    Nn = 3
    alg = vexp.Algebra(Nn, sl)
    a = rnd_cochain(3, (1, 1), sl)
    b = rnd_cochain(3, (2, 1), sl)
    br = bracket(a, 3, b, 3, sl)
    ok = True
    for trial in range(4):
        Xs = [alg.random_element(3) for _ in range(4)]
        ps = [(random.randint(-4, 4), random.randint(-4, 4)) for _ in range(3)]
        ps.append((-sum(p[0] for p in ps), -sum(p[1] for p in ps)))
        v1 = vexp.evaluate(br, 4, Xs, ps, Nn)
        v2 = vexp.bracket_explicit(a, 3, b, 3, Xs, ps, alg)
        ok = ok and (v1 == v2)
    print("%s formal bracket == explicit dual-basis bracket at N=3 (4 random points): %s" % ('sl' if sl else 'gl', ok), flush=True)
    w = rnd_cochain(3, (1, 2), sl)
    ok = True
    for trial in range(4):
        Xs = [alg.random_element(3) for _ in range(4)]
        ps = [(random.randint(-4, 4), random.randint(-4, 4)) for _ in range(3)]
        ps.append((-sum(p[0] for p in ps), -sum(p[1] for p in ps)))
        ok = ok and vexp.evaluate(d(w, 3, sl), 4, Xs, ps, Nn) == vexp.d_explicit(w, 3, Xs, ps, Nn)
    print("%s formal d == explicit CE formula at N=3: %s" % ('sl' if sl else 'gl', ok))
print("time %.1fs" % (time.time() - t0))
