"""Determine the tadpole(5,0) kernel t4 = K/pi^2 (source 1 at q, z = 0) from nilpotency at hbar^4 in a theory
with nonzero symmetrized 5-trace: GL(2) + 1 fundamental (W = kappa = 0).  Monomials: B^4 X and B^5 type
(targeted) plus a random sample.  Solve Qh2^2 + u {Q0,Qh4fan} + v {Q0,T(5,0)} = 0.  Log: test_q4_t4.log"""
import sys, itertools, random, time
from fractions import Fraction as Fr
import sympy as sp
from ht3d import Engine, mono, pstr, psum, pscale, load_kernels
from fanengine import ShapeEngine, fan_kval, tad_kval
from theories import *

nsamp = int(sys.argv[1]) if len(sys.argv) > 1 else 40
th = gl_theory(2, 1, 0)
E = Engine(th, load_kernels(4)); A = E.A
F4 = ShapeEngine(th, 4, 'fan', fan_kval(4, Fr(1, 180)))
T4 = ShapeEngine(th, 4, 'tad', tad_kval(4, Fr(1)))
Bs = [('B', a, 0) for a in range(4)]
others = [('c', a, 0) for a in range(4)] + Bs + [('p0', i, 0) for i in range(2)] + [('s0', i, 0) for i in range(2)]
rng = random.Random(5)
cands = []
for b4 in itertools.combinations_with_replacement(Bs, 4):
    for X in others:
        O = mono(A, *(list(b4) + [X]))
        if O:
            cands.append(O)
sample = rng.sample(cands, min(nsamp, len(cands)))
t0 = time.time()
rows = []
nX = nYf = nYt = 0
for O in sample:
    X = E.Q2(E.Q2(O))
    Yf = psum(E.Q0(F4(O)), F4(E.Q0(O)))
    Yt = psum(E.Q0(T4(O)), T4(E.Q0(O)))
    nX += bool(X); nYf += bool(Yf); nYt += bool(Yt)
    for m in set(X) | set(Yf) | set(Yt):
        rows.append((X.get(m, 0), Yf.get(m, 0), Yt.get(m, 0)))
    print("  %s: |Qh2^2| %d, |{Q0,F4}| %d, |{Q0,T50}| %d  [%.0fs]" % (pstr(O), len(X), len(Yf), len(Yt), time.time() - t0)); sys.stdout.flush()
u, v = sp.symbols("u v")
R = lambda x: sp.Rational(x.numerator, x.denominator) if isinstance(x, Fr) else sp.Integer(x)
eqs = [R(x) + u * R(a) + v * R(b) for x, a, b in rows if (x, a, b) != (0, 0, 0)]
print("%s: %d monomials B^4 X; Qh2^2 != 0 on %d, {Q0,F4} != 0 on %d, {Q0,T50} != 0 on %d; %d equations"
      % (th.name, len(sample), nX, nYf, nYt, len(eqs)))
sol = sp.solve(eqs, [u, v], dict=True)
print("solution:", sol)
if sol and u in sol[0] and v in sol[0]:
    t4 = sol[0][v] / sol[0][u]
    print("eps4 = %s ; t4 = K_tadpole(5,0)/pi^2 = %s = %s ; MC (ht-d2 mc_summary.log): 0.00159(1)" % (sol[0][u], t4, sp.N(t4, 10)))
