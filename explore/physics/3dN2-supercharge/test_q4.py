"""Order hbar^4: Qh2^2 + eps4 {Q0, Qh4} = 0 with Qh4 = Qh4^fan + t4 * T(5,0), on undifferentiated
monomials (kappa = 0).  Unknowns u = eps4, v = eps4 t4 solved exactly from all coefficients.
Also: path-pattern engine (fanengine) vs general engine for n = 2.  Log: test_q4_<theory>.log"""
import sys, itertools, random, time
from fractions import Fraction as Fr
from ht3d import Engine, mono, pstr, psum, pscale, load_kernels, padd
from fanengine import ShapeEngine, fan_kval, tad_kval
from theories import *

name = sys.argv[1]; maxdeg = int(sys.argv[2]); nsample = int(sys.argv[3])
TH = {"u1": lambda: u1([1]), "u1q": lambda: u1([1, 2, -3]), "su2adj": lambda: su2_adj(1),
      "su2pure": lambda: su2_adj(0), "u1pm": lambda: u1([1, -1]), "u1W": lambda: u1([1, 2, -3], W={((0, 0), (1, 0), (2, 0)): 1}), "gl2f": lambda: gl_theory(2, 1, 0), "su2adjW": lambda: su2_adj(2, W={((0, a), (1, a)): 1 for a in range(3)})}
th = TH[name]()
E = Engine(th, load_kernels(4)); A = E.A
F2 = ShapeEngine(th, 2, 'fan', fan_kval(2, Fr(1, 3)))
T2 = ShapeEngine(th, 2, 'tad', tad_kval(2, Fr(1, 3)))
F4 = ShapeEngine(th, 4, 'fan', fan_kval(4, Fr(1, 180)))
T4 = ShapeEngine(th, 4, 'tad', tad_kval(4, Fr(1)))
pool = [('c', a, 0) for a in range(th.G)] + [('B', a, 0) for a in range(th.G)]
for k, m in enumerate(th.multiplets):
    pool += [('p%d' % k, i, 0) for i in range(m[0])] + [('s%d' % k, i, 0) for i in range(m[0])]
monos = []
for deg in range(2, maxdeg + 1):
    for combo in itertools.combinations_with_replacement(pool, deg):
        O = mono(A, *combo)
        if O:
            monos.append(O)
print("theory", th.name, ": %d undifferentiated monomials of degree 2..%d" % (len(monos), maxdeg))
t0 = time.time()
# n = 2 cross-check
bad = nz = 0
for O in monos:
    if sum(len(m) for m in O) > 4:
        continue
    g = E.Q2(O)
    s = psum(F2(O), T2(O))
    nz += bool(g)
    if psum(g, pscale(s, -1)):
        bad += 1
print("n=2: general engine vs path/cycle-pattern engine: mismatches %d (Q2 != 0 on %d)  [%.0fs]" % (bad, nz, time.time() - t0))
rng = random.Random(11)
big = [O for O in monos if len(next(iter(O))) >= 4]
sample = big if len(big) <= nsample else rng.sample(big, nsample)
eqs = []   # (X, Yf, Yt) coefficient triples
nX = nYf = nYt = 0
for O in sample:
    X = E.Q2(E.Q2(O))
    Yf = psum(E.Q0(F4(O)), F4(E.Q0(O)))
    Yt = psum(E.Q0(T4(O)), T4(E.Q0(O)))
    nX += bool(X); nYf += bool(Yf); nYt += bool(Yt)
    for m in set(X) | set(Yf) | set(Yt):
        eqs.append((X.get(m, 0), Yf.get(m, 0), Yt.get(m, 0)))
print("ħ^4: %d monomials (degree >= 4) sampled; Qh2^2 != 0 on %d, {Q0,Qh4fan} != 0 on %d, {Q0,T50} != 0 on %d; %d coefficient equations  [%.0fs]"
      % (len(sample), nX, nYf, nYt, len(eqs), time.time() - t0))
# solve X + u Yf + v Yt = 0
import sympy as sp
u, v = sp.symbols("u v")
sol = None
rows = [(sp.Rational(x.numerator, x.denominator) if isinstance(x, Fr) else sp.Integer(x),
         sp.Rational(yf.numerator, yf.denominator) if isinstance(yf, Fr) else sp.Integer(yf),
         sp.Rational(yt.numerator, yt.denominator) if isinstance(yt, Fr) else sp.Integer(yt)) for x, yf, yt in eqs]
sysm = [x + u * yf + v * yt for x, yf, yt in rows if (x, yf, yt) != (0, 0, 0)]
sol = sp.solve(sysm, [u, v], dict=True)
print("solution of Qh2^2 + u {Q0,Qh4fan} + v {Q0,T(5,0)} = 0:", sol)
if sol:
    s = sol[0]
    if u in s and v in s:
        print("   eps4 = u =", s[u], "; t4 = v/u =", s[v] / s[u], "=", sp.N(s[v] / s[u], 12),
              "; MC (ht-d2 mc_summary): K_tad(5,0)/pi^2 = 0.00159(1)")
    resid = [e.subs(s) for e in sysm]
    print("   residual nonzero equations:", sum(1 for r in resid if sp.simplify(r) != 0))
