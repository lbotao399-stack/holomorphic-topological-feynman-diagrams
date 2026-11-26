"""(1) {J_T, o'} and {J_T, F_G} as cochains (explorer e1: identically zero), gl_N and sl_N; control {D, r}, r random.
(2) d{j,i_2} identity at the cochain level for sl_N: d{D,iota} = -{D, d iota} = -{D, o' - (N^2/8){D,D}}.
(3) SL2-invariant C-odd subcomplex at (3,3) for gl_N (N = 7): dim C^3, C^4, rank d3.
(4) position-space check of the derivative rule, H = 1, banana(2) at fixed Schwinger times (scipy, float):
    int d^2z [d/dz_tail e^{-|z_tail - z_head|^2/4t_0}] e^{-|z|^2/4t_1} e^{lam z}  /  int e^{...}  = -y_0(htk),
    derivative at the head (root, position w -> 0): ratio = +y_1(htk)."""
import pickle, random, sys, time
from fractions import Fraction as Fr
import numpy as np
from vlib import *
import vla

t0 = time.time()
for sl, fn in ((True, 'data/assembled_sl.pkl'), (False, 'data/assembled_gl.pkl')):
    A = pickle.load(open(fn, 'rb'))
    op = A['op']
    tag = 'sl' if sl else 'gl'
    Ts = [('D', cubic('d', 1, sl))] + ([] if sl else [('S1', cubic('s1', 1, sl)), ('S0', cubic('s0', 1, sl))])
    for nm, J in Ts:
        print("%s: {%s, o'} = 0: %s;  {%s, F_133} = 0: %s;  {%s, F_223} = 0: %s" %
              (tag, nm, is_zero(bracket(J, 3, op, 4, sl)), nm, is_zero(bracket(J, 3, A['F'][(1, 3, 3)], 4, sl)),
               nm, is_zero(bracket(J, 3, A['F'][(2, 2, 3)], 4, sl))), flush=True)
    random.seed(4)
    r = {}
    perms = all_perms(4, sl)
    monos = momentum_monos(4, (2, 2))
    for _ in range(5):
        r[(random.choice(perms), random.choice(monos))] = LP.const(random.randint(1, 5))
    r = alt(r, 4, sl)
    print("%s: control {D, r} (r random alternating 4-cochain of weight (2,2)): number of nonzero coordinates %d" %
          (tag, len(bracket(cubic('d', 1, sl), 3, r, 4, sl))), flush=True)
    if sl:
        D = cubic('d', 1, sl)
        DI = A['DI']
        lhs = d(DI, 4, sl)
        rhs = scal(-1, bracket(D, 3, lin((1, op), (-LP.N(2) * Fr(1, 8), A['DD'])), 4, sl))
        print("sl: d{D,iota} = -{D, o' - (N^2/8){D,D}}: %s;  d{D,iota} = 0: %s;  {D,{D,D}} = 0: %s" %
              (is_zero(lin((1, lhs), (-1, rhs))), is_zero(lhs), is_zero(bracket(D, 3, A['DD'], 4, sl))), flush=True)

# (3)
N = 7
SL = False
for par in (-1,):
    C3 = span_set(3, (3, 3), SL, sl2inv=True, parity=par)
    C4 = span_set(4, (3, 3), SL, sl2inv=True, parity=par)
    R3 = vla.to_rows(C3, N); R4 = vla.to_rows(C4, N)
    print("gl_N (N=7) SL2-inv C-odd (3,3): dim C^3 %d, dim C^4 %d, rank d3 %d -> C^4 = d C^3: %s" %
          (vla.rankQ(R3), vla.rankQ(R4), vla.rankQ(vla.to_rows([d(w, 3, SL) for w in C3], N)),
           vla.rankQ(R4) == vla.rankQ(vla.to_rows([d(w, 3, SL) for w in C3], N))), flush=True)

# (4)
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
import htk
from scipy import integrate
G = htk.banana(2)                     # edges (0,1), (0,1); root = 1
for (t0_, t1_) in [(0.7, 1.3), (0.25, 2.0)]:
    y = htk.brute.currents(G, [Fr(t0_).limit_denominator(100), Fr(t1_).limit_denominator(100)], [[Fr(1)]])
    y0, y1 = float(y[0][0]), float(y[1][0])          # currents per unit source (H = 1)
    for lam in [0.8 + 0.3j, -0.5 + 1.1j]:
        def base(x, yy):
            z = x + 1j * yy
            return np.exp(-abs(z) ** 2 / (4 * t0_) - abs(z) ** 2 / (4 * t1_) + lam * z)
        def tail(x, yy):        # d/dz_tail of e^{-|z - 0|^2/4t0} = -(zb/(4 t0)) e^{..}
            z = x + 1j * yy
            return -(np.conj(z) / (4 * t0_)) * base(x, yy)
        def head(x, yy):        # d/dw at w = 0 of e^{-|z - w|^2/4t1} = +((zb - wb)/(4 t1)) e^{..}
            z = x + 1j * yy
            return (np.conj(z) / (4 * t1_)) * base(x, yy)
        def cint(f):
            re = integrate.dblquad(lambda b, a: f(a, b).real, -12, 12, -12, 12, epsabs=1e-13, epsrel=1e-13)[0]
            im = integrate.dblquad(lambda b, a: f(a, b).imag, -12, 12, -12, 12, epsabs=1e-13, epsrel=1e-13)[0]
            return re + 1j * im
        I0, It, Ih = cint(base), cint(tail), cint(head)
        rt, rh = It / I0, Ih / I0
        pt, ph = -y0 * lam, +y1 * lam
        print("banana(2), t = (%.2f, %.2f), lam = %s: tail ratio %s vs -y_0 lam %s (digits %.1f); head ratio %s vs +y_1 lam %s (digits %.1f)"
              % (t0_, t1_, lam, np.round(rt, 12), np.round(pt, 12), -np.log10(abs(rt - pt) / abs(pt)),
                 np.round(rh, 12), np.round(ph, 12), -np.log10(abs(rh - ph) / abs(ph))), flush=True)
print("time %.1fs" % (time.time() - t0))
