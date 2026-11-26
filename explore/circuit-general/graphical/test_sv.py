"""Checks of the SV polylogarithm layer: x1', Bloch-Wigner, swap, z->1-z, z->1/z, derivatives (numerical)."""
import os, sys, time
import mpmath as mp
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mzv, sv
from rat import R, ONE
mp.mp.dps = 30
t = time.time()
W = sv.WMAX
x1p = sv.x1prime(W)
print('SV_WMAX = %d; x1prime: %d words (%.1fs)' % (W, len(x1p), time.time() - t))
for w in sorted(x1p, key=lambda w: (len(w), w)):
    if len(w) <= 5 and w != '1':
        print('  x1p[%s] = %s' % (w, mzv.tostr(x1p[w])))
t = time.time(); sv.subcoef(W); print('subcoef (%.1fs)' % (time.time() - t))
z = mp.mpc('0.31', '0.27')
# Bloch-Wigner
D = mp.im(mp.polylog(2, z)) + mp.arg(1 - z) * mp.log(abs(z))
v = sv.Pnum('01', z) - sv.Pnum('10', z)
print('P01-P10 - 4iD =', mp.nstr(abs(v - 4j * D), 3))
# swap: (sigma P_w)(z) = P_w(zbar)
err = 0
for w in sv.all_words(min(W, 6)):
    lhs = sum(mzv.numeric(c) * sv.Pnum(s, z) for s, c in sv.swap_word(w).items())
    err = max(err, abs(lhs - sv.Pnum(w, mp.conj(z))))
print('swap max err (|w|<=%d):' % min(W, 6), mp.nstr(err, 3))
# z -> 1-z at z0 with |z0|,|1-z0| < 1
z0 = mp.mpc('0.45', '0.21')
err = 0
for w in sv.all_words(min(W, 6)):
    lhs = sv.Pnum(w, 1 - z0, 25, N=260)
    rhs = sum(mzv.numeric(c) * sv.Pnum(v, z0, 25, N=260) for v, c in sv.trans_1mz(w).items())
    err = max(err, abs(lhs - rhs))
print('1-z max err (|w|<=%d):' % min(W, 6), mp.nstr(err, 3))
# z -> 1/z at z1 = 1.15+0.3i: P_w(1/z1) direct; P_v(z1) via 1-z transform at u = 1-z1
z1 = mp.mpc('1.15', '0.3')
u = 1 - z1
def Pz1(v):
    return sum(mzv.numeric(c) * sv.Pnum(s, u, 25, N=200) for s, c in sv.trans_1mz(v).items())
err = 0
for w in sv.all_words(min(W, 5)):
    lhs = sv.Pnum(w, 1 / z1, 25, N=260)
    rhs = sum(mzv.numeric(c) * Pz1(v) for v, c in sv.trans_inv(w).items())
    err = max(err, abs(lhs - rhs))
print('1/z max err (|w|<=%d):' % min(W, 5), mp.nstr(err, 3))
# derivatives numerically: f = x^2 y/(x-1) P_0110 + P_1001 ...
f = sv.SVF.rat(R({(2, 1): 1}, (0, 0, 1, 0, 0)), '0110') + sv.SVF.rat(ONE, '10010')
h = mp.mpf('1e-12')
def num_dz(f, z):
    # f(z,zb) with independent variation: use Wirtinger: d/dz = (d/dx - i d/dy)/2
    fx = (f.evaluate(z + h) - f.evaluate(z - h)) / (2 * h)
    fy = (f.evaluate(z + 1j * h) - f.evaluate(z - 1j * h)) / (2 * h)
    return (fx - 1j * fy) / 2, (fx + 1j * fy) / 2
mp.mp.dps = 40
a, b = num_dz(f, z)
print('dz err', mp.nstr(abs(a - f.dz().evaluate(z)), 3), ' dzb err', mp.nstr(abs(b - f.dzb().evaluate(z)), 3))
# primitives
g = f.int_z()
print('int_z check:', (g.dz() - f).iszero())
g = f.int_zb()
print('int_zb check:', (g.dzb() - f).iszero())
print('value_at_1 samples:', {w: mzv.tostr(sv.value_at_1(w)) for w in ['10', '100', '110', '1000', '0110', '10100']})
