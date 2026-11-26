"""(g) gl_N formal (kappa = tr), without SL2-invariance of I_2: i_2 = kappa_2 [p1p2]^2 f + y1 h1 + y2 h2,
h1, h2 a basis of Z^3_(2,2) = H^3_(2,2) (C^2 = 0).  Full H^4_(3,3) (dim 9): B^4 = d C^3, Z^4 = ker d4,
basis classes and coordinates (exact solves).  Classes of {J_T, h_i}, T = d, s1, s0;  C-parity of h_i and of the
classes;  classes of theta'(d), {J_T, [p1p2]^2 f} (zero by o11).
[O_3] = eps sum_{T,i} x_T y_i [{J_T, h_i}],  x = (a, b, c)."""
import time
from common import *

t0 = time.time()
sl = False
# Z^3_(2,2)
deg = (2, 2)
C3, k3 = spaces.basis(3, deg, sl)
k4 = spaces.keys_of(4, deg, sl)
img = [spaces.const_part(ce.d_ce(spaces.to_np(c), 3, sl)) for c in C3]
_, M = la.matrix_from(img, k4)
cert = la.certified_rank(M)
Z3 = []
for v in cert['kernel']:
    z = {}
    for j, cj in enumerate(v):
        if cj:
            for kk, x in C3[j].items():
                z[kk] = z.get(kk, 0) + cj * x
    Z3.append({kk: x for kk, x in z.items() if x})
print("gl_N (2,2): dim C^3 = %d, rank d3 = %d, dim Z^3 = H^3 = %d" % (len(C3), cert['rank'], len(Z3)))


def parity(w, k):
    Cw = Cop(w, k)
    if is_zero(ce.lin((1, Cw), (-1, w))):
        return +1
    if is_zero(ce.lin((1, Cw), (1, w))):
        return -1
    return 0


hs = []
for i, z in enumerate(Z3):
    zz = spaces.to_np(z)
    # split into C-even and C-odd parts
    ev = ce.scale(ce.lin((1, zz), (1, Cop(zz, 3))), Fr(1, 2))
    od = ce.scale(ce.lin((1, zz), (-1, Cop(zz, 3))), Fr(1, 2))
    print("  h%d = %s ; C-even part nonzero %s, C-odd part nonzero %s" % (i + 1, show(z, 3, maxterms=8), not is_zero(ev), not is_zero(od)))
    hs.append(zz)

# full H^4_(3,3)
deg = (3, 3)
C3, k3 = spaces.basis(3, deg, sl)
C4, k4 = spaces.basis(4, deg, sl)
k5 = spaces.keys_of(5, deg, sl)
B4 = [spaces.const_part(ce.d_ce(spaces.to_np(c), 3, sl)) for c in C3]
_, MB = la.matrix_from(B4, k4)
cB = la.certified_rank(MB)
B4 = [B4[j] for j in cB['pivcols']]
print("gl_N (3,3): dim C^3 = %d, dim C^4 = %d, rank d3 = %d (det %s)   (t=%.0f)" % (len(C3), len(C4), cB['rank'], cB['det'], time.time() - t0), flush=True)
dC4 = [spaces.const_part(ce.d_ce(spaces.to_np(c), 4, sl)) for c in C4]
_, MZ = la.matrix_from(dC4, k5)
cZ = la.certified_rank(MZ)
Z4 = []
for v in cZ['kernel']:
    z = {}
    for j, cj in enumerate(v):
        if cj:
            for kk, x in C4[j].items():
                z[kk] = z.get(kk, 0) + cj * x
    Z4.append({kk: x for kk, x in z.items() if x})
print("  rank d4 = %d (det %s), dim Z^4 = %d, dim H^4 = %d   (t=%.0f)" % (cZ['rank'], cZ['det'], len(Z4), len(Z4) - len(B4), time.time() - t0), flush=True)
chosen, cur = [], list(B4)
for i, z in enumerate(Z4):
    _, Mc = la.matrix_from(cur + [z], k4)
    if la.certified_rank(Mc)['rank'] > len(cur):
        cur.append(z)
        chosen.append(z)
print("  H^4 basis: %d Z^4 kernel vectors" % len(chosen))
_, MH = la.matrix_from(chosen + B4, k4)
cH = la.certified_rank(MH)
nH = len(chosen)


def npow_parts(c):
    pw = sorted({p for d in c.values() for p in d})
    return {p: {kk: Fr(d[p]) for kk, d in c.items() if d.get(p)} for p in pw}


def cls(c, name, verbose=True):
    out = {}
    for p, x in npow_parts(c).items():
        assert is_zero(ce.d_ce(spaces.to_np(x), 4, sl)), name
        ok, sol = la.in_span(MH, la.column(x, k4), cH)
        assert ok, name
        out[p] = sol[:nH]
    if verbose:
        print("  [%s] = %s" % (name, {p: [str(v) for v in s] for p, s in out.items() if any(s)}), flush=True)
    return out


# C on H^4
Cmat = []
for i, z in enumerate(chosen):
    r = cls(Cop(spaces.to_np(z), 4), "C basis%d" % i, verbose=False)
    Cmat.append([str(v) for v in r.get(0, [0] * nH)])
print("  C on H^4_(3,3) (rows: images of basis classes):")
for row in Cmat:
    print("    ", row)
evals = sp.Matrix([[sp.Rational(x) for x in row] for row in Cmat])
print("  trace C = %s => dim H^4_+ = %s, dim H^4_- = %s" % (evals.trace(), (nH + evals.trace()) / 2, (nH - evals.trace()) / 2))
th = load('theta.pkl')
cls(th['th']['d'], "theta'(d)")
u2f = f_bracket_sq(sl)
J = {T: J_cochain(T, sl) for T in ('d', 's1', 's0')}
for T in J:
    cls(ce.bracket(J[T], 3, u2f, 3, sl=sl), "{J_%s, [p1p2]^2 f}" % T)
R = {}
for T in J:
    for i, hh in enumerate(hs):
        R[(T, i)] = cls(ce.bracket(J[T], 3, hh, 3, sl=sl), "{J_%s, h%d}" % (T, i + 1))
save(dict(Z3=Z3, chosen=chosen, B4=B4, R=R, Cmat=Cmat), 'gl_full.pkl')
print("time %.0fs" % (time.time() - t0))
