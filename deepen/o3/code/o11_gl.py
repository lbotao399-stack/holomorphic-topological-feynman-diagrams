"""(g) gl_N, kappa = tr, formal N (trace basis; faithful for N >= 5), SL2 (U(2))-equivariant setting.
J_1 = a D + b S1 + c S0  (j = a D + b S1 + c S0 as cochains).
Part A, weight (2,2):  SL2-invariant subcomplex C^3_inv, C^4_inv, Z^4_inv; decomposition of the brackets
   (1/2){J_T,J_T}, {J_T,J_T'} = A F_0 + B F_1 + E beta  (F_223 = F_0 + N F_1, beta = d([p1p2]^2 f)), exact;
   (1/2){j,j} = lam F_223 + mu F_1 + e(a,b,c) beta;  o' = (1/48) F_223 - (N^2/40) beta;
   H^3_(2,2)^{SL2} = 0, so i_2 = kappa_2 [p1p2]^2 f,  kappa_2 = e(a,b,c) - eps nu' N^2/40  (nu' = pi nu) when
   mu = 0, lam = -eps nu'/48.
Part B, weight (3,3): SL2-invariant subcomplex, H^4_inv, basis class, C-parity; classes of theta'(d),
   {D, u2f}, {S1, u2f}, {S0, u2f}  (u2f = [p1p2]^2 f);
   [O_3] = eps kappa_2 (a[{D,u2f}] + b[{S1,u2f}] + c[{S0,u2f}]) + nu' a [theta'].
Saves data/gl_inv.pkl."""
import time
from common import *

t0 = time.time()
sl = False
obs = spaces.obstruction_cochains()
F223 = obs[(2, 2, 3)]['F']
F133 = obs[(1, 3, 3)]['F']
F0 = {k: {0: d[0]} for k, d in F223.items() if d.get(0)}
F1 = {k: {0: d[1]} for k, d in F223.items() if d.get(1)}
assert set(p for d in F223.values() for p in d) <= {0, 1}
u2f = f_bracket_sq(sl)
beta = ce.d_ce(u2f, 3, sl)
th = load('theta.pkl')
op = th['op']
J = {T: J_cochain(T, sl) for T in ('d', 's1', 's0')}


def npow_parts(c):
    pw = sorted({p for d in c.values() for p in d})
    return {p: {kk: Fr(d[p]) for kk, d in c.items() if d.get(p)} for p in pw}


def solve_in(c, basis, keys):
    """c (Npoly coefficients) = sum_p N^p sum_i x_{p,i} basis_i: returns {p: x} or None."""
    out = {}
    _, M = la.matrix_from(basis, keys)
    cert = la.certified_rank(M)
    assert cert['rank'] == len(basis)
    for p, x in npow_parts(c).items():
        ok, sol = la.in_span(M, la.column(x, keys), cert)
        if not ok:
            return None
        out[p] = sol
    return out


print("Part A: weight (2,2), formal gl_N")
deg = (2, 2)
k3, k4, k5 = spaces.keys_of(3, deg, sl), spaces.keys_of(4, deg, sl), spaces.keys_of(5, deg, sl)
S3 = inv_span(3, sl, deg)
_, M3 = la.matrix_from(S3, k3)
c3 = la.certified_rank(M3)
C3i = [S3[j] for j in c3['pivcols']]
print("  dim C^3_inv(2,2) = %d ; u2f in it: %s ; d u2f != 0: %s  =>  H^3_(2,2)^{SL2} = 0"
      % (len(C3i), la.in_span(M3, la.column(spaces.const_part(u2f), k3), c3)[0], not is_zero(beta)))
S4 = inv_span(4, sl, deg)
_, M4 = la.matrix_from(S4, k4)
c4 = la.certified_rank(M4)
C4i = [S4[j] for j in c4['pivcols']]
dC4 = [spaces.const_part(ce.d_ce(spaces.to_np(c), 4, sl)) for c in C4i]
_, M5 = la.matrix_from(dC4, k5)
c5 = la.certified_rank(M5)
print("  dim C^4_inv(2,2) = %d, rank d4|inv = %d, dim Z^4_inv(2,2) = %d" % (len(C4i), c5['rank'], len(C4i) - c5['rank']))
basisA = [spaces.const_part(F0), spaces.const_part(F1), spaces.const_part(beta)]
_, MA = la.matrix_from(basisA, k4)
print("  rank{F_0, F_1, beta} = %d" % la.certified_rank(MA)['rank'])
Tn = ['d', 's1', 's0']
dec = {}
for i, T in enumerate(Tn):
    for T2 in Tn[i:]:
        X = ce.bracket(J[T], 3, J[T2], 3, sl=sl)
        if T == T2:
            X = ce.scale(X, Fr(1, 2))
        s = solve_in(X, basisA, k4)
        dec[(T, T2)] = s
        print("  %s{J_%s, J_%s} = %s   (N-power: [F_0, F_1, beta] coefficients)" % (
            "(1/2)" if T == T2 else "", T, T2, {p: [str(v) for v in x] for p, x in s.items()} if s else ("0 (zero cochain: %s)" % is_zero(X) if s == {} else "NOT in span")), flush=True)
Nn = sp.Symbol('N')
a, b, c = sp.symbols('a b c')
x = {'d': a, 's1': b, 's0': c}
coef = [0, 0, 0]
for (T, T2), s in dec.items():
    for p, v in s.items():
        for j in range(3):
            coef[j] += x[T] * x[T2] * Nn ** p * sp.Rational(v[j].numerator, v[j].denominator)
A_, B_, E_ = [sp.expand(q) for q in coef]
lam = sp.factor(A_)
mu = sp.factor(sp.expand(B_ - Nn * A_))
print("  (1/2){j,j} = lam F_223 + mu F_1 + e beta:  lam = %s ; mu = %s ; e = %s" % (lam, mu, sp.factor(E_)))
print("  deform d6: lam == (2ac + b^2 + Nbc)/48: %s ; mu == -(2a+Nb)(2b+Nc)/48: %s" % (
    sp.expand(lam - (2 * a * c + b ** 2 + Nn * b * c) / 48) == 0,
    sp.expand(mu + (2 * a + Nn * b) * (2 * b + Nn * c) / 48) == 0))
so = solve_in(op, basisA, k4)
print("  o' = %s  (N-power: [F_0, F_1, beta])" % {p: [str(v) for v in xx] for p, xx in so.items()})

print("Part B: weight (3,3), formal gl_N, SL2-invariant subcomplex")
deg = (3, 3)
k3, k4, k5 = spaces.keys_of(3, deg, sl), spaces.keys_of(4, deg, sl), spaces.keys_of(5, deg, sl)
S3 = inv_span(3, sl, deg)
_, M3 = la.matrix_from(S3, k3)
c3 = la.certified_rank(M3)
C3i = [S3[j] for j in c3['pivcols']]
B4i_all = [spaces.const_part(ce.d_ce(spaces.to_np(cc), 3, sl)) for cc in C3i]
B4i_all = [v for v in B4i_all if v]
_, MB = la.matrix_from(B4i_all, k4)
cB = la.certified_rank(MB)
B4i = [B4i_all[j] for j in cB['pivcols']]
S4 = inv_span(4, sl, deg)
_, M4 = la.matrix_from(S4, k4)
c4 = la.certified_rank(M4)
C4i = [S4[j] for j in c4['pivcols']]
dC4 = [spaces.const_part(ce.d_ce(spaces.to_np(cc), 4, sl)) for cc in C4i]
_, M5 = la.matrix_from(dC4, k5)
c5 = la.certified_rank(M5)
Z4i = []
for v in c5['kernel']:
    z = {}
    for j, cj in enumerate(v):
        if cj:
            for kk, val in C4i[j].items():
                z[kk] = z.get(kk, 0) + cj * val
    Z4i.append({kk: val for kk, val in z.items() if val})
print("  dim C^3_inv = %d, rank d3|inv = %d, dim C^4_inv = %d, rank d4|inv = %d, dim Z^4_inv = %d, dim H^4_inv = %d   (t=%.0f)"
      % (len(C3i), cB['rank'], len(C4i), c5['rank'], len(Z4i), len(Z4i) - cB['rank'], time.time() - t0), flush=True)
h = None
for i, z in enumerate(Z4i):
    _, M = la.matrix_from(B4i + [z], k4)
    if la.certified_rank(M)['rank'] > len(B4i):
        h = z
        print("  basis class h_gl = Z^4_inv vector %d:" % i, show(h, 4, maxterms=10))
        break
_, MH = la.matrix_from([h] + B4i, k4)
cH = la.certified_rank(MH)


def cls(cc, name):
    out = {}
    for p, xx in npow_parts(cc).items():
        assert is_zero(ce.d_ce(spaces.to_np(xx), 4, sl)), name
        ok, sol = la.in_span(MH, la.column(xx, k4), cH)
        assert ok, "not in Z^4_inv: " + name
        out[p] = sol[0]
    print("  [%s] = %s [h_gl]   (N-power: coefficient)" % (name, {p: str(v) for p, v in out.items() if v}), flush=True)
    return out


Ch = spaces.const_part(Cop(spaces.to_np(h), 4))
cC = cls(spaces.to_np(Ch), "C h_gl")
r = {}
r['theta'] = cls(th['th']['d'], "theta'(d)")
for T in Tn:
    r[T] = cls(ce.bracket(J[T], 3, u2f, 3, sl=sl), "{J_%s, [p1p2]^2 f}" % T)
save(dict(dec=dec, so=so, lam=lam, mu=mu, e=E_, h=h, classes=r, C3i=C3i, C4i=C4i, B4i=B4i, Z4i=Z4i), 'gl_inv.pkl')
print("time %.0fs" % (time.time() - t0))
