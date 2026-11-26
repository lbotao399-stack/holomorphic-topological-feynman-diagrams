"""(f) Charge conjugation C: X -> -X^T (automorphism of gl_N and sl_N preserving kappa = tr).
On trace-basis k-cochains: (C w)(x_1..x_k) = w(C x_1, .., C x_k), tau -> (-1)^k tau^rev (each cycle reversed),
momenta unchanged.  C f = f, C d = -d (d_abc), C D = -D, C S1 = -S1, C S0 = -S0, C commutes with d and with {,}.
Checks (formal N, exact): C d w = d C w and C{a,b} = {Ca,Cb} on the cochains used; parities of
j (= D), i_2 (iota), o', theta', {D,iota}; C acting on H^4_(3,3) (sl_N; gl_N if data present)."""
import sys, time
from common import *

t0 = time.time()


def rev_tau(t):
    return ce.canon_tau([tuple(reversed(c)) for c in t])


def Cop(w, k):
    s = (-1) ** k
    out = {}
    for (t, m), d in w.items():
        ce.add_term(out, (rev_tau(t), m), d if isinstance(d, dict) else {0: Fr(d)}, s)
    return out


def parity(w, k):
    Cw = Cop(w, k)
    if is_zero(ce.lin((1, Cw), (-1, w))):
        return +1
    if is_zero(ce.lin((1, Cw), (1, w))):
        return -1
    return 0


th = load('theta.pkl')
i2 = load('i2_sl.pkl')
for sl in (True, False):
    nm = 'sl_N' if sl else 'gl_N'
    R = (lambda c: ce.restrict_sl(c)) if sl else (lambda c: c)
    D = J_cochain('d', sl)
    u2f = f_bracket_sq(sl)
    print("%s parities: f-cochain s_0 %+d, D %+d, [p1p2]^2 f %+d, o' %+d, theta'(d) %+d" % (
        nm, parity(ce.s3_cochain() if not sl else ce.restrict_sl(ce.s3_cochain()), 3), parity(D, 3), parity(u2f, 3),
        parity(R(th['op']), 4), parity(R(th['th']['d']), 4)))
    if not sl:
        print("  gl_N: S1 %+d, S0 %+d" % (parity(J_cochain('s1'), 3), parity(J_cochain('s0'), 3)))
    if sl:
        br = ce.bracket(D, 3, i2['iota'], 3, sl=True)
        print("  {D, iota} %+d" % parity(br, 4))
    # commutation checks
    w = R(th['th']['d'])
    print("  C d theta' == d C theta':", is_zero(ce.lin((1, Cop(ce.d_ce(w, 4, sl), 5)), (-1, ce.d_ce(Cop(w, 4), 4, sl)))))
    a = u2f
    print("  C {D, u2f} == {C D, C u2f}:", is_zero(ce.lin((1, Cop(ce.bracket(D, 3, a, 3, sl=sl), 4)),
                                                        (-1, ce.bracket(Cop(D, 3), 3, Cop(a, 3), 3, sl=sl)))))

# C on H^4_(3,3)(sl_N)
H = load('h4w3_sl.pkl')
B4, Z4, k4, Phi = H['B4'], H['Z4'], H['k4'], H['Phi']
coord = lambda c: sum(Fr(c.get(kk, 0)) * v for kk, v in Phi.items())
h = H['h'][1]
Ch = spaces.const_part(Cop(spaces.to_np(h), 4))
print("sl_N weight (3,3): [C h] = %s [h]  (h = %s)" % (coord(Ch), H['h'][0]))
print("  Z^4 basis parities:", [parity(spaces.to_np(z), 4) for z in Z4])
print("  B^4 basis parities:", [parity(spaces.to_np(b), 4) for b in B4])
print("time %.1fs" % (time.time() - t0))
