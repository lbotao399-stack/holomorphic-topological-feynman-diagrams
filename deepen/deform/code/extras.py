"""Extras (formal N): (1) gl_N quadratic map in the basis (pi^*D_sl, M1, U); (2) class of the double-trace
cochain G = alt([p1p2][p3p4] tr(12) tr(34)) on sl_N; (3) eta with (1/2){D,D}_sl - F_223/(12 N^2) = d eta;
(4) non-closed elements of C^3(1,1); (5) j_M decomposition."""
import sys, pickle
from fractions import Fraction as Fr
import sympy as sp
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/deform/code')
import ce, la, spaces

H4 = pickle.load(open((__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/deform/code/h4_data.pkl', 'rb'))
N = sp.Symbol('N')

# (1)
a, b, c, al, be, ga = sp.symbols('a b c alpha beta gamma')
lam = (2 * a * c + b ** 2 + N * b * c) / 48          # [F_223] coefficient (basis F_223, F_1)
mu = -(N * b + 2 * a) * (N * c + 2 * b) / 48           # [F_1] coefficient
# check against quad.py output: q = -b(Nb+4a)/48 [F_1] + (Nbc+2ac+b^2)/48 [F_0], [F_223] = N[F_1] + [F_0]
q1 = -b * (N * b + 4 * a) / 48
q0 = (N * b * c + 2 * a * c + b ** 2) / 48
print("(1) lam [F_223] + mu [F_1] == q1 [F_1] + q0 [F_0]:",
      sp.simplify(lam * N + mu - q1) == 0 and sp.simplify(lam - q0) == 0)
sub = {a: al, b: be - 2 * al / N, c: ga - 3 * be / N + 4 * al / N ** 2}
print("    J = alpha pi^*D_sl + beta M1 + gamma U, M1 = S1 - (3/N) S0, U = S0, pi^*D_sl = D - (2/N) S1 + (4/N^2) S0:")
print("    lam =", sp.factor(sp.expand(lam.subs(sub, simultaneous=True))), "; mu =", sp.factor(sp.expand(mu.subs(sub, simultaneous=True))))
print("    Moyal j_M = D/2 = (1/2) pi^*D_sl + (1/N) M1 + (1/N^2) U:  lam =", lam.subs({a: sp.Rational(1, 2), b: 0, c: 0}),
      ", mu =", mu.subs({a: sp.Rational(1, 2), b: 0, c: 0}))
# verify the decomposition pi^*D_sl = D - (2/N) S1 + (4/N^2) S0 formally: evaluate D - (2/N)S1 + (4/N^2)S0 on
# a central argument: must vanish (pulled back from sl_N) -- check via restriction X_0 -> identity
import evalx, random
rng = random.Random(3)
u = ce.bracket_poly(0, 1, 3)
Dg = ce.canon_raw(ce.raw_from(u, {((0, 1, 2),): {0: Fr(1)}, ((0, 2, 1),): {0: Fr(1)}}), 3)
S1 = ce.canon_raw(ce.raw_from(u, {((0,), (1, 2)): {0: Fr(1)}, ((0, 2), (1,)): {0: Fr(1)}, ((0, 1), (2,)): {0: Fr(1)}}), 3)
S0 = ce.canon_raw(ce.raw_from(u, {((0,), (1,), (2,)): {0: Fr(1)}}), 3)
for Nn in (3, 4, 5):
    pD = ce.lin((1, Dg), (Fr(-2, Nn), S1), (Fr(4, Nn * Nn), S0))
    ok = True
    for t in range(3):
        mats = [evalx.rand_mat(Nn, 4, rng=rng) for _ in range(3)]
        moms = evalx.rand_moms(3, 4, rng=rng)
        for slot in range(3):
            m2 = list(mats)
            m2[slot] = evalx.ident(Nn)
            ok &= evalx.eval_cochain(pD, m2, moms, Nn) == 0
        # equals D_sl on traceless parts
        tl = [[[x - (evalx.trace(M) / Nn if i == j else 0) for j, x in enumerate(row)] for i, row in enumerate(M)] for M in mats]
        ok &= evalx.eval_cochain(pD, mats, moms, Nn) == evalx.eval_cochain(Dg, tl, moms, Nn)
    print("    N=%d: D - (2/N) S1 + (4/N^2) S0 vanishes with a central argument and equals D(X') on traceless parts: %s" % (Nn, ok))

# (2) double-trace cochain on sl_N
d = H4['sl_N']
Phi = d['Phi']
G = ce.alt_raw(ce.raw_from(ce.poly_mul(ce.bracket_poly(0, 1, 4), ce.bracket_poly(2, 3, 4)),
                           {((0, 1), (2, 3)): {0: Fr(1)}}), 4, sl=True)
print("(2) G = alt([p1p2][p3p4] tr(12)tr(34)) on sl_N: closed", ce.d_ce(G, 4, True) == {},
      "; class in units of [F_1] (= [F_223]/N):", [sum(Fr(ce.np_eval(G.get(k, {}), 1)) * v for k, v in ph.items()) for ph in Phi])

# (3) eta
Dsl = ce.restrict_sl(Dg)
half = ce.scale(ce.bracket(Dsl, 3, Dsl, 3, sl=True), Fr(1, 2))
F = ce.restrict_sl(spaces.obstruction_cochains()[(2, 2, 3)]['F'])
# (1/2){D,D} - F/(12 N^2): Npoly arithmetic
R = ce.lin((1, half), ({-2: Fr(-1, 12)}, F))
C3 = d['C3']
for Nn in (5, 6, 7):
    RN = ce.at_N(R, Nn)
    dC3 = [spaces.const_part(ce.d_ce(spaces.to_np(cc_), 3, True)) for cc_ in C3]
    keys, M = la.matrix_from(dC3 + [RN])
    A = [row[:-1] for row in M]
    bb = [row[-1] for row in M]
    ok, x = la.in_span(A, bb)
    eta = {}
    if ok:
        for j, xj in enumerate(x):
            for k, v in C3[j].items():
                eta[k] = eta.get(k, 0) + xj * v
    print("(3) N=%d: (1/2){D,D}_sl - F_223/(12N^2) = d eta exactly: %s; eta = %s" % (
        Nn, ok, {(ce.tau_str(k[0]), k[1]): str(v) for k, v in eta.items() if v}))

# (4) non-closed elements of C^3(1,1)
for sl in (True, False):
    C3, keys3 = spaces.basis(3, (1, 1), sl)
    print("(4) %s C^3(1,1) basis and d:" % ("sl_N" if sl else "gl_N"))
    for cc_ in C3:
        dd = ce.d_ce(spaces.to_np(cc_), 3, sl)
        print("    %s  ->  d = 0: %s" % ({(ce.tau_str(k[0]), k[1]): str(v) for k, v in cc_.items()}, dd == {}))
