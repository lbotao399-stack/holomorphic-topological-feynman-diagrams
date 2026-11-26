"""Quadratic map [J] -> [(1/2){J,J}] from H^3(1,1) to H^4(2,2), formal gl_N and sl_N.
J = a D + b S1 + c S0 (gl_N), J = a D (sl_N);  D = [p1p2] tr(1{23}), S1 = [p1p2](tr1 tr23 + tr2 tr13 + tr3 tr12),
S0 = [p1p2] tr1 tr2 tr3.  Coordinates in the H^4 basis of h4.py (dual functionals Phi)."""
import sys, time, pickle
from fractions import Fraction as Fr
import sympy as sp
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/deform/code')
import ce, la, spaces

H4 = pickle.load(open((__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/deform/code/h4_data.pkl', 'rb'))
H3 = pickle.load(open((__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/deform/code/h3_data.pkl', 'rb'))
Nsym = sp.Symbol('N')


def npoly_sym(d):
    return sum(sp.Rational(v.numerator, v.denominator) * Nsym ** p for p, v in d.items())


def coords_formal(c, Phi):
    """coordinates (sympy in N) of a formal closed cochain c (Npoly coefficients)."""
    out = []
    for ph in Phi:
        s = {}
        for k, v in ph.items():
            if k in c:
                s = ce.np_add(s, ce.np_scale(c[k], v))
        out.append(sp.expand(npoly_sym(s)))
    return out


t0 = time.time()
RES = {}
for sl in (False, True):
    name = "sl_N" if sl else "gl_N"
    Js = {k.split('=')[0]: spaces.to_np(v) for k, v in H3[(name, 'GL2')].items()}
    Phi = H4[name]['Phi']
    hnames = [nm for nm, c in H4[name]['Hbasis']]
    print("==== %s: H^4 basis %s ====" % (name, hnames), flush=True)
    names = list(Js)
    Q = {}
    for i, X in enumerate(names):
        for Y in names[i:]:
            br = ce.bracket(Js[X], 3, Js[Y], 3, sl=sl)
            closed = ce.d_ce(br, 4, sl) == {}
            co = coords_formal(br, Phi)
            Q[(X, Y)] = co
            print("{%s,%s}: closed %s, weights %s, class = %s" % (X, Y, closed, sorted(ce.weight_of(br)),
                  dict(zip(hnames, [sp.factor(x) for x in co]))), flush=True)
    RES[name] = Q
    # quadratic map (1/2){J,J}
    a, b, c = sp.symbols('a b c')
    coef = {'d': a, 's1': b, 's0': c}
    qJ = [0] * len(hnames)
    for (X, Y), co in Q.items():
        w = coef[X] * coef[Y] * (sp.Rational(1, 2) if X == Y else 1)
        for t in range(len(hnames)):
            qJ[t] += w * co[t]
    qJ = [sp.factor(sp.expand(x)) for x in qJ]
    print("(1/2){J,J} class:", dict(zip(hnames, qJ)))
    RES[name + '_q'] = qJ
pickle.dump(RES, open((__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/deform/code/quad_data.pkl', 'wb'))
print("total %.0fs" % (time.time() - t0))
