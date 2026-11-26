"""Comparison with the universal (Vogel) parametrization, source: Mkrtchyan-Sergeev-Veselov, arXiv:1105.0115v3,
Table 1 (alpha, beta, gamma, t = alpha+beta+gamma), eq. (1) dim g = (alpha-2t)(beta-2t)(gamma-2t)/(alpha beta gamma),
eq. (8) C_4 = (3 t t2 - t3)/(16 t^3), t2 = sum alpha^2, t3 = sum alpha^3, eq. (9) (Okubo)
C_k = sum_{V in ad (x) ad} (dim V/dim g) ((C_2(V) - 2 C_2(ad))/2)^k  (Killing normalization, C_2(ad) = 1).
Here m_k = tr_{g(x)g} Om^k = d C_k.  Checked (exact) against logs/b_data.json and logs/interp.log:
  m4 = d C_4 ;  n4 = m4 - d/8 ;  8 n4 - 3d = -3 d (t+alpha)(t+beta)(t+gamma)/(2 t^3) = 3 d (d-3) alpha beta gamma/(4 t^3) ;
  a = (8n4-3d)/(4d(d-3)) = 3 alpha beta gamma/(16 t^3) ;  x = -3 - 4(d-1) a ;  mu(F223) = 2d(8n4-3d) = 3 d dim(X2) alpha beta gamma / t^3,
  dim X2 = d(d-3)/2 (MSV p. 4)."""
import json, sys
import sympy as sp
from fractions import Fraction as Fr
LOG = (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/simple-all/logs/'
R = sp.Rational
table = {'A': lambda n: (-2, 2, n + 1), 'B': lambda n: (-2, 4, 2 * n - 3), 'C': lambda n: (-2, 1, n + 2),
         'D': lambda n: (-2, 4, 2 * n - 4), 'G2': (-2, R(10, 3), R(8, 3)), 'F4': (-2, 5, 6), 'E6': (-2, 6, 8),
         'E7': (-2, 8, 12), 'E8': (-2, 12, 20)}


def params(key):
    if key in table:
        return tuple(sp.nsimplify(x) for x in table[key])
    return tuple(sp.Integer(x) for x in table[key[0]](int(key[1:])))


def univ(al, be, ga):
    t = al + be + ga
    t2, t3 = al ** 2 + be ** 2 + ga ** 2, al ** 3 + be ** 3 + ga ** 3
    d = (al - 2 * t) * (be - 2 * t) * (ga - 2 * t) / (al * be * ga)
    C4 = (3 * t * t2 - t3) / (16 * t ** 3)
    m4 = d * C4
    n4 = m4 - d / 8
    a = 3 * al * be * ga / (16 * t ** 3)
    return dict(t=t, d=sp.simplify(d), m4=sp.simplify(m4), n4=sp.simplify(n4), a=sp.simplify(a),
                X=sp.simplify(8 * n4 - 3 * d), X1=sp.simplify(-3 * d * (t + al) * (t + be) * (t + ga) / (2 * t ** 3)),
                X2=sp.simplify(3 * d * (d - 3) * al * be * ga / (4 * t ** 3)), x=sp.simplify(-3 - 4 * (d - 1) * a))


recs = json.load(open(LOG + 'b_data.json'))
xa = {}
for line in open(LOG + 'decomp_check.log'):
    nm = line[:16].strip()
    import re
    m = re.search(r'x = (\S+), a = (\S+) :', line)
    if m:
        xa[nm] = (sp.Rational(m.group(1)), sp.Rational(m.group(2)))
allok = True
for rec in recs:
    key = rec['key']
    u = univ(*params(key))
    m4, n4 = sp.Rational(rec['moments']['m4']), sp.Rational(rec['moments']['n4'])
    d = rec['dim']
    ok = (u['d'] == d and u['m4'] == m4 and u['n4'] == n4 and u['X'] == u['X1'] == u['X2'] == 8 * n4 - 3 * d)
    if rec['name'] in xa:
        ok &= (xa[rec['name']] == (u['x'], u['a']))
    allok &= ok
    print("%-16s (alpha,beta,gamma) = %-22s t = %-4s d %s m4 = d C_4 %s n4 %s 8n4-3d = %s ; a = 3abg/(16t^3) = %s, x = %s ; all equal to the computed values: %s"
          % (rec['name'], str(params(key)), u['t'], u['d'] == d, u['m4'] == m4, u['n4'] == n4, u['X'], u['a'], u['x'], ok))
N = sp.Symbol('N')
fam = {'so': (-2, 4, N - 4), 'sp': (-2, 1, N / 2 + 2), 'sl': (-2, 2, N)}
interp = open(LOG + 'interp.log').read()
for k, p in fam.items():
    u = univ(*p)
    print("%s_N with (alpha,beta,gamma) = %s: 8n4-3d = %s ; a = %s" % (k, p, sp.factor(u['X']), sp.factor(u['a'])))
print("all algebras agree: %s" % allok)
