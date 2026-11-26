"""Exact interpolation of n4(N) (kappa = tr_V, degree <= 8, family.py) through the 9 smallest N, check at all further N;
Killing ratio r(N) (degree 1 through 2 points, checked at all N); dimension d(N); then
8 n4_K - 3d with n4_K = n4(tr) / r^4 (Killing normalization), kk-coefficient a = (8 n4_K - 3d)/(4d(d-3)),
mu(F223) = 2 d (8 n4_K - 3d) (logs/b_analysis.log), factored."""
import json, sys
import sympy as sp

N = sp.Symbol('N')
LOG = (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/simple-all/logs/'
dims = {'so': N * (N - 1) / 2, 'sp': N * (N + 1) / 2, 'sl': N ** 2 - 1}
for kind in ('so', 'sp', 'sl'):
    data = json.load(open(LOG + 'family_%s.json' % kind))
    pts = [(sp.Integer(x['N']), sp.Rational(x['n4']), sp.Rational(x['r']), x['dim']) for x in data]
    fit = pts[:9]
    n4 = sp.expand(sp.interpolate([(a, b) for a, b, _, _ in fit], N))
    r = sp.expand(sp.interpolate([(a, c) for a, _, c, _ in pts[:2]], N))
    ok_n4 = all(n4.subs(N, a) == b for a, b, _, _ in pts)
    ok_r = all(r.subs(N, a) == c for a, _, c, _ in pts)
    ok_d = all(dims[kind].subs(N, a) == dd for a, _, _, dd in pts)
    d = dims[kind]
    n4K = n4 / r ** 4
    X = sp.factor(sp.simplify(8 * n4K - 3 * d))
    a = sp.factor(sp.simplify((8 * n4K - 3 * d) / (4 * d * (d - 3))))
    mu = sp.factor(sp.simplify(2 * d * (8 * n4K - 3 * d)))
    print("%s: N in %s" % (kind, [int(p[0]) for p in pts]))
    print("   n4(tr)(N) = %s   [interpolated through N = %s; equal at all %d computed N: %s]"
          % (sp.factor(n4), [int(p[0]) for p in fit], len(pts), ok_n4))
    print("   Killing/tr_V = %s (all N: %s) ; d = %s (all N: %s)" % (r, ok_r, d, ok_d))
    print("   8 n4_K - 3d = %s" % X)
    print("   a = (8n4_K-3d)/(4d(d-3)) = %s" % a)
    print("   mu(F223) = 2d(8n4_K - 3d) = %s" % mu)
    num = sp.numer(sp.together(8 * n4K - 3 * d))
    print("   integer roots of the numerator: %s" % sorted(set(x for x in sp.roots(sp.Poly(num, N)).keys() if x.is_integer)))
