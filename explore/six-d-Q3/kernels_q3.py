"""Three-loop kernels of Q_3: densities in the ratios (product of local factors, S2 Section 7.2),
displacement maps, zero-displacement values.  Sources at the five degree-two vertices only
(lambda = 0 at the two trivalent vertices).  Writes kernels_q3.json for the operator check."""
import json
import sys
import time
from itertools import combinations

import sympy as sp

from q3lib import GRAPHS, Reduction, LAM, Z, R, integrate_domain, outputs, tlabel, ROOT

SIGN = {'G1': 1, 'G2': -1, 'G3': 1, 'G4': 1, 'G5': -1}      # s_Gamma from check_kernels.py


def lin(vecexpr, srcs):
    """Coefficients c_u with vec = sum_u c_u lambda_u."""
    out = {}
    for u in srcs:
        c = sp.expand(sp.diff(vecexpr[0], LAM[u][0]))
        for a in range(1, 3):
            assert sp.expand(sp.diff(vecexpr[a], LAM[u][a]) - c) == 0
        if c != 0:
            out[u] = c
    return out


def brackets(vs, srcs):
    """det(v1, v2, v3) = sum_{i<j<k} c_{ijk} <ijk>."""
    cs = [lin(v, srcs) for v in vs]
    out = {}
    for u1, c1 in cs[0].items():
        for u2, c2 in cs[1].items():
            for u3, c3 in cs[2].items():
                trip = (u1, u2, u3)
                if len(set(trip)) < 3:
                    continue
                key = tuple(sorted(trip))
                perm = [key.index(x) for x in trip]
                sgn = 1 if perm in ([0, 1, 2], [1, 2, 0], [2, 0, 1]) else -1
                out[key] = sp.expand(out.get(key, 0) + sgn * c1 * c2 * c3)
    return {k: sp.factor(v) for k, v in out.items() if v != 0}


def bstr(key):
    return '<' + ''.join(tlabel(v) for v in key) + '>'


def show_lin(d):
    return ' + '.join('(%s) %s' % (sp.factor(c), bstr(k if isinstance(k, tuple) else (k,))) for k, c in sorted(d.items()))


out = {}
names = sys.argv[1:] or list(GRAPHS)
for name in names:
    t0 = time.time()
    red = Reduction(name)
    outs = outputs(name)
    zero = [v for v in range(ROOT) if v not in outs]
    sub0 = {LAM[v][a]: 0 for v in zero for a in range(3)}
    print('=' * 100)
    print(name, ' edges', [(tlabel(a), tlabel(b)) for a, b in red.edges])
    print('   spokes (letter j: edge, z_e):', [(j + 1, (tlabel(red.edges[e][0]), tlabel(red.edges[e][1])), '%+dz_%d' % (s, j + 1))
                                           for j, (e, s) in enumerate(red.spokes)])
    print('   outputs at', [tlabel(v) for v in outs], ' s_Gamma =', SIGN[name])
    # moves with reduced sources (general sources), then specialized
    facs_b = []
    for k, st in enumerate(red.steps):
        r = st['r']
        if st['kind'] != 'box':
            lv = [sp.expand(x.subs(sub0)) for x in st['lv']]
            lw = [sp.expand(x.subs(sub0)) for x in st['lw']]
            Y = [sp.expand(x.subs(sub0)) for x in st['Y']]
            pre = r[0] ** 2 if st['kind'] == 'chain' else r[0] * (1 - r[0])
            print('   move %d: %s (p,v,v\',q) = (%s,%s,%s,%s), ratios %s' % (k + 1, st['kind'], tlabel(st['p']), tlabel(st['v']),
                                                                       tlabel(st['w']), tlabel(st['q']), r))
            print('       lambda_v =', show_lin(lin(lv, outs)).replace('<', '').replace('>', ''), ' | lambda_v\' =',
                  show_lin(lin(lw, outs)).replace('<', '').replace('>', ''))
            print('       y_red(p->q) =', show_lin(lin(Y, outs)).replace('<', '').replace('>', ''))
            b = brackets([Y, lv, lw], outs)
            facs_b.append((pre, b))
            print('       factor %s * det(y_red, lambda_v, lambda_v\') = %s * [%s]' % (st['kind'], pre, show_lin(b)))
        else:
            ls = [[sp.expand(x.subs(sub0)) for x in l] for l in st['lam']]
            print('   move %d: box O-%s-%s-%s-O, ratios %s, 1 > %s > %s > %s > 0' % ((k + 1,) + tuple(tlabel(v) for v in st['verts']) + (r,) + r))
            for v, l in zip(st['verts'], ls):
                print('       lambda_%s^red =' % tlabel(v), show_lin(lin(l, outs)).replace('<', '').replace('>', ''))
            b = brackets(ls, outs)
            facs_b.append((1, b))
            print('       factor box det(lambda_a, lambda_b, lambda_c) = [%s]' % show_lin(b))
    # bracket expansion of the product; numerical (exact rational) comparison with the product of
    # the local factors of Reduction.factors() at random points; integral of each coefficient
    dom = red.domain()
    terms = {}
    for b1, c1 in facs_b[0][1].items():
        for b2, c2 in facs_b[1][1].items():
            for b3, c3 in facs_b[2][1].items():
                key = (b1, b2, b3)
                terms[key] = sp.expand(SIGN[name] * facs_b[0][0] * facs_b[1][0] * facs_b[2][0] * c1 * c2 * c3)
    import random
    rng = random.Random(5)
    Jprod = SIGN[name] * sp.prod(red.factors())
    for _ in range(3):
        lamv = {LAM[v][a]: (sp.Integer(rng.randint(-5, 5)) if v in outs else 0) for v in range(ROOT) for a in range(3)}
        rv = {r: sp.Rational(rng.randint(1, 50), 51) for r in R}
        def bnum(key):
            return sp.Matrix.hstack(*[sp.Matrix(LAM[v]) for v in key]).subs(lamv).det()
        lhs = sum(c.subs(rv) * bnum(b1) * bnum(b2) * bnum(b3) for (b1, b2, b3), c in terms.items())
        assert sp.simplify(lhs - Jprod.subs(lamv).subs(rv)) == 0, 'bracket expansion'
    ints = {key: integrate_domain(c, dom) for key, c in terms.items()}
    print('   density J_Gamma = s_Gamma * prod(factors), I_Gamma(lambda; 0) = int_D J_Gamma; bracket expansion checked at 3 points')
    print('   zero-displacement value K(lambda; 0) = sum over products of brackets:')
    agg = {}
    for (b1, b2, b3), v in ints.items():
        if v != 0:
            agg[(b1, b2, b3)] = v
    for key, v in sorted(agg.items()):
        print('       %12s  %s %s %s' % (v, bstr(key[0]), bstr(key[1]), bstr(key[2])))
    bsym = lambda key: sp.Symbol('B' + ''.join(tlabel(v) for v in key))
    K0b = sum(v * bsym(k[0]) * bsym(k[1]) * bsym(k[2]) for k, v in agg.items())
    print('   factored in bracket symbols (brackets treated as independent):', sp.factor(K0b))
    Jb = sum(c * bsym(k[0]) * bsym(k[1]) * bsym(k[2]) for k, c in terms.items())
    print('   J_Gamma in bracket symbols, factored:', sp.factor(Jb))
    # displacement
    X, expo = red.displacement()
    print('   displacement maps X_v (v outputs):')
    for v in outs:
        cs = []
        for j in range(4):
            c = sp.factor(sp.diff(X[v][0], Z[j][0]))
            cs.append(c)
        print('       X_%s = ' % tlabel(v) + ' + '.join('(%s) z_%d' % (c, j + 1) for j, c in enumerate(cs)),
              '   sum of coefficients:', sp.simplify(sum(cs)))
    out[name] = dict(sign=SIGN[name], outputs=outs,
                     terms=[[list(k[0]), list(k[1]), list(k[2]), str(c)] for k, c in terms.items()],
                     K0=[[list(k[0]), list(k[1]), list(k[2]), str(v)] for k, v in agg.items()],
                     X={str(v): [str(sp.expand(sp.diff(X[v][0], Z[j][0]))) for j in range(4)] for v in outs},
                     domain=[(k, [str(r) for r in rs]) for k, rs in dom],
                     factors=[str(f) for f in red.factors()])
    print('   [%.1fs]' % (time.time() - t0), flush=True)

if not sys.argv[1:]:
    json.dump(out, open('kernels_q3.json', 'w'), indent=1)
