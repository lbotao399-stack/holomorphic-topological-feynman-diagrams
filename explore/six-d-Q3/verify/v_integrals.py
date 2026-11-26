"""C4 (convexity), C5: zero-displacement values int_{D_Gamma} J_Gamma by sympy iterated integration
(no monomial formula), bracket expansion compared with the explorer's json (key intJ), values at the
explorer's sources, factorization test in independent bracket symbols."""
import itertools
import json
import os
import random
import re
import sys
import time
from fractions import Fraction as Fr

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
R = sp.symbols('r1:10')
r1, r2, r3, r4, r5, r6, r7, r8, r9 = R

MOVES = {'G1': ['chain', 'chain', 'box'], 'G2': ['chain', 'bubble', 'box'], 'G3': ['chain', 'bubble', 'box'],
         'G4': ['bubble', 'bubble', 'box'], 'G5': ['bubble', 'bubble', 'box']}
OUTS = {'G1': [1, 2, 4, 6, 7], 'G2': [1, 2, 4, 6, 7], 'G3': [1, 2, 5, 6, 7], 'G4': [1, 2, 4, 6, 7], 'G5': [1, 2, 5, 6, 7]}
ISIGN = {'G1': 1, 'G2': 1, 'G3': 1, 'G4': -1, 'G5': -1}


def V(**kw):
    return {int(k[1:]): sp.sympify(v) for k, v in kw.items()}


def lin(*terms):
    out = {}
    for c, vec in terms:
        for k, v in vec.items():
            out[k] = sp.expand(out.get(k, 0) + c * v)
    return {k: v for k, v in out.items() if v != 0}


def L(i):
    return {i: sp.Integer(1)}


def det_exp(a, b, c):
    """multilinear expansion: {sorted triple: coefficient}."""
    out = {}
    for i, ca in a.items():
        for j, cb in b.items():
            for k, cc in c.items():
                if len({i, j, k}) < 3:
                    continue
                key = tuple(sorted((i, j, k)))
                perm = [key.index(x) for x in (i, j, k)]
                s = 1
                for x in range(3):
                    for y in range(x + 1, 3):
                        if perm[x] > perm[y]:
                            s = -s
                out[key] = sp.expand(out.get(key, 0) + s * ca * cb * cc)
    return {k: v for k, v in out.items() if v != 0}


def J_expansion(name):
    mu = lin((1 - r2, L(1)), (1 - r3, L(2)))
    if name == 'G1':
        pref = r1 ** 2 * r4 ** 2
        dets = [det_exp(L(1), L(2), lin((r6 + (1 - r4) * (1 - r6) * r7, L(4)), ((1 - r4) * r8, L(6)), ((1 - r4) * r9, L(7)))),
                det_exp(mu, L(4), lin((r8, L(6)), (r9, L(7)))),
                det_exp(lin((1 - r5, mu), (1 - r6, L(4))), L(6), L(7))]
    elif name in ('G2', 'G4'):
        pref = r1 ** 2 * r4 * (1 - r4) if name == 'G2' else r1 * (1 - r1) * r4 * (1 - r4)
        dets = [det_exp(L(1), L(2), lin((r8, L(4)), (r9 * r5, L(6)), (r9 * r6, L(7)))),
                det_exp(lin((1 - r7, mu), (1 - r8, L(4))), L(6), L(7)),
                det_exp(mu, L(4), lin((r5, L(6)), (r6, L(7))))]
    else:
        pref = r1 ** 2 * r4 * (1 - r4) ** 2 if name == 'G3' else r1 * (1 - r1) * r4 * (1 - r4) ** 2
        dets = [det_exp(L(1), L(2), lin((r7 * (1 - r6), L(5)), (r8, L(6)), (r9, L(7)))),
                det_exp(mu, L(5), lin((r8, L(6)), (r9, L(7)))),
                det_exp(lin((1 - r5, mu), (1 - r6, L(5))), L(6), L(7))]
    out = {}
    for (k1, c1), (k2, c2), (k3, c3) in itertools.product(*[d.items() for d in dets]):
        out[(k1, k2, k3)] = sp.expand(pref * c1 * c2 * c3)
    return out


def integrate(expr, name):
    """Iterated sympy integration over D_Gamma."""
    e = expr
    for m, kind in enumerate(MOVES[name]):
        a, b, c = R[3 * m:3 * m + 3]
        if kind == 'chain':        # 0<a<1, 1>b>c>0
            e = sp.integrate(e, (c, 0, b))
            e = sp.integrate(e, (b, 0, 1))
            e = sp.integrate(e, (a, 0, 1))
        elif kind == 'bubble':
            e = sp.integrate(e, (a, 0, 1), (b, 0, 1), (c, 0, 1))
        else:                      # 1>a>b>c>0
            e = sp.integrate(e, (c, 0, b))
            e = sp.integrate(e, (b, 0, a))
            e = sp.integrate(e, (a, 0, 1))
    return sp.nsimplify(sp.expand(e))


def bsym(t):
    return sp.Symbol('b%d%d%d' % t)


def parse_intJ(s):
    m = re.match(r'\(1/(\d+)\) \* \((.*)\)$', s.strip())
    den = int(m.group(1))
    body = m.group(2)
    terms = re.findall(r'([+-]\d+) <(\d)(\d)(\d)><(\d)(\d)(\d)><(\d)(\d)(\d)>', body)
    poly = 0
    for t in terms:
        c = int(t[0])
        b = [tuple(int(x) for x in t[1:4]), tuple(int(x) for x in t[4:7]), tuple(int(x) for x in t[7:10])]
        poly += c * bsym(b[0]) * bsym(b[1]) * bsym(b[2])
    return sp.Rational(1, den) * poly, len(terms), den


def det3(a, b, c):
    return (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]))


EXPL_LAM = {  # numerics_simplex.log: lam at the outputs (in output order)
    'G1': [[0, 2, -1], [-1, 3, 3], [3, 4, 4], [3, 2, 0], [-4, 3, -4]],
    'G2': [[3, 2, -1], [1, 4, 4], [3, 1, -4], [-4, -2, 4], [2, 0, 3]],
    'G3': [[-1, -1, 4], [-1, -4, -3], [-4, -4, 2], [3, 3, -1], [-4, 3, 3]],
    'G4': [[3, -4, 1], [0, -1, -3], [-1, -3, -1], [0, -1, -3], [1, -3, -2]],
    'G5': [[-1, -2, -2], [1, 3, -1], [-4, -1, 0], [1, -1, -4], [0, 4, 2]],
}
EXPL_VAL = {'G1': Fr(52924313, 311040), 'G2': Fr(-85183, 25920), 'G3': Fr(58733, 25920),
            'G4': Fr(7757, 1440), 'G5': Fr(-18409, 51840)}


def evaluate(poly, lam):
    subs = {}
    for s in poly.free_symbols:
        i, j, k = (int(x) for x in s.name[1:])
        subs[s] = det3(lam[i], lam[j], lam[k])
    return Fr(str(sp.nsimplify(poly.subs(subs))))


if __name__ == '__main__':
    KJ = json.load(open(os.path.join(HERE, '..', 'q3_results.json')))
    res = {}
    for name in ['G1', 'G2', 'G3', 'G4', 'G5']:
        t0 = time.time()
        ex = J_expansion(name)
        mine = 0
        coeffs = {}
        for key, P in ex.items():
            v = integrate(P, name)
            coeffs[key] = v
            mine += v * bsym(key[0]) * bsym(key[1]) * bsym(key[2])
        mine = sp.expand(mine)
        theirs, nterms, den = parse_intJ(KJ[name]['intJ'])
        theirs = sp.expand(theirs)
        nz = sum(1 for v in coeffs.values() if v != 0)
        lcm = sp.ilcm(*[sp.Rational(v).q for v in coeffs.values() if v != 0])
        lam = {v: [Fr(x) for x in vec] for v, vec in zip(OUTS[name], EXPL_LAM[name])}
        val = evaluate(mine, lam)
        fac = sp.factor_list(mine)
        nonconst = [(f, k) for f, k in fac[1] if f.free_symbols]
        print('%s: %d positional bracket triples (%d nonzero), common denominator %s; explorer: %d terms, denominator %d; '
              'equal as polynomials in independent bracket symbols: %s' % (name, len(ex), nz, lcm, nterms, den,
                                                                          sp.expand(mine - theirs) == 0))
        print('    at the explorer sources: int J = %s, I = %+d * int J = %s (explorer I(lambda;0) = %s) equal %s'
              % (val, ISIGN[name], ISIGN[name] * val, EXPL_VAL[name], ISIGN[name] * val == EXPL_VAL[name]))
        print('    factor_list in independent bracket symbols: content %s, %d non-constant factor(s), degrees %s  [%.1fs]'
              % (fac[0], len(nonconst), [sp.Poly(f).total_degree() for f, k in nonconst], time.time() - t0), flush=True)
        res[name] = {'%d%d%d.%d%d%d.%d%d%d' % (k[0] + k[1] + k[2]): str(v) for k, v in coeffs.items()}
    json.dump(res, open(os.path.join(HERE, 'v_integrals.json'), 'w'), indent=1)
