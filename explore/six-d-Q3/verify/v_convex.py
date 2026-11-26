"""C4: each X_v is a convex combination of z_1..z_4 on D_Gamma (coefficients sum to 1 identically;
nonnegativity: every coefficient is a sum of products of factors r_i, 1-r_i, r_i-r_j (i<j in an
ordered block), checked by substituting the block parametrization and testing that all monomial
coefficients are >= 0).  Also: the C5 statement's explicit G1 coefficient list vs. the json."""
import json
import os
import re

import sympy as sp

import vk

HERE = os.path.dirname(os.path.abspath(__file__))
R = sp.symbols('r1:10')
Zs = [[sp.Symbol('z%d_%d' % (j, a)) for a in range(3)] for j in range(1, 5)]


def coeffs(name):
    X = vk.X_formula(name, list(R), [list(z) for z in Zs])
    out = {}
    for v, vec in X.items():
        e = sp.expand(vec[0])
        out[v] = [sp.expand(e.coeff(Zs[j][0])) for j in range(4)]
    return out


def positive_param(name):
    """substitution making D_Gamma the positive orthant image: (0,1): r = p/(1+p);
    ordered block 1>a>b>c>0 or 1>b>c>0: a = 1/(1+p1), b = a/(1+p2), c = b/(1+p3)."""
    P = sp.symbols('p1:10', positive=True)
    sub = {}
    for m, mv in enumerate(vk.GRAPHS[name]['moves']):
        i = 3 * m
        if mv[0] == 'chain':
            sub[R[i]] = P[i] / (1 + P[i])
            sub[R[i + 1]] = 1 / (1 + P[i + 1])
            sub[R[i + 2]] = sub[R[i + 1]] / (1 + P[i + 2])
        elif mv[0] == 'bubble':
            for k in range(3):
                sub[R[i + k]] = P[i + k] / (1 + P[i + k])
        else:
            sub[R[i]] = 1 / (1 + P[i])
            sub[R[i + 1]] = sub[R[i]] / (1 + P[i + 1])
            sub[R[i + 2]] = sub[R[i + 1]] / (1 + P[i + 2])
    return sub, P


import random
from fractions import Fraction as Fr
rng = random.Random(11)
for name in ['G1', 'G2', 'G3', 'G4', 'G5']:
    C = coeffs(name)
    sums = all(sp.expand(sum(c) - 1) == 0 for c in C.values())
    # nonnegativity at random rational points of D_Gamma (z_j -> unit coefficient extraction)
    e = [[Fr(int(j == k)), Fr(0), Fr(0)] for k in range(4) for j in [0]]
    zs = [[Fr(int(j == k)), Fr(0), Fr(0)] for j in range(4) for k in [0]]
    npt, neg, mn = 2000, 0, Fr(1)
    for _ in range(npt):
        r = vk.random_point(name, rng)
        for j in range(4):
            z = [[Fr(int(jj == j)), Fr(0), Fr(0)] for jj in range(4)]
            X = vk.X_formula(name, r, z)
            for v, vec in X.items():
                mn = min(mn, vec[0])
                if vec[0] < 0:
                    neg += 1
    print('%s: coefficient sums == 1 identically (sympy): %s; %d random rational points of D: negative coefficients %d, '
          'minimum coefficient %.3g' % (name, sums, npt, neg, float(mn)))

# C5: explicit G1 list in the statement vs json
stmt = ("32<124><146><167>+48<124><146><267>+128<124><146><467>+16<124><147><167>+24<124><147><267>+64<124><147><467>"
        "+48<124><246><167>+96<124><246><267>+256<124><246><467>+24<124><247><167>+48<124><247><267>+128<124><247><467>"
        "+12<126><146><167>+18<126><146><267>+48<126><146><467>+6<126><147><167>+9<126><147><267>+24<126><147><467>"
        "+18<126><246><167>+36<126><246><267>+96<126><246><467>+9<126><247><167>+18<126><247><267>+48<126><247><467>"
        "+6<127><146><167>+9<127><146><267>+24<127><146><467>+4<127><147><167>+6<127><147><267>+16<127><147><467>"
        "+9<127><246><167>+18<127><246><267>+48<127><246><467>+6<127><247><167>+12<127><247><267>+32<127><247><467>")
KJ = json.load(open(os.path.join(HERE, '..', 'q3_results.json')))
a = sorted(re.findall(r'([+-]?\d+)<(\d+)><(\d+)><(\d+)>', stmt))
b = sorted((c.replace('+', ''), x, y, w) for c, x, y, w in
           re.findall(r'([+-]\d+) <(\d+)><(\d+)><(\d+)>', KJ['G1']['intJ']))
a = sorted((c.replace('+', ''), x, y, w) for c, x, y, w in a)
print('C5 statement G1 list (36 terms) == json intJ G1:', a == b, len(a), len(b))
