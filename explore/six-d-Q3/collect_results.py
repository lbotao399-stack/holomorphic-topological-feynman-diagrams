"""Collects the formulas of Q_3 into q3_results.json (input for verification): labelled graphs,
Wick tensors T_Gamma (leg pattern), kernel tables int_{D_Gamma} J_Gamma, coefficients of Q_3."""
import json
from fractions import Fraction as Fr

import sympy as sp

import opcheck as oc
from q3lib import tlabel

KJ = json.load(open('kernels_q3.json'))
EPS = {'G1': 1, 'G2': 1, 'G3': 1, 'G4': -1, 'G5': -1}         # I_Gamma = EPS * int J_Gamma (compact_forms.py)
KAPPA = {'G1': Fr(-1, 2), 'G2': Fr(1, 2), 'G3': Fr(-1), 'G4': Fr(1, 8), 'G5': Fr(-1, 2)}   # q3check.py
out = {}
for g in ['G1', 'G2', 'G3', 'G4', 'G5']:
    G = oc.op(g)
    E = G.red.edges

    def legname(typ, x):
        if typ == 'A':
            return 'A%d' % (x + 1)
        if typ == 'D':
            return 'D%s' % tlabel(x)
        a, b = E[x]
        return 'R%s%s' % (tlabel(a), tlabel(b))
    colour = ' '.join('C_{%s}' % ' '.join(legname(t, x) for t, x in L) for L in G.legs)
    fields = ' '.join('C^{D%s}(x+X_%s)' % (tlabel(v), tlabel(v)) for v in G.outs)
    # int J = EPS * K0
    terms = []
    den = 1
    for b1, b2, b3, c in KJ[g]['K0']:
        q = EPS[g] * sp.Rational(c)
        den = sp.ilcm(den, q.q)
        terms.append((q, (b1, b2, b3)))
    num = ['%+d <%s><%s><%s>' % (int(q * den), ''.join(tlabel(v) for v in b1), ''.join(tlabel(v) for v in b2),
                                   ''.join(tlabel(v) for v in b3)) for q, (b1, b2, b3) in sorted(terms, key=lambda x: x[1])]
    coeffQ3 = KAPPA[g] * EPS[g]           # coefficient of int J T delta^4 in Q_3 / (i / 512 pi^9)
    out[g] = dict(edges=[(tlabel(a), tlabel(b)) for a, b in E],
                  spokes=[('%s->%s' % (tlabel(E[e][0]), tlabel(E[e][1])), 'z_e = %+d z_%d' % (s, j + 1))
                          for j, (e, s) in enumerate(G.red.spokes)],
                  outputs=[tlabel(v) for v in G.outs], aut=oc.AUT[g], eps_I=EPS[g], kappa=str(KAPPA[g]),
                  T=colour + ' * ' + fields, intJ='(1/%d) * (%s)' % (den, ' '.join(num)),
                  coefficient_in_Q3_over_i_512pi9=str(coeffQ3),
                  domain=KJ[g]['domain'])
    print(g, '|Aut| =', oc.AUT[g], ' I = %+d int J' % EPS[g], ' kappa =', KAPPA[g], ' coefficient of int J T:', coeffQ3)
    print('   T_%s =' % g, out[g]['T'])
    print('   int J_%s = 1/%d (%d terms)' % (g, den, len(num)))
json.dump(out, open('q3_results.json', 'w'), indent=1)
