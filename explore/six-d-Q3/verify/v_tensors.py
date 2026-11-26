"""C6: Wick tensors T_Gamma regenerated from the labelled graphs with the leg pattern of S2 (1.4),
(1.6), (1.7) (in-edges in edge order, output, out-edges; spoke j -> A_j in edge order; tree edge (a,b) ->
R_ab; outputs in vertex order), compared with the explorer's strings; S2's T_A, T_B and the box
reproduced; arithmetic of the Q_3 coefficients; counting of letters and derivatives."""
import json
import os
import re
from fractions import Fraction as Fr

from vops import GRAPH_SPECS

HERE = os.path.dirname(os.path.abspath(__file__))
KJ = json.load(open(os.path.join(HERE, '..', 'q3_results.json')))


def tensor_string(name):
    g = GRAPH_SPECS[name]
    E = g['edges']
    spokes = g['spokes']
    parts = []
    for v in g['verts']:
        legs = []
        for leg in g['legs'][v]:
            if leg[0] == 'out':
                legs.append('D%d' % leg[1])
            elif leg[1] in spokes:
                legs.append('A%d' % (spokes.index(leg[1]) + 1))
            else:
                a, b = E[leg[1]]
                legs.append('R%s%s' % (a, b))
        parts.append('C_{%s}' % ' '.join(legs))
    outs = ' '.join('C^{D%d}(x+X_%d)' % (v, v) for v in g['outputs'])
    return ' '.join(parts) + ' * ' + outs


for name in ['G1', 'G2', 'G3', 'G4', 'G5']:
    mine = tensor_string(name)
    theirs = KJ[name]['T']
    print('%s: regenerated %s\n     equal to explorer string: %s' % (name, mine, mine == theirs))
for name in ['box', 'A', 'B']:
    print('%s: %s' % (name, tensor_string(name)))

# coefficient arithmetic
aut = {'G1': 2, 'G2': 2, 'G3': 1, 'G4': 8, 'G5': 2}
kappa = {'G1': Fr(-1, 2), 'G2': Fr(1, 2), 'G3': Fr(-1), 'G4': Fr(1, 8), 'G5': Fr(-1, 2)}
isign = {'G1': 1, 'G2': 1, 'G3': 1, 'G4': -1, 'G5': -1}
stated = {'G1': Fr(-1), 'G2': Fr(1), 'G3': Fr(-2), 'G4': Fr(-1, 4), 'G5': Fr(1)}
for g in aut:
    c = 2 * kappa[g] * isign[g]   # (i/512 pi^9) kappa I = (i/1024 pi^9) (2 kappa sgn_I) int J
    print('%s: |kappa| = %s, 1/|Aut| = %s; coefficient of int_D J T in units i/(1024 pi^9): %s (stated %s) %s'
          % (g, abs(kappa[g]), Fr(1, aut[g]), c, stated[g], c == stated[g]))
print('(i/(8 pi^3))^3 = i/(512 pi^9): 8^3 =', 8 ** 3)
# letters and derivatives: each bracket eps^{abc} d_a d_b d_c adds one d_1, d_2, d_3; three brackets
for g in aut:
    spec = GRAPH_SPECS[g]
    print('%s: inputs %d, outputs %d, brackets per term 3 -> derivatives (3,3,3)' % (g, len(spec['spokes']), len(spec['outputs'])))
