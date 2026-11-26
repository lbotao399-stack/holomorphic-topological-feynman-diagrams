"""C8 detail: on O = bbb, Q_Gamma Q_0 with only one of the two components of hQ0 b
(f c b, or f g1 g2) kept: each part is nonzero, the sum vanishes."""
import random
from fractions import Fraction as Fr
import numpy as np
import vops
from vops import Op, pair, so_n
dim, f = so_n(5)
rng = random.Random(3)
O = [(('O', k), 'b') for k in range(3)]
cols = {('O', k): [rng.randint(-2, 2) for _ in range(dim)] for k in range(3)}
F = (['c'] * 5, np.array([[rng.randint(-2, 2) for _ in range(dim)] for _ in range(5)], dtype=float),
     [[rng.randint(-4, 4) for _ in range(3)] for _ in range(5)])
full = dict(vops.Q0_RULE)
for g in ['G1', 'G2', 'G3', 'G4', 'G5']:
    vals = []
    for part in ([full['b'][0]], [full['b'][1]], full['b']):
        vops.Q0_RULE['b'] = part
        vals.append(pair([Op('q0', 'q'), Op(g, 'd')], O, cols, {}, F, dim, f, {}))
    vops.Q0_RULE['b'] = full['b']
    print('%s Q0 (bbb): part f c b -> %s, part f g1 g2 -> %s, both -> %s' % (g, vals[0], vals[1], vals[2]))
for X in ['A', 'B']:
    for order in (['box', X], [X, 'box']):
        ops = [Op(order[0], 'u'), Op(order[1], 'd')]
        hist = sum(1 for _ in vops.histories(ops, O))
        print('%s acting after %s on bbb: %d histories (all types), value %s' % (order[1], order[0], hist,
              pair(ops, O, cols, {}, F, dim, f, {})))
