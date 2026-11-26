"""Three-chain move for theta(3,3,3) on C^4 (no two path lengths sum to 5): local density as the 8x8 determinant
    L = [ y_red^a delta_{i i'} | (w_i - delta_{ij}) Q^a_{j,l} ],  rows (i, a), i = 2, 3, a = 1..4;
    columns i' = 2, 3 and (j, l), j = 1, 2, 3, l = 2, 3   (CONVENTIONS Sec. 5, local formula),
    Q_{j,2} = lam_{j,1}, Q_{j,3} = lam_{j,1} + lam_{j,2}, y_red = - sum_{j,k} (1 - sigma_{j,k}) lam_{j,k},
    sigma_{j,1} = x_{j,2} + x_{j,3}, sigma_{j,2} = x_{j,3}, w_j = t_red / R_j  (t_red = 1).
Compared with the own brute-force density D(t)|det[t; dt/dxi]| at random rational points (sign fixed once)."""
import sys, random
from fractions import Fraction as Fr
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/msym-D8/code')
from graphs_c4 import theta_O
import brute_c4 as bf
import sympy as sp

G, spokes, outputs, paths = theta_O(3, 3, 3)
cot = bf.cotree_of(G.nv, G.edges)
rng = random.Random(9)
names = ['w1_2', 'w1_3', 'x1_1_2', 'x1_1_3', 'x1_2_2', 'x1_2_3', 'x1_3_2', 'x1_3_3']
syms = [sp.Symbol(n) for n in names]


def section(xi):
    w = [1 - xi['w1_2'] - xi['w1_3'], xi['w1_2'], xi['w1_3']]
    t = [None] * 9
    for i, (P, _) in enumerate(paths):
        xs = [1 - xi['x1_%d_2' % (i + 1)] - xi['x1_%d_3' % (i + 1)], xi['x1_%d_2' % (i + 1)], xi['x1_%d_3' % (i + 1)]]
        for j, e in enumerate(P):
            t[e] = xs[j] / w[i]
    return t


texpr = section(dict(zip(names, syms)))
Jm = sp.Matrix([texpr] + [[sp.diff(x, s) for x in texpr] for s in syms])
ratios = []
for it in range(6):
    while True:
        a = [Fr(rng.randint(1, 60)) for _ in range(3)]; sa = sum(a)
        w = [x / sa for x in a]
        xs = []
        for _ in range(3):
            b = [Fr(rng.randint(1, 60)) for _ in range(3)]; sb = sum(b)
            xs.append([x / sb for x in b])
        break
    xi = {'w1_2': w[1], 'w1_3': w[2]}
    for j in range(3):
        xi['x1_%d_2' % (j + 1)] = xs[j][1]; xi['x1_%d_3' % (j + 1)] = xs[j][2]
    lam = {v: [rng.randint(-9, 9) for _ in range(4)] for v in outputs}
    t = section(xi)
    J = Jm.subs({s: sp.Rational(xi[n].numerator, xi[n].denominator) for s, n in zip(syms, names)}).det()
    Dt, y = bf.density_D(G.nv, G.edges, t, lam, cot)
    g = Dt * abs(Fr(int(sp.Rational(J).p), int(sp.Rational(J).q)))
    ints = [p_[1] for p_ in paths]
    sig = [[xs[j][1] + xs[j][2], xs[j][2]] for j in range(3)]
    yred = [-sum((1 - sig[j][k]) * lam[ints[j][k]][a] for j in range(3) for k in range(2)) for a in range(4)]
    Q = {}
    for j in range(3):
        Q[(j, 2)] = lam[ints[j][0]]
        Q[(j, 3)] = [lam[ints[j][0]][a] + lam[ints[j][1]][a] for a in range(4)]
    cols = [('y', 2), ('y', 3)] + [(j, l) for j in range(3) for l in (2, 3)]
    rows = []
    for i in (1, 2):                 # chains 2, 3 (0-based 1, 2)
        for a in range(4):
            row = []
            for c in cols:
                if c[0] == 'y':
                    row.append(yred[a] if c[1] - 1 == i else Fr(0))
                else:
                    j, l = c
                    row.append((w[i] - (1 if i == j else 0)) * Q[(j, l)][a])
            rows.append(row)
    dL = bf.fdet(rows)
    ratios.append(g / dL if dL else None)
    print('point %d: brute g = %s, det L = %s, ratio %s' % (it, g, dL, ratios[-1]))
print('g = (%s) det L at all points:' % ratios[0], all(r == ratios[0] for r in ratios))
