"""Symbolic check of the compact forms of the report: I_Gamma = eps_Gamma int_{D_Gamma} J_Gamma with J_Gamma
written as r-prefactor * det * det * det, nested displacement maps, automorphism groups (text labels)."""
import itertools
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher
import sympy as sp
from q3lib import GRAPHS, Reduction, LAM, Z, R, ROOT, tlabel

r1, r2, r3, r4, r5, r6, r7, r8, r9 = R
l = {k: sp.Matrix(LAM[k - 1]) for k in range(1, 8)}
z1, z2, z3, z4 = [sp.Matrix(Z[j]) for j in range(4)]
det = lambda a, b, c: sp.Matrix.hstack(a, b, c).det()
SIGN = {'G1': 1, 'G2': -1, 'G3': 1, 'G4': 1, 'G5': -1}

# compact densities (sources at the outputs only) and eps_Gamma
mu = (1 - r2) * l[1] + (1 - r3) * l[2]
FORMS = {
    'G1': (1, r1 ** 2 * r4 ** 2
           * det(l[1], l[2], (r6 + (1 - r4) * (1 - r6) * r7) * l[4] + (1 - r4) * (r8 * l[6] + r9 * l[7]))
           * det(mu, l[4], r8 * l[6] + r9 * l[7])
           * det((1 - r5) * mu + (1 - r6) * l[4], l[6], l[7])),
    'G2': (1, r1 ** 2 * r4 * (1 - r4)
           * det(l[1], l[2], r8 * l[4] + r9 * (r5 * l[6] + r6 * l[7]))
           * det((1 - r7) * mu + (1 - r8) * l[4], l[6], l[7])
           * det(mu, l[4], r5 * l[6] + r6 * l[7])),
    'G3': (1, r1 ** 2 * r4 * (1 - r4) ** 2
           * det(l[1], l[2], r7 * (1 - r6) * l[5] + r8 * l[6] + r9 * l[7])
           * det(mu, l[5], r8 * l[6] + r9 * l[7])
           * det((1 - r5) * mu + (1 - r6) * l[5], l[6], l[7])),
    'G4': (-1, r1 * (1 - r1) * r4 * (1 - r4)
           * det(l[1], l[2], r8 * l[4] + r9 * (r5 * l[6] + r6 * l[7]))
           * det((1 - r7) * mu + (1 - r8) * l[4], l[6], l[7])
           * det(mu, l[4], r5 * l[6] + r6 * l[7])),
    'G5': (-1, r1 * (1 - r1) * r4 * (1 - r4) ** 2
           * det(l[1], l[2], r7 * (1 - r6) * l[5] + r8 * l[6] + r9 * l[7])
           * det(mu, l[5], r8 * l[6] + r9 * l[7])
           * det((1 - r5) * mu + (1 - r6) * l[5], l[6], l[7])),
}

# nested displacement maps
W1 = (1 - r1) * z1 + r1 * z2


def maps(name):
    if name == 'G1':
        W2 = (1 - r4) * W1 + r4 * z3
        P5 = z4 + r7 * (W2 - z4)
        P3 = r5 * W1 + (1 - r5) * P5
        return {1: r2 * z1 + (1 - r2) * P3, 2: r3 * z1 + (1 - r3) * P3, 4: r6 * W1 + (1 - r6) * P5,
                6: z4 + r8 * (W2 - z4), 7: z4 + r9 * (W2 - z4)}
    if name in ('G2', 'G4'):
        W2 = (1 - r4) * z3 + r4 * z4
        P3 = r7 * W1 + (1 - r7) * W2
        P5 = r9 * W1 + (1 - r9) * W2
        X2 = r3 * z1 + (1 - r3) * P3 if name == 'G2' else r3 * z2 + (1 - r3) * P3
        return {1: r2 * z1 + (1 - r2) * P3, 2: X2, 4: r8 * W1 + (1 - r8) * W2,
                6: r5 * P5 + (1 - r5) * z3, 7: r6 * P5 + (1 - r6) * z4}
    W2 = (1 - r4) * W1 + r4 * z3
    P4 = z4 + r7 * (W2 - z4)
    P3 = r5 * W1 + (1 - r5) * P4
    X2 = r3 * z1 + (1 - r3) * P3 if name == 'G3' else r3 * z2 + (1 - r3) * P3
    return {1: r2 * z1 + (1 - r2) * P3, 2: X2, 5: r6 * z3 + (1 - r6) * P4,
            6: z4 + r8 * (W2 - z4), 7: z4 + r9 * (W2 - z4)}


for name in GRAPHS:
    red = Reduction(name)
    outs = [v + 1 for v in range(ROOT) if sum((a == v) + (b == v) for a, b in red.edges) == 2]
    zero = {LAM[v][a]: 0 for v in range(ROOT) if v + 1 not in outs for a in range(3)}
    dens = sp.expand(SIGN[name] * sp.prod(red.factors()).subs(zero))
    eps, Jc = FORMS[name]
    ok_J = sp.expand(dens - eps * Jc) == 0
    X, _ = red.displacement()
    M = maps(name)
    ok_X = all(sp.expand(X[v - 1][a] - M[v][a]) == 0 for v in outs for a in range(3))
    # automorphisms fixing O, in text labels
    G = nx.Graph()
    G.add_edges_from(red.edges)
    nx.set_node_attributes(G, {v: v == ROOT for v in G.nodes()}, 'O')
    auts = list(GraphMatcher(G, G, node_match=lambda a, b: a['O'] == b['O']).isomorphisms_iter())
    cyc = []
    for m in auts:
        seen, cs = set(), []
        for v in sorted(m):
            if v in seen or m[v] == v:
                continue
            c, x = [], v
            while x not in seen:
                seen.add(x); c.append(tlabel(x)); x = m[x]
            cs.append('(' + ' '.join(c) + ')')
        cyc.append(''.join(cs) or 'id')
    print('%s: I = %+d int J (compact form) %s; nested X_v %s; |Aut| = %d: %s'
          % (name, eps, ok_J, ok_X, len(auts), ', '.join(cyc)))
