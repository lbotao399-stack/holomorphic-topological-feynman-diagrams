"""Direct contraction check of T_1, T_3 for sl2, sl3 (trace form): T_n(x) / (tr x^2)^2 at random x."""
import json, random
from fractions import Fraction as Fr
import numpy as np
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher
import colour as C
rng = random.Random(71)
recs = json.load(open('kernels_n3.json'))
def aut(nv, edges):
    g = nx.Graph(); g.add_edges_from(edges)
    return sum(1 for _ in GraphMatcher(g, g).isomorphisms_iter())
for alg in (C.sl(2), C.sl(3)):
    for _ in range(2):
        q = 0
        while q == 0:
            x = [rng.randint(-3, 3) for _ in range(alg.dim)]
            M = sum((np.array(alg.B[e], dtype=object) * x[e] for e in range(alg.dim)), np.zeros_like(alg.B[0]))
            q = Fr(int(np.trace(M.dot(M)))) ** 2
        t1 = Fr(-1, 6) * Fr(C.einsum_W(alg, 4, [(0, 1), (1, 2), (2, 3), (3, 0)], [0, 1, 2, 3], [x] * 4), alg.kinv_scale ** 4) / 8
        t3, t3z = Fr(0), Fr(0)
        for r in recs:
            c = {k: Fr(v) for k, v in r['bare'].items()}
            if not c:
                continue
            E = [tuple(e) for e in r['edges']]
            w = Fr(C.einsum_W(alg, r['nv'], E, r['legs'] + [r['nv'] - 1], [x] * 4), alg.kinv_scale ** len(E)) / aut(r['nv'], E)
            t3 += c.get('1', 0) * w; t3z += c.get('z3', 0) * w
        print('%s: T_1/(tr x^2)^2 = %s, T_3/(tr x^2)^2 = %s, zeta3 part %s' % (alg.name, t1 / q, t3 / q, t3z), flush=True)
