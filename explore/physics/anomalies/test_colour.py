"""ribbon (formal gl(t)) vs explicit einsum contraction for gl2, gl3 at random integer matrices; box, both thetas,
all nine three-loop bulk graphs."""
import json, random
import numpy as np
from fractions import Fraction as Fr
import colour as C
rng = random.Random(2)
G3 = json.load(open('kernels_n3.json'))
graphs = [('box', 4, [(0, 1), (1, 2), (2, 3), (3, 0)], [0, 1, 2, 3]),
          ('theta133', 6, [(0, 5), (0, 1), (1, 2), (2, 5), (0, 3), (3, 4), (4, 5)], [1, 2, 3, 4]),
          ('theta223', 6, [(0, 1), (1, 5), (0, 2), (2, 5), (0, 3), (3, 4), (4, 5)], [1, 2, 3, 4])]
for r in G3:
    graphs.append(('n3#%d' % r['idx'], r['nv'], [tuple(e) for e in r['edges']], r['legs'] + [r['root']]))
for N in (2, 3):
    alg = C.gl(N)
    for name, nv, edges, ext in graphs:
        rib = C.ribbon(nv, edges, ext)
        oks = []
        for _ in range(2):
            mats = [np.array([[rng.randint(-3, 3) for _ in range(N)] for _ in range(N)], dtype=object) for _ in range(4)]
            xs = [[int(m[a, b]) for a in range(N) for b in range(N)] for m in mats]  # coords in basis E_ab: x^{ab} with tr(E_ab E_cd)
            # coordinates: m = sum m[a,b] E_ab
            v1 = C.ribbon_eval(rib, mats, N)
            v2 = Fr(C.einsum_W(alg, nv, edges, ext, xs), alg.kinv_scale ** len(edges))
            oks.append(v1 == v2)
        print('gl%d %-9s ribbon == einsum: %s  (%d ribbon terms)' % (N, name, oks, len(rib)), flush=True)
