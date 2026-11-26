"""Sym parts of the colour weights of all four-leg cubic gauge graphs with <= 2 loops (tree, box, both thetas):
formal gl(t) (single x) and so8 (fit in tr x^4, (tr x^2)^2, Pf x) -- the colour structures available to a
gauge sub-blob with four external legs inside a three-loop graph containing closed-string (Kodaira-Spencer) lines."""
import random
import numpy as np
import sympy as sp
from fractions import Fraction as Fr
import colour as C
from assemble_classical import pf
rng = random.Random(61)
graphs = [('tree', 6, None, None),
          ('box', 4, [(0, 1), (1, 2), (2, 3), (3, 0)], [0, 1, 2, 3]),
          ('theta133', 6, [(0, 5), (0, 1), (1, 2), (2, 5), (0, 3), (3, 4), (4, 5)], [1, 2, 3, 4]),
          ('theta223', 6, [(0, 1), (1, 5), (0, 2), (2, 5), (0, 3), (3, 4), (4, 5)], [1, 2, 3, 4])]
for name, nv, edges, ext in graphs[1:]:
    print('%-9s formal gl(t) W(x,x,x,x):' % name, C.ribbon_single(C.ribbon(nv, edges, ext)))
alg = C.so(8)
for name, nv, edges, ext in graphs[1:]:
    rows, vals = [], []
    for _ in range(5):
        x = [rng.randint(-2, 2) for _ in range(alg.dim)]
        M = sum((np.array(alg.B[e], dtype=object) * x[e] for e in range(alg.dim)), np.zeros_like(alg.B[0]))
        M2 = M.dot(M)
        rows.append([int(np.trace(M2.dot(M2))), int(np.trace(M2)) ** 2, int(pf(M.tolist()))])
        vals.append(Fr(C.einsum_W(alg, nv, edges, ext, [x] * 4), alg.kinv_scale ** len(edges)))
    A = sp.Matrix(rows)
    sol = A[:3, :].LUsolve(sp.Matrix(vals[:3]))
    ok = all(sum(A[i, j] * sol[j] for j in range(3)) == vals[i] for i in range(5))
    print('so8 %-9s Sym W in (tr x^4, (tr x^2)^2, Pf): %s  fit consistent: %s' % (name, list(sol), ok), flush=True)
