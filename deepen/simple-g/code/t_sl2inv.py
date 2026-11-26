"""SL2(C) acting on (z1,z2): on momentum polynomials of n independent 2-vectors q_i (after conservation) of
torus weight (2,2) (degree 2 in the q_i^1, degree 2 in the q_i^2) the invariants = kernel of the raising operator
E = sum_i q_i^1 d/dq_i^2 (weight-0 vectors killed by E span the trivial isotypic part).  Exact (sympy):
dim ker E for n = 2 (3-cochains), 3 (4-cochains), 4 (5-cochains), and comparison with [q0q1]^2 / the 6 quadratics QB /
the 20 bracket monomials MON20."""
import sys, itertools
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/simple-g/code')
import sympy as sp
import closedness as cl
for n in (2, 3, 4):
    q = [[sp.Symbol('q%d_%d' % (i, c)) for c in range(2)] for i in range(n)]
    g1 = [q[i][0] for i in range(n)]; g2 = [q[i][1] for i in range(n)]
    m1 = [e for e in itertools.product(range(3), repeat=n) if sum(e) == 2]
    mons = [sp.Mul(*[g1[i] ** a[i] * g2[i] ** b[i] for i in range(n)]) for a in m1 for b in m1]
    allg = g1 + g2
    imgs = [sp.Poly(sp.expand(sum(q[i][0] * sp.diff(m, q[i][1]) for i in range(n))), *allg).as_dict() if True else None for m in mons]
    keys = sorted(set().union(*imgs))
    E = sp.Matrix([[im.get(k, 0) for im in imgs] for k in keys])
    ker = E.nullspace()
    B = lambda a, b: q[a][0] * q[b][1] - q[a][1] * q[b][0]
    if n == 2:
        cand = [B(0, 1) ** 2]
    elif n == 3:
        u, v, w = B(0, 1), B(0, 2), B(1, 2)
        cand = [u*u, v*v, w*w, u*v, u*w, v*w]
    else:
        bs = [B(a, b) for (a, b) in cl.PAIRS]
        cand = [bs[i] * bs[j] for (i, j) in cl.MON20]
    pc = [sp.Poly(sp.expand(c), *allg).as_dict() for c in cand]
    mk = sorted(set().union(*pc))
    rk = sp.Matrix([[p.get(k, 0) for k in mk] for p in pc]).rank()
    Ec = all(sp.expand(sum(q[i][0] * sp.diff(c, q[i][1]) for i in range(n))) == 0 for c in cand)
    print("n = %d: weight-(2,2) monomials %d, dim ker E = %d ; candidates: %d, rank %d, killed by E: %s"
          % (n, len(mons), len(ker), len(cand), rk, Ec))
