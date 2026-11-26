"""(d) Kernel table: K_G^(v)(lambda;0)/pi on sum p = 0 as polynomials in the brackets [ij] = [p_i p_j]
(1-based external labels i = 1..4 = ext[0..3] = vertices 1,2,3,4; p4 = -p1-p2-p3; [12],[13],[23] algebraically
independent).  Insertion [q_h1 q_h2], (h1,h2) = first two half-edges of L'_v."""
import itertools
from common import *

KD = load('kernels.pkl')
p = PSYM
sub3 = {p[3][c]: -(p[0][c] + p[1][c] + p[2][c]) for c in range(2)}
B = lambda i, j: p[i][0] * p[j][1] - p[i][1] * p[j][0]
b12, b13, b23 = sp.symbols('b12 b13 b23')
gens = [b12, b13, b23]
vals = [B(0, 1), B(0, 2), B(1, 2)]


def in_brackets(K):
    K = sp.expand(K.subs(sub3, simultaneous=True))
    w = sp.Poly(K, *[p[i][c] for i in range(3) for c in range(2)]).total_degree() // 2
    monos = [m for m in itertools.product(range(w + 1), repeat=3) if sum(m) == w]
    cs = sp.symbols('c0:%d' % len(monos))
    ans = sum(cc * vals[0] ** m[0] * vals[1] ** m[1] * vals[2] ** m[2] for cc, m in zip(cs, monos))
    eqs = sp.Poly(sp.expand(ans - K), *[p[i][c] for i in range(3) for c in range(2)]).coeffs()
    sol = sp.solve(eqs, cs, dict=True)
    assert len(sol) == 1
    expr = sum(sol[0].get(cc, cc) * gens[0] ** m[0] * gens[1] ** m[1] * gens[2] ** m[2] for cc, m in zip(cs, monos))
    return sp.factor(expr)


print("brackets: b12 = [p1p2], b13 = [p1p3], b23 = [p2p3];  [p3p4] = [p1p3] + [p2p3] = b13 + b23 on sum p = 0")
for pqr, dat in KD.items():
    print("theta%s edges %s" % (pqr, dat['edges']))
    print("  K/pi (no insertion) =", in_brackets(dat['K0']))
    for v in range(6):
        print("  v=%d  L'_v = %s :  K^(v)/pi = %s" % (v, dat['local'][v], in_brackets(dat['Kv'][v])))
