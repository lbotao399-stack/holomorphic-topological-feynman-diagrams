"""Numerical cross-check (engine oracle: direct lattice integration of the brute-force density,
no reduction) of the four-leg kernels of hCS on C^3 at momentum conservation:
 box: K = -Delta/6;  theta(1,3,3): K = -Delta^2/432;  theta(2,2,3): K = 0;  Delta = det(p1,p2,p3)."""
import sys, time, random
from fractions import Fraction as Fr
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
import htk
from htk import theta, polygon

rng = random.Random(11)


def det3(a, b, c):
    return (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]))


p = [[rng.randint(-4, 4) for _ in range(3)] for _ in range(3)]
p.append([-(p[0][i] + p[1][i] + p[2][i]) for i in range(3)])
D = det3(p[0], p[1], p[2])
print("momenta", p, "Delta =", D)
box = polygon(4)
lam = [p[0], p[1], p[2]]
t0 = time.time()
ex = htk.kernel_exact(box, lam, T=0)
val, err = htk.kernel_numeric(box, lam, T=0)
print("box: exact %s, -Delta/6 = %s, oracle %.10f +- %.1e (%.0fs)" % (ex, Fr(-D, 6), val, err, time.time() - t0), flush=True)
for pqr, pred in [((1, 3, 3), Fr(-D * D, 432)), ((2, 2, 3), Fr(0))]:
    G, interior = theta(*pqr)
    ext = [v for P in interior for v in P]
    lamd = [[0, 0, 0] for _ in range(G.nv - 1)]
    for i, v in enumerate(ext):
        lamd[v] = p[i]
    t0 = time.time()
    ex = htk.kernel_exact(G, lamd, T=0)
    val, err = htk.kernel_numeric(G, lamd, T=0)
    print("theta%s: exact %s, prediction %s, oracle %.10f +- %.1e, |oracle - exact| = %.1e (%.0fs)" % (
        pqr, ex, pred, val, err, abs(val - float(ex)), time.time() - t0), flush=True)
