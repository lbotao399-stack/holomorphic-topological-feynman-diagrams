"""Direct numerical cross-check (no reduction formulas):
K = int_{RP^5} S_R(R) d^5 log R,   S_R(R) = int_x S(t(x,R)) / prod_e x_e dx,   t_{E,j} = R_E x_{E,j},
S = D prod t (htk.oracle.scaled_density_batch, brute force), x-integral by Gauss product rules (exact for
polynomials in x; agreement of m=3,4 checked).  Trapezoidal rule on h Z^5 (+ shift h/2), log R_6 = 0,
flood fill from the origin through points with |S_R| > tol * max."""
import sys, os, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'code'))
from htk import oracle
from k4num import setup, x_rule, SR

def run(which, h, tol, m=3, shift=0.0, maxpts=4_000_000, chunk=4000):
    name, edges = __import__('graphs').K4SUB[which]
    G, lam, lam_all, cubic, chains, Dc, others = setup(edges)
    st = oracle._Setup(G)
    xn, xw = x_rule(chains, m)
    def F(keys):
        s = np.array(keys, dtype=float) * h + shift
        R = np.ones((len(keys), 6)); R[:, :5] = np.exp(s)
        return SR(R, G, lam, chains, xn, xw, st)
    seen = {}
    front = [tuple([0] * 5)]
    seen[front[0]] = None
    total = 0.0; vmax = 0.0; below = 0.0; n = 0
    while front and n < maxpts:
        nxt = []
        for i in range(0, len(front), chunk):
            blk = front[i:i + chunk]
            vals = F(blk)
            n += len(blk)
            for k, v in zip(blk, vals):
                if not np.isfinite(v):
                    v = 0.0
                seen[k] = v
                total += v
                av = abs(v)
                vmax = max(vmax, av)
                if av > tol * vmax:
                    for d in range(5):
                        for sgn in (-1, 1):
                            kk = list(k); kk[d] += sgn; kk = tuple(kk)
                            if kk not in seen:
                                seen[kk] = None
                                nxt.append(kk)
                else:
                    below += av
        front = nxt
    return total * h ** 5, below * h ** 5, n

if __name__ == "__main__":
    which = int(sys.argv[1]); h = float(sys.argv[2]); tol = float(sys.argv[3])
    t0 = time.time()
    v0, b0, n0 = run(which, h, tol)
    v1, b1, n1 = run(which, h, tol, shift=h / 2)
    print("case %d h=%g tol=%g: K(shift 0) = %.10e  K(shift h/2) = %.10e  points %d %d  sub-threshold mass %.2e  time %.0fs"
          % (which, h, tol, v0, v1, n0, n1, b0, time.time() - t0))
