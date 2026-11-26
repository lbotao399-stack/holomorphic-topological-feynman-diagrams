"""K33 numerics (b3, b4). usage: python3 run_k33.py <h> <tol> <m> <p>
configs: [b4] lam0=e1, lam1=e2, lam2=(1/3,2/5), lam3=-(sum), lam4=lam5=0;
         [b3a] sources {0,1,3,5}: lam0=e1, lam1=e2, lam3=(1/3,2/5), lam5=-(sum);
         [b3b] sources {0,1,2}: lam0=e1, lam1=e2, lam2=-(e1+e2)."""
import sys, time
import numpy as np
from k33 import K33Int
from lattice2 import lattice
h, tol, m, p = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
lams = [[[1, 0], [0, 1], [1 / 3, 2 / 5], [-4 / 3, -7 / 5], [0, 0]],
        [[1, 0], [0, 1], [0, 0], [1 / 3, 2 / 5], [0, 0]],
        [[1, 0], [0, 1], [-1, -1], [0, 0], [0, 0]]]
I = K33Int(lams, m, p)
t0 = time.time()
v0, b0, n0 = lattice(I, 5, h, tol, 0.0, chunk=200)
v1, b1, n1 = lattice(I, 5, h, tol, h / 2, chunk=200)
print("K33 h=%.2f tol=%.0e m=%d p=%d points %d %.0fs" % (h, tol, m, p, n0 + n1, time.time() - t0))
for k, lab in enumerate(["b4 src{0,1,2,3}", "b3 src{0,1,3,5}", "b3 src{0,1,2}"]):
    print("  %s: shift0 %.8e  shift h/2 %.8e  boundary-layer mass %.1e" % (lab, v0[k], v1[k], b0))
