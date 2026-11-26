"""Own-coordinate QMC (no engine input) of K for theta(1,2,2) at H = 1 (banana construction of b21.py),
versus -pi l1 l2/3 (exact, v_a5_theta) and the a5 closed form +pi l1 l2/3."""
import sys, math
import numpy as np
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/ht-families/verify")
from v_b2 import banana
from v_b1 import own_qmc
nv, edges, ints = banana((1, 2, 2))
print("edges", edges, "interiors", ints)
for lam in ([2, 3, -5], [-1, 4, 7]):
    mean, err, nex = own_qmc(nv, edges, np.array([[float(x)] for x in lam]), 15, 16, 7)
    l1, l2 = lam[1], lam[2]
    print("lambda %s: own QMC K = %.10g +- %.2g ; -pi l1 l2/3 = %.10g ; a5 closed form +pi l1 l2/3 = %.10g" % (
        lam, mean, err, -math.pi * l1 * l2 / 3, math.pi * l1 * l2 / 3))
