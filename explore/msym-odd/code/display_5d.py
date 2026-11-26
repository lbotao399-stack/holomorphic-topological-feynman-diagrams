"""Explicit X_v maps and densities of the 5d Q_2 kernels (format of S2 eqs. (1.6)-(1.9))."""
import sympy as sp
import htk_path  # noqa
from ops5 import KernelData

for kind, names in (('A', ['1', '2', '4', '5']), ('B', ['u', 'w', 'a', 'b'])):
    K = KernelData(kind)
    TK = K.TK
    z = sp.symbols('z1 z2 z3')
    print("graph (%s): paths (edges) %s, n_P = %s, alpha_P = %s, s = %+d" % (kind, TK.paths, TK.nP, [str(a) for a in TK.alpha], K.s))
    print("   D_%s(lam) =" % kind, sp.factor(K.D))
    for v, nm in enumerate(names):
        X = sum(K.cz[v][j] * z[j] for j in range(3))
        print("   X_%s = %s" % (nm, sp.collect(sp.expand(X), z)))
    # coefficient check: each X_v is an affine combination of z1, z2, z3 (coefficients sum to 1)
    print("   sum of coefficients of z in X_v:", [sp.simplify(sum(K.cz[v])) for v in range(4)])
