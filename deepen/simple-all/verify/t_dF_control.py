"""Control for va_points.dF: the same code with the QB monomials permuted (uu <-> uv) must give d != 0."""
import va_points, vpoint, numpy as np
orig = vpoint.Engine.F
def Fperm(self, X):
    F1, F2 = orig(self, X)
    perm = (3, 1, 2, 0, 4, 5)
    return [F1[k] for k in perm], [F2[k] for k in perm]
for key in ('so5', 'G2', 'E6'):
    print('unmodified:', va_points.dF(key, 1048573))
    vpoint.Engine.F = Fperm
    print('control (coefficients uu<->uv permuted):', va_points.dF(key, 1048573))
    vpoint.Engine.F = orig
