# Exact check: S3 eq-rho2 at the dipole Laplacian M = [[x1+x3, x3],[x3, x2+x3]]
# equals (-1/2) * BHP eq. (5.5) (i=1) and pi * BH Def. 2 / eq. (1.5).
from sympy import symbols, sqrt, simplify, pi, Rational
x1, x2, x3 = symbols('x1 x2 x3', positive=True)
# 2-forms in basis (dx1dx2, dx1dx3, dx2dx3)
def wedge(u, v):  # u, v: 1-forms as (c1,c2,c3)
    return (u[0]*v[1]-u[1]*v[0], u[0]*v[2]-u[2]*v[0], u[1]*v[2]-u[2]*v[1])
M11, M12, M22 = x1+x3, x3, x2+x3
dM11, dM12, dM22 = (1, 0, 1), (0, 0, 1), (0, 1, 1)
det = M11*M22 - M12**2
num = [M11*a - M12*b + M22*c for a, b, c in zip(wedge(dM12, dM22), wedge(dM11, dM22), wedge(dM11, dM12))]
rho2 = [-n/(2*det**Rational(3, 2)) for n in num]                    # S3 eq-rho2
Psi = x1*x2 + x1*x3 + x2*x3
Omega3 = (x3, -x2, x1)                                              # BHP (5.1), m=3: x1 dx2dx3 - x2 dx1dx3 + x3 dx1dx2
phiBHP = [-o/Psi**Rational(3, 2) for o in Omega3]                   # BHP (5.5), i=1
phiBH = [o/(2*pi*Psi**Rational(3, 2)) for o in Omega3]              # BH (1.5) as read
print("det M == Psi:", simplify(det - Psi) == 0)
print("rho2 == (-1/2) phi_BHP:", all(simplify(r + phi/2) == 0 for r, phi in zip(rho2, phiBHP)))
print("rho2 == pi phi_BH:", all(simplify(r - pi*phi) == 0 for r, phi in zip(rho2, phiBH)))
