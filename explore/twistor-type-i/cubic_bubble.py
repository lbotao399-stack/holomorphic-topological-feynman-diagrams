"""Two independent reductions of the (3,0,0) three-point WZW momentum master."""
from itertools import permutations
import sympy as s


def on_shell_pair(i, j):
    return (0, 1, -1)[(j - i) % 3]


def main():
    z = s.symbols('z')
    x, y = s.symbols('x y')
    # Shift q to its mean -x p1 + y p3. Normalize the three cyclic
    # Kaehler pairings to one only in the polynomial check.
    mean = [-x * on_shell_pair(i, 1) + y * on_shell_pair(i, 3) for i in (1, 2, 3)]
    product = mean[0] * (mean[1] - 1) * mean[2]
    assert s.expand(product + x*y*(1-x-y)) == 0
    # The Dirichlet integral gives -z/Gamma(5+z) after removing the
    # common Gamma(2+z-d/2). Recompute by a large-q tensor expansion.
    coefficient = 0
    for derivative_legs, vertex_sign in (((1, 2, 3), 1), ((1, 3), -1)):
        h = len(derivative_legs)
        for r in range(h+1):
            # r insertions from (q+p1)^2, the others from (q-p3)^2.
            routing = (1,)*r + (3,)*(h-r)
            pairings = sum(s.prod(on_shell_pair(i, j) for i, j in zip(derivative_legs, ordering))
                           for ordering in permutations(routing))
            # The angular moment and the massive radial integral cancel
            # their dimension-dependent rising factorials.
            coefficient += vertex_sign * (-1)**r * pairings / s.gamma(2+z+h)
    assert s.simplify(s.expand_func(coefficient * s.gamma(5+z)) + z) == 0
    # At z=epsilon, d=4-2epsilon: epsilon*Gamma(2epsilon) -> 1/2.
    remaining_integral = -s.Rational(1, 2) / s.factorial(4)
    bubble_residue = -s.Rational(1, 4*3)
    master = bubble_residue * remaining_integral
    assert master == s.Rational(1, 576)
    print('Three-point (3,0,0) momentum pole in units of omega^3/(4*pi)^4:', master)


if __name__ == '__main__':
    main()
