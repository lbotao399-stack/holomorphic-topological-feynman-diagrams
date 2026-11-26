"""Ordered Poisson vertex kernel in the Li–Liu–Chang fan orientation."""
from collections import defaultdict
from fractions import Fraction
from math import comb, factorial, prod


def bracket(p, q):
    return p[0] * q[1] - p[1] * q[0]


def dot(p, z):
    return p[0] * z[0] + p[1] * z[1]


def zero_displacement(n, p, q):
    if n < 1:
        raise ValueError("the quantum fan starts at one loop")
    denominator = 2**n * factorial(n) * factorial(n + 1) * factorial(2*n + 1)
    return Fraction(bracket(p, q)**(2*n + 1), denominator)


def rising(k, degree):
    return prod(range(k, k + degree))


def shifted(power, offset):
    """Coefficients of (x + offset)**power in increasing powers of x."""
    return [comb(power, i) * offset**(power-i) for i in range(power+1)]


def taylor(n, p, q, path, spokes, order):
    """Return coefficients through order in exp(-u sum_e y_e.z_e), including the kernel.

    path and spokes follow the orientation in the paper. Mixed jet coefficients
    follow by polarization. All entries may be integers or rational numbers.
    """
    if len(path) != n or len(spokes) != n+1 or order < 0:
        raise ValueError("expected n path displacements, n+1 spokes and nonnegative order")
    kernel = zero_displacement(n, p, q)
    if not kernel:
        return [Fraction(0)] * (order + 1)

    offsets = []
    for i, z in enumerate(path):
        value = tuple(z[j] + spokes[i+1][j] - spokes[i][j] for j in (0, 1))
        offsets.append((dot(p, value), dot(q, value)))

    # A monomial (k,i,j) is s**k * (p.D)**i * (q.D)**j.  Here s*p is the
    # accumulated source, and D the displacement transferred along the path.
    layers = [defaultdict(Fraction) for _ in range(order+1)]
    xp, xq = offsets[-1]
    for m in range(order+1):
        for a in range(m+1):
            b = m-a
            weight = Fraction((-1)**a * rising(2, a) * rising(2*n, b),
                              factorial(a)*factorial(b)*rising(2*n+3, m))
            for i, left in enumerate(shifted(a, xp)):
                for j, right in enumerate(shifted(b, xq)):
                    layers[m][a, i, j] += weight * left * right

    for vertex in range(n-1, 0, -1):
        xp, xq = offsets[vertex-1]
        new = [defaultdict(Fraction) for _ in layers]
        alpha = 2*(n-vertex)+2
        for m, coefficients in enumerate(layers):
            for (k, i, j), coefficient in coefficients.items():
                if not coefficient:
                    continue
                for extra in range(order-m+1):
                    odd_moment = Fraction(alpha, alpha+k+extra)
                    even_moment = Fraction(rising(3, extra)*rising(2*vertex, i+j),
                                           rising(3+2*vertex, extra+i+j))
                    weight = coefficient * odd_moment * even_moment * Fraction((-1)**extra, factorial(extra))
                    for a, left in enumerate(shifted(i+extra, xp)):
                        for b, right in enumerate(shifted(j, xq)):
                            new[m+extra][k+extra, a, b] += weight * left * right
        layers = new

    # Set s=1 and D=0, then restore the phase of the two endpoint spokes.
    reduced = [sum(c for (k, i, j), c in layer.items() if i == j == 0) for layer in layers]
    constant = -dot(p, spokes[0]) - dot(q, spokes[-1])
    return [kernel * sum(reduced[m-l] * Fraction(constant**l, factorial(l)) for l in range(m+1))
            for m in range(order+1)]
