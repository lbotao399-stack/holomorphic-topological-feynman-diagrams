"""Cartan projection of the two contour orders in 2201.02595 (7.3.5)."""
from itertools import product

from sympy import Matrix, Rational, Symbol, diag, eye, factor, simplify, zeros


def comm(x, y):
    return x*y - y*x


def verify(entries):
    size = len(entries)
    x = diag(*entries)
    assert x.trace() == 0
    units = eye(size)
    basis = [units[:, i] * units[:, j].T
             for i, j in product(range(size), repeat=2) if i != j]
    for i in range(size-1):
        h = zeros(size)
        h[i, i], h[i+1, i+1] = 1, -1
        basis.append(h)
    metric = Matrix([[(u*v).trace() for v in basis] for u in basis])
    inverse = metric.inv()

    # A contour around the displaced J first produces [[X,t_j],X].
    # The second contraction attaches this to the zero-mode [X,t_i].
    colour = sum(inverse[i, j] * (x*comm(comm(x, u), comm(comm(x, v), x))).trace()
                 for i, u in enumerate(basis) for j, v in enumerate(basis) if inverse[i, j])
    fourth = sum((a-b)**4 for a, b in product(entries, repeat=2))
    assert fourth != 0 and colour == -fourth
    squared = (x*x).trace()**2
    lam_squared = fourth/squared
    period = Symbol('period', nonzero=True)  # 2 pi i
    coupling_squared = lam_squared / (12 * period**3)
    coefficient = factor(3*coupling_squared*squared / (2*colour))
    inner_then_outer = 2*coefficient*colour - coupling_squared*squared
    outer_channels = 2*coupling_squared*squared
    assert simplify(inner_then_outer-outer_channels) == 0
    assert coefficient == -Rational(1, 8)/period**3
    print(entries, 'Tr_ad(X^4) =', fourth, 'one ordered colour =', colour,
          'C =', coefficient)


if __name__ == '__main__':
    verify((1,-1))
    verify((2,1,-3))
