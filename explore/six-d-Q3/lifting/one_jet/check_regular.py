"""Exact jet-slice reconstruction and interpolation on regular diagonal matrices."""

from itertools import product
from math import comb, factorial
import sympy as sp


class MatrixJets:
    def __init__(self, rank, order):
        self.rank = rank
        self.order = order
        self.zero = (0, 0, 0)
        self.degrees = sorted((k for k in product(range(order+1), repeat=3) if sum(k) <= order),
                              key=lambda k: (sum(k), k))
        self.identity = {self.zero: sp.eye(rank)}

    def multiply(self, left, right):
        result = {}
        for a, x in left.items():
            for b, y in right.items():
                degree = tuple(i+j for i, j in zip(a, b))
                if sum(degree) <= self.order:
                    result[degree] = result.get(degree, sp.zeros(self.rank)) + x*y
        return result

    def exponential(self, positive, sign=1):
        result = dict(self.identity)
        term = self.identity
        for i in range(1, self.order+1):
            term = self.multiply(term, positive)
            for k, m in term.items():
                result[k] = result.get(k, sp.zeros(self.rank)) + m * sp.Rational(sign**i, factorial(i))
        return result

    def conjugate(self, parameter, matrix, sign):
        return self.multiply(self.multiply(self.exponential(parameter, sign), matrix),
                             self.exponential(parameter, -sign))


def verify(rank, order):
    ring = MatrixJets(rank, order)
    diagonal, known = {}, {}
    for number, degree in enumerate(ring.degrees):
        if number == 0:
            roots = [i+1 for i in range(rank-1)]
            roots.append(-sum(roots))
            diagonal[degree] = sp.diag(*roots)
            continue
        entries = [(number+1)**(i+1) for i in range(rank-1)]
        diagonal[degree] = sp.diag(*entries, -sum(entries))
        known[degree] = sp.Matrix(rank, rank, lambda i,j:
                                  0 if i == j else ((3*number+2*i+j) % 7)-3)
    matrix = ring.conjugate(known, diagonal, 1)
    reconstructed = {}
    for total in range(1, order+1):
        residual = ring.conjugate(reconstructed, matrix, -1)
        for degree in ring.degrees:
            if sum(degree) == total:
                reconstructed[degree] = sp.Matrix(rank, rank, lambda i,j:
                    0 if i == j else -residual[degree][i,j] / (roots[i]-roots[j]))
    assert reconstructed == known and ring.conjugate(reconstructed, matrix, -1) == diagonal

    for i, j in product(range(rank), repeat=2):
        if i == j:
            continue
        mat = sp.zeros(len(ring.degrees))
        for m, target in enumerate(ring.degrees):
            for n, source in enumerate(ring.degrees):
                delta = tuple(a-b for a,b in zip(target, source))
                if min(delta) >= 0:
                    mat[m,n] = diagonal[delta][i,i] - diagonal[delta][j,j]
        assert mat.is_lower and mat.det() == (roots[i]-roots[j])**comb(order+3, 3)

    trace_rows = sp.Matrix([[roots[i]**k for i in range(rank)] for k in range(rank)])
    assert trace_rows.det() == sp.prod(roots[j]-roots[i] for i in range(rank) for j in range(i+1,rank))
    for m in range(len(ring.degrees)):
        values = sp.Matrix([(-1)**(i+m)*(i+2)**(m%3) for i in range(rank)])
        assert trace_rows.inv()*(trace_rows*values) == values


def verify_weyl():
    # B2 acts by signed permutations. G2 acts in dihedral coordinates
    # by interchange and opposite sixth-root rotations. z is central.
    a, b, z = sp.symbols("a b z")
    cases = (
        ((z, a*a+b*b, a*a*b*b), 4*a*b*(a*a-b*b)),
        ((z, a*b, a**6+b**6), 6*(b**6-a**6)),
    )
    for generators, reflecting_product in cases:
        jac = sp.Matrix(generators).jacobian((z, a, b))
        assert sp.expand(jac.det() - reflecting_product) == 0
        at_regular = jac.subs({a: 2, b: 3, z: 7})
        direction = sp.Matrix([2, -1, 3])
        assert at_regular.inv() * (at_regular * direction) == direction


if __name__ == "__main__":
    verify_weyl()
    print("B2 and G2 with one central factor: Weyl Jacobians exact")
    for rank, order in ((2,3), (3,2), (4,2)):
        verify(rank, order)
        print(f"sl{rank}, jet order {order}: slice, root units and Cartan traces exact")
