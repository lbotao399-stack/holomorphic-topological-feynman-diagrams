"""Normal-symbol Ward matrix for the trilinear one-loop tangent vector."""
from itertools import combinations, product
from sympy import LeviCivita, Matrix, symbols


def triple(p, i, j, k):
    return Matrix([p[i], p[j], p[k]]).det()


def derivative(p, v, repeated, kept):
    others = sorted(set(range(4)) - {repeated, kept})
    i, j = others
    shifted = list(p)
    shifted[repeated] = v
    return (triple(shifted, kept, repeated, i)*triple(p, kept, repeated, j)
            + triple(p, kept, repeated, i)*triple(shifted, kept, repeated, j))


def check():
    pairs = list(combinations(range(4), 2))
    tangent = [(1, 0), (0, 1), (2, 3), (3, 4)]
    normal = [(1,2,3), (2,-1,1), (4,0,1), (1,3,-2)]
    rows = []
    for j in range(4):
        basis = []
        for v in ((1,0,0), (0,1,0), (0,0,1)):
            numbers = [derivative(normal, v, i, j) for i in range(4) if i != j]
            assert sum(numbers) == 0
            basis.append(numbers)
            for a in range(2):
                rows.append([tangent[i][a]*derivative(normal, v, i, j) if j in pair else 0
                             for pair in pairs for i in [pair[0] if pair[1] == j else pair[1]]])
        assert Matrix(basis).rank() == 2
    assert Matrix(rows).rank() == 6
    w, r = symbols('w:3'), symbols('r:3')
    for s, i, j in product(range(3), repeat=3):
        residual = sum((LeviCivita(s,j,k)*w[i]-LeviCivita(s,i,k)*w[j]-LeviCivita(i,j,k)*w[s])*r[k]
                       for k in range(3))
        unit = sum(x*y for x, y in zip(w, r))
        assert (residual+LeviCivita(s,i,j)*unit).expand() == 0
    print('The strict shear Ward kernel has dimension zero; all 27 regular-source residuals agree.')


if __name__ == '__main__':
    check()
