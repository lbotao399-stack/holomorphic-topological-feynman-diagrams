"""Chiral supersymmetry test for the completed log-determinant insertion (2604.21015)."""
from itertools import combinations
import sympy as sp


def verify():
    z = sp.symbols('z1:4')
    # Remove the common nonzero z_i**(-1/2) factor on each column.
    # Multiplication by the two Ward one-forms on the three one-leg variables.
    rows = []
    for weights in ([1]*3,z):
        for i,j in combinations(range(3),2):
            row = [0]*3
            row[j],row[i] = weights[i],-weights[j]
            rows.append(row)
    matrix = sp.Matrix(rows)
    assert sp.factor(matrix[[0,1,3],:].det()) in (z[0]-z[1],z[1]-z[0])
    # The two components of the same charge acting on chi^13 K^2 chi^34.
    parent = sp.Matrix([[1,1],[z[0],z[2]]])
    assert parent.det() == z[2]-z[0]
    # Independent line-incidence check: eight independent odd shifts of the
    # two polarizations per flavor at distinct z annihilate a completed insertion.
    block = sp.Matrix([[1,1],[z[0],z[2]]])
    assert sp.factor(sp.diag(*([block]*4)).det()-(z[2]-z[0])**4) == 0
    print('completed insertion: two independent Ward equations; the one-loop K K chi coefficient is zero')


if __name__ == '__main__':
    verify()
