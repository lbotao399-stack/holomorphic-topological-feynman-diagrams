"""Coadjoint cocycles and the trace pairing of constant-field BF sections."""
from itertools import product
from sympy import Matrix, eye, zeros


def check(n, centre=False):
    units = eye(n)
    basis = [units[:,i]*units[:,j].T for i,j in product(range(n),repeat=2) if i != j]
    for i in range(n-1):
        h = zeros(n)
        h[i,i], h[i+1,i+1] = 1,-1
        basis.append(h)
    if centre:
        basis.append(eye(n))
    dim = len(basis)
    metric = Matrix([[(x*y).trace() for y in basis] for x in basis])
    duals = [sum((coef*y for coef,y in zip(row,basis)),zeros(n)) for row in metric.inv().tolist()]
    brackets = [[[((x*y-y*x)*z).trace() for z in duals] for y in basis] for x in basis]
    # phi([x,y])(z) = -phi(y)([x,z])+phi(x)([y,z]).
    rows = []
    for i,j,k in product(range(dim),repeat=3):
        if i >= j:
            continue
        row = [0]*(dim**2)
        for h in range(dim):
            row[h*dim+k] += brackets[i][j][h]
            row[j*dim+h] += brackets[i][k][h]
            row[i*dim+h] -= brackets[j][k][h]
        rows.append(row)
    cocycles = Matrix(rows).nullspace()
    pullbacks = []
    for col in cocycles:
        form = col.reshape(dim,dim)
        sym = form + form.T
        assert all(sym[i,j] == 0 for i,j in product(range(dim),repeat=2)
                   if not centre or i != dim-1 or j != dim-1)
        pullbacks.append(sym.reshape(dim**2,1))
    rank = Matrix.hstack(*pullbacks).rank()
    assert len(cocycles) == n*n-1+int(centre) and rank == int(centre)
    print(('gl' if centre else 'sl')+str(n), 'cocycles',len(cocycles),'image of symmetric part',rank)


if __name__ == '__main__':
    check(2)
    check(3)
    check(2,centre=True)
