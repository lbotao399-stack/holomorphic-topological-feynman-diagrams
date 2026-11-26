"""Equivariant polynomial sewing inverse for a fivefold reciprocal potential."""
from fractions import Fraction
from random import Random

from flint import fmpz_mat


def det(rows):
    return int(fmpz_mat(rows).det())


def minor_and_derivatives(p,slots):
    rows = [p[k] for k in slots]
    value = det(rows)
    # This divergence differentiates slot k in the direction of p_1.
    deriv = [0]*len(p)
    for j,k in enumerate(slots):
        deriv[k] = det(rows[:j]+[p[0]]+rows[j+1:])
    return value,deriv


def sew(weights,p):
    return [[sum(w*q[i] for w,q in zip(row,p)) for i in range(5)]
            for row in weights]


def check(p,m):
    n = len(p)
    d,dd = minor_and_derivatives(p,list(range(5)))
    assert d
    for long_hook in (False,True):
        if long_hook and (n<6 or not m):
            continue
        weight = Fraction(m+4,2) if long_hook else m+5
        f = [0]*n
        derivative = [[0]*n for _ in p]
        if long_hook:
            for i in range(6):
                minors,changes = minor_and_derivatives(p,[j for j in range(6) if j!=i])
                f[i] = (-1)**i*minors*d**(m-1)
                for j in range(n):
                    derivative[i][j] = (-1)**i*changes[j]*d**(m-1)
                    if m>1:
                        derivative[i][j] += (-1)**i*minors*(m-1)*d**(m-2)*dd[j]
            # The six-row determinant relation is the first term in sewing
            # the highest vector with a six-box column.
            assert all(sum(p[i][c]*f[i] for i in range(n))==0 for c in range(5))
            assert f[5]==-d**m
            if m>=2:
                cofactors = [value//d**(m-1) for value in f]
                syzygy = [[x*y*d**(m-2) for x in cofactors] for y in cofactors]
                assert syzygy[5][5]==d**m and not any(any(row) for row in sew(syzygy,p))
        else:
            f[0] = d**m
            if m:
                derivative[0] = [m*d**(m-1)*v for v in dd]

        potential = [[p[0][c]*f[i] for c in range(5)] for i in range(n)]
        diverged = [[Fraction(derivative[i][j]+derivative[j][i]
                      +5*((j==0)*f[i]+(i==0)*f[j]),2) for j in range(n)] for i in range(n)]
        expected = [[weight*x for x in a] for a in potential]
        assert sew(diverged,p)==expected
        alpha,beta = m+5,Fraction(m+4,2)
        inverse = [[(alpha+beta-weight)*x/(alpha*beta) for x in row] for row in diverged]
        assert sew(inverse,p)==potential

        # Lift to N redundant flavors, choose zero last row, and compare
        # every row difference after imposing the canonical zero diagonal.
        lifted = [row+[0] for row in inverse]+[[0]*(n+1)]
        balanced = p+[tuple(-sum(c) for c in zip(*p))]
        zero_diagonal = [[x-Fraction(lifted[i][i]+lifted[j][j],2) for j,x in enumerate(row)]
                         for i,row in enumerate(lifted)]
        rows = sew(zero_diagonal,balanced)
        restored = [[v-last for v,last in zip(row,rows[-1])] for row in rows[:-1]]
        assert restored==potential


if __name__=='__main__':
    rng = Random(1269)
    for size in (5,6,11):
        p = [tuple(rng.randrange(-5,6) for _ in range(5)) for _ in range(size)]
        for m in (0,1,4,8):
            check(p,m)
    print('Five-, six- and eleven-flavor components, m=0,1,4,8, inverse and zero-diagonal lift passed.')
