"""Scalar-source projections of a transverse fivefold bulk covector potential."""
from itertools import combinations
from random import Random

from soft_consistency import interior, polarized_weights, vector


def pairing(p,v):
    return sum(x*y for x,y in zip(p,v))


def test_potential(p):
    # Symmetric Green matrix: the constant five-vector identity is exact
    # independently of the individual numerical transfer coefficients.
    green = [[pairing(x,y)**2 for y in p] for x in p]
    return [tuple(sum(g*y[k] for g,y in zip(row,p)) for k in range(5)) for row in green]


def rank_four(p,a,k,i,j,v,w):
    left = tuple(x-y for x,y in zip(a[i],a[k]))
    right = tuple(x-y for x,y in zip(a[j],a[k]))
    return pairing(left,v)*pairing(p[k],w)-pairing(right,w)*pairing(p[k],v)


def generic_plane(rng,n):
    momenta = [tuple(rng.randrange(-3,4) for _ in range(5)) for _ in range(n-1)]
    return [*momenta,tuple(-sum(c) for c in zip(*momenta))]


def tangent(p):
    # Two small rational-independent transverse coefficients suffice for the fixture.
    return (p[1],-p[0],p[3],-p[2],0)


def source_projection(p,i,j,v,w,source_size):
    scalars = [k for k in range(len(p)) if k not in (i,j)]
    value = 0
    for group in combinations(scalars,source_size):
        other = [k for k in range(len(p)) if k not in group]
        q = [*(p[k] for k in other),tuple(sum(p[k][d] for k in group) for d in range(5))]
        a = test_potential(q)
        value += rank_four(q,a,len(q)-1,other.index(i),other.index(j),v,w)
    return value


def variational_check():
    import sympy as sp
    from marked_consistency import du, derivative, u, variational_derivative

    # A one-form J=(du_2)^2 dz^1 gives the bivector coefficient L=du wedge J.
    # Its derivative is nonzero, while the single-trivector condition kills it.
    one = [du[1]**2,0*u,0*u,0*u,0*u]
    two = {(i,j):du[i]*one[j]-du[j]*one[i] for i,j in combinations(range(5),2)}
    three = {(i,j,k):derivative(two[j,k],i)-derivative(two[i,k],j)+derivative(two[i,j],k)
             for i,j,k in combinations(range(5),3)}
    assert any(three.values())
    for indices in combinations(range(5),4):
        coefficient = sum((-1)**j*du[i]*three[indices[:j]+indices[j+1:]]
                          for j,i in enumerate(indices))
        assert sp.expand(coefficient)==0

    # Start with alpha(h)=integral (D_1 u)^2 h. Its exterior derivative
    # after integration by parts is K[h]=-2(D_1^2 u)h-4(D_1 u)D_1 h.
    second = derivative(du[0],0)
    primitive = sp.expand((-2*u*second+4*derivative(u*du[0],0))/3)
    exact_difference = variational_derivative(-u*du[0]**2/3)
    assert sp.expand(primitive-du[0]**2-exact_difference)==0


def check():
    rng = Random(117)
    p = generic_plane(rng,12)
    a = test_potential(p)
    for i,j in combinations(range(5),2):
        assert sum(x[i]*y[j]-x[j]*y[i] for x,y in zip(p,a))==0
    # Keep the two vectors in index order and put the rank-four slot first;
    # moving an even four-vector through other slots has no Koszul sign.
    total = 0
    for i,j in combinations(range(12),2):
        v,w = tangent(p[i]),tangent(p[j])
        for k in range(12):
            if k in (i,j):
                continue
            blade = interior(p[k],((31,1),))
            coeffs = polarized_weights((blade,vector(v),vector(w)))
            direct = sum(x*(a[k],a[i],a[j])[slot][c] for slot,c,x in coeffs)
            assert direct==rank_four(p,a,k,i,j,v,w)
            total += abs(direct)
        assert source_projection(p,i,j,v,w,1)==0
    assert total>0
    # A nonlinear source has no constant-Hamiltonian cancellation in general.
    extended = generic_plane(rng,13)
    nonlinear = source_projection(extended,0,1,tangent(extended[0]),tangent(extended[1]),2)
    assert (total,nonlinear)==(531781712,6544)
    variational_check()
    print('Twelve-wave rank-four detector, total absolute response:',total)
    print('Linear scalar insertion vanishes; quadratic-source control:',nonlinear)
    print('Bivector Koszul identity and variational scaling primitive checked.')


if __name__=='__main__':
    check()
