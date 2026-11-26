"""All-genus normal pair tensors on five projected waves with two fewer planar derivatives available."""
from fractions import Fraction
from functools import cache
from itertools import combinations
from math import prod

from flint import fmpz_mpoly_ctx


GENERATORS = ((1,2),(1,3),(1,4),(2,4))


def signed_triple(z,i,j,k):
    if len({i,j,k})<3:
        return 0
    remaining = [s for s in range(5) if s not in (i,j,k)]
    a,b = remaining
    signs = sum(s>t for s,t in combinations((i,j,k),2)) + a+b-1
    return (-1)**signs*(z[a]-z[b])


def generators(centered,pair):
    i,j = pair
    a,b = centered[i],centered[j]
    out = []
    for r,s in GENERATORS:
        out.append(sum(((a**r-b**r)*w**s+(a**s-b**s)*w**r)*signed_triple(centered,i,j,k)
                       for k,w in enumerate(centered)))
    return out


@cache
def powers(degree):
    return [(a,b,c,d) for a in range(degree//2+1) for b in range(degree//3+1)
            for c in range(degree//4+1) for d in range(degree//5+1)
            if 2*a+3*b+4*c+5*d == degree]


def tensor_row(z,pair,degree):
    mean = Fraction(sum(z),5)
    z = tuple(t-mean for t in z)
    elementary = [sum(prod(c) for c in combinations(z,k)) for k in range(2,6)]
    return [v*prod(w**e for w,e in zip(elementary,exponents))
            for (r,s),v in zip(GENERATORS,generators(z,pair))
            for exponents in powers(degree-1-r-s)]


def dimensions(genus):
    if genus<2:
        raise ValueError("The two planar transfer periods require genus at least two")
    return tuple(sum(len(powers(2*genus-d-r-s)) for r,s in GENERATORS) for d in range(3))


def check():
    from quartic_seven_plane_space import PAIRS,as_matrix,normal_row,FRAMES
    # Polynomial row sums and momentum-relation contractions for the four explicit generators.
    ctx = fmpz_mpoly_ctx.get(('x0','x1','x2','x3'))
    z = (*ctx.gens(),-sum(ctx.gens()))
    rows = {ij:generators(z,ij) for ij in PAIRS}
    for i in range(5):
        assert all(not sum(col) for col in zip(*(v for ij,v in rows.items() if i in ij)))
    assert all(not sum(z[i]*z[j]*v[d] for (i,j),v in rows.items()) for d in range(4))
    # The Hessian of the degree-ten Vandermonde supplies the full alternating target.
    five = fmpz_mpoly_ctx.get(('w0','w1','w2','w3','w4'))
    w = five.gens()
    delta = prod(w[i]-w[j] for i,j in PAIRS)
    diagonal = [delta.derivative(i).derivative(i) for i in range(5)]
    assert not sum(diagonal)
    out = {ij:3*delta.derivative(ij[0]).derivative(ij[1])+diagonal[ij[0]]+diagonal[ij[1]]
           for ij in PAIRS}
    assert all(not sum(v for ij,v in out.items() if i in ij) for i in range(5))
    assert 2*sum(w[i]*w[j]*v for (i,j),v in out.items()) == 170*delta
    special = (-2,-1,0,1,2)
    assert as_matrix([generators(special,PAIRS[k]) for k in (0,1,2,4)]).det() == 82944
    counts = [((2 if n==2 else 3)*dims[0]+dims[1]+dims[2])
              for n in range(2,12) for dims in [dimensions(n)]]
    assert counts == [3,12,25,44,71,108,155,215,287,373]
    # Same polynomial spaces as the direct quartic Ward equations, also on a held-out frame.
    for d in range(3):
        table,old = [],[]
        for frame in FRAMES:
            for pair in PAIRS:
                row = tensor_row(frame,pair,9-d)
                if as_matrix(table+[row]).rank()>len(table):
                    table.append(row)
                    old.append(normal_row(frame,pair,d))
        change = as_matrix(table).inv()*as_matrix(old)
        assert change.det()
        frame = (1,2,4,8,16)
        for pair in PAIRS:
            assert as_matrix([tensor_row(frame,pair,9-d)])*change == as_matrix([normal_row(frame,pair,d)])
    print('Four free normal generators; seven-planar ranks 6,4,3. Genus counts',counts)


if __name__=='__main__':
    check()
