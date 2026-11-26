"""Genus-two fivefold vertex with a moving level-one scalar on a hard mark.

The four other marked momenta are a*u, b*u, c*v and d*v. Pulling the hard
cotangent line back from Mbar(2,4) gives the threefold simplex period with
one psi class. The difference contributes the four primary fusions. A
rational tail containing two or more forgotten points has positive dimension
and its Jacobian is trivial, so it has no contribution with one psi insertion.
"""
from fractions import Fraction as F
from math import lcm

from flint import fmpq_mat
from bulk_eight import BULK
from marked_boundaries import quartics
from worldsheet_primary import fourfold_primary


def diagram_terms(hard,u,v,strengths):
    """Return the polarized ancestor followed by fusions with a*u,b*u,c*v,d*v."""
    if len(hard)!=4 or len(strengths)!=4 or any(len(p)!=5 for p in (*hard,u,v)):
        raise ValueError('Expected four hard fivefold covectors, two directions and four strengths')
    a,b,c,d = map(F,strengths)
    if any(h+(a+b)*i+(c+d)*j for h,i,j in zip(map(sum,zip(*hard)),u,v)):
        raise ValueError('All eight momenta must sum to zero')
    # On each forgotten pair the ordered two-point Abel map contributes
    # 2! times the product of its squared strengths. Base-dependent shifts
    # of the pair do not change the top-degree integral over its two marks.
    frame = fmpq_mat([[str(x) for x in p] for p in (*hard[:3],u,v)])
    delta = F(str(frame.det()))
    # Threefold graph periods 1/20160 + 29/544320 = 1/9720.
    ancestor = 4*(a*b*c*d)**2/F(9720)
    # Each fusion leaves the opposite pair saturated. The remaining fourfold
    # simplex contributes its fourth determinant power divided by 62208.
    tails = ((c*d)**2*b**4, (c*d)**2*a**4, (a*b)**2*d**4, (a*b)**2*c**4)
    return (delta**4*ancestor,*(delta**4*t/31104 for t in tails))


def primary(waves):
    scale = lcm(*(F(x).denominator for p in waves for x in p))
    integral = [tuple(int(scale*F(x)) for x in p) for p in waves]
    return -sum(w*coef for w,coef in zip(quartics(integral),BULK))/scale**20


def check_fixture(hard,u,v,strengths):
    probes = [tuple(s*x for x in direction) for s,direction in zip(strengths,(u,u,v,v))]
    expected = diagram_terms(hard,u,v,strengths)
    for mark in range(4):
        for j,probe in enumerate(probes):
            fused = list(hard)
            fused[mark] = tuple(x+y for x,y in zip(fused[mark],probe))
            fused.extend(p for i,p in enumerate(probes) if i!=j)
            assert primary(fused)==expected[j+1]
    return sum(expected)


def check():
    assert F(1,20160)+F(29,544320)==F(1,9720)
    basis = [tuple(int(i==j) for i in range(5)) for j in range(5)]
    u,v = basis[3:]
    simplex = [*basis[:3],(-1,-1,-1,0,0)]
    unit = (1,-1,1,-1)
    assert diagram_terms(simplex,u,v,unit)==(F(1,2430),*[F(1,31104)]*4)
    assert check_fixture(simplex,u,v,unit)==F(7,12960)
    # Verify the surviving fourfold normalization independently of its fivefold lift.
    projected = [tuple(x[:4]) for x in basis[:4]]+[(-1,-1,-1,-1)]
    assert fourfold_primary(2,projected)==(F(1,62208),0)

    # Unequal rational strengths, arbitrary longitudinal hard components, a
    # zero primary (string equation), and dependent axes.
    u = (1,-2,0,1,3)
    v = (-1,2,1,1,0)
    first = [(1,3,0,2,-1),(2,1,1,-1,0),(F(3,2),4,0,1,2)]
    strengths = (2,F(-3,5),4,F(1,3))
    for second,weights in ((v,strengths),(v,(0,*strengths[1:])),(u,strengths)):
        a,b,c,d = weights
        balance = tuple(-sum(col)-(a+b)*x-(c+d)*y for col,x,y in zip(zip(*first),u,second))
        total = check_fixture([*first,balance],u,second,weights)
        assert (total>0)==(second!=u)

    # Fully polarized output on root e^1, vector e^2 with polarization e_1
    # and scalar level-one e^3. All other fields are primary.
    assert sum(x*y for x,y in zip(basis[1],basis[0]))==0
    print('Fivefold labeled vertex:',F(7,12960))
    print('Elementary scalar in the A_1 source convention, times (2 pi i)^10:',-F(7,12960))


if __name__=='__main__':
    check()
