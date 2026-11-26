"""Exact sign of the genus-three trace on two distinct real collinear pairs."""
from collections import Counter
from fractions import Fraction
from functools import cache
from itertools import combinations
from math import comb

from flint import fmpz_mat, fmpz_mpoly_ctx

from two_pair_scalar import certificate


def chamber_coefficients(trace):
    # Pair strengths are m+h,m-h and n+k,n-k. U,V are the first two
    # coordinates of the fourth unpaired momentum after the shift below.
    ctx = fmpz_mpoly_ctx.get(('m','h','U','n','k','V','gap1','gap2','gap3','gap4'))
    m,h,u,n,k,v,*gaps = ctx.gens()
    gaps[0] += gaps[3]
    total = sum(gaps)
    partials = [sum(gaps[:j]) for j in (1,2,3)]
    center = sum(partials)+total
    # The mean of the four highest ordered projected cofactor coordinates
    # is center/(4*total). Translate each fiber coordinate by its pair sum
    # times that mean. Multiplication by (2*total)^4 clears these fractions.
    @cache
    def pairs(i,j,b,d):
        return ((m-h)**i*(m+h)**(2-i-b)*(n-k)**j*(n+k)**(2-j-d)
                *(2*total*u-m*center)**b*(2*total*v-n*center)**d
                *(2*total)**(4-b-d))

    @cache
    def projected(i,j,k):
        return ((-partials[0])**i*(-partials[1])**j*(-partials[2])**k
                *total**(6-i-j-k))

    s,t,*_ = trace.context().gens()
    reduced = trace//(s*s*t*t)
    assert reduced*s*s*t*t == trace
    polynomial = ctx.constant(0)
    for (i,j,b,d,x,y,z),a in reduced.terms():
        assert i+b <= 2 and j+d <= 2 and x+y+z <= 6
        polynomial += a*pairs(i,j,b,d)*projected(x,y,z)
    blocks = {}
    for powers,value in polynomial.terms():
        assert sum(powers[:3]) == sum(powers[3:6]) == 2
        assert powers[1]%2 == powers[4]%2 == 0 and sum(powers[6:]) == 10
        blocks.setdefault(powers[6:],{})[powers[:6]] = int(value)
    return blocks


def gram_matrix(block):
    # Four times the Gram matrix on (m,h,U) tensor (n,k,V).
    gram = fmpz_mat(9,9)
    for i in range(9):
        for j in range(i,9):
            powers = [0]*6
            for slot in (i,j):
                powers[slot//3] += 1
                powers[3+slot%3] += 1
            multiplicity = (1 if i==j else 2)*(2 if i//3!=j//3 and i%3!=j%3 else 1)
            gram[i,j] = gram[j,i] = (4//multiplicity)*block.get(tuple(powers),0)
    if gram[0,0] == 0:
        # (mn)(UV) and (mV)(Un) are the same scalar monomial.
        gram[2,6] += gram[0,8]
        gram[6,2] = gram[2,6]
        gram[0,8] = gram[8,0] = 0
    return gram


def check(channels):
    partial = sum(Fraction(1,j**3) for j in range(1,257))
    assert Fraction(1202056,10**6) < partial+Fraction(1,2*257**2)
    assert partial+Fraction(1,2*256**2) < Fraction(1202057,10**6)
    rational, zeta = [chamber_coefficients(p) for p,_ in channels]
    assert set(rational) == set(zeta) and len(rational) == comb(13,3)
    for endpoint in (1202056,1202057):
        ranks = Counter()
        for powers in rational:
            first,second = rational[powers],zeta[powers]
            block = {k:-10**6*channels[1][1]*first.get(k,0)
                     -endpoint*channels[0][1]*second.get(k,0)
                     for k in first.keys()|second.keys()}
            gram = gram_matrix(block)
            groups = ((0,2,6,8),(1,7),(3,5),(4,))
            for group in groups:
                assert all(gram[i,j]==0 for i in group for j in range(9) if j not in group)
                for size in range(1,len(group)+1):
                    for sub in combinations(group,size):
                        assert fmpz_mat([[gram[i,j] for j in sub] for i in sub]).det()>=0
            rank = gram.rank()
            ranks[rank] += 1
            if 10 in powers:
                if powers[0]==10:
                    assert rank==5 and [i for i in range(9) if not gram[i,i]]==[0,1,3,4]
                else:
                    assert rank==9
        assert ranks=={9:282,5:4}
    print('Both rational endpoints: 282 positive definite and four rank-five positive semidefinite forms.')


if __name__=='__main__':
    check(certificate())
