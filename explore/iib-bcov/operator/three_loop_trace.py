"""Nine-field trace boundary on two triples of collinear bulk momenta."""
from itertools import combinations
from math import prod
from flint import fmpz_mpoly_ctx


def determinant(rows):
    if len(rows)==1:
        return rows[0][0]
    # Expand the sparsest row in a matrix of linear forms.
    pivot = min(range(len(rows)),key=lambda i: sum(bool(x) for x in rows[i]))
    rest = rows[:pivot]+rows[pivot+1:]
    return sum((-1)**(pivot+j)*x*determinant([r[:j]+r[j+1:] for r in rest])
               for j,x in enumerate(rows[pivot]) if x)


def trace(rows):
    # Unordered triples of five-slot sets, covering nine labels, with
    # pair-intersection sizes {2,2,3}. Every minor is squared.
    nonzero = []
    for subset in combinations(range(9),5):
        d = determinant([rows[i] for i in subset])
        if d:
            nonzero.append((sum(1<<i for i in subset),d*d))
    out = 0
    for pos, (a,x) in enumerate(nonzero):
        for i in range(pos):
            b,y = nonzero[i]
            ab = (a&b).bit_count()
            if ab not in (2,3):
                continue
            for c,z in nonzero[:i]:
                if a|b|c == 511 and sorted((ab,(b&c).bit_count(),(a&c).bit_count()))==[2,2,3]:
                    out += x*y*z
    return out


def verify():
    ctx = fmpz_mpoly_ctx.get(('a1','a2','a3','b1','b2','b3'))
    factors = ctx.gens()
    axes = (0,0,0,1,1,1,2,3,4)
    p = [[s*int(i==j) for j in range(5)] for s,i in zip((*factors,1,1,1),axes)]
    p.append([-sum(column) for column in zip(*p)])
    potential = [[0]*5 for _ in range(10)]
    for i,j in combinations(range(10),2):
        fused = [x+y for x,y in zip(p[i],p[j])]
        value = trace([q for h,q in enumerate(p) if h not in (i,j)]+[fused])
        for k in range(5):
            potential[i][k] += p[j][k]*value
            potential[j][k] += p[i][k]*value
    square = prod(factors)**2
    for i in range(3):
        potential[i][1] -= 96*sum(factors[3:])*square
        potential[i+3][0] -= 96*sum(factors[:3])*square
    common = [potential[7 if k==2 else 6][k] for k in range(5)]
    for row,momentum in zip(potential,p):
        delta = [x-y for x,y in zip(row,common)]
        for k,h in combinations(range(5),2):
            assert delta[k]*momentum[h] == delta[h]*momentum[k]
    print('100 coordinate polynomials vanish over Z[a1,a2,a3,b1,b2,b3].')
    print('Unit-strength transverse boundary pairing:', 96*3)


if __name__=='__main__':
    verify()
