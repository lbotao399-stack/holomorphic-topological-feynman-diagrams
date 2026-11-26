"""Two nine-slot trace primitives for six arbitrary planar external covectors."""
from collections import Counter
from itertools import combinations
from flint import fmpz_mpoly_ctx


def bracket(p,q):
    return p[0]*q[1]-p[1]*q[0]


def matchings(p):
    if not p:
        return 1
    return len(p)//2*sum(bracket(p[0],p[j])**2*matchings(p[1:j]+p[j+1:])
                         for j in range(1,len(p)))


def incidence():
    # A nonzero five-row minor uses three or four of the four retained rows.
    minor = []
    pairs = list(combinations(range(5),2))
    for rows in combinations(range(9),5):
        plane = tuple(i for i in rows if i < 5)
        if len(plane) in (1,2):
            index = 10+plane[0] if len(plane)==1 else pairs.index(plane)
            minor.append((sum(1<<i for i in rows),index))
    out = [Counter(), Counter()]
    for (a,i),(b,j),(c,k) in combinations(minor,3):
        if a|b|c != 511:
            continue
        counts = sorted(((a&b).bit_count(),(b&c).bit_count(),(a&c).bit_count()))
        if counts in ([2,2,3],[2,3,3]):
            out[counts==[2,3,3]][tuple(sorted((i,j,k)))] += 1
    return out


def trace(rows,coefficients):
    total = [sum(x) for x in zip(*rows)]
    minors = [bracket(p,q)**2 for p,q in combinations(rows,2)]
    minors.extend(bracket(q,total)**2 for q in rows)
    return [sum(weight*minors[i]*minors[j]*minors[k] for (i,j,k),weight in data.items())
            for data in coefficients]


def verify():
    ctx = fmpz_mpoly_ctx.get(tuple(f'x{i}' for i in range(8)))
    x = ctx.gens()
    p = [(1,0),(0,1)]+[(x[i],x[i+1]) for i in range(0,8,2)]
    # GL(2) normalizes the first two independent rows. S_6 then reduces all
    # longitudinality tests to the second component at the first planar slot.
    coeff = incidence()
    base = [trace(p[:j]+p[j+1:],coeff) for j in range(6)]
    raw = [0,0]
    for j in range(6):
        if j:
            fused = tuple(a+b for a,b in zip(p[0],p[j]))
            weights = trace(p[1:j]+p[j+1:]+[fused],coeff)
            for k in (0,1):
                raw[k] += p[j][1]*weights[k]
        for k in (0,1):
            raw[k] -= p[j][1]*(base[0][k]+base[j][k])
    z = matchings(p)
    first = z*sum(q[1] for q in p)
    # Six times the second transfer covector cancels the n(n-1) denominator.
    second = sum(q[1]*(z-9*bracket(p[0],q)**2*matchings(p[1:j]+p[j+1:]))
                 for j,q in enumerate(p) if j)
    assert raw[1]-raw[0] == 8*first
    assert 4*raw[0]-raw[1] == 28*second
    print('Both trace boundaries verified as polynomials in eight free planar coordinates.')


if __name__=='__main__':
    verify()
