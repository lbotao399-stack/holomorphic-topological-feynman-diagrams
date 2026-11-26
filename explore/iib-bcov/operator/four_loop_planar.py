"""Eight planar external waves in the four-loop transverse primary bulk equation."""
from fractions import Fraction
from functools import cache
from itertools import combinations
from math import factorial
from pathlib import Path
import sys

from flint import fmpz_mat, fmpz_mpoly_ctx

from planar_trace import bracket, matchings
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from planar_bulk import response


def restricted_symbols(plane):
    """Symmetric bracket sums on an odd number of planar slots."""
    count = len(plane)
    pairs = [[bracket(p,q) for q in plane] for p in plane]
    full = (1<<count)-1

    @cache
    def match(mask):
        if not mask:
            return 1
        i = (mask&-mask).bit_length()-1
        return sum(pairs[i][j]**2*match(mask^(1<<i)^(1<<j))
                   for j in range(i+1,count) if mask & (1<<j))

    total = tuple(map(sum,zip(*plane)))
    first = second = third = 0
    for i in range(count):
        second += bracket(plane[i],total)**2*match(full^(1<<i))
        for j,k in combinations([h for h in range(count) if h!=i],2):
            rest = match(full^(1<<i)^(1<<j)^(1<<k))
            first += pairs[i][j]**2*pairs[i][k]**2*rest
            third += pairs[j][k]**2*pairs[j][i]*pairs[k][i]*rest
    return first,second,third


def boundary(plane):
    # Read the second covector component at slot e1 relative to any retained slot.
    n = len(plane)//2
    assert len(plane)==2*n and plane[0]==(1,0)
    base = [restricted_symbols(plane[:j]+plane[j+1:]) for j in range(2*n)]
    out = [0,0,0]
    for j in range(1,2*n):
        fused = tuple(x+y for x,y in zip(plane[0],plane[j]))
        weight = restricted_symbols([fused,*plane[1:j],*plane[j+1:]])
        for h in range(3):
            out[h] += plane[j][1]*(weight[h]-base[0][h]-base[j][h])
    total_mass = matchings(plane)
    first = total_mass*sum(q[1] for q in plane)
    second = sum(q[1]*(total_mass-n*n*bracket(plane[0],q)**2*matchings(plane[1:j]+plane[j+1:]))
                 for j,q in enumerate(plane[1:],1))
    return out,first,second


def check_lower_vertices():
    plane = [(1,0,0,0,0)]*4+[(0,1,0,0,0)]*3+[(-4,-3,0,0,0)]
    retained = [(0,0,1,0,0),(0,0,0,1,0),(0,0,0,0,1),(0,0,-1,-1,-1)]
    waves = [*plane,*retained]
    for genus,loops in ((2,2),(3,1)):
        for internal in range(2,4):
            for ids in combinations(range(12),2*genus+3-internal):
                if 8-sum(i<8 for i in ids)<=2*loops:
                    assert fmpz_mat([waves[i] for i in ids]).rank()+internal-1<5
    for ids in combinations(range(12),6):
        other = [i for i in range(12) if i not in ids]
        assert min(fmpz_mat([waves[i] for i in side]).rank() for side in (ids,other))<5


def check_rational_scalar_waves():
    probes = [
        ([(1,0), (1,0), (2,0), (0,2), (0,2), (1,2), (-2,-1)], (168192,516096,-62208)),
        ([(2,0), (2,0), (0,-2), (0,1), (1,1), (2,2), (1,-1)], (244224,1299456,62464)),
    ]
    responses = []
    for plane, expected in probes:
        a, b, d = restricted_symbols(plane)
        assert (a, b, d) == expected and 22*a - 5*b + 18*d == 0
        responses.append((Fraction(b,576), Fraction(8*a-b,9)))
    assert responses == [(896,92160), (2256,72704)]
    (k1,j1), (k2,j2) = responses
    assert k1*j2 - k2*j1 == -142770176
    assert (141*k1 - 56*k2, 141*j1 - 56*j2) == (0,8923136)


def check():
    check_rational_scalar_waves()
    ctx = fmpz_mpoly_ctx.get(tuple('x'+str(i) for i in range(12)))
    x = ctx.gens()
    covectors = [(1,0),(0,1),*zip(x[::2],x[1::2])]
    (a,b,c),f,g = boundary(covectors)
    assert 32*a==f+3*g and 4*b==f and 22*a-5*b+18*c==0
    assert len(f)==987 and len(g)==907
    # Two transverse tests separate the two responses. Three seven-wave tests
    # distinguish the three scalar restrictions before taking their boundary.
    fixtures = [
        [(1,0),(0,1),(1,1),(0,2),(2,1),(-1,3),(2,2),(1,3)],
        [(1,0),(0,1),(1,2),(1,-1),(-1,1),(2,-2),(3,2),(2,-1)],
        [(1,0),(0,1),(2,1),(-1,2),(0,3),(2,1),(-1,-2),(-1,1)]]
    assert fmpz_mat([restricted_symbols(p[:-1]) for p in fixtures]).det()
    assert fmpz_mat([boundary(p)[1:] for p in fixtures[:2]]).det()
    for n in range(2,9):
        p = [(1,0)]*n+[(0,1)]*(n-1)+[(-n,1-n)]
        for k in (Fraction(0),Fraction(7,3)):
            assert response(tuple(p),k,Fraction(1))[n][0]==-factorial(n)**2*(n-1)
    check_lower_vertices()
    print('Quartic planar boundary has rank two, with one kernel among three scalar restrictions.')
    print('The four-loop transverse projection is -1728 J. Lower-genus insertions vanish.')


if __name__=='__main__':
    check()
