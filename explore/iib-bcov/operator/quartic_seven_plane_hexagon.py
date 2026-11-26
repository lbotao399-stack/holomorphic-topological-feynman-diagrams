"""Complete genus-three vertex insertion into the four-loop hexagon on seven planar waves."""
from fractions import Fraction as F
from functools import cache
from itertools import combinations, permutations, product
from math import factorial, prod
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory
import json
import sys

from flint import fmpq_mat, fmpz_mat
from four_loop_planar import restricted_symbols
from three_loop_soft_jet import DOUBLE_SOFT, orbit_integers
from fusion_transfer_box import shift
from planar_trace import matchings
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from marked_theta import green_moments
from bulk_eight import BULK
from marked_boundaries import quartics


def det(p):
    if len(p)==2:
        return p[0][0]*p[1][1]-p[0][1]*p[1][0]
    a,b,c=p
    return a[0]*det([b[1:],c[1:]])-a[1]*det([(b[0],b[2]),(c[0],c[2])])+a[2]*det([b[:2],c[:2]])


def actual(exe,points):
    data=json.loads(Path(__file__).with_name('three_loop_coefficients.json').read_text())
    assert tuple(map(tuple,data['basis_words']))==DOUBLE_SOFT
    coefficients=[list(map(F,x)) for x in data['rational_and_zeta']]
    return [tuple(sum(q*c[k] for q,c in zip(row,coefficients)) for k in range(2))
            for row in orbit_integers(exe,points)]


def lifted(plane,retained):
    rows=[tuple(p)+(0,0,0) for p in plane]+list(retained)
    return rows+[tuple(-sum(c) for c in zip(*rows))]


def restrict(exe):
    u=(1,0);v=(0,1)
    five=[[u,u,v,v,(1,1)],[u,u,u,v,v],[u,v,(1,1),(1,-1),(2,1)],[u,v,(2,-1),(1,2),(-3,1)]]
    axes=[(0,0,1,0,0),(0,0,0,1,0),(0,0,0,0,1)]
    mixed=[(2,1,1,2,0),(-1,2,1,0,2),(1,0,0,1,2)]
    probes=[lifted(s,axes) for s in five]+[lifted(five[-1],mixed)]
    values=actual(exe,probes)
    m=fmpq_mat([list(restricted_symbols(s)) for s in five[:3]])
    assert m.det()
    out=[]
    for col in range(2):
        b=fmpq_mat([[str(r[col])] for r in values[:3]])
        x=m.inv()*b
        out.append(tuple(F(str(x[k,0])) for k in range(3)))
        check=sum(a*c for a,c in zip(out[-1],restricted_symbols(five[-1])))
        assert check==values[3][col]
        assert check*det([p[2:] for p in mixed])**6==values[4][col]
    return tuple(out)


@cache
def green6(c,i,j):
    # Integrate (sum c_i s_i)^6 against the grounded circle Green function.
    value=F()
    for a in range(7):
        for b in range(7-a):
            k=(a,b,6-a-b)
            coeff=factorial(6)//(factorial(k[0])*factorial(k[1])*factorial(k[2]))
            value+=coeff*c[0]**k[0]*c[1]**k[1]*c[2]**k[2]*green_moments(k,i,j)[1]
    return value


def chamber_green(k,i,j):
    """Independently integrate the Green kernel on the six ordered simplices."""
    answer = F()
    for order in permutations(range(3)):
        def monomial(extra):
            powers = [k[h] + extra.count(h) for h in order]
            return F(1, prod(n + sum(powers[:n]) for n in range(1,4)))
        earlier = min((i,j), key=order.index)
        answer += monomial((earlier,)) - monomial((i,j))
    return answer


def planar_factors(plane,coefficients):
    ans=[F(),F()]
    for i,j in combinations(range(7),2):
        d=det([plane[i],plane[j]])
        if d:
            a=restricted_symbols([p for k,p in enumerate(plane) if k not in (i,j)])
            for ch,coeff in enumerate(coefficients):
                ans[ch]+=d*d*sum(x*y for x,y in zip(a,coeff))
    return tuple(ans)


def theta_planar_factor(plane):
    axes = [(0,0,1,0,0),(0,0,0,1,0),(0,0,0,0,1)]
    total = F()
    for stored in combinations(range(7),3):
        quantum = lifted([plane[i] for i in stored], axes)
        symbol = sum(c*v for c,v in zip(BULK,quartics(quantum)))
        ordinary = [p for i,p in enumerate(plane) if i not in stored]
        # matchings includes 2! for the two unordered pairs.
        total -= symbol*matchings(ordinary)
    return total


def check_planar_placements():
    plane = [(2,1),(1,-3),(0,1),(3,2)]
    paths = [(1,0),(1,1),(0,1)]
    lengths = (3,5,7)
    value = 0
    for assigned in product(range(3),repeat=4):
        rows = [[x*y for x in paths[h] for y in point] for point,h in zip(plane,assigned)]
        value += prod(lengths[h] for h in assigned)*int(fmpz_mat(rows).det())**2
    tree = sum(x*y for x,y in combinations(lengths,2))
    assert value == tree**2*matchings(plane)


def check_paired_transfers(coefficients,difference):
    # The final normal wave absorbs the negative planar sum. Its normal projection is fixed.
    hexagon = {}
    theta = {}
    for t in (-3,-1,0,1,2,3):
        plane = [(1,0)]*3+[(0,1)]*3+[(t,t)]
        theta[t] = theta_planar_factor(plane)
        hexagon[t] = planar_factors(plane,coefficients)
        assert theta[t] == F(t*t*(18+4*t+5*t*t),4320)
        assert hexagon[t] == (
            -F(t*t*(9133295639+2195940786*t+2961866876*t*t),23538138624000),
            F(t*t*(2890121+695214*t+938306*t*t),9081072000),
        )
    assert 27*theta[-1] == 19*theta[1]
    required = tuple(difference*(19*y-27*x) for x,y in zip(hexagon[-1],hexagon[1]))
    assert required == (F(9538548955331,37149840000),-F(3033930713,14332500))
    print('bare paired transfer required by local trace',required)


def normal_component(normal,i,j,direction=(0,0,1)):
    # Difference of two nonplanar potentials in the specified normal direction.
    result=F()
    for group in combinations(range(5),2):
        singles=[k for k in range(5) if k not in group]
        tr=[normal[k] for k in singles]
        d=det(tr)
        c=tuple(det((normal[group[0]],normal[group[1]],v)) for v in tr)
        for a,k in enumerate(singles):
            sign=(k==i)-(k==j)
            if sign and d:
                result-=F(sign*d*d,2)*sum(sum(x*y for x,y in zip(q,direction))*green6(c,a,b) for b,q in enumerate(tr))
    return result


def box_normal(normal):
    axes=(1,0,0)
    edges=[(i,j) for i in (0,1) for j in (2,3,4)]
    total=F()
    for mask in range(64):
        p=normal
        for k,(i,j) in enumerate(edges):
            if mask>>k&1:
                p=shift(p,i,j,axes)
        total+=(-1)**mask.bit_count()*normal_component(p,0,1)
    return total


def check():
    for a in range(7):
        for b in range(7-a):
            k = (a,b,6-a-b)
            for i in range(3):
                for j in range(3):
                    assert chamber_green(k,i,j) == green_moments(k,i,j)[1]
    with TemporaryDirectory() as tmp:
        exe=Path(tmp)/'nine'
        run(['c++','-std=c++17','-O3',str(Path(__file__).with_name('nine_sextics.cpp')),'-o',str(exe)],check=True)
        coeff=restrict(exe)
    assert coeff == (
        (F(-841479311,211843247616000),F(2104429,8826801984000),F(-70595227,31384184832000)),
        (F(53359,16345929600),F(-21479,108972864000),F(40357,21794572800)),
    )
    plane=((1,0),)*3+((0,1),)*3+((-3,-3),)
    normal=((1,2,0),(3,-1,0),(-2,1,1),(0,2,2),(-2,-4,-3))
    factors=planar_factors(plane,coeff)
    value=normal_component(normal,0,1)
    diff=box_normal(normal)
    assert factors == (F(-834350719,74724249600),F(1321319,144144000))
    assert (value,diff) == (F(653198,15),F(-49747896,35))
    assert tuple(x*value for x in factors) == (F(-2994484730491,6158592000),F(4742213891,11880000))
    assert tuple(x*diff for x in factors) == (F(52408071712547,3302208000),F(-82996010347,6370000))
    check_planar_placements()
    check_paired_transfers(coeff,diff)
    print('C3 planar coefficients',coeff,'planar factors',factors)
    print('orthogonal value, sixth difference',value,diff)
    print('one-loop insertion',tuple(x*value for x in factors),'sixth difference',tuple(x*diff for x in factors))


if __name__=='__main__':
    check()
