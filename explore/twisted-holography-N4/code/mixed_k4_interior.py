"""All six active edge constraints of the K4 Newton boundary, at zero branch powers.

Let U be the Kirchhoff polynomial and write N_e=|U j_e|^2 for the squared
numerator of the minimizing current. With l_root=1, the signed integrand is

 N_root det d(U,N_others)/d(l_root,l_others) / (3 * 2^6 * U^13).

Replacing the first Jacobian column by its Euler sum gives (3U,6N_others);
elimination of that column differentiates N_e/U^2. The numerator has degree
33 in the six lengths. Each accepted exponent vector m has a_i=1+m_i with
sum a_i=39, 0<a_i<13 and every triangle sum above 13. These conditions
are tested by membership in the full table of convergent tetrahedron periods.
Degree-six interpolation in the charges' dot products proves the same support
bound for every input. The sign includes the divergence factor (-1)^6.
"""
from collections import defaultdict
from fractions import Fraction as F
from functools import cache
from itertools import combinations, combinations_with_replacement
import gzip
import json
from math import factorial, lcm, prod
from pathlib import Path
import sys

from flint import fmpq, fmpq_mat, fmpz_mpoly_ctx

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'explore/nonSP/verify/k4_periods'))
from core import EDGES, PERMS

CTX=fmpz_mpoly_ctx.get(tuple(f'l{i}' for i in range(6)))


@cache
def trees():
    """Tree supply currents for three unit sources grounded at vertex 3."""
    lengths=CTX.gens()
    u=CTX.constant(0)
    flow=[[CTX.constant(0) for _ in range(3)] for _ in EDGES]
    for tree in combinations(range(6),3):
        unseen={0,1,2}
        toward={}
        reached={3}
        while unseen:
            progress=False
            for e in tree:
                a,b=EDGES[e]
                if (a in reached)==(b in reached):
                    continue
                far,near=(a,b) if a in unseen else (b,a)
                toward[far]=(near,e,1 if far==a else -1)
                reached.add(far)
                unseen.remove(far)
                progress=True
            if not progress:
                break
        if unseen:
            continue
        term=prod(lengths[i] for i in range(6) if i not in tree)
        u+=term
        for s in range(3):
            v=s
            while v!=3:
                v,e,sign=toward[v]
                flow[e][s]+=sign*term
    assert len(u)==16
    return u,flow


def det3(rows,cols):
    x,y,z=cols
    a,b,c=rows
    return a[x]*(b[y]*c[z]-b[z]*c[y])-a[y]*(b[x]*c[z]-b[z]*c[x])+a[z]*(b[x]*c[y]-b[y]*c[x])


@cache
def period_table():
    bank=ROOT/'explore/iib-bcov/operator/quartic_seven_plane_tetra.weights.gz'
    orbits={}
    values={}
    with gzip.open(bank,'rt') as f:
        for line in f:
            ex,a,b=json.loads(line)
            key=tuple(ex)
            factor=F(prod(factorial(i-1) for i in key),factorial(12))
            values[key]=(factor*F(a),factor*F(b))
            for p in PERMS:
                orbits[tuple(key[i] for i in p)]=key
    assert len(values)==3872
    return orbits,values


def numerator(gram,root):
    u,flow=trees()
    n=[sum((v[i]*v[j]*gram[i][j] for i in range(3) for j in range(3) if gram[i][j]),
           CTX.constant(0)) for v in flow]
    if not n[root]:
        return CTX.constant(0)
    # Relabel both the differentiated rows and columns so that the root is first.
    ids=(root,*(i for i in range(6) if i!=root))
    rows=[[u.derivative(i) for i in ids]]
    rows.extend([[n[e].derivative(i) for i in ids] for e in ids[1:]])
    det=CTX.constant(0)
    for chosen in combinations(range(6),3):
        complement=tuple(i for i in range(6) if i not in chosen)
        det+=(-1)**(sum(chosen)-3)*det3(rows[:3],chosen)*det3(rows[3:],complement)
    return n[root]*det


def gram_interior(gram,root=0,check=False):
    """Return rational and zeta(3) coefficients for the Gram matrix of charges 0,1,2."""
    if not 0<=root<6 or len(gram)!=3 or any(len(row)!=3 for row in gram):
        raise ValueError('Use a 3-by-3 charge Gram matrix and a root edge from 0 to 5.')
    gram=tuple(tuple(map(F,row)) for row in gram)
    if any(gram[i][j]!=gram[j][i] for i in range(3) for j in range(i)):
        raise ValueError('The charge Gram matrix must be symmetric.')
    scale=lcm(*(v.denominator for row in gram for v in row))
    integers=[[int(scale*a) for a in row] for row in gram]
    poly=numerator(integers,root)
    denominator=3*2**6*scale**6
    if check:
        u=trees()[0]
        for lengths in ((1,1,1,1,1,1),(2,3,5,7,11,13)):
            # For a general root length, the affine determinant has an extra factor l_root.
            measured=F(int(poly(*lengths))*lengths[root],denominator*int(u(*lengths))**13)
            assert measured==grounded(gram,root,lengths),(measured,grounded(gram,root,lengths))
    orbits,periods=period_table()
    totals=defaultdict(int)
    for ex,c in poly.terms():
        powers=tuple(int(i)+1 for i in ex)
        assert sum(powers)==39
        if powers not in orbits:
            raise ArithmeticError(f'Nonconvergent tetrahedral monomial {powers}')
        totals[orbits[powers]]+=int(c)
    return tuple(sum((c*periods[p][j] for p,c in totals.items()),F())/denominator for j in range(2))


# Monomial coefficients for x_ij=−P_i·P_j, in EDGES order. The rational
# column uses the stabilizer of edge 01; the zeta(3) column uses S_4.
# The six digits are exponents. Every distinct monomial in an orbit occurs once.
RATIONAL_ROWS = """
000024:-10395 000033:-116655 000042:-176715 000114:-20790 000123:-370755
000132:-911295 000141:-353430 000222:-1479555 000231:-1119195 000240:-176715
000330:-324555 001014:-117810 001023:4648208 001032:3640124 001041:-436590
001104:-117810 001113:9180916 001122:21186988 001131:6619588 001140:-436590
001203:4650518 001212:21184678 001221:19363474 001230:3416054 001302:3635504
001311:5824948 001320:3048764 001401:-859320 001410:-859320 002022:10843650
002031:3877592 002040:-98175 002112:26007879 002121:23521395 002130:3903002
002202:10449031 002211:25691328 002220:8774743 002301:6219229 002310:5701789
002400:-757680 003030:-125895 003111:14916840 003120:3651212 003210:5992849
003300:-1783960 004110:-734580 011112:61342204 011121:67648138 011130:14473320
011220:26245249 012120:26669365 012210:32301214 100014:-143220 100023:4223168
100032:3503834 100041:-436590 100113:8303116 100122:20103598 100131:6277708
100140:-436590 100221:18797524 100230:3210464 101013:17836680 101022:22104903
101031:6678724 101040:-196350 101103:17827440 101112:68496370 101121:48618724
101130:6627904 101202:27891817 101211:61682794 101220:19615647 101301:13614000
101310:13249020 101400:-734580 102021:17832868 102030:-781935 102111:57319666
102120:17608798 102201:29487634 102210:19071792 102300:4639189 103110:5292724
111111:145516320 111120:54691810 200004:-133980 200013:6736669 200022:7715839
200031:3237722 200040:-98175 200112:22730832 200121:18346071 200130:3161492
200220:7323210 201012:19204848 201021:15055786 201030:-989835 201102:28888882
201111:45893482 201120:15016516 201201:18819061 201210:13278855 201300:2376092
202020:-2020095 202110:14265766 202200:4120555 211110:29708140 300003:-2213620
300012:3571045 300021:1740380 300030:-333795 300111:8402640 300120:1800440
301011:2099380 301020:-1474935 301101:8153160 301110:2797000 301200:1835090
400002:-1131900 400011:-1381380 400020:-236775 400110:-838530 401010:-473550
401100:-838530
"""
ZETA_ROWS = """
001023:-163 001113:-326 001122:-815 001203:-163 001212:-815
001302:-163 001311:-326 002022:-372 002112:-915 002202:-380
002211:-1038 002301:-278 003111:-636 003300:32 011112:-2156
011121:-2669 011130:-636 011220:-1118 012210:-1364 111111:-6708
"""


@cache
def polynomial():
    columns=[]
    for text,group,d in ((RATIONAL_ROWS,[p for p in PERMS if p[0]==0],50295168000),
                         (ZETA_ROWS,PERMS,1940400)):
        result={}
        for entry in text.split():
            label,v=entry.split(':')
            ex=tuple(map(int,label))
            assert len(ex)==6 and sum(ex)==6 and ex==min(tuple(ex[k] for k in p) for p in group)
            for key in {tuple(ex[k] for k in p) for p in group}:
                assert key not in result
                result[key]=F(int(v),d)
        columns.append(result)
    assert tuple(map(len,columns))==(426,296)
    return columns


def edge_value(weights,root=0):
    """Evaluate both invariant polynomials over the coefficient field of the weights."""
    if len(weights)!=6 or not 0<=root<6:
        raise ValueError('Supply six ordered edge weights and a root in range(6).')
    permutation=next(p for p in PERMS if p[0]==root)
    axes=[weights[j] for j in permutation]
    return tuple(sum(c*prod(x**j for x,j in zip(axes,k) if j) for k,c in col.items())
                 for col in polynomial())


def interior(charges,root=0,check=False):
    """Signed all-active contribution at k=0; edges 01,02,03,12,13,23.

    Returns the exact rational and zeta(3) columns. Polynomial evaluation
    applies to arbitrary Euclidean source vectors and extends bilinearly.
    The optional independent graph integral accepts rational Gram entries.
    """
    if len(charges)!=4 or len({len(v) for v in charges})!=1 or any(map(sum,zip(*charges))):
        raise ValueError('Supply four conserved vectors of equal dimension.')
    weights=tuple(-sum(a*b for a,b in zip(charges[i],charges[j])) for i,j in EDGES)
    val=edge_value(weights,root)
    if check:
        gram=tuple(tuple(sum(F(a)*F(b) for a,b in zip(x,y)) for y in charges[:3]) for x in charges[:3])
        assert val==gram_interior(gram,root,True)
    return val


def rebuild():
    # Degree six in six dot products. The simplex cardinal functions for
    # |a|=6 are product_i product_{s<a_i}(x_i-s*sum(x)/6)/a_i!.
    # Hence agreement on these integral points proves equality over Q.
    group=[p for p in PERMS if p[0]==0]
    rows={min(tuple(r.count(i) for i in p) for p in group)
          for r in combinations_with_replacement(range(6),6)}
    for weights in sorted(rows):
        gram=[[0]*3 for _ in range(3)]
        for (i,j),c in zip(EDGES,weights):
            gram[i][i]+=c
            if j!=3:
                gram[j][j]+=c
                gram[i][j]-=c
                gram[j][i]-=c
        assert edge_value(weights)==gram_interior(gram)
    assert len(rows)==138
    print('Certified all 462 degree-six interpolation nodes on 138 rooted edge orbits.')


def grounded(gram,root,lengths):
    """Independent differential of grounded potentials with respect to resistances."""
    t=fmpq_mat([[int(i==j)*fmpq(1,lengths[i]) for j in range(6)] for i in range(6)])
    inc=fmpq_mat([[int(v==a)-int(v==b) for a,b in EDGES] for v in range(3)])
    g=fmpq_mat([[fmpq(str(x)) for x in row] for row in gram])
    inverse=(inc*t*inc.transpose()).inv()
    current=t*inc.transpose()*inverse
    dot=current*g*current.transpose()
    projection=t-current*inc*t
    rest=[i for i in range(6) if i!=root]
    jac=fmpq_mat([[-projection[i,j]*dot[i,j] for j in rest] for i in rest])
    return F(str(dot[root,root]*jac.det()/2))


if __name__=='__main__':
    from argparse import ArgumentParser
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--rebuild',action='store_true',help='Certify the full charge polynomial from the tetrahedron periods')
    if parser.parse_args().rebuild:
        rebuild()
    unit=((1,0,0),(0,1,0),(0,0,1),(-1,-1,-1))
    rank_two=((1,0),(0,1),(1,2),(-2,-3))
    for root in range(6):
        expected=F(6646939 if 3 in EDGES[root] else 11084911,16765056000),-F(3,4312)
        assert interior(unit,root,True)==expected
        assert interior(rank_two,root,True)[1]==F(2,1617)
        assert interior(((1,),(2,),(3,),(-6,)),root,True)==(0,0)
    print('Root 01, orthogonal charges:',*interior(unit))
    print('Root 01, two-dimensional charges:',*interior(rank_two))
