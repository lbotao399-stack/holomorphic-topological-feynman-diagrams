"""The primary genus-three differential of the elementary scalar."""
from fractions import Fraction as F
from itertools import combinations
from math import lcm, prod
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory
import json
import sys

from flint import fmpq_mat, fmpz_mat, fmpz_mpoly_ctx

ROOT = Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'operator'))
from three_loop_soft_jet import DOUBLE_SOFT, orbit_integers


def elementary_three_loop(scalars, vector, polarization, *, exe=None):
    """Reduced Q3(phi) as the two coefficients of 1 and zeta(3), on rational primary waves."""
    waves = [tuple(map(F,p)) for p in (*scalars,vector)]
    polarization = tuple(map(F,polarization))
    if len(waves)!=8 or any(len(p)!=5 for p in (*waves,polarization)):
        raise ValueError("Expected seven scalar waves, one vector wave and its polarization on the fivefold.")
    if sum(x*y for x,y in zip(waves[-1],polarization)):
        raise ValueError("The vector wave must be transverse.")
    root = tuple(-sum(column) for column in zip(*waves))
    pairing = sum(x*y for x,y in zip(root,polarization))
    if not pairing:
        return F(0),F(0)
    scale = lcm(*(x.denominator for p in waves for x in p))
    points = [tuple(int(x*scale) for x in p) for p in (*waves,root)]
    with TemporaryDirectory() as tmp:
        if exe is None:
            exe = Path(tmp)/'nine'
            run(['c++','-std=c++17','-O3',str(ROOT/'operator'/'nine_sextics.cpp'),
                 '-o',str(exe)],check=True)
        row, = orbit_integers(exe,[points])
    data = json.loads((ROOT/'operator'/'three_loop_coefficients.json').read_text())
    assert tuple(map(tuple,data['basis_words'])) == DOUBLE_SOFT
    return tuple(pairing*sum(x*F(q) for x,q in zip(row,column))/scale**30
                 for column in zip(*data['rational_and_zeta']))


def momenta(strengths):
    p=[tuple(x if i==j else 0 for i in range(5))
       for j,group in enumerate(strengths) for x in group]
    return [*p,tuple(-sum(xs) for xs in zip(*p))]


def coefficients_from_graphs(rows):
    equations,periods=[],[]
    for name in ('three_loop_constraints.json','three_loop_four_pair.json'):
        data=json.loads((ROOT/'operator'/name).read_text())
        equations.extend(data['rows'])
        periods.extend(data['reduced_periods'])
    graph=fmpq_mat(equations).transpose()
    augmented=fmpq_mat([*equations,*rows]).transpose()
    reduced,rank=augmented.rref()
    assert rank==len(equations)==graph.rank()
    weights=fmpq_mat([[reduced[i,len(equations)+j] for i in range(len(equations))]
                      for j in range(len(rows))])
    assert weights*graph.transpose()==fmpq_mat(rows)
    return weights*fmpq_mat(periods)


def strength_polynomials(strengths):
    paired=[p for p in strengths if len(p)==2]
    common=prod(x**6 for p in strengths if len(p)==1 for x in p)
    common*=prod(x*x for p in strengths if len(p)>1 for x in p)
    if any(len(p)>3 for p in strengths):
        return []
    # Symmetrize the residual quadratic on each two-letter line. Three-letter
    # lines have their homogeneous degree six exhausted by the common factor.
    out=[]
    for k in range(len(paired)+1):
        values=0
        for use in combinations(range(len(paired)),k):
            values+=prod(a*b if i in use else (a+b)**2 for i,(a,b) in enumerate(paired))
        out.append(common*values)
    return out


def sextics(edge):
    """Three symmetric degree-six polynomials in the six projected cofactor coordinates."""
    matching=sum(x*y*z for (a,x),(b,y),(c,z) in combinations(edge.items(),3)
                 if len(set(a+b+c))==6)
    disjoint=sum(x**2*y for a,x in edge.items() for b,y in edge.items()
                 if len(set(a+b))==4)
    triangle=sum(edge[i,j]*edge[j,k]*edge[i,k] for i,j,k in combinations(range(6),3))
    return [matching,disjoint,triangle]


def triple_collinear():
    base=[tuple(int(i==j) for j in range(5)) for i in range(5)]
    extra=((1,2,3,4),(1,0,-1,2),(2,-1,3,0),(2,2,1,1),(1,2,0,1))
    graphs,design=[],[]
    for w in extra:
        p=[base[0]]*3+base[1:]+[(0,*w)]
        p.append(tuple(-sum(column) for column in zip(*p)))
        graphs.append(p)
        edge={e:int(fmpz_mat([base[0]]+[q for h,q in enumerate(p[3:]) if h not in e]).det())**2
              for e in combinations(range(6),2)}
        design.append(sextics(edge))
    with TemporaryDirectory() as tmp:
        exe=Path(tmp)/'nine'
        run(['c++','-std=c++17','-O3',str(ROOT/'operator'/'nine_sextics.cpp'),
             '-o',str(exe)],check=True)
        rows=orbit_integers(exe,graphs)
    interpolation=fmpq_mat(design[:3])
    assert fmpq_mat(design)*interpolation.solve(fmpq_mat(rows[:3]))==fmpq_mat(rows)
    coefficients=interpolation.solve(coefficients_from_graphs(rows[:3]))
    return [[str(x) for x in row] for row in coefficients.tolist()]


def chamber_sign(coefficients):
    ctx=fmpz_mpoly_ctx.get(tuple('abcde'))
    gaps=list(ctx.gens())
    # Order the six cofactor coordinates, then split by the two exterior gaps.
    # Reflection of the ordered six coordinates exchanges the two subchambers.
    gaps[0]+=gaps[-1]
    edge={pair:sum(gaps[pair[0]:pair[1]])**2 for pair in combinations(range(6),2)}
    monomials=[poly.to_dict() for poly in sextics(edge)]
    support=set.union(*(set(poly) for poly in monomials))
    lower=sum(F(1,i**3) for i in range(1,5))+F(1,50)
    upper=F(2089,1728)
    channels=[[F(x) for x in row] for row in coefficients]
    for ex in support:
        values=[int(poly.get(ex,0)) for poly in monomials]
        a,b=[sum(x*y for x,y in zip(values,col)) for col in zip(*channels)]
        assert a+b*lower<0 and a+b*upper<0,(ex,a,b)
    return support


def check_mark_changes(strengths):
    # The final wave is the vector on the fifth coordinate line. The test
    # polarization is minus the fourth coordinate vector.
    external=momenta(strengths)[:-1]
    scalars,vector=external[:-1],external[-1]
    assert vector[3]==0 and any(q[3] for q in scalars)
    for i,q in enumerate(scalars):
        if q[3]:
            fused=scalars[:]
            fused[i]=tuple(x+y for x,y in zip(q,vector))
            assert fmpz_mat(fused).rank()<5
    if len(external)==8:
        count=0
        for live in combinations(range(7),4):
            if fmpz_mat([vector,*[scalars[j] for j in live]]).det():
                spectators=[q for j,q in enumerate(scalars) if j not in live]
                assert fmpz_mat(spectators).rank()+2<5
                count+=1
        assert count==6


def certificate():
    families=([(1,1,1),(1,1),(1,),(1,),(1,)],
              [(1,1),(1,1),(1,1),(1,),(1,)])
    probes=[]
    for family in families:
        pairs=[i for i,p in enumerate(family) if len(p)==2]
        for k in range(len(pairs)+1):
            probes.append([(1,-1) if i in pairs[:k] else p for i,p in enumerate(family)])
    checks=[[(2,1,-1),(1,4),(2,),(-1,),(3,)],
            [(1,-2),(3,1),(4,1),(1,),(-2,)],
            [(1,2,3,4),(1,),(1,),(1,),(1,)]]
    with TemporaryDirectory() as tmp:
        exe=Path(tmp)/'nine'
        run(['c++','-std=c++17','-O3',str(ROOT/'operator'/'nine_sextics.cpp'),
             '-o',str(exe)],check=True)
        rows=orbit_integers(exe,[momenta(p) for p in probes+checks])
        # Exercise the general vertex away from the collinear strata, including
        # the extra power of momentum in the root-polarization contraction.
        frame=[tuple(int(i==j) for i in range(5)) for j in range(5)]
        scalars=frame[:4]+[(1,0,0,1,1),(0,1,0,0,-1),(0,0,1,-1,0)]
        polarization=(0,0,0,-1,0)
        general=(F(-597710892167297,43216022513664000),F(36860229443,3334569638400))
        for scale in (F(1),F(1,2)):
            scaled=[tuple(scale*x for x in q) for q in [*scalars,frame[-1]]]
            assert elementary_three_loop(scaled[:-1],scaled[-1],polarization,exe=exe)==tuple(
                scale**31*x for x in general)
    evaluated=coefficients_from_graphs(rows)
    for rational,zeta in evaluated.tolist()[:6]:
        # Integral comparison after four terms gives zeta(3) < 2089/1728.
        assert zeta>0 and 1728*rational+2089*zeta<0
    output=[]
    for start,end,extra in ((0,2,6),(2,6,7)):
        matrix=fmpq_mat([strength_polynomials(p) for p in probes[start:end]])
        target=fmpq_mat(evaluated.tolist()[start:end])
        coefficients=matrix.solve(target)
        # Verify the polynomial of each of the 18 basis traces separately at
        # another source, independently of the integrated values.
        check=fmpq_mat([strength_polynomials(checks[extra-6])])
        assert check*matrix.solve(fmpq_mat(rows[start:end]))==fmpq_mat([rows[extra]])
        assert check*coefficients==fmpq_mat([evaluated.tolist()[extra]])
        output.append([[str(x) for x in r] for r in coefficients.tolist()])
    assert not any(rows[-1])
    return output


if __name__=='__main__':
    calculated=certificate()
    stored=json.loads(Path(__file__).with_suffix('.json').read_text())
    assert calculated==stored['coefficients']
    three=triple_collinear()
    assert three==stored['three_collinear_sextics']
    support=chamber_sign(three)
    assert len(support)==205
    assert {ex.index(6) for ex in support if 6 in ex}=={1,2,3,4}
    check_mark_changes([(1,1),(1,),(1,),(1,),(1,)])
    check_mark_changes([(1,1,1),(1,1),(1,),(1,),(1,)])
    print("Six negative corner periods, 205 negative chamber coefficients, marked tests at two and three loops.")
