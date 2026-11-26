"""First soft derivative and manifestly doubly soft nine-primary trace symbols."""
from itertools import combinations
from math import factorial,isqrt
from pathlib import Path
from random import Random
from subprocess import run
from tempfile import TemporaryDirectory

from genus_two_probes import det
from quartic_unisolvence import FIXTURES,pivot_rows


SUBSETS=list(combinations(range(9),5))


DOUBLE_SOFT = ((59, 59, 81, 81, 113, 113),
 (5, 5, 15, 96, 110, 119),
 (21, 21, 52, 62, 110, 110),
 (10, 37, 75, 79, 123, 123),
 (9, 44, 47, 75, 78, 105),
 (8, 47, 68, 78, 88, 91),
 (45, 59, 59, 60, 81, 96),
 (2, 28, 58, 77, 92, 92),
 (0, 0, 51, 51, 112, 112),
 (13, 27, 66, 77, 92, 115),
 (0, 3, 16, 68, 121, 121),
 (11, 19, 30, 30, 68, 110),
 (25, 41, 55, 74, 119, 119),
 (11, 14, 35, 48, 85, 114),
 (3, 25, 25, 36, 63, 112),
 (49, 49, 53, 53, 70, 94),
 (50, 50, 77, 77, 121, 121),
 (23, 23, 50, 54, 122, 122))

RESTRICTION_WORDS = ((23, 23, 23, 23, 23, 23),
 (2, 2, 2, 2, 62, 62),
 (3, 3, 16, 65, 65, 82),
 (31, 37, 37, 41, 113, 115),
 (0, 8, 8, 28, 102, 102),
 (49, 49, 74, 74, 74, 74),
 (0, 12, 12, 12, 74, 74),
 (2, 31, 32, 82, 82, 109),
 (7, 30, 53, 71, 82, 125),
 (60, 60, 72, 78, 78, 117),
 (8, 8, 79, 98, 102, 102),
 (116, 116, 116, 116, 121, 121),
 (1, 51, 60, 61, 64, 98),
 (5, 17, 38, 38, 40, 90),
 (62, 62, 62, 97, 97, 97),
 (2, 2, 22, 22, 56, 124),
 (37, 37, 60, 106, 106, 106),
 (0, 37, 46, 56, 92, 116),
 (15, 65, 65, 105, 105, 122),
 (55, 68, 84, 92, 95, 121),
 (27, 27, 36, 39, 39, 41),
 (2, 2, 2, 2, 10, 34),
 (78, 94, 95, 107, 107, 108),
 (0, 0, 6, 26, 47, 105),
 (29, 29, 35, 35, 38, 38),
 (5, 5, 97, 97, 121, 121),
 (1, 19, 36, 83, 96, 111),
 (2, 2, 76, 76, 84, 84),
 (85, 85, 100, 100, 104, 104),
 (53, 59, 59, 61, 61, 104),
 (44, 45, 48, 54, 69, 112),
 (28, 28, 29, 60, 94, 95))


def native(exe,ps,words,linear):
    arrays=[(len(words),len(ps),int(linear)),*words]
    for p in ps:
        arrays.append([det([p[i] for i in sub])%1009 for sub in SUBSETS])
    stdin='\n'.join(' '.join(map(str,row)) for row in arrays)
    result=run([exe],input=stdin,text=True,capture_output=True,check=True)
    return [list(map(int,line.split())) for line in result.stdout.splitlines()]


def orbit_integers(exe,ps,words=DOUBLE_SOFT,leading_soft=0):
    """For k=leading_soft retain degree two in each of the final k slots, 0<=k<=3."""
    if leading_soft not in (0,1,2,3):
        raise ValueError("At most three of nine momenta can be soft at total rank five.")
    # CRT with a deterministic coefficient bound allows arbitrary integral momentum rows.
    minors=[[det([p[h] for h in s]) for s in SUBSETS] for p in ps]
    bound=factorial(9)*max((abs(v) for m in minors for v in m),default=0)**6
    out=[[0]*len(words) for p in ps]
    product=1
    if not bound:
        return out
    for prime in range(999983,2,-2):
        if any(prime%d==0 for d in range(3,isqrt(prime)+1,2)):
            continue
        records=[(len(words),len(ps),(3,6,7,9)[leading_soft],prime),*words,*([v%prime for v in row] for row in minors)]
        inp='\n'.join(' '.join(map(str,r)) for r in records)
        result=run([exe],input=inp,text=True,capture_output=True,check=True)
        inverse=pow(product,-1,prime)
        for row,values in zip(out,result.stdout.splitlines()):
            for h,residue in enumerate(map(int,values.split())):
                row[h] += product*((residue-row[h])*inverse%prime)
        product*=prime
        if product>2*bound:
            return [[x if x<=product//2 else x-product for x in row] for row in out]
    raise ArithmeticError('Coefficient bound exceeds the product of the available CRT primes.')


def exact_boundary(exe,p):
    """Eighteen 10-by-5 covector arrays, in DOUBLE_SOFT order, for any conserved integral fixture."""
    if len(p)!=10 or any(len(q)!=5 for q in p) or any(map(sum,zip(*p))):
        raise ValueError('Ten conserved five-component momenta are required.')
    pairs=list(combinations(range(10),2))
    q=[[v for h,v in enumerate(p) if h not in (i,j)]
       +[tuple(a+b for a,b in zip(p[i],p[j]))] for i,j in pairs]
    vals=orbit_integers(exe,q)
    arrays=[]
    for col in zip(*vals):
        a=[[0]*5 for _ in p]
        for (i,j),v in zip(pairs,col):
            for k in range(5):
                a[i][k]+=v*p[j][k]
                a[j][k]+=v*p[i][k]
        arrays.append(a)
    return arrays



def soft_boundary(exe,p,spurion,words):
    pairs=list(combinations(range(9),2))
    arrays=[(len(words),len(pairs),2),*words]
    for i,j in pairs:
        q=[v for h,v in enumerate(p) if h not in (i,j)]
        q.extend([tuple(x+y for x,y in zip(p[i],p[j])),(0,)*5])
        z=q[:]
        z[-1]=spurion
        z[-2]=tuple(x-y for x,y in zip(q[-2],spurion))
        base=[det([q[h] for h in s]) for s in SUBSETS]
        unit=[det([z[h] for h in s]) for s in SUBSETS]
        arrays.extend([base,[(b-a)%1009 for a,b in zip(base,unit)]])
    inp='\n'.join(' '.join(map(str,r)) for r in arrays)
    result=run([exe],input=inp,text=True,capture_output=True,check=True)
    values=[list(map(int,line.split())) for line in result.stdout.splitlines()]
    return physical_transfers(p,pairs,values)


def physical_transfers(p,pairs,values):
    cols=[]
    for x in zip(*values):
        a=[[0]*5 for pp in p]
        for (i,j),coef in zip(pairs,x):
            for h in range(5):
                a[i][h]+=coef*p[j][h]
                a[j][h]+=coef*p[i][h]
        col=[]
        for i,j in combinations(range(len(p)-1),2):
            for h in range(5):
                if p[i][h]:continue
                for k in range(5):
                    if not p[j][k]:
                        col.append(((a[i][h]-a[-1][h])*p[-1][k]
                                    -(a[j][k]-a[-1][k])*p[-1][h])%1009)
        cols.append(col)
    return list(zip(*cols))


def constant_points():
    rng=Random(43286)
    points=[]
    for _ in range(30):
        p=[tuple(int(i==k) for i in range(5)) for k in range(5)]
        p.extend(tuple(rng.randrange(-3,4) for j in range(5)) for i in range(2))
        p.append(tuple(-sum(z) for z in zip(*p)))
        # Advance the deterministic sequence past the companion soft covector.
        tuple(rng.randrange(-3,4) for i in range(5))
        p.append((0,)*5)
        points.append(p)
    return points


def random_points(count):
    rng=Random(76455)
    out=[]
    for j in range(count):
        p=[tuple(int(i==k) for i in range(5)) for k in range(5)]
        p += [tuple(rng.randrange(-3,4) for i in range(5)) for k in range(3)]
        p.append(tuple(-sum(z) for z in zip(*p)))
        out.append(p)
    return out


def choose(rows):
    ids=pivot_rows(rows)
    m=[rows[j] for j in ids]
    cols=pivot_rows(zip(*m))
    value=det([[r[k] for k in cols] for r in m])%1009
    return ids,cols,value


def check(exe):
    for w in DOUBLE_SOFT:
        assert all(sum(i in SUBSETS[j] for j in w)>=2 for i in range(9))
    mat=native(exe,random_points(18),DOUBLE_SOFT,False)
    assert det(mat)%1009==840
    rows=native(exe,constant_points(),RESTRICTION_WORDS,False)
    rows+=soft_boundary(exe,FIXTURES[0],(1,-1,2,-2,3),RESTRICTION_WORDS)
    ids=tuple(range(16))+(30,31,32,33,34,38,42,46,47,48,49,62,63,64,65,78)
    assert det([rows[j] for j in ids])%1009==868
    from quintic_unisolvence import FIXTURES as ten_wave
    p=ten_wave[0]
    pairs=list(combinations(range(10),2))
    q=[[v for h,v in enumerate(p) if h not in (i,j)]
       +[tuple(a+b for a,b in zip(p[i],p[j]))] for i,j in pairs]
    exact=orbit_integers(exe,q)
    mod=native(exe,q,DOUBLE_SOFT,False)
    assert [[z%1009 for z in row] for row in exact]==mod
    rows=physical_transfers(p,pairs,exact)
    ids=(0,1,2,3,4,8,12,16,17,18,19,32,33,34,35,48,49,50)
    assert det([rows[j] for j in ids])%1009==904
    planar=(((1,0),(1,0),(1,0),(0,1),(0,1),(0,1)),
            ((1,0),(0,1),(1,1),(1,-1),(1,0),(0,1)))
    for u,mult,coeffs in zip(planar,(576,192),((-4,-2,4,1),(-12,-16,12,5))):
        p=[(*v,0,0,0) for v in u]+[(0,0,1,0,0),(0,0,0,1,0),(0,0,0,0,1)]
        p.append(tuple(-sum(z) for z in zip(*p)))
        expected=[0]*18
        for j,x in zip((2,7,8,17),coeffs):
            expected[j]=mult*x
        assert [a[0][1]-a[6][1] for a in exact_boundary(exe,p)]==expected
    print('The two exact planar coefficient rows passed.')
    print('The ten-wave physical boundary is injective; minor 904 modulo 1009.')
    print('18 double-soft symbols; independent evaluation minor 840 modulo 1009.')
    print('32 constant and linear constraints; independent minor 868 modulo 1009.')


if __name__=='__main__':
    with TemporaryDirectory() as tmp:
        exe=Path(tmp)/'nine-sextic'
        run(['c++','-O3','-std=c++17',str(Path(__file__).with_name('nine_sextics.cpp')),
             '-o',str(exe)],check=True)
        check(exe)
