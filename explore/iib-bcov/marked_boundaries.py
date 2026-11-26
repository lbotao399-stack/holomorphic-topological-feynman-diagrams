"""Test genus-two primary scalar counterterms against two integrated-cubic projections."""
from fractions import Fraction as F
from itertools import combinations
from flint import fmpz_mat
from integrated_cubic import run
from marked_theta import MOMENTA


def quartics(q):
    assert len(q)==7 and not any(map(sum,zip(*q)))
    # Two relations among seven momenta: the constant row and this cofactor row.
    c=[(-1)**i*int(fmpz_mat(q[:i]+q[i+1:6]).det()) for i in range(6)]+[0]
    assert not any(sum(a*b for a,b in zip(c,column)) for column in zip(*q))
    centered=[7*x-sum(c) for x in c]
    return sum(x*x for x in centered)**2,sum(x**4 for x in centered)


def images(p,charge):
    out=[0,0]
    for i in range(7):
        fused=tuple(a+b for a,b in zip(p[i],p[-1]))
        q=(*p[:i],fused,*p[i+1:7])
        values=quartics(q)
        for j in range(6):
            assert values==quartics((*q[:j],q[j+1],q[j],*q[j+2:]))
        for j in range(2):
            out[j]+=charge[i]*values[j]
    return tuple(out)



def bulk_contacts(p,charge):
    # Tree bracket with the integrated scalar cubic. Its two spectators fuse.
    values=[0,0]
    for i,j in combinations(range(7),2):
        remaining=[q for k,q in enumerate(p[:7]) if k not in (i,j)]
        base=int(fmpz_mat(remaining).det())
        if base:
            replaced=[int(fmpz_mat(remaining[:k]+[p[7]]+remaining[k+1:]).det()) for k in range(5)]
            pair=charge[i]+charge[j]
            values[0]+=pair*base**4
            values[1]+=pair*base**2*sum(x*x for x in replaced)
    return tuple(values)

def probe(b):
    p=list(MOMENTA[:7])
    p[3]=(0,b,0,0,0)
    p.append(tuple(-sum(col) for col in zip(*p)))
    # Polarization e_2-(1+b)e_5, transverse to the last input.
    charge=tuple(row[1]-(1+b)*row[4] for row in p)
    assert charge[-1]==0
    shapes=((1,5,5),(2,4,5),(3,3,5),(3,4,4))
    integrated=sum(sum(run(shape,tuple(p),charge)) for shape in shapes)
    return images(tuple(p),charge),integrated,bulk_contacts(tuple(p),charge)


if __name__=='__main__':
    rows=[probe(b) for b in (1,2)]
    print(rows)
    assert rows==[((19208,9604),F(1,1215),(4,44)),((172872,86436),F(2,405),(36,360))]
    assert [x-9*y for x,y in zip(rows[1][0],rows[0][0])]==[0,0]
    assert rows[1][1]-9*rows[0][1]==-F(1,405)
    for _,value,(a,b) in rows:
        assert value+F(2,3645)*a-F(1,14580)*b==0
