"""The nine-active-edge quotient with two successive charged series reductions."""
from collections import defaultdict
from fractions import Fraction as F
from functools import cache
from math import factorial, lcm, prod
from pathlib import Path
import sys

from flint import fmpz_mpoly_ctx
import sympy as sp

from genus_five_k4_period import CTX, folded_weights, moment
from genus_five_nine_faces import loads
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'nonSP/verify/k4_periods'))
from current_minor import numerator
from kernel import G, tree_flows

EDGES=((0,1),(0,5),(1,2),(1,5),(2,3),(2,4),(3,4),(3,5),(4,5))
OUTSIDE=fmpz_mpoly_ctx.get(('l23','l24','l34','l35','l45'))


def groups(poly):
    data=defaultdict(dict)
    for ex,c in poly.terms():
        k=tuple(map(int,ex))
        data[k[:4]][k[4:]]=c
    return {k:OUTSIDE.from_dict(row) for k,row in data.items()}


@cache
def beta(a,b):
    if min(a,b)<=0:raise ArithmeticError((a,b))
    return F(factorial(a-1)*factorial(b-1),factorial(a+b-1))


@cache
def path_integral(k):
    a,b,c,d=(i+1 for i in k)
    merged=a+b+d-14
    return beta(a,b)*beta(14-a-b,14-d)*beta(merged,c)


def core_value(core,scale):
    """Use a short convergent kite cut and an independently normalized reciprocal cut."""
    bank=folded_weights(core,14)
    total=[F(),F()]
    for p,c in bank.items():
        val=moment(p)
        assert val==moment(tuple(14-i for i in reversed(p)))
        for i,w in enumerate(val):total[i]-=c*w/(192*scale)
    return len(bank),tuple(total)


def integrate(marks):
    graph=G(6,EDGES)
    charge=loads(marks,6)
    trees=tree_flows(graph,charge)
    raw=numerator(graph,charge,next(iter(trees.values())))
    ctx=raw.context()
    x=ctx.gens()
    tree={m:prod(x[i] for i in m) for m in trees}
    current=[[sum((v[e][d]*tree[m] for m,v in trees.items()),ctx.constant(0))
              for d in range(2)] for e in range(9)]
    b=sp.Matrix([[int(v==i)-int(v==j) for i,j in EDGES] for v in range(5)])
    chords=tuple(j for j in range(9) if j not in b.rref()[1])
    assert chords==(3,6,7,8)
    root=sum(v*v for v in current[2])
    insert=sum((2*current[3][d]+current[7][d]+current[8][d])**2+
               2*current[7][d]**2+2*current[8][d]**2 for d in range(2))
    for lengths in ((2,3,1,5,7,11,4,6,3),(3,2,1,4,5,6,7,8,9)):
        c=sp.diag(*(sp.Rational(1,i) for i in lengths))
        tr=c*b.T*(b*c*b.T).inv()
        j=tr*sp.Matrix(charge[:-1])
        h=-(c-tr*b*c).multiply_elementwise(j*j.T)
        u=int(sum(tree.values())(*lengths))
        assert h.minor_submatrix(2,2).det()==F(int(raw(*lengths))**2,u**10)
        observed=j.row(2).dot(j.row(2))*(
            (2*j.row(3)+j.row(7)+j.row(8)).norm()**2+2*j.row(7).dot(j.row(7))+2*j.row(8).dot(j.row(8)))
        assert observed==F(int(root(*lengths))*int(insert(*lengths)),u**4)
    # In the projective measure the squared current minor contributes one
    # root resistance. Combine only the five unaffected lengths in FLINT.
    left=groups(x[2]*raw**2)
    right=groups(root*insert)
    sums={tuple(a+b for a,b in zip(p,q)) for p in left for q in right}
    scale=lcm(*(path_integral(k).denominator for k in sums))
    weights={k:int(path_integral(k)*scale) for k in sums}
    reduced=defaultdict(lambda:OUTSIDE.constant(0))
    for p,A in left.items():
        for q,B in right.items():
            k=tuple(a+b for a,b in zip(p,q))
            power=sum(k)+3-14
            reduced[power]+=weights[k]*A*B
    terms={}
    for k,poly in reduced.items():
        for ex,c in poly.terms():
            exp=tuple(map(int,ex))
            terms[(*exp[:2],k,*exp[2:])]=c
    core=CTX.from_dict(terms)
    count,total=core_value(core,scale)
    print(marks,'convergent core monomials and orbits',len(core),count,*total,flush=True)
    expected={(3,4,0):(F(-2030121825437,111965087975989248000000),F(597739,43196407398144000)),
              (0,3,4):(F(-7425850380601,895720703807913984000000),F(2127563,345571259185152000))}
    assert total==expected[marks]
    return total


if __name__=='__main__':
    for marks in ((3,4,0),(0,3,4)):integrate(marks)
