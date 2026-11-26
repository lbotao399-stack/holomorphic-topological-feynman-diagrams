"""Exact eight-primary BCOV theta response on arbitrary integral conserved momenta."""
from fractions import Fraction as F
from functools import cache
from pathlib import Path
import sys
from itertools import combinations,combinations_with_replacement
from math import factorial,lcm,prod
from collections import Counter
from flint import fmpz_mat
from marked_theta import partitions,conductance,green_moments,MOMENTA,square
from marked_boundaries import quartics

SHAPES=((1,5,5),(2,4,5),(3,3,5),(3,4,4))
SLOTS=range(8)


@cache
def template(shape):
    groups=next(partitions(shape))
    paths={i:j for j,g in enumerate(groups) for i in g}
    # Two slots on the same path have dependent images in the determinant.
    pairs=[(i,j) for i,j in combinations(SLOTS,2) if paths[i]!=paths[j]]
    products=list(combinations_with_replacement(range(len(pairs)),2))
    matrix=[]
    scale=1
    deg=[6-v for v in shape]
    mass=conductance(deg)/(2*prod(factorial(k) for k in Counter(shape).values()))
    for h,j in products:
        word=pairs[h]+pairs[j]
        k=tuple(word.count(i) for i in SLOTS)
        row=[]
        for a,b in combinations_with_replacement(SLOTS,2):
            xy,br=green_moments(k,a,b)
            value=mass*(1 if h==j else 2)*(xy+(F(9,deg[paths[a]])*br if paths[a]==paths[b] else 0))
            row.append(value)
            scale=lcm(scale,value.denominator)
        matrix.append(row)
    moments=fmpz_mat([[int(x*scale) for x in row] for row in matrix])
    # Two five-row blocks; the first and second L columns go in different blocks.
    left,middle,right=groups
    splits=[]
    for chosen in combinations(middle,4-len(left)):
        top=left+chosen
        bottom=tuple(i for i in middle if i not in chosen)+right
        split_sign=(-1)**(sum([0,*[2+i for i in top]])-10)
        splits.append((split_sign,top,bottom))
    return pairs,products,moments,scale,splits


def determinants(p):
    table={}
    for ids in combinations(SLOTS,5):
        table[ids]=int(fmpz_mat([p[i] for i in ids]).det())
    out={}
    for group in combinations(SLOTS,4):
        row=[0]*8
        for i in SLOTS:
            if i not in group:
                row[i]=(-1)**sum(j<i for j in group)*table[tuple(sorted((*group,i)))]
        out[group]=row
    return out


def coeffs(order,splits,pairs,minors):
    q=[0]*len(pairs)
    for sign,top,bottom in splits:
        # Each path is in increasing label order. Restore sorting across paths.
        def linear(slots):
            idx=[order[i] for i in slots]
            parity=(-1)**sum(a>b for a,b in combinations(idx,2))
            return parity,minors[tuple(sorted(idx))]
        t,left=linear(top)
        b,right=linear(bottom)
        scale=sign*t*b
        for h,(i,j) in enumerate(pairs):
            a,b=order[i],order[j]
            q[h]+=scale*(left[a]*right[b]+left[b]*right[a])
    return q


def potentials(p,check=False,shapes=SHAPES):
    if any(map(sum,zip(*p))):
        raise ValueError('Eight conserved five-component momenta required.')
    return grounded_potentials(p,check,shapes)


def grounded_potentials(p,check=False,shapes=SHAPES):
    """Ground the balancing charge at the second trivalent endpoint."""
    if len(p)!=8 or any(len(q)!=5 for q in p):
        raise ValueError('Eight five-component momenta required.')
    minors=determinants(p)
    out=[[F() for _ in SLOTS] for _ in SLOTS]
    for shape in shapes:
        pairs,products,green,den,splits=template(shape)
        orders=[]
        values=[]
        for groups in partitions(shape):
            order=sum(groups,())
            q=coeffs(order,splits,pairs,minors)
            if not any(q):
                continue
            if check and not orders:
                sparse=square(groups,p)
                for times in ((2,1,4,3,7,3,2,1),(3,5,2,7,4,1,8,6)):
                    test=sum(c*times[order[i]]*times[order[j]] for c,(i,j) in zip(q,pairs))
                    assert test**2==sum(c*prod(s**e for s,e in zip(times,k)) for k,c in sparse.items())
            values.append([q[i]*q[j] for i,j in products])
            orders.append(order)
        if not values:
            continue
        weighted=(fmpz_mat(values)*green).tolist()
        total=[[0]*8 for _ in SLOTS]
        for order,row in zip(orders,weighted):
            for (i,j),value in zip(combinations_with_replacement(SLOTS,2),row):
                a,b=sorted((order[i],order[j]))
                total[a][b]+=int(value)
        for i in SLOTS:
            for j in SLOTS:
                a,b=sorted((i,j))
                out[i][j]+=F(total[a][b],den)
    return [tuple(sum(row[j]*p[j][k] for j in SLOTS) for k in range(5)) for row in out]


def primitives(p):
    out=[[[0]*5 for _ in SLOTS] for _ in range(2)]
    for i,j in combinations(SLOTS,2):
        fused=tuple(x+y for x,y in zip(p[i],p[j]))
        q=tuple(p[h] for h in SLOTS if h not in (i,j))+(fused,)
        for a,coef in enumerate(quartics(q)):
            for k in range(5):
                out[a][i][k]+=coef*p[j][k]
                out[a][j][k]+=coef*p[i][k]
    return out


def transverse_rows(p,a):
    # Eliminate arbitrary shifts of each potential by its source covector
    # and the single common translation of all eight potentials.
    from flint import fmpq_mat
    gauges=[[int(i==h)*p[i][j] for i in SLOTS for j in range(5)] for h in SLOTS]
    gauges.extend([[int(j==h) for i in SLOTS for j in range(5)] for h in range(5)])
    basis,rank=fmpq_mat(gauges).rref()
    flat=[F(x) for row in a for x in row]
    for i in range(rank):
        row=[F(str(basis[i,j])) for j in range(40)]
        pivot=next(j for j,z in enumerate(row) if z)
        factor=flat[pivot]
        flat=[x-factor*y for x,y in zip(flat,row)]
    return flat


BULK=(F(-31,2987228160),F(1,71124480))


def residual(p,fields):
    basis=primitives(p)
    return [[fields[i][k]-sum(BULK[a]*basis[a][i][k] for a in range(2))
             for k in range(5)] for i in SLOTS]


def cubic_graph(p,fields,vector):
    return tuple(10*fields[vector][k]-sum(row[k] for row in fields) for k in range(5))


def cubic_contact(p,vector):
    result=[F()]*5
    for i,j in combinations((i for i in SLOTS if i!=vector),2):
        fused=tuple(x+y for x,y in zip(p[i],p[j]))
        q=tuple(p[k] for k in SLOTS if k not in (i,j))+(fused,)
        weight=sum(x*y for x,y in zip(BULK,quartics(q)))
        for k in range(5):
            result[k]+=fused[k]*weight
    return result


if __name__=='__main__':
    from flint import fmpq_mat
    sys.path.insert(0,str(Path(__file__).resolve().parent/'operator'))
    from genus_two_probes import MOMENTA as probe,SECOND_MOMENTA,evaluate_potential,evaluate_strong
    ps=primitives(probe)
    graph=potentials(probe,check=True)
    arrays=[evaluate_potential(x) for x in (*ps,graph)]
    solution,rank=fmpq_mat([[str(a) for a in row] for row in zip(*arrays)]).rref()
    assert rank==2 and [F(str(solution[i,2])) for i in range(2)]==list(BULK)
    assert not any(evaluate_potential(residual(probe,graph)))
    second=potentials(SECOND_MOMENTA,check=True)
    errors=residual(probe,graph),residual(SECOND_MOMENTA,second)
    assert not any(evaluate_strong(*errors))
    for p,fields in ((probe,graph),(SECOND_MOMENTA,second)):
        ps=primitives(p)
        for i in SLOTS:
            raw=cubic_graph(p,fields,i)
            bulk=cubic_contact(p,i)
            marked=[9*sum(BULK[h]*ps[h][i][k] for h in range(2)) for k in range(5)]
            remainder=[a+b-c for a,b,c in zip(raw,bulk,marked)]
            assert all(remainder[j]*p[i][k]==remainder[k]*p[i][j]
                       for j,k in combinations(range(5),2))
    print('Fourteen bulk and twenty joint polynomial identity coordinates vanish.')
    print('Seven-primary trace coefficients for the four theta sum:',*BULK)
    each=[evaluate_potential(potentials(probe,shapes=(shape,))) for shape in SHAPES]
    print('Single-shape ranks with the two universal contacts:',
          [fmpq_mat([[str(x) for x in row] for row in zip(*arrays[:2],v)]).rref()[1] for v in each])
    print('Rank of the four individual theta responses:',
          fmpq_mat([[str(x) for x in row] for row in each]).rref()[1])
    physical=potentials(MOMENTA,check=True)
    charge=(1,0,-2,0,0)
    raw=sum(x*y for x,y in zip(cubic_graph(MOMENTA,physical,7),charge))
    assert raw==F(1,1215)
    print('Independent integrated-cubic test:',raw)
