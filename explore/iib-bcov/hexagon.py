"""Primary BCOV hexagon on C^5, with the closed polyvector propagator."""
from fractions import Fraction
from itertools import permutations
from math import factorial,prod
from flint import fmpz_mat


def dot(x,y):
    return sum(a*b for a,b in zip(x,y))


def det(rows):
    return int(fmpz_mat(rows).det())


def currents(order,momenta,time):
    # Return y and grounded potentials multiplied by the total resistance.
    prefix=[0]*5
    sums=[]
    for label in order:
        prefix=[a+b for a,b in zip(prefix,momenta[label])]
        sums.append(prefix)
    assert not any(prefix)
    total=sum(time)
    offset=[dot(time,col) for col in zip(*sums)]
    y=[[total*a-b for a,b in zip(p,offset)] for p in sums]
    phi=[[0]*5]
    for t,current in zip(time,y):
        phi.append([a-t*b for a,b in zip(phi[-1],current)])
    assert not any(phi.pop())
    return y,phi


def odd_integral(order,y,polarizations):
    # Literal 30-variable top exterior coefficient of the six propagators and six fields.
    rows=[]
    scale=1
    for i,current in enumerate(y):
        pivot=next(j for j,x in enumerate(current) if x)
        scale*=(-1)**pivot*current[pivot]**3
        for j in range(5):
            if j==pivot:
                continue
            row=[0]*30
            for v,s in ((i,1),((i+1)%6,-1)):
                row[5*v+j]=s*current[pivot]
                row[5*v+pivot]=-s*current[j]
            rows.append(row)
    for label,u in enumerate(polarizations):
        row=[0]*30
        start=5*order.index(label)
        row[start:start+5]=u
        rows.append(row)
    return Fraction(det(rows),scale)


def averaged_green():
    # Positions 1,...,5 are measured from the grounded vertex at position zero.
    g=[[Fraction(min(i,j)*(6-max(i,j)),42) for j in range(6)] for i in range(6)]
    trace=sum(g[i][i] for i in range(6))
    diagonal=trace/6
    off_diagonal=(sum(map(sum,g))-trace)/30
    assert diagonal-off_diagonal==Fraction(1,12)
    return (diagonal-off_diagonal)/factorial(5)


def run():
    unit=[tuple(int(i==j) for j in range(5)) for i in range(5)]
    momenta=unit+[(-1,)*5]
    transverse=unit[1:]+[unit[0],(1,-1,0,0,0)]
    assert not any(dot(p,u) for p,u in zip(momenta,transverse))
    longitudinal=transverse[:-1]+[unit[0]]
    time=(1,2,3,4,5,6)
    for polarizations in (transverse,longitudinal):
        raw=0
        nonzero=0
        for order in permutations(range(6)):
            y,phi=currents(order,momenta,time)
            minor=det([momenta[j] for j in order[:5]])
            klein=det([list(u)+[dot(u,phi[order.index(i)])]
                       for i,u in enumerate(polarizations)])
            value=odd_integral(order,y,polarizations)
            # Currents and potentials in the test were both multiplied by sum(time).
            assert value==minor*sum(time)**4*klein, (order,value,klein)
            nonzero+=bool(value)
            raw+=value*minor
        expected=bool(dot(momenta[-1],polarizations[-1]))
        assert nonzero and bool(raw)==expected
        print('ordered nonzero',nonzero,'of 720, polarized exterior sum',raw)
    print('integrated diagonal minus off-diagonal resistance response',averaged_green())


if __name__=='__main__':
    run()
