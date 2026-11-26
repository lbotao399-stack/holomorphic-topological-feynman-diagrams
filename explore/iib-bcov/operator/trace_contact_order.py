"""Order of the nine-primary trace along conserved real momentum arcs."""
from functools import cache
from math import inf

from flint import fmpq, fmpq_poly
from three_loop_soft_jet import SUBSETS


def covering_order(orders):
    """Return the trace order and a minimizing triple of zero-based five-row sets."""
    if len(orders)!=len(SUBSETS) or any(n!=inf and (int(n)!=n or n<0) for n in orders):
        raise ValueError("Expected 126 nonnegative integer orders or infinity.")
    masks = [(sum(1<<i for i in subset),n,subset) for subset,n in zip(SUBSETS,orders) if n!=inf]
    best = {0:(0,())}
    for _ in range(3):
        following = {}
        for covered,(cost,witness) in best.items():
            for mask,order,subset in masks:
                key = covered|mask
                new = cost+order
                if key not in following or new<following[key][0]:
                    following[key] = (new,(*witness,subset))
        best = following
    result,witness = best.get((1<<9)-1,(inf,()))
    return 2*result,witness


def polynomial(value):
    if isinstance(value,fmpq_poly):
        return value
    coeffs = value if isinstance(value,(tuple,list)) else (value,)
    return fmpq_poly([fmpq(str(x)) for x in coeffs])


def polynomial_order(points):
    """Each momentum entry is rational or a finite coefficient list in ascending powers."""
    if len(points)!=9 or any(len(row)!=5 for row in points):
        raise ValueError("Expected nine five-component momentum curves.")
    points = [[polynomial(v) for v in row] for row in points]
    if any(sum(column) for column in zip(*points)):
        raise ValueError("Momentum must be conserved as a polynomial identity.")

    @cache
    def minor(rows,columns):
        if not rows:
            return fmpq_poly(1)
        return sum(((-1)**i*points[rows[0]][j]*minor(rows[1:],columns[:i]+columns[i+1:])
                    for i,j in enumerate(columns) if points[rows[0]][j]),fmpq_poly())

    orders = [next((i for i,c in enumerate(minor(s,tuple(range(5)))) if c),inf) for s in SUBSETS]
    return covering_order(orders)


def check():
    t = fmpq_poly([0,1])
    frame = [tuple(int(i==j) for i in range(5)) for j in range(5)]
    for bend in (0,1):
        # Three hard labels stay on line 1. A fourth soft label has the same
        # velocity and can first leave that line at an arbitrary higher order.
        for power in (2,3,5):
            hard = [frame[0]]*3+frame[1:]
            soft = (t,bend*t**power,0,0,0)
            eight = [*hard,soft]
            points = [*eight,tuple(-sum(column) for column in zip(*eight))]
            order,witness = polynomial_order(points)
            assert order==(2*power if bend else inf)
            assert bool(witness)==bool(bend)
            # This diagonal rank degeneration adds 3 to every finite minor order.
            scaled = [tuple(x*factor for x,factor in zip(p,(1,1,1,t,t**2))) for p in points]
            assert polynomial_order(scaled)[0]==order+18
    print('Finite, infinite and rank-degenerate contact orders verified.')


if __name__=='__main__':
    check()
