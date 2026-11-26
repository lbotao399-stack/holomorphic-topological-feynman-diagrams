"""Exact equal-recoil soft expansion of the seven-primary genus-two trace."""
from fractions import Fraction as F
from math import lcm

from flint import fmpq_mat
from bulk_eight import BULK
from marked_boundaries import quartics


def recoil(hard, velocity):
    """Coefficients of t^2, t^3, t^4 with hard inputs shifted by -t*velocity/6."""
    if len(hard)!=6 or any(len(p)!=5 for p in (*hard,velocity)) or any(map(sum,zip(*hard))):
        raise ValueError("Six conserved hard covectors and one soft velocity are required.")
    frame = fmpq_mat([[str(x) for x in p] for p in hard[:-1]])
    determinant = F(str(frame.det()))
    if not determinant:
        raise ValueError("The hard six-wave configuration must span five dimensions.")
    column = frame.transpose().solve(fmpq_mat([[str(x)] for x in velocity]))
    coordinates = [F(str(x)) for x in column.entries()]+[F(0)]
    mean = sum(coordinates)/6
    second,third,fourth = [sum((x-mean)**k for x in coordinates) for k in (2,3,4)]
    return tuple(-determinant**4*x/1244160 for x in
                 (48*second,-24*third,31*second**2-42*fourth))


def check():
    hard = [tuple(F(1+i,2) if j==i else F(int(j==i-1),3) for j in range(5))
            for i in range(5)]
    hard.append(tuple(-sum(column) for column in zip(*hard)))
    for v in (hard[0],hard[-1],(2,-3,4,7,-1),(0,0,0,0,0)):
        coefficients = recoil(hard,v)
        for t in (F(-3),F(0),F(2,5),F(6),F(9)):
            soft = tuple(t*x for x in v)
            waves = [tuple(x-y/6 for x,y in zip(p,soft)) for p in hard]+[soft]
            scale = lcm(*(x.denominator for p in waves for x in p))
            integer = [tuple(int(x*scale) for x in p) for p in waves]
            direct = sum(x*c for x,c in zip(quartics(integer),BULK))/scale**20
            assert sum(c*t**k for k,c in enumerate(coefficients,2))==direct
        first,middle,last = coefficients
        assert middle**2<=4*first*last
        assert (middle**2==4*first*last)==(v in (hard[0],hard[-1]) or not any(v))
    print('Seven-primary recoil polynomial and equality directions reproduced on rational inputs.')


if __name__=='__main__':
    check()
