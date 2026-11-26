"""All bare primary kernels for positive-degree scalar colors on the compactified threefold."""
from collections import Counter
from fractions import Fraction
from math import factorial, lcm, prod
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'kodaira-spencer'))
from middle_point import rooted_box
from marked_theta import conductance, green_moments, partitions, square
from torus_quartic import quartic, shape_integral


def theta_integral(shape, waves, polarization):
    """Ordered positions on three parallel paths; the marked vertex is one trivalent end."""
    degrees = tuple(4-k for k in shape)
    mass = conductance(degrees)
    bridge = [conductance(tuple(x-(j==i) for j,x in enumerate(degrees))) for i in range(3)]
    charges = [sum(x*y for x,y in zip(p,polarization)) for p in waves]
    total = Fraction(0)
    for groups in partitions(shape):
        path = {i:j for j,group in enumerate(groups) for i in group}
        for powers, coefficient in square(groups, waves).items():
            for i,charge in enumerate(charges):
                if not charge:
                    continue
                xy, same_path = green_moments(powers, 3, i)
                total += coefficient*charge*(mass*xy + (bridge[path[i]]*same_path if path[i]==path[3] else 0))
    return total/prod(factorial(v) for v in Counter(shape).values())


def kernels(waves, polarization):
    """Scalar geometric coefficients on labeled waves, with the vector last, before ring contractions."""
    if len(waves) < 3 or len(polarization) != 3 or any(len(p)!=3 for p in waves):
        raise ValueError('Expected at least three three-dimensional waves and one polarization.')
    waves = tuple(tuple(Fraction(x) for x in p) for p in waves)
    polarization = tuple(Fraction(x) for x in polarization)
    if sum(x*y for x,y in zip(waves[-1], polarization)):
        raise ValueError('The vector polarization must be transverse.')
    if len(waves) == 3:
        return {(4,): rooted_box(waves, polarization)}
    if len(waves) == 5:
        return {shape:value/24 for shape,value in quartic(waves, polarization).items()}
    if len(waves) > 5:
        return {}
    scale = lcm(*(x.denominator for p in waves for x in p))
    integer_waves = tuple(tuple(int(scale*x) for x in p) for p in waves)
    return {shape:theta_integral(shape, integer_waves, polarization)/scale**13
            for shape in ((1,3,3),(2,2,3))}


if __name__ == '__main__':
    for waves, vector in (
        (((1,0,0), (1,0,0), (0,1,0), (0,0,1)), (1,1,0)),
        (((1,1,0), (0,1,-1), (2,0,1), (-3,-2,0)), (2,-3,1)),
        (((-2,1,1), (1,-3,5), (3,0,1), (0,2,1)), (0,1,-2)),
    ):
        values = kernels(waves, vector)
        for shape, value in values.items():
            original = shape_integral(tuple(k-1 for k in shape), waves, vector, ((0,1),)*3)
            assert original == value, (shape, original, value)
        print('Two-loop labeled threefold kernel:', *values.values(), 'sum', sum(values.values()))
