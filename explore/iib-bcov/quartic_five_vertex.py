"""Eight double residues for the first gamma master on positive K3,3 propagator powers."""
from fractions import Fraction as F
from math import prod

from quartic_nonplanar_cut import projection
from quartic_planar_vertex import node_residue


def double_pole(powers, shift):
    if min(powers) < 1:
        raise ValueError('Node charts use positive propagator powers.')
    target = tuple(p-1 for p in powers)
    positive = sum(node_residue(target, shift, axis) for axis in range(4))
    negative = sum(node_residue(target, shift, axis) for axis in range(4, 8))
    scale = F(4)**(1-2*shift)*F(8**shift, prod((2*j+1)*(2*j+2)*(2*j+3) for j in range(shift)))
    return -6*projection(powers, shift)+scale*(positive+2*negative)/8


def check():
    from json import loads
    from pathlib import Path
    from random import Random
    from quartic_five_pole import physical_poles
    rows = loads(Path(__file__).with_name('quartic_planar_terminal.json').read_text())
    for row in rows:
        p = row['uncut']
        if min(p) > 0:
            assert double_pole(p, 1) == F(row['top,planar_double,elementary_double'][2])
    rng = Random(8560)
    for trial in range(48):
        shift = trial % 4
        p = [1]*8
        for i in rng.sample(range(8), 1+trial % 5):
            p[i] += 1
        if trial % 7 == 0:
            p[trial % 8] += 1
        assert physical_poles(p, shift) == (0, double_pole(p, shift))
    print('All 48 dotted-node tests at four dimensions/shifts agree with the independent five-cut projector.')


if __name__ == '__main__':
    check()
