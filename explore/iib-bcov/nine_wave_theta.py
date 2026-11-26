"""Bare two-loop Ward potentials on nine waves with one external level-one scalar."""
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import lcm

from bulk_eight import BULK, grounded_potentials, potentials
from genus_one_graph_ward import plus, separate_divergences
from longitudinal_theta import closed_potentials
from marked_boundaries import quartics
from two_loop_moving_descendant import descendant_vertices


@lru_cache(maxsize=96)
def endpoint_potentials(waves):
    """Eight bivalent marks with balancing charge at the grounded branch vertex."""
    denom = lcm(*(Fraction(x).denominator for p in waves for x in p))
    integral = tuple(tuple(int(x*denom) for x in p) for p in waves)
    rows = grounded_potentials(integral, check=True)
    return tuple(tuple(x/denom**21 for x in p) for p in rows)


def bare_theta(waves, descendant, direct=False):
    if len(waves) != 9 or any(len(p) != 5 for p in waves) or any(plus(waves)):
        raise ValueError('Supply nine conserved five-dimensional waves.')
    labels = tuple(i for i in range(9) if i != descendant)
    ordinary = tuple(waves[i] for i in labels)
    out = [[Fraction(0)]*5 for p in waves]
    # The quartic can replace either trivalent endpoint. Both give its zero grounded row.
    for i, p in zip(labels, endpoint_potentials(ordinary)):
        out[i] = [2*x for x in p]
    # Otherwise a bivalent site has the descendant and one ordinary external field.
    for partner in labels:
        groups = tuple((i, descendant) if i == partner else (i,) for i in labels)
        charges = tuple(plus(waves[i] for i in group) for group in groups)
        rows = potentials(charges, check=True) if direct else closed_potentials(charges)
        for group, row in zip(groups, rows):
            for i in group:
                out[i] = [x+y for x, y in zip(out[i], row)]
    return out


def tree_potentials(waves, descendant):
    """One cubic cut, or the quartic cut carrying the level-one wave."""
    out = [[Fraction(0)]*5 for p in waves]
    labels = tuple(i for i in range(9) if i != descendant)
    for i,j in combinations(labels, 2):
        group = (i,j)
        keep = [k for k in range(9) if k not in group]
        charge = plus(waves[k] for k in group)
        fused = tuple(waves[k] for k in keep) + (charge,)
        value = descendant_vertices(fused)[keep.index(descendant)]
        for x in group:
            out[x] = [a+value*b for a,b in zip(out[x],charge)]
        # Moving the marked wave into this classical end gives the quartic cut.
        group = (i,j,descendant)
        charge = plus(waves[k] for k in group)
        fused = tuple(waves[k] for k in range(9) if k not in group) + (charge,)
        denom = lcm(*(Fraction(x).denominator for p in fused for x in p))
        integral = tuple(tuple(int(x*denom) for x in p) for p in fused)
        value = -sum(c*y for c,y in zip(BULK, quartics(integral)))/denom**20
        for x in group:
            out[x] = [a+value*b for a,b in zip(out[x],charge)]
    return out


def check_theta():
    from two_loop_moving_descendant import BASIS
    ordinary = (*BASIS, (1,2,-1,3,2), (-2,1,3,1,-2))
    ordinary += (tuple(-x for x in plus(ordinary)),)
    zero = (0,)*5
    soft = bare_theta((*ordinary,zero),8)
    raw = potentials(ordinary)
    assert all(soft[i][j]-soft[0][j] == 10*(raw[i][j]-raw[0][j])
               for i in range(8) for j in range(5))
    waves = (*ordinary[:-1], (1,1,0,0,0), tuple(a-b for a,b in zip(ordinary[-1],(1,1,0,0,0))))
    reduced, direct = bare_theta(waves,7), bare_theta(waves,7,direct=True)
    assert all(reduced[i][j]-reduced[0][j] == direct[i][j]-direct[0][j]
               for i in range(9) for j in range(5))
    print('The two endpoint and eight paired-site contributions recover the ten-vertex dilaton factor.')
    print('All nine bare potential differences match the four uncompressed theta shapes.')


if __name__ == '__main__':
    check_theta()
