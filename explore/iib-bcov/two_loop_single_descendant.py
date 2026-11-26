"""Genus-two fivefold elementary vertices with one scalar descendant, levels zero to four."""
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import comb, lcm

from flint import fmpz_mat
from bulk_eight import BULK
from genus_one_graph_ward import plus
from longitudinal_theta import longitudinal
from marked_boundaries import quartics
from nine_wave_theta import endpoint_potentials
from two_loop_elliptic_insertions import single_mark
from two_loop_moving_descendant import descendant_vertices, fused_symbols


def primary(waves):
    return -sum(c*y for c,y in zip(BULK, quartics(waves)))


def complement(items, subset):
    return tuple(p for i,p in enumerate(items) if i not in subset)


def cut_difference(waves, first, second, project):
    """Eight-site theta potential difference after contraction with a test vector."""
    values = fused_symbols(waves)
    def row(i):
        return (longitudinal(waves[:i]+waves[i+1:]) * project(waves[i])
                + sum(value * project(waves[b if a == i else a])
                      for (a,b), value in values.items() if i in (a,b)))
    return row(first)-row(second)


@lru_cache(maxsize=32768)
def _single(marked, others):
    level = len(others)-6
    waves = (marked,*others)
    if level == 0:
        return primary(waves)
    for i,p in enumerate(others):
        if not any(p):
            return _single(marked,others[:i]+others[i+1:])
    if fmpz_mat(waves).rank() < 5:
        return Fraction(0)
    if not any(marked):
        # The path shifts two momenta in one common direction. Four five-minors
        # bound its degree by four. Avoid also producing a zero primary.
        axis = next(i for i in range(5) if any(others[0][j] for j in range(5) if j!=i))
        total = Fraction(0)
        for t in range(1,6):
            moved = tuple(t*int(i==axis) for i in range(5))
            recoil = tuple(a-b for a,b in zip(others[0],moved))
            total += (-1)**(t-1)*comb(5,t)*_single(moved,tuple(sorted((recoil,*others[1:]))))
        return total
    if level == 1:
        return descendant_vertices(waves)[0]
    # The preceding Ward equation has the marked slot at level h.
    h = level-1
    partner,x,y,divisor = next((r,x,y,det) for r,p in enumerate(others)
                               for x,y in combinations(range(5),2)
                               if (det := marked[x]*p[y]-marked[y]*p[x]))
    test = others[partner]
    def project(p):
        return p[x]*test[y]-p[y]*test[x]
    labels = tuple(i for i in range(len(others)) if i!=partner)
    result = Fraction(0)
    # The marked bivalent theta site contains h other primary waves.
    for group in combinations(labels,h):
        charge = plus((marked,*(others[i] for i in group)))
        keep = tuple(i for i in range(len(others)) if i not in group)
        sites = (charge,*(others[i] for i in keep))
        result += cut_difference(sites,0,1+keep.index(partner),project)
    # The marked trivalent site contains h-1 other primary waves. Integrate with
    # that endpoint grounded; the two choices of endpoint give the same term.
    for group in combinations(labels,h-1):
        keep = tuple(i for i in range(len(others)) if i not in group)
        sites = tuple(others[i] for i in keep)
        result -= 2*project(endpoint_potentials(sites)[keep.index(partner)])
    # In a cubic cut disjoint from the distinguished slot, only a pair containing
    # the test slot contributes to this difference of grounded potentials.
    for i in labels:
        q = plus((others[i],test))
        if project(q):
            reduced = tuple(sorted((*complement(others,(i,partner)),q)))
            result -= project(q)*_single(marked,reduced)
    # The other classical end has the descendant and h+1 ordinary fields.
    for group in combinations(labels,h+1):
        q = plus((marked,*(others[i] for i in group)))
        if project(q):
            result += project(q)*primary((*complement(others,group),q))
    circle = single_mark(waves,0)
    result += project(circle[0])-project(circle[1+partner])
    return result/divisor


def single_vertex(momenta, marked=0, *, full_source=False):
    """Worldsheet vertex with one marked scalar of level n-7; full sources admit any level."""
    if len(momenta)<7 or not 0<=marked<len(momenta):
        raise ValueError('Expected at least seven waves and one valid marked slot.')
    if not full_source and len(momenta)>11:
        raise ValueError('The minimal scalar field has descendant level at most four.')
    if any(len(p)!=5 for p in momenta) or any(plus(momenta)):
        raise ValueError('Expected conserved five-dimensional momenta.')
    denom = lcm(*(Fraction(x).denominator for p in momenta for x in p))
    waves = tuple(tuple(int(denom*Fraction(x)) for x in p) for p in momenta)
    return _single(waves[marked],tuple(sorted(waves[:marked]+waves[marked+1:])))/denom**20


def top_fixture():
    retained = [(1,-1,0,0,0),(1,0,0,0,0),(-2,1,0,0,0)]
    extra = [tuple(sign*int(i==j) for j in range(5)) for i in range(1,5) for sign in (1,-1)]
    return retained, extra


if __name__ == '__main__':
    from ancestor_projection import fusions
    from two_loop_moving_descendant import BASIS
    waves = [*BASIS[:3],(-1,-1,-1,-1,-1),(0,0,0,1,1),BASIS[3],(0,0,0,-1,0),BASIS[4],(0,0,0,0,-1)]
    value = single_vertex(waves,4)
    assert value == Fraction(649,77760)
    print('Two-loop moving level two:', value, flush=True)
    retained, extra = top_fixture()
    totals = [sum(single_vertex(q) for q,k in fusions(retained,extra,4,j)) for j in range(5)]
    pulled_back = [sum((-1)**(l-j)*comb(l,j)*totals[l] for l in range(j,5)) for j in range(5)]
    assert pulled_back == [Fraction(91,360), Fraction(347,540), Fraction(3227,6480),
                           Fraction(137,1215), Fraction(179,38880)]
    assert sum(pulled_back) == totals[0] == Fraction(19579,12960)
    print('Two-loop moving level four:', totals[0], 'ancestor:', pulled_back[0], flush=True)
