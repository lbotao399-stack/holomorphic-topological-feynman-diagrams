"""Primary two-loop equation for the integrated quartic in minimal BCOV on C^5."""
from fractions import Fraction as F
from functools import cache
from itertools import combinations, combinations_with_replacement, permutations
from math import factorial, lcm
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from bulk_eight import BULK, grounded_potentials
from marked_boundaries import quartics
from flint import fmpz_mat


def det(p):
    return int(fmpz_mat(p).det())


def add(*rows):
    return tuple(map(sum, zip(*rows)))


def dot(a,b):
    return sum(x*y for x,y in zip(a,b))


def quartic_eight(p):
    # Four S_8-orbits of products of squared five-row determinants.
    assert len(p) == 8 and not any(add(*p))
    slots = tuple(combinations(range(8), 5))
    minors = [det([p[i] for i in s])**2 for s in slots]
    values = [0] * 4
    for (i, a), (j, b) in combinations_with_replacement(enumerate(slots), 2):
        values[len(set(a) & set(b))-2] += minors[i] * minors[j]
    return values


@cache
def tail_moment(positions):
    # T_i is the (6-i)-th ascending order statistic of five uniform variables.
    out = F(1)
    for h, position in enumerate(sorted(positions, reverse=True)):
        out *= F(6 - position + h, 6 + h)
    return out


@cache
def hexagon_moments():
    k = {}
    for j in range(5):
        for ell in range(5):
            for a,b in combinations_with_replacement(range(5),2):
                # Entries depend just on the equality pattern; use restricted-growth labels.
                pattern = (j, ell, a, b)
                labels = {v: i for i,v in enumerate(dict.fromkeys(pattern))}
                key = tuple(labels[v] for v in pattern)
                if key in k:
                    continue
                m = max(key)+1
                x = F()
                for order in permutations(range(1,6),m):
                    v,l,s,t = (order[i] for i in key)
                    x += tail_moment((max(v,l),s,t)) - tail_moment((v,l,s,t))
                k[key] = x * F(factorial(5-m),factorial(5))
    scale = lcm(*(val.denominator for val in k.values()))
    out = {}
    for j in range(5):
        for ell in range(5):
            for a,b in combinations_with_replacement(range(5),2):
                pattern = (j,ell,a,b)
                labels = {v:i for i,v in enumerate(dict.fromkeys(pattern))}
                key = tuple(labels[v] for v in pattern)
                out[j,ell,a,b] = int(k[key]*scale)*(1 if a==b else 2)
    return scale,out


def one_loop_on_primary_mark(p):
    """The 6-scalar genus-one marked primitive inserted in a hexagon."""
    scale, moments = hexagon_moments()
    values = [[0]*5 for _ in p]
    for active in combinations(range(9),5):
        rows = [p[i] for i in active]
        square = det(rows)**2
        if not square:
            continue
        spectators = [p[i] for i in range(9) if i not in active]
        minors = [det([*spectators,z]) for z in rows]
        for j, slot in enumerate(active):
            for ell, wave in enumerate(rows):
                weight = square*sum(moments[j,ell,a,b]*minors[a]*minors[b]
                                    for a,b in combinations_with_replacement(range(5),2))
                for k in range(5):
                    values[slot][k] += weight*wave[k]
    return [[F(x,8*scale) for x in row] for row in values]


@cache
def averaged_green(j, ell, tail):
    pattern = (j,ell,*tail)
    labels = {v:i for i,v in enumerate(dict.fromkeys(pattern))}
    key = tuple(labels[v] for v in pattern)
    return green_pattern(key)


@cache
def green_pattern(key):
    m = max(key)+1
    total = F()
    for order in permutations(range(1,6),m):
        j,ell,*tail = (order[i] for i in key)
        total += tail_moment((max(j,ell),*tail)) - tail_moment((j,ell,*tail))
    return total*F(factorial(5-m),factorial(5))


def general_one_loop(p):
    # The entire squared-minor symbol at the marked vertex is at most quadratic in the tails.
    values = [[F()]*5 for _ in p]
    for active in combinations(range(len(p)),5):
        rows = [p[i] for i in active]
        square = det(rows)**2
        if not square:
            continue
        spectator = [p[i] for i in range(len(p)) if i not in active]
        poly = {(): sum(det(q)**2 for q in combinations(spectator,5))}
        for size in (3,4):
            for group in combinations(spectator,size):
                extra = [add(*rows)] if size==3 else []
                coef = [det([*group,v,*extra]) for v in rows]
                if size==4:
                    poly[()] += sum(coef)**2
                    for i,c in enumerate(coef):
                        poly[i,] = poly.get((i,),0) - 2*sum(coef)*c
                factor = 1 if size==3 else 2
                for i,j in combinations_with_replacement(range(5),2):
                    poly[i,j] = poly.get((i,j),0) + factor*(1 if i==j else 2)*coef[i]*coef[j]
        for j,slot in enumerate(active):
            for ell,wave in enumerate(rows):
                weight = F(square,48)*sum(v*averaged_green(j,ell,t) for t,v in poly.items() if v)
                for k in range(5):
                    values[slot][k] += weight*wave[k]
    return values


def overlap_symbols(p):
    # Each ordered five-subset pair contributes according to its overlap; divide off the ordering.
    subsets = list(combinations(range(len(p)),5))
    vals = [det([p[i] for i in subset])**2 for subset in subsets]
    masks = [sum(1<<i for i in subset) for subset in subsets]
    out = [0]*6
    for i,v in enumerate(vals):
        if not v:
            continue
        mask = masks[i]
        for j in range(i,len(vals)):
            out[(mask&masks[j]).bit_count()] += v*vals[j]
    return out


def general_fused_marked(p):
    out = [[[0]*5 for _ in p] for h in range(6)]
    for i,j in combinations(range(len(p)),2):
        q = [p[h] for h in range(len(p)) if h not in (i,j)] + [add(p[i],p[j])]
        for h,val in enumerate(overlap_symbols(q)):
            for k in range(5):
                out[h][i][k] += val*p[j][k]
                out[h][j][k] += val*p[i][k]
    return out


def fused_marked(p):
    # Advect each scalar by the external vector; its momentum label is the other fused slot.
    out = [[[0]*5 for _ in p] for h in range(4)]
    for i, j in combinations(range(9),2):
        q = [p[h] for h in range(9) if h not in (i,j)] + [add(p[i],p[j])]
        for h, val in enumerate(quartic_eight(q)):
            for k in range(5):
                out[h][i][k] += val*p[j][k]
                out[h][j][k] += val*p[i][k]
    return out


def tree_bulk(p):
    out = [[F()]*5 for _ in p]
    for group in combinations(range(len(p)),len(p)-6):
        rem = [i for i in range(len(p)) if i not in group]
        spectator = add(*(p[i] for i in group))
        weight = dot(BULK, quartics([p[i] for i in rem] + [spectator]))
        for i in rem:
            for k in range(5):
                out[i][k] += weight*spectator[k]
    return out


@cache
def sorted_response(rows):
    # Labels of equal momentum have equal integrated potential by slot permutation invariance.
    phi = grounded_potentials(rows)
    return dict(zip(rows,phi))


def theta_response(rows):
    potential = sorted_response(tuple(sorted(rows)))
    return [potential[p] for p in rows]


def theta_terms(p):
    root_two = [[F()]*5 for _ in p]
    root_three = [[F()]*5 for _ in p]
    for rem in combinations(range(len(p)),7):
        spectator = tuple(-x for x in add(*(p[h] for h in rem)))
        phi = theta_response([p[h] for h in rem] + [spectator])
        for index,h in enumerate(rem):
            for k in range(5):
                root_two[h][k] += phi[index][k]-phi[-1][k]
    for rem in combinations(range(len(p)),8):
        phi = theta_response([p[h] for h in rem])
        for index,h in enumerate(rem):
            for k in range(5):
                root_three[h][k] += 2*phi[index][k]
    return root_two,root_three


def transverse(p, arrays):
    rows = []
    for i, v in enumerate(p):
        first = next(j for j,x in enumerate(v) if x)
        for j in range(5):
            if j!=first:
                rows.append([a[i][j]*v[first]-a[i][first]*v[j] for a in arrays])
    return rows


def terms(p):
    assert not any(add(*p))
    one_loop = general_one_loop(p)
    if len(p) == 9:
        assert one_loop == one_loop_on_primary_mark(p)
    return [*theta_terms(p), one_loop, tree_bulk(p)]


QUARTIC_MARK = (F(43,207360),F(1,69120),-F(1,311040),F(1,124416))


if __name__ == '__main__':
    from flint import fmpq_mat
    from quartic_unisolvence import FIXTURES
    all_rows = []
    for p in FIXTURES:
        arrays = [*fused_marked(p),*terms(p)]
        rows = transverse(p,arrays)
        assert all(dot(QUARTIC_MARK,r[:4])+sum(r[4:])==0 for r in rows)
        all_rows.extend(r[:4] for r in rows)
    assert fmpq_mat([[str(x) for x in r] for r in all_rows]).rref()[1] == 4
    print('All 108 transverse quartic equations vanish; the scalar boundary has rank four.')
    print('Signed marked primitive:',*QUARTIC_MARK)
