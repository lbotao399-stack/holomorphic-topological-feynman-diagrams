"""Labeled genus-two flat fivefold descendant vertices from the graph Ward equation.

The positive-geometric-loop graphs are the four theta shapes and one elliptic
vertex on a hexagon. Two loops meeting only at a four-valent vertex have an
independent resistance scale on each circuit, so their projective heat form
vanishes. The one-edge trees join genus zero to genus two or genus one to
genus one. Every genus-two vertex on a tree has smaller arity.
"""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from math import comb, lcm, prod

from flint import fmpz_mat
from genus_one_graph_ward import circle_groups, classical, plus
from longitudinal_theta import closed_potentials
from nine_wave_theta import endpoint_potentials
from theta_endpoint_sources import theta_potentials
from two_loop_elliptic_insertions import elliptic_tree, hexagon
from two_loop_single_descendant import primary, single_vertex

ZERO = (0,)*5


def branches(labels,levels,maximum):
    yield ()
    for size in range(1,maximum+1):
        for group in combinations(labels,size):
            if sum(levels[i] for i in group)==size:
                yield group


def theta_groups(levels):
    """Eight unordered bivalent groups and two unordered branch decorations."""
    labels = tuple(range(len(levels)))
    spare = len(labels)-8
    for left in branches(labels,levels,spare):
        other = tuple(i for i in labels if i not in left)
        for right in branches(other,levels,spare-len(left)):
            if right<left:
                continue
            mid = tuple(i for i in other if i not in right)
            for ordinary in circle_groups(mid,levels,8):
                yield ordinary,left,right


@lru_cache(maxsize=256)
def site_rows(waves,left,right):
    # Reversal exchanges the charged endpoints and translates all potentials
    # by the same covector. The tables already divide by the reversal order.
    return tuple(map(tuple,theta_potentials(waves,left,right)))


def theta_rows(momenta,levels):
    result = [[F() for _ in range(5)] for _ in levels]
    for ordinary,left,right in theta_groups(levels):
        groups = (*ordinary,left,right)
        charges = tuple(plus(momenta[i] for i in group) if group else ZERO for group in groups)
        # A zero total bivalent charge leaves its marked-position column zero.
        if any(not any(q) for q in charges[:8]):
            continue
        if not left and right:
            # Ground the occupied branch, so all external rows are supplied
            # by the eight-site one-endpoint integration.
            rows = (*endpoint_potentials(charges[:8]),ZERO,ZERO)
        elif not left and not right:
            rows = (*closed_potentials(charges[:8]),ZERO,ZERO)
        else:
            rows = site_rows(charges[:8],charges[8],charges[9])
        weight = (1 if left==right else 2)*prod(classical(levels,g) for g in groups)
        for group,row in zip(groups,rows):
            for i in group:
                for j,x in enumerate(row):
                    result[i][j] += weight*x
    return result


def tree_groups(levels):
    labels = tuple(range(len(levels)))
    for size in range(2,len(labels)-5):
        for group in combinations(labels,size):
            if sum(levels[i] for i in group)==size-2:
                yield group


def tree_slots(momenta,levels,group,charge):
    slots = [(k,p) for i,(k,p) in enumerate(zip(levels,momenta)) if i not in group]
    slots.append((0,charge))
    return tuple(sorted(slots))


def tree_rows(momenta,levels):
    """The genus-zero end of a tree attached to the local genus-two vertex."""
    result = [[F() for _ in range(5)] for _ in levels]
    for group in tree_groups(levels):
        q = plus(momenta[i] for i in group)
        if not any(q):
            continue
        weight = classical(levels,group)*integer_vertex(tree_slots(momenta,levels,group,q))
        for i in group:
            for j,x in enumerate(q):
                result[i][j] += weight*x
    return result


def divergences(momenta,rows):
    # At a zero wave the divergence vanishes. Its raised value requires a
    # separate polynomial evaluation; it does not affect the other slots.
    for i,j in combinations(range(len(momenta)),2):
        for x,y in combinations(range(5),2):
            p,q = momenta[i],momenta[j]
            divisor = p[x]*q[y]-p[y]*q[x]
            if not divisor:
                continue
            value = ((rows[i][x]-rows[j][x])*q[y]-(rows[i][y]-rows[j][y])*q[x])/divisor
            common = tuple(a-value*b for a,b in zip(rows[i],p))
            result = []
            for wave,row in zip(momenta,rows):
                axis = next((a for a,x in enumerate(wave) if x),None)
                coeff = (row[axis]-common[axis])/wave[axis] if axis is not None else None
                assert all(r-c==(coeff or 0)*b for r,c,b in zip(row,common,wave))
                result.append(coeff)
            return tuple(result)
    return (F(),)*len(momenta)


@lru_cache(maxsize=128)
def recovered(momenta,levels):
    if sum(levels)!=len(momenta)-8:
        raise ValueError('The n-wave genus-two Ward source has total level n-8')
    rows = (theta_rows(momenta,levels),tree_rows(momenta,levels),
            hexagon(momenta,levels),elliptic_tree(momenta,levels))
    return divergences(momenta,[plus(parts) for parts in zip(*rows)])


def projected_ward_difference(momenta,levels,first,second,polarization):
    """Return (Ward row[first] - Ward row[second]) on one polarization."""
    if len(momenta)!=len(levels) or sum(levels)!=len(levels)-8 or len(polarization)!=5:
        raise ValueError('A fivefold n-wave Ward source has total level n-8.')
    components = [(j,x) for j,x in enumerate(polarization) if x]
    if first==second or not components:
        return F()
    tree_weights = Counter()
    for group in tree_groups(levels):
        sign = int(first in group)-int(second in group)
        if not sign:
            continue
        q = plus(momenta[i] for i in group)
        projected = sum(x*q[j] for j,x in components)
        if projected:
            slots = tree_slots(momenta,levels,group,q)
            tree_weights[slots] += sign*projected*classical(levels,group)
    value = sum((w*integer_vertex(slots) for slots,w in tree_weights.items() if w),F())
    for diagram in (theta_rows,hexagon,elliptic_tree):
        rows = diagram(momenta,levels)
        value += sum(x*(rows[first][j]-rows[second][j]) for j,x in components)
    return value


def raised_divergence(momenta,levels,raised):
    # A nonparallel anchor with high level removes more descendants from the
    # contributing tree cuts. A two-coordinate polarization kills its own row
    # after subtraction of the common potential.
    p = momenta[raised]
    for anchor in sorted(range(len(levels)),key=levels.__getitem__,reverse=True):
        q = momenta[anchor]
        for x,y in combinations(range(5),2):
            determinant = p[x]*q[y]-p[y]*q[x]
            if determinant:
                v = [0]*5
                v[x],v[y] = q[y],-q[x]
                return projected_ward_difference(momenta,levels,raised,anchor,v)/determinant
    return F()


@lru_cache(maxsize=2048)
def integer_vertex(slots):
    levels = tuple(k for k,_ in slots)
    momenta = tuple(p for _,p in slots)
    if fmpz_mat(momenta).rank()<5:
        return F()
    if not any(levels):
        return primary(momenta)
    for i,(k,p) in enumerate(slots):
        if any(p) or k>1:
            continue
        rest = slots[:i]+slots[i+1:]
        if k==1:
            return (len(slots)+1)*integer_vertex(rest)
        return sum(integer_vertex(tuple(sorted((*rest[:j],(s-1,q),*rest[j+1:]))))
                   for j,(s,q) in enumerate(rest) if s)
    positive = [i for i,k in enumerate(levels) if k]
    if len(positive)==1:
        return single_vertex(momenta,positive[0])
    # Lower a nonzero marked wave, prioritizing the largest level; partition
    # symmetry permits sorting all slots at every lower-arity tree cut.
    moving = [i for i in positive if any(momenta[i])]
    if moving:
        i = max(moving,key=levels.__getitem__)
        lower = tuple(k-int(j==i) for j,k in enumerate(levels))
        return raised_divergence(momenta,lower,i)
    # Along a balanced two-slot deformation every five-minor is affine, and
    # the genus-two trace has four five-minors. Five samples recover t=0.
    i = positive[0]
    j = next(j for j,p in enumerate(momenta) if any(p))
    axis = next(a for a in range(5) if any(x for b,x in enumerate(momenta[j]) if a!=b))
    result = F()
    for t in range(1,6):
        step = tuple(t*int(a==axis) for a in range(5))
        shifted = list(slots)
        shifted[i] = levels[i],step
        shifted[j] = levels[j],tuple(x-y for x,y in zip(momenta[j],step))
        result += (-1)**(t-1)*comb(5,t)*integer_vertex(tuple(sorted(shifted)))
    return result


def vertex(momenta,levels):
    """Polarized genus-two worldsheet action, arbitrary conserved rational fivefold waves."""
    if (len(momenta)!=len(levels) or len(momenta)<7 or
        any(len(p)!=5 for p in momenta) or any(plus(momenta)) or
        any(type(k) is not int or not 0<=k<=4 for k in levels) or sum(levels)!=len(levels)-7):
        raise ValueError('The n-wave genus-two fivefold vertex requires total level n-7 in levels 0 to 4')
    scale = lcm(*(F(x).denominator for p in momenta for x in p))
    slots = tuple(sorted((k,tuple(int(F(x)*scale) for x in p)) for k,p in zip(levels,momenta)))
    return integer_vertex(slots)/scale**20


def check():
    from level_two_scalar_reduction import momenta, vertex as saturated_vertex
    from two_loop_moving_descendant import BASIS, descendant_vertices

    hard,extra = momenta(1,1,(1,-1,1,-1))
    nine = hard+extra
    systems = []
    for i in range(9):
        lower = tuple(int(j==i) for j in range(9))
        system = recovered(nine,lower)
        assert system[i]==single_vertex(nine,i)
        systems.append(system)
    assert all(systems[i][j]==systems[j][i] for i,j in combinations(range(9),2))
    assert systems[4][4]==saturated_vertex(1,1,False)==F(649,77760)
    assert systems[4][2]==saturated_vertex(1,1,True)==F(421,38880)

    # Test string and dilaton on the unsimplified source, including the
    # unknown divergence of the zero wave itself where applicable.
    eight = (*BASIS,(1,2,3,4,5),(-2,1,3,-1,2))
    eight += (tuple(-x for x in plus(eight)),)
    levels = (1,)+(0,)*8
    assert recovered(eight,(0,)*8)==descendant_vertices(eight)
    const_primary = recovered((*eight,ZERO),levels)
    below = descendant_vertices(eight)
    assert const_primary[0]==below[0] and const_primary[-1] is None
    assert all(const_primary[i]==below[0]+below[i] for i in range(1,8))
    const_dilaton = recovered((*eight,ZERO),(0,)*8+(1,))
    assert const_dilaton[:-1]==tuple(10*x for x in below)
    assert vertex((*eight,ZERO),(0,)*8+(2,))==single_vertex((*eight,ZERO),8)

    ten = list(nine)+[BASIS[0]]
    ten[3] = tuple(x-y for x,y in zip(ten[3],ten[-1]))
    ten = tuple(ten)
    mixed = tuple(int(i in (2,4)) for i in range(10))
    level_two = tuple(2*int(i==4) for i in range(10))
    assert recovered(ten,mixed)[4]==recovered(ten,level_two)[2]==F(8449,38880)
    assert vertex(ten,tuple(a+int(i==4) for i,a in enumerate(mixed)))==F(8449,38880)
    # A second ten-wave fixture has a nonzero tree joining the two elliptic
    # action vertices; the single level-two Ward source has no such tree.
    generic = (*BASIS,(1,1,1,1,1),(2,1,-1,3,2),(0,1,2,1,3),(-2,3,1,-1,2))
    generic += (tuple(-x for x in plus(generic)),)
    mixed = (1,1)+(0,)*8
    assert elliptic_tree(generic,mixed)[1][0]==F(13613,72)
    assert recovered(generic,mixed)[0]==recovered(generic,(2,)+(0,)*9)[1]==F(2015618689,19440)
    full = recovered(generic,mixed)
    assert all(raised_divergence(generic,mixed,i)==full[i] for i in (0,1,3))
    v = (1,2,-1,0,3)
    assert projected_ward_difference(generic,mixed,2,4,v)==sum(
        x*(full[2]*generic[2][j]-full[4]*generic[4][j]) for j,x in enumerate(v))
    for raised in (0,1,2):
        lower = tuple(int(i<3 and i!=raised) for i in range(10))
        assert recovered(generic,lower)[raised]==F(13441892659,38880)
    print('Nine-wave systems:',9*9,'raised values and 36 cross-slot identities')
    print('Nine-wave moving level two and two level-one waves:',systems[4][4],systems[4][2])
    print('Ten-wave moving levels two and one on two fixtures:',F(8449,38880),F(2015618689,19440))
    print('Three moving first descendants at ten waves:',F(13441892659,38880))


if __name__=='__main__':
    check()
