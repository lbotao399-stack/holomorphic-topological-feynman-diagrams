"""Genus-two differential of a normal-ordered scalar point polynomial in minimal fivefold BCOV.

The three kernels multiply the ordinary derivatives f', f'' and f''' of the original
point observable, with common physical factor lambda^4/(2*pi*i)^10. External scalar
waves may have levels zero to four. The one external primary vector is transverse.
"""
from fractions import Fraction
from itertools import combinations
from math import lcm, prod

from genus_one_graph_ward import circle_groups, classical, plus
from genus_two_graph_ward import ZERO, branches, site_rows, theta_groups, vertex
from longitudinal_theta import closed_potentials
from nine_wave_theta import endpoint_potentials
from two_loop_elliptic_insertions import _groupings, integrated_vertex


def voltage(rows, groups, vector, point, polarization):
    vector_site = next(i for i,g in enumerate(groups) if vector in g)
    point_site = next(i for i,g in enumerate(groups) if point in g)
    return sum((a-b)*v for a,b,v in zip(rows[vector_site],rows[point_site],polarization))


def square_theta(p, levels, vector, point, polarization):
    result = Fraction(0)
    for ordinary,left,right in theta_groups(levels):
        if (point,) not in ordinary:
            continue
        groups = (*ordinary,left,right)
        charges = tuple(plus(p[i] for i in group) if group else ZERO for group in groups)
        if any(not any(q) for q in charges[:8]):
            continue
        if not left and not right:
            rows = (*closed_potentials(charges[:8]),ZERO,ZERO)
        elif not left:
            rows = (*endpoint_potentials(charges[:8]),ZERO,ZERO)
        else:
            rows = site_rows(charges[:8],charges[8],charges[9])
        weight = (1 if left==right else 2)*prod(classical(levels,g) for g in groups)
        result += weight*voltage(rows,groups,vector,point,polarization)
    return result


def square_hexagon(p, levels, vector, point, polarization):
    result = Fraction(0)
    labels = tuple(range(len(levels)))
    for bulk,ordinary in _groupings(labels,levels):
        if (point,) not in ordinary:
            continue
        single = tuple(plus(p[i] for i in group) for group in ordinary)
        if any(not any(q) for q in single):
            continue
        root = tuple(p[i] for i in bulk)
        rows = (*integrated_vertex(root,tuple(levels[i] for i in bulk),single),ZERO)
        weight = prod(classical(levels,group) for group in ordinary)
        result += weight*voltage(rows,(*ordinary,bulk),vector,point,polarization)
    return result


def cube_theta(p, levels, vector, point, polarization):
    result = Fraction(0)
    external = tuple(range(point))
    for other in branches(external,levels,len(external)-8):
        keep = tuple(i for i in external if i not in other)
        tip = plus(p[i] for i in other) if other else ZERO
        for ordinary in circle_groups(keep,levels,8):
            charges = tuple(plus(p[i] for i in group) for group in ordinary)
            if any(not any(q) for q in charges):
                continue
            groups = (*ordinary,other,(point,))
            if not other:
                rows = (*endpoint_potentials(charges),ZERO,ZERO)
            else:
                rows = site_rows(charges,tip,p[point])
            # Fixing the operator at one branch removes the usual exchange of theta endpoints.
            weight = 2*prod(classical(levels,group) for group in (*ordinary,other))
            result += weight*voltage(rows,groups,vector,point,polarization)
    return result


def kernel_terms(scalars, levels, vector, polarization):
    """Return (action tree, classical theta, elliptic hexagon) for the unique active order."""
    if len(scalars)!=len(levels) or any(len(p)!=5 for p in (*scalars,vector,polarization)):
        raise ValueError('Supply fivefold momenta, levels of all scalar waves, and a vector polarization.')
    if any(type(k) is not int or not 0<=k<=4 for k in levels):
        raise ValueError('Fivefold scalar descendant levels range from zero to four.')
    if sum(x*v for x,v in zip(vector,polarization)):
        raise ValueError('The external primary vector must be transverse.')
    order = len(scalars)-4-sum(levels)
    if order not in (1,2,3):
        return (Fraction(0),)*3
    # Include the vector in the integer momentum scale; the scalar polarization stays unchanged.
    waves = (*scalars,vector)
    denominator = lcm(*(Fraction(x).denominator for p in waves for x in p))
    p = tuple(tuple(int(denominator*Fraction(x)) for x in row) for row in waves)
    point = len(p)
    vslot = point-1
    root = tuple(-x for x in plus(p))
    p += (root,)
    k = (*levels,0,0)
    tree = theta = circle = Fraction(0)
    if order==1:
        tree = -sum(x*v for x,v in zip(root,polarization))*vertex(p,k)
    elif order==2:
        theta = square_theta(p,k,vslot,point,polarization)
        circle = square_hexagon(p,k,vslot,point,polarization)
    else:
        theta = cube_theta(p,k,vslot,point,polarization)
    return tuple(x/denominator**21 for x in (tree,theta,circle))


def kernel(scalars, levels, vector, polarization):
    """Active coefficient of f^(k) for k=number of scalar waves - 4 - total descendant level."""
    return sum(kernel_terms(scalars,levels,vector,polarization),Fraction(0))


def polynomial_coefficient(scalars, levels, vector, polarization):
    """Full polarized value on phi^r/r! at zero background, including undifferentiated spectators.

    With r=len(scalars)-4-sum(levels), multiplication by f^(r)(0) gives the value for f.
    """
    degree = len(scalars)-4-sum(levels)
    primary = tuple(i for i,k in enumerate(levels) if k==0)
    out = Fraction(0)
    for order in range(1,min(degree,3)+1):
        for spectators in combinations(primary,degree-order):
            keep = tuple(i for i in range(len(scalars)) if i not in spectators)
            out += kernel(tuple(scalars[i] for i in keep),tuple(levels[i] for i in keep),
                          vector,polarization)
    return out


def one_descendant_check(primary, d, polarization, order):
    """Replay the separate marked topologies with the original integrated heat-kernel tables."""
    from bulk_eight import grounded_potentials, potentials
    from theta_endpoint_sources import theta_potentials
    from two_loop_descendant_hexagon import group_potentials
    end = tuple(-x for x in plus((*primary,d)))
    vector = len(primary)-1
    def paired(a,b=ZERO):
        return sum((x-y)*v for x,y,v in zip(a,b,polarization))
    theta = circle = Fraction(0)
    for i in range(len(primary)):
        merged = tuple(plus((p,d)) if i==j else p for j,p in enumerate(primary))
        if order==2:
            rows = potentials((*merged,end))
            theta += paired(rows[vector],rows[-1])
        else:
            theta += 2*paired(grounded_potentials(merged)[vector])
    if order==2:
        rows = grounded_potentials((*primary,end))
        theta += 2*paired(rows[vector],rows[-1])
        for group in combinations(range(len(primary)),3):
            outside = [i for i in range(len(primary)) if i not in group]
            rows = group_potentials((d,*(primary[i] for i in group)),(*(primary[i] for i in outside),end))
            source = rows[outside.index(vector)] if vector in outside else ZERO
            circle += paired(source,rows[-1])
    else:
        # Reverse the endpoints and ground the descendant, to check the marked-root gauge.
        rows = theta_potentials(primary,end,d)
        theta += 2*paired(rows[vector],rows[8])
    return theta,circle


def check():
    from longitudinal_theta import primary_coefficient
    from two_loop_moving_descendant import BASIS as b
    descendant = (0,-1,1,0,1)
    other = (1,0,0,-1,1)
    fixtures = (
        ((0,0,1,2,3),(0,0,0,-1,0),Fraction(-1,31104),Fraction(-79,12960)),
        ((0,0,1,1,2,3),(1,1,1,1,0),Fraction(1,405),Fraction(625,2592)),
        ((0,0,1,1,2,2,3),(1,1,1,1,0),Fraction(7,720),Fraction(139,270)))
    for order,(indices,polar,old,moving) in enumerate(fixtures,1):
        scalar = tuple(b[i] for i in indices)
        levels = (0,)*len(scalar)
        assert kernel(scalar,levels,b[4],polar)==primary_coefficient(scalar,b[4],polar)==old
        assert kernel((*scalar,ZERO),(*levels,1),b[4],polar)==9*old
        terms = kernel_terms((*scalar,descendant),(*levels,1),b[4],polar)
        assert sum(terms)==moving
        print('Order',order,'moving A1 tree, theta, elliptic hexagon:',*terms,flush=True)
        if order==1:
            continue
        assert terms[1:]==one_descendant_check((*scalar,b[4]),descendant,polar,order)
        full = polynomial_coefficient((*scalar,descendant),(*levels,1),b[4],polar)
        assert full=={2:Fraction(1745,5184),3:Fraction(20251,9720)}[order]
        double = kernel((*scalar,descendant,other),(*levels,1,1),b[4],polar)
        assert double=={2:Fraction(5318,243),3:Fraction(455621,9720)}[order]
        assert double==kernel((other,*scalar,descendant),(1,*levels,1),b[4],polar)
        assert kernel((*scalar,descendant,ZERO),(*levels,1,1),b[4],polar)==10*moving
        assert kernel((*scalar,descendant,ZERO),(*levels,2,0),b[4],polar)==moving
        print('Order',order,'full polarized monomial',full,'two moving descendants',double,flush=True)
    print('Rooted topology, primary, moving-permutation, string and dilaton tests passed.')


if __name__=='__main__':
    check()
