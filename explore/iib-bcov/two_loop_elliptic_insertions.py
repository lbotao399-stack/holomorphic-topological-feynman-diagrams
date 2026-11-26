"""Elliptic vertices in the all-arity genus-two fivefold Ward source."""
from fractions import Fraction
from itertools import combinations
from math import factorial, lcm, prod

from elliptic_descendants import blocks, elliptic_psi, vertex
from genus_one_graph_ward import circle_groups, classical
from two_loop_descendant_hexagon import add_rows, det


def mark_partitions(levels):
    """Ordinary blocks and moving-block psi weight in the genus-one vertex."""
    labels = tuple(range(len(levels)))
    for partition in blocks(labels,6,levels,1):
        ordinary = []
        moving = None
        weight = Fraction(1)
        for group in partition:
            k = tuple(sorted(levels[i] for i in group))
            degree = sum(k)
            if degree==len(group)-1:
                ordinary.append(group)
                weight *= Fraction(factorial(degree),prod(factorial(v) for v in k))
            elif degree==len(group) and moving is None:
                moving = k
            else:
                break
        else:
            if moving is not None:
                yield tuple(ordinary),weight*elliptic_psi(moving)


def single_quadratic(root, root_levels, singletons):
    """Cauchy-Binet and a rank-one update give the constant and quadratic current moments."""
    total = add_rows(singletons)
    primary = [p for p,k in zip(root,root_levels) if k==0]
    # Twice the Gram matrix of the exterior primaries and the two centered currents.
    gram = [[2*sum(p[i]*p[j] for p in primary)+total[i]*total[j] for j in range(5)]
            for i in range(5)]
    # The cofactor is polynomial and includes singular Gram matrices.
    cofactor = [[0]*5 for _ in range(5)]
    for i in range(5):
        for j in range(i+1):
            minor = tuple(tuple(x for l,x in enumerate(row) if l!=j) for h,row in enumerate(gram) if h!=i)
            cofactor[i][j] = cofactor[j][i] = (-1)**(i+j)*det(minor)
    contracted = [[sum(x*y for x,y in zip(p,col)) for col in zip(*cofactor)] for p in singletons]
    quadratic = [[Fraction(sum(x*y for x,y in zip(row,q)),192) for q in singletons] for row in contracted]
    return Fraction(det(gram),768),quadratic


def integrated_vertex(root, root_levels, singletons):
    """Quadratic Green moments of the elliptic vertex at the rooted circle point."""
    volume = det(singletons)
    if not volume:
        return [(Fraction(0),)*5]*5
    if sum(k>0 for k in root_levels)==1:
        constant,quadratic = single_quadratic(root,root_levels,singletons)
    else:
        m = len(root)
        total = add_rows(singletons)
        # Center the two internal currents at total/2, with opposite displacements.
        doubled = [tuple(2*x for x in p) for p in root]+[total,total]
        constant = Fraction(0)
        quadratic = [[Fraction(0)]*5 for _ in range(5)]
        for ordinary,weight in mark_partitions((*root_levels,0,0)):
            rows = [add_rows(doubled[i] for i in group) for group in ordinary]
            constant += weight*Fraction(det(rows)**2,1024)
            linear = [Fraction(0)]*5
            for i,group in enumerate(ordinary):
                direction = int(m in group)-int(m+1 in group)
                if direction:
                    for j,p in enumerate(singletons):
                        replacement = (*rows[:i],p,*rows[i+1:])
                        linear[j] += Fraction(direction*det(replacement),16)
            for i,a in enumerate(linear):
                for j,b in enumerate(linear):
                    quadratic[i][j] += weight*a*b
    size = sum(quadratic[i][i] for i in range(5))
    kernel = [[(constant/6+size/72-quadratic[i][i]/180 if i==j else
                constant/12+size/144-(quadratic[i][i]+quadratic[j][j]-quadratic[i][j])/360)
               for j in range(5)] for i in range(5)]
    return [tuple(Fraction(volume**2,2)*sum(t*x for t,x in zip(row,coordinate))
                  for coordinate in zip(*singletons)) for row in kernel]


def _integer_input(momenta, levels):
    if len(momenta)!=len(levels) or any(len(p)!=5 for p in momenta) or not momenta or any(add_rows(momenta)):
        raise ValueError('Expected conserved fivefold momenta and one level per wave.')
    if any(type(k) is not int or k<0 for k in levels) or sum(levels)!=len(levels)-8:
        raise ValueError('The genus-two Ward source at n waves has total level n-8.')
    scale = lcm(*(Fraction(x).denominator for p in momenta for x in p))
    return tuple(tuple(int(scale*Fraction(x)) for x in p) for p in momenta),scale


def _groupings(labels,levels):
    positive = [i for i in labels if levels[i]]
    if len(positive)==1:
        for singles in combinations((i for i in labels if i!=positive[0]),5):
            yield tuple(i for i in labels if i not in singles),tuple((i,) for i in singles)
        return
    for size in range(4,len(labels)-4):
        for group in combinations(labels,size):
            if sum(levels[i] for i in group)!=size-3:
                continue
            others = tuple(i for i in labels if i not in group)
            for ordinary in circle_groups(others,levels,5):
                yield group,ordinary


def hexagon(momenta, levels):
    """The single positive-genus bulk insertion in the one-loop hexagon at any arity."""
    p,scale = _integer_input(momenta,levels)
    labels = tuple(range(len(levels)))
    potential = [[Fraction(0)]*5 for _ in labels]
    for group,ordinary in _groupings(labels,levels):
        q = tuple(add_rows(p[i] for i in block) for block in ordinary)
        classical_weight = prod(classical(levels,block) for block in ordinary)
        root = tuple(p[i] for i in group)
        rows = integrated_vertex(root,tuple(levels[i] for i in group),q)
        for block,row in zip(ordinary,rows):
            for i in block:
                for j,x in enumerate(row):
                    potential[i][j] += classical_weight*x
    return [tuple(x/scale**21 for x in row) for row in potential]


def single_mark(momenta, descendant, level=None):
    """Single A_k source with primary waves elsewhere. n=8+k, including k=0."""
    if not 0<=descendant<len(momenta):
        raise ValueError('The distinguished slot must be a wave label.')
    if level is None:
        level=len(momenta)-8
    return hexagon(momenta,tuple(level if i==descendant else 0 for i in range(len(momenta))))


def elliptic_tree(momenta, levels):
    """Two elliptic vertices joined by one edge, grounding the component containing wave zero."""
    p,scale = _integer_input(momenta,levels)
    n = len(levels)
    result = [[Fraction(0)]*5 for _ in levels]
    if sum(k>0 for k in levels)<2:
        return [tuple(row) for row in result]
    for size in range(5,n-4):
        for group in combinations(range(1,n),size):
            if sum(levels[i] for i in group)!=size-4:
                continue
            other = tuple(i for i in range(n) if i not in group)
            charge = add_rows(p[i] for i in group)
            left = vertex([*(p[i] for i in group),tuple(-x for x in charge)],
                          [*(levels[i] for i in group),0])
            if not left:
                continue
            right = vertex([*(p[i] for i in other),charge],[*(levels[i] for i in other),0])
            for i in group:
                for j,x in enumerate(charge):
                    result[i][j] += left*right*x
    return [tuple(x/scale**21 for x in row) for row in result]


def independent_group_test(root, levels, singletons):
    """Interpolate the complete elliptic vertex in the two currents, then use raw cube moments."""
    from marked_theta import green_moments
    center = tuple(Fraction(x,2) for x in add_rows(singletons))
    def at(displacement):
        first = tuple(x+y for x,y in zip(center,displacement))
        second = tuple(x-y for x,y in zip(center,displacement))
        return vertex((*root,first,second),(*levels,0,0))
    zero = (0,)*5
    constant = at(zero)
    axis = [at(p)-constant for p in singletons]
    terms = {(i,i):axis[i] for i in range(5)}
    for i,j in combinations(range(5),2):
        terms[i,j] = at(add_rows((singletons[i],singletons[j])))-constant-axis[i]-axis[j]
    matrix = [[Fraction(0)]*5 for _ in range(5)]
    for i in range(5):
        for j in range(5):
            flat = green_moments(zero,i,j)[1]
            result = constant*flat
            for (k,l),coef in terms.items():
                e = tuple(int(a==k) for a in range(5))
                f = tuple(int(a==l) for a in range(5))
                both = add_rows((e,f))
                moment = (green_moments(both,i,j)[1]
                          -(green_moments(e,i,j)[1]+green_moments(f,i,j)[1])/2+flat/4)
                result += coef*moment
            matrix[i][j] = result
    expected = [tuple(Fraction(det(singletons)**2,2)*sum(a*b for a,b in zip(row,col))
                      for col in zip(*singletons)) for row in matrix]
    assert integrated_vertex(root,levels,singletons)==expected


def check():
    from random import Random
    from two_loop_descendant_hexagon import potentials
    random = Random(78511)
    q = tuple(tuple(int(i==j) for j in range(5)) for i in range(5))
    p = []
    for _ in range(6):
        p.append(tuple(random.randrange(-2,3) for i in range(5)))
    def root(count):
        return (*p[:count-1],tuple(-x for x in add_rows((*p[:count-1],*q))))
    for levels in ((1,0,0,0),(0,0,0,0,2),(0,1,0,1,0),(0,1,0,0,2,0)):
        independent_group_test(root(len(levels)),levels,q)
    nine = (*root(4),*q)
    base = single_mark(nine,0)
    assert base==potentials(nine,0)
    assert any(any(row) for row in base)
    extended = (*nine,(0,)*5)
    constant_dilaton = (1,)+(0,)*8+(1,)
    assert hexagon(extended,constant_dilaton)[:9]==[tuple(11*x for x in row) for row in base]
    assert not any(any(row) for row in elliptic_tree(extended,constant_dilaton))
    assert single_mark(extended,0)[:9]==base
    mixed_levels = (1,1)+(0,)*8
    mixed = hexagon(extended,mixed_levels)
    second = single_mark(nine,1)
    assert mixed[:9]==[add_rows((a,b)) for a,b in zip(base,second)]
    ten = (*root(5),*q)
    # Separate the two first descendants or keep both on one elliptic vertex.
    diag = elliptic_tree(ten,mixed_levels)
    assert any(any(row) for row in diag)
    assert not any(diag[0])
    eleven = (*ten,(0,)*5)
    plus_dilaton = (*mixed_levels,1)
    tree = elliptic_tree(eleven,plus_dilaton)
    assert tree[:10]==[tuple(12*x for x in row) for row in diag]
    # The tree gauge is set by the root at label zero; compare differences under relabeling.
    order=(9,2,1,6,5,7,0,4,3,8)
    perm = elliptic_tree(tuple(ten[i] for i in order),tuple(mixed_levels[i] for i in order))
    for j,i in enumerate(order):
        assert perm[j]==tuple(a-b for a,b in zip(diag[i],diag[order[0]]))
    assert not any(any(row) for row in elliptic_tree(ten,(2,)+(0,)*9))
    print('All quadratics matched directly interpolated elliptic vertices in the two internal currents.')
    print('Nine-wave graph, dilaton factor 11, both string limits and tree dilaton factor 12 passed.')


if __name__=='__main__':
    check()
