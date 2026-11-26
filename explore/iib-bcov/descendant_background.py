"""Evaluate primary-background and descendant wave profiles against an ordinary vertex bank."""
from collections import Counter
from fractions import Fraction
from itertools import combinations_with_replacement, product
from math import comb, factorial, prod

from elliptic_descendants import vertex as elliptic
from genus_two_graph_ward import vertex as two_loop
from scalar_ancestor_response import scalar_power
from betti_flow_polytope import weak


def terms(counts, levels):
    """Distribute indistinguishable transverse waves among retained descendant fields."""
    def visit(i, remaining, attached, weight):
        if i == len(levels):
            yield remaining, attached, weight / prod(map(factorial,remaining))
            return
        for entry in product(*(range(min(m,levels[i])+1) for m in remaining)):
            spent = sum(entry)
            if spent <= levels[i]:
                after = tuple(m-s for m,s in zip(remaining,entry))
                yield from visit(i+1,after,(*attached,entry),
                                 weight*(-1)**spent/prod(map(factorial,entry)))
    yield from visit(0, counts, (), Fraction(1))


def coefficient(vertex, charges, levels, counts, directions=None, shifts=None):
    """Background power coefficient, with primary root last; return the action normalization."""
    d = len(charges[0])
    if directions is None:
        directions = tuple(tuple(int(i==j) for i in range(5-d)) for j in range(5-d))
    if shifts is None:
        shifts = ((0,)*(5-d),)*(len(charges)-1)
    if (len(charges)!=len(levels) or levels[-1] or any(len(p)!=d for p in charges)
            or any(map(sum,zip(*charges))) or any(n<0 for n in (*counts,*levels))
            or len(directions)!=len(counts) or any(len(v)!=5-d for v in directions)
            or len(shifts)!=len(charges)-1 or any(len(v)!=5-d for v in shifts)):
        raise ValueError('Supply conserved projected wave momenta, legal powers, and a primary root.')
    def transfer(multiplicities):
        return tuple(sum(n*v[i] for n,v in zip(multiplicities,directions)) for i in range(5-d))
    totals = Counter()
    for remaining, attached, weight in terms(counts,levels[:-1]):
        hard = tuple((k-sum(f),(*p,*(x+y for x,y in zip(t,transfer(f)))))
                     for p,k,f,t in zip(charges,levels,attached,shifts))
        free = tuple((0, (0,)*d+v) for v,n in zip(directions,remaining) for _ in range(n))
        recoil = tuple(-x-sum(col) for x,col in zip(transfer(counts),zip(*shifts)))
        root = ((0,(*charges[-1],*recoil)),)
        totals[tuple(sorted(hard+free+root))] += weight
    return sum(weight*vertex(tuple(p for k,p in slots),tuple(k for k,p in slots))
               for slots,weight in totals.items() if weight)


def check_hessian():
    from flint import fmpq, fmpq_mpoly_ctx, fmpz_mat
    from itertools import combinations
    for genus,vertex,charges,levels,directions in (
            (1,elliptic,((1,0),(0,1),(1,-1),(-2,0)),(1,1,0,0),
             ((1,0,0),(0,1,0),(0,0,1),(1,-1,2))),
            (2,two_loop,((1,0,0),(0,1,0),(0,0,1),(1,-1,0),(-2,0,-1)),(1,1,0,0,0),
             ((1,0),(0,1),(1,-2)))):
        dimension = len(directions[0])
        ctx = fmpq_mpoly_ctx.get(tuple('u'+str(i) for i in range(len(directions))))
        axes = ctx.gens()
        determinant = sum(int(fmpz_mat([directions[i] for i in group]).det())**2
                          *prod(axes[i] for i in group)
                          for group in combinations(range(len(directions)),dimension))
        expected = (determinant**genus).to_dict()
        ground = scalar_power(genus,len(charges[0])*genus,charges,levels)
        shifts = ((1,)*dimension,)+((0,)*dimension,)*(len(charges)-2)
        response = ctx.constant(0)
        for counts in weak(genus*dimension,len(directions)):
            direct = coefficient(vertex,charges,levels,counts,directions,shifts)
            assert direct == ground*Fraction(str(expected.get(counts,0))), (counts,direct)
            # Root contraction against the first transverse coordinate vector.
            root = -1-sum(n*v[0] for n,v in zip(counts,directions))
            response += fmpq(str(direct*root))*prod(t**n for t,n in zip(axes,counts))
        power = determinant**genus
        derivative = sum(v[0]*t*power.derivative(i) for i,(v,t) in enumerate(zip(directions,axes)))
        assert response == fmpq(str(ground))*(-power-derivative)
        print('Fivefold graph check with oblique transverse waves:',genus,determinant,flush=True)


def check():
    for g,vertex,charges,levels in (
            (1,elliptic,((1,0),(0,1),(1,-1),(-2,0)),(1,1,0,0)),
            (2,two_loop,((1,0,0),(0,1,0),(0,0,1),(1,-1,0),(-2,0,-1)),(2,0,0,0,0)),
            (2,two_loop,((1,0,0),(0,1,0),(0,0,1),(1,-1,0),(-2,0,-1)),(1,1,0,0,0))):
        d=len(charges[0])
        ground = scalar_power(g,d*g,charges,levels)
        for counts in weak(g*(5-d),5-d):
            measured=coefficient(vertex,charges,levels,counts)
            assert measured == (ground if all(m==g for m in counts) else 0), (counts,measured,ground)
        print('Complete transverse-amplitude polynomial:',g,levels,ground,flush=True)


def soft_root_lowerings(levels):
    for j, level in enumerate(levels):
        if level:
            yield tuple(x - int(i == j) for i, x in enumerate(levels))


def soft_root_ranks():
    from flint import fmpz_mat

    triples = [tuple(reversed(row)) for row in combinations_with_replacement(range(5), 3)]
    relations = {1: (-2, 1), 2: (3, -3, 1), 3: (2, -2, 1, 0),
                 5: (12, -6, -9, 6, -2)}
    for genus in range(1, 11):
        inputs = sorted((row for row in triples if sum(row) == genus + 1), reverse=True)
        vertices = sorted(row for row in triples if sum(row) == genus)
        columns = [Counter(tuple(sorted(row, reverse=True)) for row in soft_root_lowerings(levels))
                   for levels in inputs]
        matrix = fmpz_mat([[counts[v] for counts in columns] for v in vertices])
        assert matrix.rank() == len(inputs) - int(genus in relations)
        if genus in relations:
            coeff = relations[genus]
            assert len(coeff) == len(inputs)
            assert all(sum(c * counts[v] for c, counts in zip(coeff, columns)) == 0 for v in vertices)
        print('Minimal soft roots at genus', genus, 'inputs', len(inputs), 'rank', matrix.rank())


def soft_root_check():
    # The primary momentum is t*e1. The signed root contraction already
    # contains t, so its linear coefficient uses the remaining vertex at zero.
    hard = ((1, 0), (0, 1), (-1, -1))
    # Shifting one scalar along a single recoil direction has degree at most
    # four in its amplitude: change planar basis to put the shift on one axis.
    for levels, expected in [((3, 0, 0), Fraction(1, 3456)),
                             ((2, 1, 0), Fraction(1, 1080)),
                             ((1, 1, 1), Fraction(11, 5760))]:
        direct = coefficient(two_loop, (*hard, (0, 0)), (*levels, 0), (2, 2, 2))
        assert direct == expected == sum(scalar_power(2, 4, hard, row) for row in soft_root_lowerings(levels))
        for p, q, recoil in [(hard[0], hard[1], (1, 0)), ((2, 1), (1, 2), (-1, 1))]:
            area = p[0]*q[1] - p[1]*q[0]
            root = recoil[0]*q[1] - recoil[1]*q[0]
            samples = []
            for t in range(1, 6):
                soft = tuple(t*x for x in recoil)
                moving = tuple(-x-y-z for x, y, z in zip(p, q, soft))
                samples.append(scalar_power(2, 4, (p, q, moving, soft), (*levels, 0)))
            at_zero = sum((-1)**(k-1)*comb(5, k)*v for k, v in enumerate(samples, 1))
            assert -root*at_zero == -root*area**4*direct
        print('Genus-two signed unit soft root, levels', levels, ':', -direct)
    # Complete genus-five three-mark sums, grouped in edge ranges 0-8, 9, 10.
    subtotal = (Fraction(24942350470871679106273, 3221948086174512414720000000),
                -Fraction(153185788543248886281611, 6443896172349024829440000000),
                Fraction(2295582284289643992637, 143197692718867218432000000))
    assert 3*sum(subtotal) == Fraction(51787, 965667225600000)
    vertex = {(2, 2, 1): sum(subtotal), (4, 1, 0): Fraction(14951, 6083703521280000),
              (3, 1, 1): Fraction(5783, 506975293440000), (3, 2, 0): Fraction(463, 86910050304000)}
    for levels, expected in [((4, 2, 0), Fraction(15787, 2027901173760000)),
                             ((4, 1, 1), Fraction(49649, 3041851760640000)),
                             ((3, 3, 0), Fraction(463, 43455025152000)),
                             ((3, 2, 1), Fraction(191417, 5530639564800000)),
                             ((2, 2, 2), Fraction(51787, 965667225600000))]:
        measured = sum(vertex[tuple(sorted(row, reverse=True))] for row in soft_root_lowerings(levels))
        assert measured == expected
        print('Genus-five signed unit soft root, levels', levels, ':', measured)
    soft_root_ranks()


if __name__=='__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--soft-root', action='store_true', help='replay the soft primary at genera two and five')
    if parser.parse_args().soft_root:
        soft_root_check()
    else:
        check()
        check_hessian()
