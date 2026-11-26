"""Generic Betti rank and leading soft coefficients on rational wave configurations."""
from fractions import Fraction
from itertools import combinations
from random import Random

from flint import fmpq_mat

from primary_rank import partition_witness, rank_of


def subsets(size):
    for length in range(size+1):
        yield from combinations(range(size), length)


def independent_count(rows, genus):
    """Greedy basis of the union of g linear matroids; partition_witness augments exactly."""
    basis = []
    for row in rows:
        if partition_witness([*basis, row], genus)[0] is not None:
            basis.append(row)
    return len(basis)


def integral(rows):
    from math import lcm
    result = []
    for row in rows:
        p = tuple(map(Fraction, row))
        multiple = lcm(*(x.denominator for x in p))
        result.append(tuple(int(multiple*x) for x in p))
    return result


def betti_rank(genus, hard):
    if genus < 2 or len(hard) < 2 or any(map(sum, zip(*hard))):
        raise ValueError("Expected conserved marked momenta at genus at least two")
    rows = integral(hard)
    return min(genus*rank_of(rows, range(len(rows))), 3*genus-3+independent_count(rows, genus))


def quotient_rows(hard, soft):
    """Evaluate soft waves on an exact basis of the annihilator of the hard span."""
    dimension = len(hard[0])
    reduced, rank = fmpq_mat([list(row) for row in hard]).rref()
    pivots = [next(i for i in range(dimension) if reduced[j,i]) for j in range(rank)]
    free = [i for i in range(dimension) if i not in pivots]
    return integral([[Fraction(p[i])-sum(Fraction(str(reduced[j,i]))*p[b]
                                         for j,b in enumerate(pivots)) for i in free] for p in soft])


def sharp_soft(genus, hard, soft):
    """Test strict positivity of the first allowed coefficient, of degree 2 per soft slot."""
    if genus < 2 or not soft or len(hard) < 2:
        raise ValueError("Expected soft slots and at least two hard slots at genus at least two")
    dimension = len(hard[0])
    if (len(hard)+len(soft) != (dimension-3)*genus+3 or
            any(len(p) != dimension for p in [*hard,*soft]) or any(map(sum,zip(*hard)))):
        raise ValueError("Incompatible primary arity, momentum dimension, or conservation")
    hard, soft = integral(hard), integral(soft)
    if partition_witness([*hard,*soft], genus)[0] is None:
        return False
    codimension = dimension-rank_of(hard, range(len(hard)))
    if not codimension:
        return True
    return independent_count(quotient_rows(hard, soft), genus) == genus*codimension


def exhaustive_check(genus, hard, soft):
    """Direct dimension of each relative column space via the hard projected rank."""
    hard = integral(hard)
    soft = integral(soft)
    dimension = len(hard[0])
    target = genus*dimension-len(soft)
    if betti_rank(genus, hard) != target:
        return False
    for selected in subsets(len(soft)):
        directions = [soft[i] for i in selected]
        d = rank_of(soft, selected)
        total = [*directions,*hard]
        def projected(indices):
            return rank_of(total, [*range(len(directions)),*(i+len(directions) for i in indices)])-d
        v = min(len(hard)-len(t)+genus*projected(t) for t in subsets(len(hard)))
        r = min(genus*projected(range(len(hard))), 3*genus-3+v)
        if genus*d+r-target < len(selected):
            return False
    return True


def check():
    rng = Random(91137)
    count = 0
    for dimension, genus in ((4,2),(5,2),(5,3),(5,4),(6,2)):
        basis = [tuple(int(i==j) for j in range(dimension)) for i in range(dimension)]
        n = (dimension-3)*genus+3
        for _ in range(9):
            k = rng.randrange(1, n-2)
            d = rng.randrange(2,dimension+1)
            hard = [rng.choice(basis[:d]) for i in range(n-k-1)]
            hard.append(tuple(-sum(col) for col in zip(*hard)))
            soft = [rng.choice(basis) for i in range(k)]
            assert sharp_soft(genus,hard,soft) == exhaustive_check(genus,hard,soft)
            count += 1
    basis = [tuple(int(i==j) for j in range(5)) for i in range(5)]
    for k in (1,2,3):
        hard = [*basis,*[basis[0]]*(3-k)]
        hard.append(tuple(-sum(col) for col in zip(*hard)))
        assert sharp_soft(3,hard,[(1,1,0,0,0)]*k)
        assert not sharp_soft(3,hard,[basis[0]]*k)
    # At genus two, four conserved hard waves span a three-plane. Three soft
    # labels cannot supply two bases of the two-dimensional quotient, even
    # though the seven labels split into two independent sets.
    hard = [*basis[:3],(-1,-1,-1,0,0)]
    soft = [*basis[3:],(0,0,0,1,1)]
    assert partition_witness([*hard,*soft],2)[0] is not None
    assert not sharp_soft(2,hard,soft)
    # With a hard hyperplane and two independent copies of one quotient basis
    # the leading order is attained; making one soft velocity tangent kills it.
    hard = [*basis[:4],(-1,-1,-1,-1,0)]
    assert sharp_soft(2,hard,[basis[4],(1,0,0,0,1)])
    assert not sharp_soft(2,hard,[basis[4],basis[0]])
    print("Soft leading coefficient: projected ranks and partition tests agree on",count,"frames.")
    print("Low-genus sharpness, tangent velocity and quotient-base obstructions verified.")


if __name__ == '__main__':
    check()
