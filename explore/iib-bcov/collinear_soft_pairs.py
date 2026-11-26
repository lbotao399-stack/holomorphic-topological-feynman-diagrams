"""Fiber moment identity for g+l collinear marks and l fourfold soft directions."""
from itertools import combinations
from math import factorial
from random import Random

from flint import fmpz_mat


def mixed_determinant(matrices):
    """Squarefree coefficient in det(sum x_i M_i) for r matrices of size r."""
    count = len(matrices)
    result = 0
    for length in range(1,count+1):
        for slots in combinations(range(count),length):
            total = fmpz_mat(count,count)
            for i in slots:
                total += matrices[i]
            result += (-1)**(count-length)*total.det()
    return int(result)


def matchings(indices, pairs):
    """Unordered disjoint pairs, with unpaired positions allowed."""
    if not pairs:
        yield ()
        return
    if len(indices)<2*pairs:
        return
    first,*rest = indices
    yield from matchings(rest,pairs)
    for i,other in enumerate(rest):
        for tail in matchings(rest[:i]+rest[i+1:],pairs-1):
            yield ((first,other),*tail)


def verify(genus, excess, slopes, projection):
    count = genus+excess
    full = []
    normal = []
    for slope in slopes:
        image = fmpz_mat([[sum(slope[j]*projection[i,j*genus+k] for j in range(4))
                           for k in range(genus)] for i in range(excess)])
        columns = fmpz_mat([[int(i==j) for j in range(genus)] for i in range(genus)]+image.tolist())
        full.append(columns*columns.transpose())
        normal.append(image)
    value = mixed_determinant(full)
    lower = 0
    for pairs in matchings(list(range(count)),excess):
        differences = [normal[i]-normal[j] for i,j in pairs]
        lower += mixed_determinant([d*d.transpose() for d in differences])
    assert value == factorial(genus-excess)*lower
    return value


def check():
    rng = Random(41983)
    for genus, excess in ((3,1),(4,2),(5,2),(5,3),(6,3)):
        count = genus+excess
        projection = fmpz_mat([[rng.randrange(-3,4) for j in range(4*genus)] for i in range(excess)])
        slopes = [tuple(rng.randrange(-2,3) for j in range(4)) for i in range(count)]
        value = verify(genus,excess,slopes,projection)
        assert value>0
        for multiplicity in (genus,genus+1):
            repeated = [(0,0,0,0)]*multiplicity+[(1,0,0,0)]*(count-multiplicity)
            assert bool(verify(genus,excess,repeated,projection)) == (multiplicity==genus)
        print("genus",genus,"pairs",excess,"unreduced mixed determinant",value)


if __name__=='__main__':
    check()
