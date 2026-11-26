"""Stable equivariant matrix CE cohomology by the number of external central inputs."""
from random import Random

from sympy import Matrix
from invariant_ranks import balance, differential, dimension, independent_indices


def symmetric_multiplicity(w):
    return sum(2*i+3*j+4*k == w for k in range(w//4+1)
               for j in range(w//3+1) for i in range(w//2+1))


def prediction(w, central):
    if central == 0:
        return int(w == 2)
    if central == 2:
        return int(w % 2 == 0)
    if central == 4:
        return symmetric_multiplicity(w if w % 2 else w-6)
    return 0


def verify(w, central):
    three, _ = dimension(3, w, 5, central)
    four, samples = dimension(4, w, 5, central)
    r3 = len(independent_indices(Matrix.vstack(*(differential(3, three, p, central)
                                                  for p in samples)))) if three else 0
    h3 = len(three) - r3  # SL(2)-invariant C^2 vanishes at positive weight.
    assert h3 == int(central == 3 and w % 2 == 1)
    h4 = prediction(w, central)
    upper = len(four) - r3 - h4
    assert upper >= 0
    if central == 4:
        assert upper == 0 and not three
    r4 = 0
    rng = Random(2009)
    matrix = Matrix.zeros(0, len(four))
    for _ in range(50):
        if r4 == upper:
            break
        p = balance([(rng.randrange(-5, 6), rng.randrange(-5, 6)) for _ in range(4)])
        matrix = matrix.col_join(differential(4, four, p, central))
        r4 = len(independent_indices(matrix))
    assert r4 == upper
    return h3, len(four)-r3-r4


if __name__ == '__main__':
    for w in range(2, 9):
        result = [verify(w, central) for central in range(5)]
        print(w, 'H3', *(r[0] for r in result), 'H4', *(r[1] for r in result), flush=True)
