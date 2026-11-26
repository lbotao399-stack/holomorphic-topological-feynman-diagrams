"""Five- and six-mark Abel moments of the fourfold primary with g-1 line waves."""
from fractions import Fraction as F
from itertools import combinations
from math import factorial
from random import Random

from flint import fmpz_mat
from worldsheet_primary import fourfold_primary


def check_fiber():
    rng = Random(862)
    for genus in range(3,10):
        for remaining in (2,3):
            rows = [[rng.randrange(-5,6) for _ in range(remaining)] for _ in range(genus)]
            hard = fmpz_mat(rows)
            unit = [[int(i==j) for j in range(genus)] for i in range(genus)]
            total = 0
            for soft in combinations(range(genus),genus-remaining):
                extended = fmpz_mat([row+[unit[i][j] for j in soft] for i,row in enumerate(rows)])
                total += factorial(len(soft))*int(extended.det())**2
            assert total == factorial(genus-remaining)*(hard.transpose()*hard).det() > 0


def genus_three(a,b):
    waves = [(a,0,0,0),(b,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1),(-a-b,-1,-1,-1)]
    return fourfold_primary(3,waves)


if __name__ == '__main__':
    check_fiber()
    plus,minus = genus_three(1,1),genus_three(1,-1)
    five = tuple((x+y)/4 for x,y in zip(plus,minus))
    assert five == (F(740466719,423686495232000),-F(469153,326918592000))
    assert minus == (F(841479311,211843247616000),-F(53359,16345929600))
    for a,b in ((1,3),(3,-2),(1,F(1,2)),(-4,-1)):
        value = tuple((a*b)**2*((a+b)**2*x-a*b*y) for x,y in zip(five,minus))
        assert value == genus_three(a,b)
    print('Genus-three five-mark moment',*five)
    print('Genus-three six-mark moment',*minus)
