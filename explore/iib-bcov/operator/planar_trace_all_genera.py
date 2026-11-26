"""Polynomial-in-genus replay of the two planar trace boundaries."""
from collections import Counter
from fractions import Fraction
from itertools import permutations
from math import factorial, prod
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory

from flint import fmpz_mpoly_ctx, fmpz_poly

from four_loop_planar import boundary, restricted_symbols
from bulk_eight import BULK, quartics
from elementary_scalar import elementary_three_loop


def assignment_count(omissions):
    """Ways to distribute omitted transverse labels over the distinct pair minors."""
    permutations_count = factorial(4) // prod(factorial(k) for k in Counter(omissions).values())
    arrangements = factorial(sum(omissions)) // prod(factorial(k) for k in omissions)
    return permutations_count * arrangements


def elementary_plane(n, t):
    return [(1,0)]*(n-1) + [(0,1)]*(n-2) + [(t,1),(1-n-t,1-n)]


def check_integrated_vertices():
    vector = (0,0,1,0,0)
    root = (0,0,0,1,0)
    fixed = [(0,0,0,0,1), (0,0,-1,-1,-1)]
    for t in (0,1,2):
        waves = [(*row,0,0,0) for row in elementary_plane(2,t)]
        assert sum(c*m for c,m in zip(BULK, quartics(waves + fixed + [vector,root]))) == -Fraction(1,17280)
    with TemporaryDirectory() as tmp:
        exe = Path(tmp)/'nine'
        run(['c++','-std=c++17','-O3',str(Path(__file__).with_name('nine_sextics.cpp')),
             '-o',str(exe)],check=True)
        values = []
        for t in (0,1,-1,2):
            waves = [(*row,0,0,0) for row in elementary_plane(3,t)]
            values.append(elementary_three_loop(waves+fixed, vector, root, exe=exe))
    j3 = (-Fraction(2078212919,635529742848000), Fraction(131371,49037788800))
    other = (-Fraction(49188779,2715939072000), Fraction(7783,523908000))
    assert 30*(other[0]*j3[1]-other[1]*j3[0]) == Fraction(171285139,96092097118617600000)
    for (at_zero,at_one,at_minus_one,at_two),period in zip(zip(*values),j3):
        assert at_one == at_minus_one and at_two == at_zero+4*(at_one-at_zero)
        assert at_zero - Fraction(7,9)*(at_one-at_zero) == 30*period
    print('Independent full genus-two and genus-three elementary vertices reproduce the period.')


def check():
    # After polarization against an arbitrary common measure with second
    # moment equal to the identity, the defects have degree at most three
    # in n-1. The README counts the remaining indistinguishable pairs.
    for n in range(2, 6):
        coords = fmpz_mpoly_ctx.get(tuple(f'x{i}' for i in range(4*n-4))).gens()
        plane = [(1,0), (0,1), *zip(coords[::2], coords[1::2])]
        (a,b,d), f, g = boundary(plane)
        assert n*factorial(n)*a == (n-1)*f + (2*n+1)*g
        assert factorial(n)*b == 2*(n-1)*f
        assert n*factorial(n)*d == (n-1)*f - (3*n-1)*g
        print('genus', n, 'matching and transfer monomials', len(f), len(g))

    # Take e^4 as the elementary point root and e_4 as the polarization
    # of the e^3 vector. The same four evaluations extract the physical
    # bulk source in the original point coordinates.
    for n in range(2, 8):
        first, second, recoil = (1,0), (0,1), (-n,1-n)
        probes = [
            [first]*(n-1) + [second]*(n-2) + [recoil, (1,1)],
            [first]*(n-1) + [second]*(n-1) + [recoil],
            [first]*n + [second]*(n-2) + [(-n,2-n)],
            [first]*n + [second]*(n-1),
        ]
        table = [restricted_symbols(p) for p in probes]
        row = [n*(a-b-c+d) for a,b,c,d in zip(*table)]
        assert row[0] == -factorial(n)*(n-1)**2*(2*n+1)
        assert row[1] == 0
        if n == 4:
            assert table == [(3258,0,-1134),(3492,216,-1728),(288,0,0),(36,216,0)]
            assert 22*row[0] + 18*row[2] == 0

    # This single elementary-field measurement kills both the unknown
    # quartic boundary kernel and the first resistance period.
    t = fmpz_poly([0,1])
    a,b,d = restricted_symbols([(1,0)]*3 + [(0,1)]*2 + [(t,1),(-3-t,-3)])
    assert a == 18*(64*t*t + 117) and b == 0
    assert 22*a - 5*b + 18*d == 576*(44*t*t + 45)
    assert 11*8*a - 9*9072 == 2304*(44*t*t + 45)

    # Two real inputs independently measure K_4 and J_4 before fixing the free trace.
    a,b,d = restricted_symbols([(1,0)]*3 + [(0,1)]*3 + [(t,t)])
    assert a == 18*t*t*(t*t+4) and b == 36*t*t*(t+3)**2
    assert d == 18*t*t*(2*t-1)
    relation = 3*t*t-6*t-5
    assert 22*a - 5*b + 18*d == 72*t*t*relation
    assert b - (2400*t+1600) == (12*t*t+96*t+320)*relation
    assert a - (408*t+290) == (6*t*t+12*t+58)*relation
    assert Fraction(25,9)*96 - Fraction(25,6)*80 == -Fraction(200,3)

    for n in range(2, 8):
        plane = elementary_plane(n, t)
        a,b,d = restricted_symbols(plane)
        factor = factorial(n-1)*(n-1)**2
        assert 2*(n-1)*a == factor*(2*(n-1)**2*(n*n-n+1) + n**3*(n-2)*t*t)
        assert b == 0 and d == -factor*(n-1)*(2*n-1)
        paired = [(1,0)]*(n-1) + [(0,1)]*(n-1) + [(t,t)]
        a,b,d = restricted_symbols(paired)
        factor = (n-1)*factorial(n-1)*t*t
        assert (a,b,d) == (factor*(t*t+n), 2*factor*(t+n-1)**2, factor*(2*t-1))
        assert 2*(3*n-1)*a - 5*b + 2*(2*n+1)*d == 2*factor*(n-2)*(3*t*t-6*t+3-2*n)
    check_integrated_vertices()

    # Directly enumerate the two quartic assignments used in the separate
    # four-loop replay, including the per-label permutation degeneracy.
    assert assignment_count((1,1,1,1)) == len(set(permutations(range(4))))
    assert assignment_count((0,1,1,1)) == 4 * len(set(permutations(range(3))))
    for n in range(4, 13):
        for pairs in (n, n-1):
            omissions = tuple(pairs//4 + int(i < pairs%4) for i in range(4))
            assert min(n-k for k in omissions) >= 3
            assert assignment_count(omissions) > 0


if __name__ == '__main__':
    check()
