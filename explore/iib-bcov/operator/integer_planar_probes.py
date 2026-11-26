"""Two positive integer planar elementary-field probes separating K_n and J_n."""
from fractions import Fraction
from math import factorial, isqrt, lcm, prod


def paired_offsets(n, t):
    """Rational offsets of magnitude below one with the prescribed squared norm."""
    from sympy.solvers.diophantine.diophantine import sum_of_four_squares

    count = 2 * ((n - 1) // 2)
    norm = Fraction((n - 1) * (n - 2) * (2*n - 3*(t-1)**2), n * (3*n-1))
    squares = sum_of_four_squares(norm.numerator * norm.denominator)
    seed = [Fraction(int(x), norm.denominator) for x in squares] + [Fraction(0)] * (count - 4)
    if max(seed) < 1:
        return seed

    # Rational reflections approaching the constant vector of norm squared 'norm'.
    # norm/count < 1, so sufficiently accurate dyadic approximations work.
    scale = 1
    while True:
        entry = Fraction(isqrt(scale**2 * norm.numerator // (count * norm.denominator)), scale)
        direction = [x - entry for x in seed]
        multiplier = 2 * sum(x*y for x,y in zip(seed,direction)) / sum(y*y for y in direction)
        rotated = [x - multiplier*y for x,y in zip(seed,direction)]
        if all(abs(x) < 1 for x in rotated):
            assert sum(x*x for x in rotated) == norm
            return rotated
        scale *= 2


def integer_waves(n, t):
    """Return the two axial strength lists and the oblique strength in the planar frame."""
    if n < 2 or t not in (1, 2):
        raise ValueError("Expected genus >= 2 and probe 1 or 2.")
    if n == 2:
        return [1], [1], t
    if n == 4:
        if t == 1:
            return [1,16,16], [3,12,18], 11
        return [5,29,32], [11,23,32], 44
    if n == 3:
        offsets = [Fraction(1,2), Fraction(1,2) if t == 1 else Fraction(0)]
    else:
        offsets = paired_offsets(n, t)
    split = len(offsets)//2
    groups = []
    for half in (offsets[:split], offsets[split:]):
        strengths = [1+x for v in half for x in (v,-v)]
        if (n-1) % 2:
            strengths.append(Fraction(1))
        groups.append(strengths)
    scale = lcm(*(x.denominator for group in groups for x in group))
    return [int(scale*x) for x in groups[0]], [int(scale*x) for x in groups[1]], scale*t


def trace_symbols(first, second, oblique):
    """Exact three bracket polynomials, with all labeled multiplicities restored."""
    m = len(first)
    assert len(second) == m
    s, h = sum(first), sum(second)
    factor = m * factorial(m) * (prod(first) * prod(second) * oblique)**2
    a = oblique**2 + Fraction(m+1, 2*m) * sum(x*x for x in first+second)
    b = (s+oblique)**2 + (h+oblique)**2 + Fraction((s-h)**2, m)
    d = Fraction(oblique*(s+h), m) - Fraction(s*h, m*m)
    return tuple(factor*x for x in (a,b,d))


def check():
    from four_loop_planar import restricted_symbols

    for n in range(2, 25):
        rows = []
        for t in (1,2):
            first, second, oblique = integer_waves(n,t)
            assert len(first) == len(second) == n-1 and min(first+second) > 0
            a,b,d = trace_symbols(first,second,oblique)
            assert 2*(3*n-1)*a - 5*b + 2*(2*n+1)*d == 0
            assert d*2*(n-1+t)**2 == b*(2*t-1)
            if n <= 7:
                waves = [(x,0) for x in first] + [(0,y) for y in second] + [(oblique,oblique)]
                assert restricted_symbols(waves) == (a,b,d)
            rows.append((a,b))
        assert rows[0][0]*rows[1][1] != rows[1][0]*rows[0][1]
    print("Two positive integer probes at genera 2 through 24; direct labeled matchings through 7.")


if __name__ == "__main__":
    check()
