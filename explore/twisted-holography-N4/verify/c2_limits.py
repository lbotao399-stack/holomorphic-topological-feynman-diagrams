"""C2 large-m behaviour from the reduced core algorithm (validated against gluing for m <= 9, faces for m <= 11)."""
import sys
from fractions import Fraction
from math import comb, factorial
import mpmath as mp
from c2_counts import core_statistics, P_reduced, B_formula, elem_sym

H = int(sys.argv[1]) if len(sys.argv) > 1 else 3
stats = {h: core_statistics(h) for h in range(H + 1)}


def P(m, h):
    return P_reduced(m, stats[h], closed=True)


def B_cores(m, h):
    return sum(c * comb(m, s) for s, (c, L, D) in stats[h].items())


def interp(points):
    """exact polynomial through (x, y) points -> list of coefficients, low to high."""
    n = len(points)
    coeffs = [Fraction(0)] * n
    for i, (xi, yi) in enumerate(points):
        # basis polynomial
        basis = [Fraction(1)]
        denom = Fraction(1)
        for j, (xj, _) in enumerate(points):
            if j == i:
                continue
            basis = [Fraction(0)] + basis
            for t in range(len(basis) - 1):
                basis[t] -= xj * basis[t + 1]
            denom *= (xi - xj)
        for t in range(n):
            coeffs[t] += yi * basis[t] / denom
    return coeffs


def evalpoly(c, x):
    return sum(ci * x**i for i, ci in enumerate(c))


one_handle = lambda n: comb(n + 2, 4) + n * (n - 1) + comb(n, 3) + 2 + (-1)**n * (n - 1)

print('closed one-handle formula vs reduced P_{n+1,1}, n=2..300:',
      all(P(n + 1, 1) == one_handle(n) for n in range(2, 301)))
for h in range(H + 1):
    print(f'h={h}: unmarked from cores == e_(2h+1)(1..m)/C(m+1,2) for m=1..150:',
          all(B_cores(m, h) == B_formula(m, h) for m in range(1, 151)))

for h in range(1, H + 1):
    deg = 4 * h
    start = 4 * h
    res = {}
    for parity in (0, 1):
        ms = [m for m in range(start, start + 2 * (deg + 1) + 2) if m % 2 == parity][:deg + 1]
        cP = interp([(Fraction(m), Fraction(P(m, h))) for m in ms])
        cD = interp([(Fraction(m), Fraction(P(m, h)) - B_formula(m, h)) for m in ms])
        check = all(evalpoly(cP, m) == P(m, h) for m in range(start, start + 80) if m % 2 == parity)
        checkD = all(evalpoly(cD, m) == P(m, h) - B_formula(m, h) for m in range(start, start + 80) if m % 2 == parity)
        smallest = min(m for m in range(1, start + 1) if all(evalpoly(cP, mm) == P(mm, h)
                                                          for mm in range(m, start + 1) if mm % 2 == parity))
        res[parity] = (cP, cD, check, checkD, smallest)
        print(f'h={h} parity={parity}: degree<={deg} fit verified on m<{start+80}: P {check}, P-B {checkD};'
              f' quasi-polynomial valid from m>={smallest}')
        print(f'   leading coeff of P (m^{deg}) = {cP[deg]}, claim 1/(2^{2*h}(2h+1)!) = {Fraction(1, 2**(2*h)*factorial(2*h+1))}')
        print(f'   coeff of P-B: m^{deg} -> {cD[deg]}, m^{deg-1} -> {cD[deg-1]}, claim 4h/(2^{{2h}}(2h+1)!) = '
              f'{Fraction(4*h, 2**(2*h)*factorial(2*h+1))}')

mp.mp.dps = 30


def ratio_sum(n, N):
    m = n + 1
    return sum(mp.mpf(P(m, h)) / mp.mpf(N)**(2 * h) for h in range(H + 1))


def sub_sum(n, N):
    m = n + 1
    return m * sum((mp.mpf(P(m, h)) - mp.mpf(B_formula(m, h).numerator) / B_formula(m, h).denominator)
                   / mp.mpf(N)**(2 * h) for h in range(H + 1))


for x in (mp.mpf(1), mp.mpf(2)):
    y = x / 2
    lim1 = mp.sinh(y) / y
    lim2 = 2 * (mp.cosh(y) - mp.sinh(y) / y)
    print(f'x={x}: sinh(x/2)/(x/2) = {mp.nstr(lim1, 15)}, 2(cosh(x/2)-sinh(x/2)/(x/2)) = {mp.nstr(lim2, 15)}')
    for N in (10**3, 10**4, 10**5, 10**6, 10**7, 10**8):
        n = int(mp.nint(mp.sqrt(x * N)))
        r = ratio_sum(n, N)
        s = sub_sum(n, N)
        print(f'   N={N:>10} n={n:>6} n^2/N={mp.nstr(mp.mpf(n)**2/N, 8)}  sum_h P/N^2h = {mp.nstr(r, 12)}'
              f' (diff {mp.nstr(r - lim1, 3)}, (n+1)*diff {mp.nstr((n+1)*(r - lim1), 6)})'
              f'   m*sum(P-B)/N^2h = {mp.nstr(s, 10)} (diff {mp.nstr(s - lim2, 3)})', flush=True)
