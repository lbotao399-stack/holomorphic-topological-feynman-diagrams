"""dim V^g for V = S^2, Lambda^2, S^3, Lambda^3, T^3, S^4, Lambda^4, T^4 of the adjoint representation,
by  dim V^g = sum_{w in W} sgn(w) mult_V(rho - w rho)  (Weyl character formula), own root systems from Cartan
matrices (positive roots by root strings), plethysm formulas for S^k, Lambda^k, exact integer arithmetic."""
import itertools
from fractions import Fraction as Fr
from collections import defaultdict

CARTAN = {
    'A1': [[2]],
    'A2': [[2, -1], [-1, 2]],
    'A3': [[2, -1, 0], [-1, 2, -1], [0, -1, 2]],
    'A4': [[2, -1, 0, 0], [-1, 2, -1, 0], [0, -1, 2, -1], [0, 0, -1, 2]],
    'B2': [[2, -2], [-1, 2]],          # A_ij = <alpha_i, alpha_j^vee>; alpha_2 short
    'B3': [[2, -1, 0], [-1, 2, -2], [0, -1, 2]],
    'C3': [[2, -1, 0], [-1, 2, -1], [0, -2, 2]],
    'D4': [[2, -1, 0, 0], [-1, 2, -1, -1], [0, -1, 2, 0], [0, -1, 0, 2]],
    'G2': [[2, -3], [-1, 2]],
    'F4': [[2, -1, 0, 0], [-1, 2, -2, 0], [0, -1, 2, -1], [0, 0, -1, 2]],
}
NAMES = {'A1': 'sl2', 'A2': 'sl3', 'A3': 'sl4 = so6', 'A4': 'sl5', 'B2': 'so5 = sp4', 'B3': 'so7', 'C3': 'sp6',
         'D4': 'so8', 'G2': 'G2', 'F4': 'F4'}


def pair(lam, i, A):
    """<lam, alpha_i^vee> for lam in root coordinates: sum_j lam_j A[j][i]."""
    return sum(lam[j] * A[j][i] for j in range(len(A)))


def positive_roots(A):
    r = len(A)
    simple = [tuple(int(i == j) for j in range(r)) for i in range(r)]
    roots = set(simple)
    frontier = list(simple)
    while frontier:
        new = []
        for b in frontier:
            for i in range(r):
                # p = max k with b - k alpha_i a root
                p = 0
                while True:
                    c = tuple(b[j] - (p + 1) * (j == i) for j in range(r))
                    if c in roots:
                        p += 1
                    else:
                        break
                q = p - pair(b, i, A)
                if q > 0:
                    c = tuple(b[j] + (j == i) for j in range(r))
                    if c not in roots:
                        roots.add(c)
                        new.append(c)
        frontier = new
    return sorted(roots)


def weyl_rho_points(A, pos):
    """{(rho - w rho, sgn w)} via BFS on the orbit of rho (regular), length = BFS depth."""
    r = len(A)
    rho = [Fr(0)] * r
    for a in pos:
        for j in range(r):
            rho[j] += Fr(a[j], 2)
    rho = tuple(rho)
    seen = {rho: 0}
    frontier = [rho]
    while frontier:
        new = []
        for lam in frontier:
            for i in range(r):
                c = pair(lam, i, A)
                mu = tuple(lam[j] - c * (j == i) for j in range(r))
                if mu not in seen:
                    seen[mu] = seen[lam] + 1
                    new.append(mu)
        frontier = new
    out = []
    for lam, L in seen.items():
        d = tuple(rho[j] - lam[j] for j in range(r))
        assert all(x.denominator == 1 for x in d)
        out.append((tuple(int(x) for x in d), (-1) ** L))
    return out


def mul(a, b):
    out = defaultdict(int)
    for k1, v1 in a.items():
        for k2, v2 in b.items():
            out[tuple(x + y for x, y in zip(k1, k2))] += v1 * v2
    return out


def scale(a, k):
    return {tuple(k * x for x in key): v for key, v in a.items()}


def lin(*terms):
    out = defaultdict(Fr)
    for c, ch in terms:
        for k, v in ch.items():
            out[k] += c * v
    return out


def main():
    print("%-4s %-10s %5s %5s | S2 L2 S3 L3 T3 S4 L4 T4" % ('type', 'algebra', '#pos', '|W|'))
    for t, A in CARTAN.items():
        r = len(A)
        pos = positive_roots(A)
        chi = defaultdict(int)
        for a in pos:
            chi[a] += 1
            chi[tuple(-x for x in a)] += 1
        chi[(0,) * r] += r
        pts = weyl_rho_points(A, pos)
        c2 = scale(chi, 2); c3 = scale(chi, 3); c4 = scale(chi, 4)
        x2 = mul(chi, chi)
        x3 = mul(x2, chi)
        x4 = mul(x2, x2)
        c2x = mul(c2, chi); c2x2 = mul(c2, x2); c2c2 = mul(c2, c2); c3x = mul(c3, chi)
        V = {
            'S2': lin((Fr(1, 2), x2), (Fr(1, 2), c2)),
            'L2': lin((Fr(1, 2), x2), (Fr(-1, 2), c2)),
            'S3': lin((Fr(1, 6), x3), (Fr(1, 2), c2x), (Fr(1, 3), c3)),
            'L3': lin((Fr(1, 6), x3), (Fr(-1, 2), c2x), (Fr(1, 3), c3)),
            'T3': lin((Fr(1), x3)),
            'S4': lin((Fr(1, 24), x4), (Fr(1, 4), c2x2), (Fr(1, 8), c2c2), (Fr(1, 3), c3x), (Fr(1, 4), c4)),
            'L4': lin((Fr(1, 24), x4), (Fr(-1, 4), c2x2), (Fr(1, 8), c2c2), (Fr(1, 3), c3x), (Fr(-1, 4), c4)),
            'T4': lin((Fr(1), x4)),
        }
        dims = {}
        for k, ch in V.items():
            d = sum(sg * ch.get(p, 0) for p, sg in pts)
            assert d.denominator == 1
            dims[k] = int(d)
        # sanity: dim g = 2 #pos + r
        print("%-4s %-10s %5d %5d | %2d %2d %2d %2d %2d %2d %2d %2d   (dim %d)" % (
            t, NAMES[t], len(pos), len(pts), dims['S2'], dims['L2'], dims['S3'], dims['L3'], dims['T3'],
            dims['S4'], dims['L4'], dims['T4'], 2 * len(pos) + r))


if __name__ == '__main__':
    main()
