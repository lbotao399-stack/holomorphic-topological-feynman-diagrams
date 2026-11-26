"""Endpoint current wedge in the ordered 11d two-loop theta kernel (S3, eq-insertion)."""
from fractions import Fraction
from itertools import combinations, product
from pathlib import Path
import random
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "verify"))
sys.path.insert(0, str(HERE.parents[2] / "code"))
from htk.simplex import dirichlet
from vcommon import D_theta, PQR, theta_graph


def endpoint_component(prefix, alpha, powers, i, j):
    """Normalized component after tilting each simplex by a single monomial.

    prefix is Q^P_j in path order A,B,C. powers contains the w exponents
    followed by x^A, x^B, x^C exponents. The base measure is S3 eq-insertion.
    Multiply by the scalar integral of the given monomial for its raw moment.
    """
    weights = [a + k for a, k in zip(alpha, powers[0])]
    means = []
    for path, exponents in zip(prefix, powers[1:]):
        mass = sum(exponents) + len(path)
        means.append([sum((k + 1) * q[slot] for k, q in zip(exponents, path))
                      / Fraction(mass) for slot in (i, j)])
    total = Fraction(0)
    for p, q, r in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
        total += weights[p] * (means[q][0] * means[r][1] - means[q][1] * means[r][0])
    return total / sum(weights)


def mass(pqr, powers):
    """S3 scalar monomial integral divided by pi, including the 1/2, before D_theta."""
    w, *x = powers
    value, pi_power = dirichlet(tuple(Fraction(11, 2) - n + k for n, k in zip(pqr, w)))
    assert pi_power == 1
    for exponents in x:
        factor, pi_power = dirichlet(tuple(exponents))
        assert pi_power == 0
        value *= factor
    return value / 2


def edge_terms(prefix, path, slot):
    # Expand ell_P directly as w_P * sum mhat_R - mhat_P at lambda_N = 0.
    for r, seq in enumerate(prefix):
        for k, vec in enumerate(seq):
            if vec[slot]:
                yield vec[slot], ((0, path), (r + 1, k))
                if r == path:
                    yield -vec[slot], ((r + 1, k),)


def direct_component(prefix, pqr, powers, i, j):
    # Independent direct product of the two currents, preserving w_A*w_B terms.
    result = Fraction(0)
    for u, v, sign in ((i, j, 1), (j, i, -1)):
        for (a, slots_a), (b, slots_b) in product(edge_terms(prefix, 0, u),
                                                   edge_terms(prefix, 1, v)):
            shifted = [row[:] for row in powers]
            for row, col in slots_a + slots_b:
                shifted[row][col] += 1
            result += sign * a * b * mass(pqr, shifted)
    return result


def check():
    rng = random.Random(211103049)
    for pqr in PQR:
        nv, _, paths = theta_graph(*pqr)
        momenta = {v: [rng.randint(-3, 3) for _ in range(5)] for v in range(1, nv - 2)}
        momenta[nv - 2] = [-sum(vec[i] for vec in momenta.values()) for i in range(5)]
        if pqr == (1, 6, 6):
            for k in range(5):
                momenta[k + 1] = [int(i == k) for i in range(5)]
                momenta[k + 6] = [-int(i == 4-k) for i in range(5)]
        assert all(sum(vec[i] for vec in momenta.values()) == 0 for i in range(5))
        det = D_theta(momenta, paths)
        assert det != 0
        alpha = [Fraction(13, 2) - n for n in pqr]
        prefix = []
        for verts, _ in paths:
            seq = [[0] * 5]
            for v in verts:
                seq.append([x + y for x, y in zip(seq[-1], momenta[v])])
            prefix.append(seq)
        zero = [[0] * n for n in (3, *pqr)]
        for powers in (zero, [[rng.randint(0, 2) for _ in row] for row in zero]):
            for i, j in combinations(range(5), 2):
                expected = mass(pqr, powers) * endpoint_component(prefix, alpha, powers, i, j)
                assert direct_component(prefix, pqr, powers, i, j) == expected
        if pqr == (1, 6, 6):
            value = det * direct_component(prefix, pqr, zero, 0, 1)
            assert det == -1 and value == Fraction(1, 1123200)
            print(f"theta_166, D_theta={det}, integrated component [1,2] / pi = {value}")
        print(f"{pqr}: 10 exterior components, 2 monomial weights passed")


if __name__ == "__main__":
    check()
