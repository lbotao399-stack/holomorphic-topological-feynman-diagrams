"""Longitudinal Gaussian, inverse transverse coordinates and native one-loop residue."""

from pathlib import Path
from random import Random
import sys
import sympy as sp

from coefficients import coefficient, sub

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from chain import gauge_box


def run():
    rng = Random(12701)
    counts = [0, 0]
    for m in range(1, 10):
        size = m + 1
        edges = [(j, j + 1) for j in range(m)] + [(j, size + j) for j in range(size)]
        incidence = sp.zeros(len(edges), 2 * size)
        for e, (left, right) in enumerate(edges):
            incidence[e, left] = 1
            incidence[e, right] = -1
        resistance = [
            sp.Rational(rng.randrange(1, 19), rng.randrange(1, 9)) for _ in edges
        ]
        for ground in (0, m, size, 2 * size - 1):
            tangent = incidence.copy()
            tangent.col_del(ground)
            assert tangent.det() in (-1, 1)
            weights = sp.diag(*(1 / t for t in resistance))
            assert (tangent.T * weights * tangent).det() == 1 / sp.prod(resistance)
            # A translation column is invisible to each edge difference.
            volume = incidence.col_join(sp.eye(2 * size)[m, :])
            assert abs(volume.det()) == 1
            counts[0] += 1
        transverse = incidence[:, :size]
        response = (
            sp.diag(*(1 / t for t in resistance))
            * transverse
            * (transverse.T * weights * transverse).inv()
        )
        alpha = list(response[m:, 0])
        beta = list(response[m:, m])
        assert sum(alpha) == sum(beta) == 1
        first = second = 0
        for j in range(m):
            first += alpha[j]
            second += beta[j]
            denominator = alpha[j + 1] * second + beta[j + 1] * (1 - first)
            assert (
                alpha[j] * second + beta[j] * (1 - first)
            ) / denominator == resistance[m + j + 1] / resistance[m + j]
            assert (
                alpha[j] * beta[j + 1] - alpha[j + 1] * beta[j]
            ) / denominator == resistance[j] / resistance[m + j]
        counts[1] += m
    cases = 0
    for extra in range(5):
        for a in range(extra + 1):
            total = (a + 1, extra - a + 1)
            for i in range(total[0] + 1):
                for j in range(total[1] + 1):
                    r = (i, j)
                    t = sub(total, r)
                    for x in range(total[0]):
                        for y in range(total[1]):
                            k = (x, y)
                            l = sub(sub(total, (1, 1)), k)
                            direct = 16 * sum(gauge_box((*r, *t), (*k, *l)))
                            assert direct == coefficient(r, t, (k, l))
                            cases += 1
    print(
        *counts,
        "tangent groundings and inverse path edges;",
        cases,
        "native fixed-line residue checks",
    )


if __name__ == "__main__":
    run()
