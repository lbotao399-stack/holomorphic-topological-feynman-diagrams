"""Check the holomorphic gamma Ward identity at K=0 on the five 11d theta graphs."""

import random
import sys
from collections import Counter
from fractions import Fraction
from itertools import product
from math import comb
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "verify"))
from vcommon import D_theta, H, PQR, gauss_det, kirchhoff, theta_graph
from v_integrand import vertex_tensor


def check_site(vertex, kind, tails, heads, edges, current, external):
    """Compare the ordered Wick words with the two reduced tensors, in every component."""
    momentum = external[vertex]["lam"]
    gamma = external[vertex]["vg"]
    longitudinal = {vertex: dict(external[vertex], vg=momentum)}
    args = (vertex, kind, tails, heads, edges, current)
    _, measured = vertex_tensor(*args, external)
    _, pure_gradient = vertex_tensor(*args, longitudinal)
    if kind == "b":
        expected = {(i, j): -(momentum[i] * gamma[j] - momentum[j] * gamma[i])
                    for i, j in measured}
    else:
        signed = [[x if edges[e][1] == vertex else -x for x in current[e]] for e in heads]
        assert all(x + y + z == 0 for x, y, z in zip(*signed, momentum))
        dual = [[int(i == j) for j in range(H)] for i in range(H)]
        expected = {(i, j): gauss_det([dual[i], dual[j], gamma, signed[0], momentum])
                    for i, j in measured}
    assert measured == expected, (vertex, kind, "polarized tensor")
    assert not any(pure_gradient.values()), (vertex, kind, "longitudinal tensor")
    return len(measured), sum(value != 0 for value in measured.values())


def check_graph(lengths, rng):
    n, edges, paths = theta_graph(*lengths)
    interior = range(1, n - 1)

    def random_vector():
        return [rng.randrange(-4, 5) for _ in range(H)]

    external = {
        v: dict(u=random_vector(), vg=random_vector(),
                lam=random_vector() if v in interior else [0] * H)
        for v in range(n)
    }
    # The trivalent vertices have no external legs; all ten interior momenta sum to zero.
    external[n - 2]["lam"] = [-sum(external[v]["lam"][i] for v in range(1, n - 2))
                              for i in range(H)]
    momenta = {v: external[v]["lam"] for v in range(n)}
    assert D_theta(momenta, paths) != 0
    resistance = [Fraction(rng.randrange(1, 7), rng.randrange(1, 6)) for _ in edges]
    current = kirchhoff(n, edges, resistance, momenta)
    counts = Counter()
    occurrences = 0
    distinct = {}
    for reverse in product((False, True), repeat=len(edges)):
        heads = [[] for _ in range(n)]
        tails = [[] for _ in range(n)]
        for e, ((v, w), flipped) in enumerate(zip(edges, reverse)):
            if flipped:
                v, w = w, v
            tails[v].append(e)
            heads[w].append(e)
        if any(len(heads[v]) not in (1, 3) for v in (0, n - 1)):
            continue
        cubic_gamma = [len(heads[v]) >= 2 for v in range(n)]
        m = sum(cubic_gamma)
        counts[m] += 1
        gamma_sites = [v for v in interior if len(heads[v]) in (0, 2)]
        assert len(gamma_sites) == 2 * m - 1
        assert len(interior) - len(gamma_sites) == 11 - 2 * m
        internal_derivatives = 0
        for v in range(n):
            if cubic_gamma[v]:
                internal_derivatives += 1 if v in interior else 2
            elif heads[v]:
                internal_derivatives += 1
        assert internal_derivatives == 13 - m
        occurrences += len(gamma_sites)
        for v in gamma_sites:
            kind = "c" if cubic_gamma[v] else "b"
            key = v, kind
            if key not in distinct:
                distinct[key] = check_site(v, kind, tails[v], heads[v], edges, current, external)
    assert counts == {m: 4 * comb(10, 2 * m - 1) for m in range(1, 6)}
    assert len(distinct) == 20 and occurrences == 10240
    components = sum(tested for tested, nonzero in distinct.values())
    nonzero = sum(nonzero for tested, nonzero in distinct.values())
    print(f"{lengths}: flows by m={tuple(counts[m] for m in range(1, 6))}; "
          f"gamma legs={occurrences}; components={components}; generic nonzeros={nonzero}")
    return sum(counts.values()), occurrences, components, nonzero


if __name__ == "__main__":
    rng = random.Random(20261121)
    totals = tuple(map(sum, zip(*(check_graph(lengths, rng) for lengths in PQR))))
    assert totals[:3] == (10240, 51200, 2500)
    print(f"PASS: {totals[0]} flows; {totals[1]} gamma legs; {totals[2]} reduced tensor "
          f"equalities and {totals[2]} longitudinal zeros; {totals[3]} generic nonzeros.")
