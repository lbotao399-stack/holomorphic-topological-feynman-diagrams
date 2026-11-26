"""Exact scalar theta_166 Wick and determinant for a fixed mixed-current labeling."""

from collections import Counter, defaultdict
from fractions import Fraction
from functools import cache
from itertools import permutations
from math import comb, factorial, prod

N = 5
IDENTITY = {0: 1}
MASK = (1 << 20) - 1


def collect(out, poly, w=1, shift=0):
    for key, c in poly.items():
        out[shift + key] = out.get(shift + key, 0) + w * c


def multiply(a, b):
    out = {}
    for key, c in a.items():
        collect(out, b, c, key)
    return {k: c for k, c in out.items() if c}


def epsilon(indices):
    if len(set(indices)) != len(indices):
        return 0
    return (-1) ** sum(a > b for i, a in enumerate(indices) for b in indices[i + 1 :])


def graph_factors(
    data, order, paths=(1, 2), endpoints=True, lie=False, shape=(1, 6, 6)
):
    momenta, polar = data
    gamma = order.index(9)
    edges, vertices = [], []
    edge_start = vertex_start = 0
    for length in shape:
        edges.append(list(range(edge_start, edge_start + length)))
        vertices.append(list(range(vertex_start, vertex_start + length - 1)))
        edge_start += length
        vertex_start += length - 1
    if edge_start != 13 or vertex_start != 10:
        raise ValueError("invalid two-loop shape")
    marked = next(h for h, slots in enumerate(vertices) if gamma in slots)
    currents = [None] * 13
    for h in range(3):
        offsets = [0] * N
        for pos, edge in enumerate(edges[h]):
            currents[edge] = []
            for i in range(N):
                v = {0: offsets[i]} if offsets[i] else {}
                if h != 2:
                    v[1 << (4 * i)] = 1 if h == 1 else -1
                if h != 1:
                    v[1 << (20 + 4 * i)] = 1 if h == 2 else -1
                currents[edge].append(v)
            if pos < len(vertices[h]):
                p = momenta[order[vertices[h][pos]]]
                offsets = [x + y for x, y in zip(offsets, p)]

    def local(h, labels, pos, forward):
        v = labels[pos]
        p, u = momenta[order[v]], polar[order[v]]
        edge = edges[h][pos]
        before, after = forward[edge], forward[edge + 1]
        mat = [[{} for j in range(N)] for i in range(N)]
        if before != after:
            if before:
                for i, j, g, k, l in permutations(range(N)):
                    collect(
                        mat[i][j],
                        currents[edge][k],
                        epsilon((i, j, g, k, l)) * u[g] * p[l],
                    )
            else:
                for i in range(N):
                    for j in range(N):
                        mat[i][j][0] = u[i] * p[j] - p[i] * u[j]
        else:
            curr = currents[edge + int(not before)]
            sign = 1 if before else -1
            trace = {}
            for i in range(N):
                collect(trace, curr[i], u[i])
            for i in range(N):
                for j in range(N):
                    if lie:
                        mat[i][j][0] = -u[i if before else j] * p[j if before else i]
                    else:
                        collect(
                            mat[i][j],
                            curr[j if before else i],
                            sign * u[i if before else j],
                        )
                    if i == j:
                        collect(mat[i][j], trace, -sign)
        return [[{k: x for k, x in t.items() if x} for t in row] for row in mat]

    def endpoint(right, forward):
        site = [es[-1 if right else 0] for es in edges]
        selected = [i for i, e in enumerate(site) if forward[e] == right]
        answer = defaultdict(dict)
        if len(selected) == 3:
            for i, j, k, a, b in permutations(range(N)):
                collect(
                    answer[i, j, k],
                    multiply(currents[site[1]][a], currents[site[2]][b]),
                    -epsilon((i, j, k, a, b)),
                )
        else:
            (h,) = selected
            j, k = [i for i in range(3) if i != h]
            for a in range(N):
                for b in range(N):
                    if a == b:
                        continue
                    ix = [0] * 3
                    ix[j], ix[k], ix[h] = a, b, a
                    collect(answer[tuple(ix)], currents[site[h]][b], 1 if right else -1)
                    ix[h] = b
                    collect(answer[tuple(ix)], currents[site[h]][a], -1 if right else 1)
        return answer

    for flow in range(4):
        forward, tails = [], []
        for h in range(3):
            nodes = [0, *(v + 1 for v in vertices[h]), 11]
            seen = False
            for left, right in zip(nodes, nodes[1:]):
                if h == marked:
                    f = seen if flow < 2 else not seen
                elif flow < 2:
                    f = flow == 1
                else:
                    f = (h == (marked + 1) % 3) == (flow == 2)
                forward.append(f)
                tails.append(left if f else right)
                if right == gamma + 1:
                    seen = True
        sign = (-1) ** sum(a > b for i, a in enumerate(tails) for b in tails[i + 1 :])
        chains = {}
        for h in paths:
            chain = [[IDENTITY if i == j else {} for j in range(N)] for i in range(N)]
            for pos in range(len(vertices[h])):
                mat = local(h, vertices[h], pos, forward)
                new = [[{} for j in range(N)] for i in range(N)]
                for i in range(N):
                    for k in range(N):
                        for j in range(N):
                            collect(new[i][j], multiply(chain[i][k], mat[k][j]))
                chain = new
            chains[h] = chain
        ends = (
            (endpoint(False, forward), endpoint(True, forward)) if endpoints else None
        )
        yield sign, chains, ends


def polynomials(data, order, lie=False):
    raw = {}
    for sign, chains, (left, right) in graph_factors(data, order, lie=lie):
        for (i, j, k), x in left.items():
            for (h, s, t), y in right.items():
                if h != i:
                    continue
                collect(
                    raw,
                    multiply(
                        multiply(x, chains[1][j][s]), multiply(y, chains[2][k][t])
                    ),
                    sign,
                )
    return {k: v for k, v in raw.items() if v}


@cache
def components(key):
    return tuple((key >> (4 * i)) & 15 for i in range(N))


@cache
def splits(key):
    if not key:
        return [(0, 0, 1)]
    j = next(i for i, n in enumerate(components(key)) if n)
    n = components(key)[j]
    unit = 1 << (4 * j)
    return [
        (k + x * unit, d + x, c * comb(n, x))
        for k, d, c in splits(key - n * unit)
        for x in range(n + 1)
    ]


def resistance(src):
    prefix = (0,) * N
    gaps = []
    for p in src:
        prefix = tuple(a + b for a, b in zip(prefix, p))
        gaps.append(prefix)

    @cache
    def homogeneous(gap, key):
        if not key:
            return 1
        if not gap:
            return 0
        return homogeneous(gap - 1, key) + sum(
            value * homogeneous(gap, key - (1 << (4 * i)))
            for i, (n, value) in enumerate(zip(components(key), gaps[gap - 1]))
            if n and value
        )

    @cache
    def moment(key):
        powers = components(key)
        return Fraction(
            factorial(len(gaps))
            * prod(map(factorial, powers))
            * homogeneous(len(gaps), key),
            factorial(len(gaps) + sum(powers)),
        )

    return moment


@cache
def conduct(a, b, c, d):
    return sum(
        Fraction(
            (-1) ** (a - i + c - j)
            * comb(a, i)
            * comb(c, j)
            * prod(range(1, 2 * (i + b), 2))
            * prod(range(1, 2 * (j + d), 2)),
            prod(range(13, 13 + 2 * (i + b + j + d), 2)),
        )
        for i in range(a + 1)
        for j in range(c + 1)
    )


def expectation(numerator, momenta, order):
    b_mean = resistance([momenta[j] for j in order[:5]])
    c_mean = resistance([momenta[j] for j in order[5:]])
    # Collect common split monomials before performing rational operations.
    terms = Counter()
    for key, coeff in numerator.items():
        b, c = key & MASK, key >> 20
        db, dc = sum(components(b)), sum(components(c))
        for ub, ib, cb in splits(b):
            for uc, ic, cc in splits(c):
                terms[b - ub + c - uc, ub + uc, db - ib, ib, ic, dc - ic] += (
                    coeff * cb * cc
                )
    return sum(
        v * b_mean(k) * c_mean(l) * conduct(a, b, c, d)
        for (k, l, a, b, c, d), v in terms.items()
        if v
    )


def fixture():
    p = [
        (1, 0, 0, 0, 0),
        (0, 2, 0, 0, 0),
        (3, 0, 0, 0, 0),
        (0, 1, 1, 0, 0),
        (0, 0, 1, 0, 0),
        (0, 0, 1, 1, 0),
        (0, 0, 0, 1, 0),
        (0, 0, 0, 1, 1),
        (0, 0, 0, 0, 1),
    ]
    p.append(tuple(-sum(column) for column in zip(*p)))
    u = [
        (0, 1, 0, 0, 0),
        (-1, 0, 0, 0, 0),
        (0, 1, 0, 0, 0),
        (-1, 0, 0, 0, 0),
        (0, 0, 0, 1, 0),
        (0, 0, -1, 1, 0),
        (0, 0, 0, 0, 1),
        (0, 0, 0, -1, 1),
        (0, 0, 1, 0, 0),
        (0, 0, 1, 0, 0),
    ]
    return p, u


if __name__ == "__main__":
    data = fixture()
    order = [int(j) for j in "0146723589"]
    raw = polynomials(data, order)
    result = expectation(raw, data[0], order)
    print(
        "Wick monomials",
        len(raw),
        "exact normalized Schwinger mean",
        result,
        "determinant-weighted",
        18 * result,
    )

    def at_prime(v, p):
        return v.numerator * pow(v.denominator, -1, p) % p

    assert at_prime(result, 1009) == 750
    assert at_prime(7**3 * 18 * result, 7) == 6
