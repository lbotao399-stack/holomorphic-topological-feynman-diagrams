"""Independent fixed-length integration via a grounded electrical network and a polynomial adjugate."""

from argparse import ArgumentParser
from itertools import combinations, permutations
from math import factorial, prod
from gzip import open as gzip_open
from pathlib import Path
import networkx as nx
from flint import fmpq, fmpq_mat, fmpq_mpoly_ctx
from four_loop_arakelov import cubic_cores, graph_shapes, skeleton_chains


def adjugate(rows, ctx):
    """Return determinant and adjugate, up to the same row pivot sign."""
    n = len(rows)
    a = [r + [ctx.constant(int(i == j)) for j in range(n)] for i, r in enumerate(rows)]
    old = ctx.constant(1)
    for k in range(n):
        idx = next((i for i in range(k, n) if a[i][k]), None)
        if idx is None:
            return ctx.constant(0), None
        a[k], a[idx] = a[idx], a[k]
        pivot = a[k][k]
        for i, row in enumerate(a):
            if i == k:
                continue
            m = row[k]
            for j in range(2 * n):
                if j != k:
                    row[j] = (row[j] * pivot - m * a[k][j]) / old
            row[k] = ctx.constant(0)
        old = pivot
    return old, [row[n:] for row in a]


def cycle_matrix(edges):
    net = nx.Graph(edges)
    ids = {frozenset(e): i for i, e in enumerate(edges)}
    result = []
    for cycle in nx.cycle_basis(net, 0):
        row = [0] * len(edges)
        for v, w in zip(cycle, cycle[1:] + cycle[:1]):
            e = ids[frozenset((v, w))]
            row[e] = 1 if edges[e] == (v, w) else -1
        result.append(row)
    return result


def solve_shape(edges, times):
    nodes = nx.Graph(edges)
    branch, chains = skeleton_chains(edges, len(nodes))
    ordered = []
    for vs, es in chains:
        i, j = map(branch.index, (vs[0], vs[-1]))
        if i > j:
            i, j, vs, es = j, i, vs[::-1], es[::-1]
        ordered.append(((i, j), es, vs))
    ordered.sort(key=lambda row: (row[0], tuple(row[1])))
    core = [r[0] for r in ordered]
    ctx = fmpq_mpoly_ctx.get(("s0", "s1", "s2", "s3"))
    zero, one = ctx.constant(0), ctx.constant(1)
    segments, groups, full_edges, paths = [], [], [], []
    offset = 0
    for length, (_, _, vs) in zip(times, ordered):
        group = list(range(offset, offset + len(vs) - 2))
        offset += len(group)
        fractions = [ctx.gen(k) for k in group]
        paths.append(list(range(len(full_edges), len(full_edges) + len(vs) - 1)))
        segments.extend([length * f for f in fractions + [one - sum(fractions, zero)]])
        groups.append(group)
        full_edges.extend(zip(vs, vs[1:]))
    assert offset == 4
    net = nx.Graph(full_edges)
    marks = [v for v, d in nodes.degree if d == 2]
    E, nv = len(full_edges), len(net)
    C = cycle_matrix(full_edges)
    g = len(C)
    # If the metric of a core edge is fixed, the cycle Gram matrix is constant.
    gram = [
        [
            sum((c[e] * d[e] * length for e, length in enumerate(segments)), zero)
            for d in C
        ]
        for c in C
    ]
    assert all(entry.is_constant() for row in gram for entry in row)
    gram = fmpq_mat(
        [[row[j].coefficient(0) if row[j] else 0 for j in range(g)] for row in gram]
    )
    tree_polynomial = gram.det()
    inverse = gram.inv().tolist()
    path_flow = []
    ids = {frozenset(e): i for i, e in enumerate(full_edges)}
    for v in range(nv):
        path = nx.shortest_path(net, v, 0)
        row = [0] * E
        for left, right in zip(path, path[1:]):
            e = ids[frozenset((left, right))]
            row[e] = 1 if full_edges[e] == (left, right) else -1
        path_flow.append(row)
    currents = []
    for row in path_flow:
        drops = [
            sum((c[e] * row[e] * s for e, s in enumerate(segments)), zero) for c in C
        ]
        corr = [sum((b * v for b, v in zip(a, drops)), zero) for a in inverse]
        currents.append(
            [row[e] - sum((c[e] * v for c, v in zip(C, corr)), zero) for e in range(E)]
        )
    green = [
        [
            sum((p * s * j for p, s, j in zip(pp, segments, flow)), zero)
            for pp in path_flow
        ]
        for flow in currents
    ]
    assert all(green[v][w] == green[w][v] for v, w in combinations(range(nv), 2))
    R = [
        [green[v][v] + green[w][w] - 2 * green[v][w] for w in range(nv)]
        for v in range(nv)
    ]
    core_res = [
        [R[v][w].coefficient(0) if R[v][w] else fmpq(0) for w in branch] for v in branch
    ]
    assert all(R[v][w].is_constant() for v, w in combinations(branch, 2))
    # Zhang's tau = one quarter the integral of the squared derivative of resistance to vertex zero.
    tau = fmpq(0)
    for (v, w), length in zip(core, times):
        slope = (core_res[0][w] - core_res[0][v]) / length
        quadratic = 1 - core_res[v][w] / length
        tau += length * (slope**2 + quadratic**2 / 3) / 4
    theta = 2 * sum(core_res[v][w] for v, w in combinations(range(len(branch)), 2))
    epsilon = (4 * g - 4) * tau / g + theta / (2 * g)
    phi = (5 * g - 2) * tau / g + theta / (4 * g) - fmpq(sum(times), 4)

    # Sum the signed maximal minors by appending a unit row before differencing
    # the successive columns on each path. This gives sum(lengths) times the
    # oriented projective cofactor, even when every core edge carries a mark.
    columns = [(path[i], path[i - 1]) for path in paths for i in range(1, len(path))]
    columns += [(path[0], None) for path in paths]
    assert len(columns) == E

    def block(flow):
        return [
            [
                c[e] * flow[e] - (c[prev] * flow[prev] if prev is not None else 0)
                for e, prev in columns
            ]
            for c in C
        ]

    origin = currents[marks[-1]]
    blocks = [block([x - y for x, y in zip(flow, origin)]) for flow in currents]
    matrix = [r for v in marks[:3] for r in blocks[v]]
    matrix.append([ctx.constant(int(prev is None)) for _, prev in columns])
    det, adj = adjugate(matrix, ctx)
    if not det:
        return fmpq(0)
    for i, row in enumerate(matrix):
        for j in range(len(row)):
            assert sum((a * adj[k][j] for k, a in enumerate(row)), zero) == (
                det if i == j else zero
            )
    # Jacobi's formula differentiates the unsuppressed electrical determinant at a vertex.
    total = zero
    for v in range(nv):
        ds = [
            sum(
                (
                    adj[j][i * g + k] * blocks[v][k][j]
                    for k in range(g)
                    for j in range(E)
                ),
                zero,
            )
            for i in range(3)
        ]
        ds.append(g * det - sum(ds, zero))
        for i, h, k in permutations(range(4), 3):
            u, b, c = marks[i], marks[h], marks[k]
            voltage_twice = R[v][c] - R[v][u] + R[b][u] - R[b][c]
            voltage_derivative_twice = R[v][c] + R[b][v] - R[b][c]
            total += ds[i] * voltage_twice + det * voltage_derivative_twice / 2
    integrand = det * (total - 24 * det * (epsilon + phi))
    aut = sum(1 for _ in nx.isomorphism.GraphMatcher(nodes, nodes).isomorphisms_iter())
    factor = (
        tree_polynomial ** (2 * g + 1)
        * prod(t ** len(a) for t, a in zip(times, groups))
        / (aut * sum(times) ** 2)
    )
    value = fmpq(0)
    for ex, c in integrand.terms():
        powers = tuple(map(int, ex))
        value += (
            c
            * prod(factorial(j) for j in powers)
            / prod(factorial(len(a) + sum(powers[k] for k in a)) for a in groups)
        )
    return factor * value


def folded_evaluation(filename):
    """Undo the cut factors and average the polynomial over nine single doubled edges."""
    result = fmpq(0)
    factorials = [factorial(i) for i in range(12)]
    with gzip_open(filename, "rt") as file:
        for line in file:
            *powers, weight, unused = line.split()
            exponents = [int(p) - 1 for p in powers]
            assert len(exponents) == 9 and sum(exponents) == 39 and unused == "0"
            inverse_cut = fmpq(factorials[11], prod(factorials[j] for j in exponents))
            result += fmpq(weight) * inverse_cut * sum(2**j for j in exponents) / 9
    return result


def check(folded=None, selected=None):
    from fractions import Fraction
    from hashlib import sha256

    if selected is None:
        from four_loop_arakelov import period_polynomial
        from square_three import tight_graphs

        for genus, expected in (
            (2, fmpq(797, 2)),
            (3, fmpq(1078033588583880388, 2625)),
        ):
            total = fmpq(0)
            for graph in tight_graphs(3, genus):
                if max(dict(nx.Graph(graph.edges).degree()).values()) != 3:
                    continue
                core, reduced = period_polynomial(graph.edges)
                times = list(range(1, len(core) + 1))
                direct = solve_shape(graph.edges, times)
                assert direct == reduced(*times)
                total += direct
            assert total == expected

    sums = []
    folded_count = 0
    raw_count = 0
    for core in cubic_cores():
        diagrams = list(graph_shapes(core))
        if not diagrams:
            continue
        name = {3: "two-double", 8: "one-double", 28: "prism", 11: "bipartite"}[
            len(diagrams)
        ]
        values = []
        for i, edges in enumerate(diagrams):
            if selected and (name, i) != selected:
                continue
            result = solve_shape(edges, list(range(1, 10)))
            values.append(str(result))
            if selected:
                print(name, i, "at (1,...,9):", result)
            if folded is None:
                continue
            base = folded / f"{name}-{i:02d}"
            log = base.with_suffix(".log")
            if log.exists():
                for row in log.read_text().splitlines():
                    if row.startswith("Core-length control "):
                        assert fmpq(row.split()[-1]) == result, (
                            name,
                            i,
                            "core-length control",
                        )
                        raw_count += 1
            weights = base.with_suffix(".tsv.gz")
            if not weights.exists() or not base.with_suffix(".done").exists():
                continue
            direct = (
                sum(
                    solve_shape(edges, [1 + int(k == j) for k in range(9)])
                    for j in range(9)
                )
                / 9
            )
            actual = folded_evaluation(weights)
            assert direct == actual, (
                name,
                i,
                "automorphism-folded check",
                str(direct),
                str(actual),
            )
            folded_count += 1
        if values:
            sums.append(str(sum(map(Fraction, values))))
    if not selected:
        replay_id = sha256(" ".join(sums).encode()).hexdigest()
        print("Four core sums at distinct lengths:", *sums)
        assert (
            replay_id
            == "b491c8f83a592937813dcd6b4aed3210ba554ed0a5f53cad67f64f64aed1f875"
        )
    if folded is not None:
        print(
            "Matched raw polynomials:",
            raw_count,
            "matched folded polynomials:",
            folded_count,
        )


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument(
        "--folded", type=Path, help="Directory of completed orbit-weight files and logs"
    )
    parser.add_argument("--shape", nargs=2, metavar=("CORE", "INDEX"))
    args = parser.parse_args()
    check(args.folded, (args.shape[0], int(args.shape[1])) if args.shape else None)
