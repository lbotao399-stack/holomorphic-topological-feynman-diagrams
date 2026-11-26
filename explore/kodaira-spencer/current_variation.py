"""Differentiate the projective current numerator at fixed graph and external sources."""

from fractions import Fraction as F
from math import lcm
from flint import fmpq_mat, fmpz_mpoly_ctx


def adjugate(matrix, ctx):
    """Fraction-free Gauss–Jordan elimination with sparsity-preserving pivots."""
    size = len(matrix)
    zero, one = ctx.constant(0), ctx.constant(1)
    a = [
        [*row, *[one if i == j else zero for j in range(size)]]
        for i, row in enumerate(matrix)
    ]
    sign, previous = 1, one
    order = list(range(size))
    for k in range(size):
        rows = [sum(bool(a[i][j]) for j in range(k, size)) for i in range(k, size)]
        cols = [sum(bool(a[i][j]) for i in range(k, size)) for j in range(k, size)]
        key = min(
            (len(a[i][j]), rows[i - k] * cols[j - k], i, j)
            for i in range(k, size)
            for j in range(k, size)
            if a[i][j]
        )
        i, j = key[-2:]
        if i != k:
            a[i], a[k] = a[k], a[i]
            sign = -sign
        if j != k:
            for row in a:
                row[j], row[k] = row[k], row[j]
            order[j], order[k] = order[k], order[j]
            sign = -sign
        pivot = a[k][k]
        for i in range(size):
            if i == k:
                continue
            source = a[i][k]
            for j in [*range(k + 1, size), *range(size, 2 * size)]:
                a[i][j] = (pivot * a[i][j] - source * a[k][j]) // previous
            a[i][k] = zero
        previous = pivot
    out = [None] * size
    for i, j in enumerate(order):
        out[j] = [sign * x for x in a[i][size:]]
    return sign * previous, out


def gradients(graph, flow, variations):
    """Differentiate the same Schur complement used in current_minor.numerator.

    The constant block B has one row per edge. After pivoting its columns,
    E=B_p^{-1} B_rest. The polynomial block has rows W in the same order.
    H=W_rest-E.transpose()*W_p, with independent integer row scalings.
    A single adjugate of H gives all conserved source derivatives.
    """
    n, m, h = graph.n, graph.N, len(flow[0])
    dim = n * (n + h)
    cycles = list(map(list, zip(*graph.C)))
    extended = [c + y for c, y in zip(cycles, flow)]
    ctx = fmpz_mpoly_ctx.get(tuple(f"t{i}" for i in range(m)))
    t = ctx.gens()
    slot = lambda r, k: n * r + k if k < n else n * n + (k - n) * n + r
    mat = [[0] * dim for _ in range(m)]
    for e, b in enumerate(extended):
        for r in range(n):
            for k, v in enumerate(b):
                mat[e][slot(r, k)] = b[r] * v
    reduced, rank = fmpq_mat(mat).rref()
    assert rank == m
    pivots = [next(j for j in range(dim) if reduced[i, j]) for i in range(m)]
    rest = [j for j in range(dim) if j not in pivots]
    bp = fmpq_mat([[row[j] for j in pivots] for row in mat])
    inv = bp.inv()
    sign = (-1) ** (sum(j - i for i, j in enumerate(pivots)) + (n * n - 1) * m)
    prefactor = sign * F(str(bp.det()))
    # Scale each null row to an integer vector, undoing the scales in the prefactor.
    ec = fmpq_mat([[reduced[i, j] for j in rest] for i in range(m)])
    projection = []
    for p, j in enumerate(rest):
        scale = lcm(*(int(ec[i, p].denominator) for i in range(m)))
        v = [0] * dim
        v[j] = scale
        for i, k in enumerate(pivots):
            v[k] = -int(scale * ec[i, p])
        projection.append((v, scale))
        prefactor /= scale
    rs = [(r, s) for r in range(n) for s in range(n) if (r, s) != (n - 1, n - 1)]

    def variable_block(src):
        block = [[ctx.constant(0) for pair in rs] for j in range(dim)]
        for col, (r, s) in enumerate(rs):
            for k in range(n + h):
                block[slot(r, k)][col] = sum(
                    (y[k] * c[s] * te for y, c, te in zip(src, cycles, t)),
                    ctx.constant(0),
                )
        return block

    w = variable_block(extended)

    def apply(rows, block):
        return [
            [
                sum((a * p for a, p in zip(row, column) if a and p), ctx.constant(0))
                for column in zip(*block)
            ]
            for row in rows
        ]

    H = apply([v for v, _ in projection], w)
    det, adj = adjugate(H, ctx)
    nbase = (det * prefactor.numerator) // prefactor.denominator
    assert nbase * prefactor.denominator == det * prefactor.numerator
    # First contract the adjugate with the linear coefficient of each length.
    base_wp = [w[j] for j in pivots]
    adj_wp = [
        [
            sum((a * b for a, b in zip(rr, col) if a and b), ctx.constant(0))
            for col in zip(*adj)
        ]
        for rr in base_wp
    ]
    out = []
    for shift in variations:
        dext = [[0] * n + q for q in shift]
        db = [[0] * dim for _ in range(m)]
        for e, y in enumerate(shift):
            for r in range(n):
                for k, x in enumerate(y):
                    db[e][slot(r, n + k)] = cycles[e][r] * x
        dBP = fmpq_mat([[row[j] for j in pivots] for row in db])
        de0 = inv * dBP
        trace = sum(de0[i, i] for i in range(m))
        # dE = B_p^{-1}(dB_rest - dB_p E).
        de = inv * fmpq_mat([[row[j] for j in rest] for row in db]) - de0 * ec
        dw = apply([v for v, _ in projection], variable_block(dext))
        derivative = sum(
            (
                adj[j][i] * dw[i][j]
                for j in range(len(rs))
                for i in range(len(rs))
                if dw[i][j] and adj[j][i]
            ),
            ctx.constant(0),
        )
        extra = [(F(str(trace)), det)]
        for i in range(m):
            for j, (_, scale) in enumerate(projection):
                if de[i, j] and adj_wp[i][j]:
                    extra.append((-F(str(scale * de[i, j])), adj_wp[i][j]))
        scale = lcm(*(v.denominator for v, _ in extra))
        derivative = scale * derivative + sum(
            (int(c * scale) * p for c, p in extra), ctx.constant(0)
        )
        denom = prefactor.denominator * scale
        der = (derivative * prefactor.numerator) // denom
        assert der * denom == derivative * prefactor.numerator
        out.append(der)
    return nbase, out
