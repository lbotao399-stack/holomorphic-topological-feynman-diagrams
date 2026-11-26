"""Cohomology of the ordinary transverse scalar sector and its marked divergence."""

from collections import defaultdict
from sys import stdin
from mixed_probe_jets import pivots, matrix
from relative_dimensions import dimensions
from mixed_vertical import hdim
from mixed_six_conserving import WEIGHTS


def verify(lines):
    rows = defaultdict(list)
    for line in lines:
        kind, *r = map(int, line.split())
        rows[kind].append(r)
    tangent, scalar = dimensions(3, 32, 9, -3)
    assert scalar == [0, 0, 0, 2, 16, 53, 74, 44, 11, 1]
    assert tangent == [0, 0, 0, 4, 31, 75, 70, 25, 3, 0]
    for prime in (251, 257):
        assert [hdim(6, q, (-2, -2, -2), prime) for q in range(3, 7)] == [
            (2, 0, 2, 0),
            (5, 2, 3, 0),
            (4, 3, 1, 0),
            (1, 1, 0, 0),
        ]
    degrees = rows[0][0]
    first = [j for j, v in enumerate(degrees) if v == 1]
    assert len(first) == 70 and len(pivots(rows[1])) == 28
    closed = matrix([[r[j] for j in first] for r in rows[1]]).nullspace().to_Matrix()
    assert closed.rows == 48
    assert len(pivots(r[:53] for r in rows[2])) == 39
    assert len(pivots(r[53:] for r in rows[2])) == 74
    assert len(pivots(rows[3])) == 34

    waves = rows[4]
    incoming = [r[:53] for r in waves]
    cols = pivots(incoming)
    assert len(cols) == 39
    values = matrix([[r[203 + j] for j in first] for r in waves]).to_Matrix() @ closed.T
    for j in range(values.cols):
        v = list(map(int, values.col(j)))
        if sum(x * y for x, y in zip(WEIGHTS, v)):
            eq = [r + [y] for r, y in zip(incoming, v)]
            if pivots(eq) == cols + [53]:
                break
    else:
        raise AssertionError("no closed marked representative")
    mask = pivots(zip(*([r[i] for i in cols + [53]] for r in eq)))
    square = matrix([[eq[k][i] for k in mask] for i in cols + [53]])
    num, den = square.solve_den(matrix([[0]] * 39 + [[1]]))
    weights = list(map(int, num.to_Matrix()))
    combined = [sum(x * waves[k][i] for x, k in zip(weights, mask)) for i in range(308)]
    assert not any(combined[:203])
    assert sum(x * v[k] for x, k in zip(weights, mask)) == den and den
    print("ordinary scalar cochain dimensions", 2, 16, 53, 74, 44, 11, 1)
    print("middle incoming and outgoing ranks", 39, 34)
    print("ordinary and marked boundary columns killed", 53, 150)
    print(
        "unmarked six-class detected by divergence of first-derivative representative"
    )


if __name__ == "__main__":
    verify(stdin)
