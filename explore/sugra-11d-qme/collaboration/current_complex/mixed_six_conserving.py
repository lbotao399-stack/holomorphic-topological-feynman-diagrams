"""Integer polarizations of the scalar six-cycle at one conserved momentum geometry."""

from collections import defaultdict
from pathlib import Path
from sys import stdin


# Indexed by the six polarization bits in conserving_fields.
WEIGHTS = (
    7,
    -1,
    3,
    -8,
    1,
    -2,
    7,
    9,
    -7,
    14,
    -4,
    7,
    -7,
    8,
    -3,
    -5,
    -1,
    4,
    -13,
    -6,
    -3,
    -2,
    2,
    0,
    16,
    -15,
    -7,
    13,
    -13,
    9,
    1,
    -10,
    -3,
    -10,
    3,
    9,
    -3,
    0,
    -9,
    1,
    -12,
    -3,
    -2,
    1,
    14,
    -16,
    -1,
    10,
    -11,
    -4,
    0,
    9,
    4,
    -3,
    -11,
    11,
    0,
    2,
    -6,
    -9,
    7,
    -16,
    12,
    5,
)
CLASS_RATIO = -2217817151632288990651635154940307443560736538601897342476896593773383364911585160693312


def verify(lines):
    from mixed_probe_jets import matrix, pivots

    rows = defaultdict(list)
    for line in lines:
        tag, *row = map(int, line.split())
        rows[tag].append(row)
    assert [len(rows[j]) for j in range(5)] == [45, 64, 64, 64, 64]
    old = [0] * 255
    for line in (
        Path(__file__).with_name("mixed_six_weights.tsv").read_text().splitlines()
    ):
        w, selector, mask = map(int, line.split())
        old = [x + w * y for x, y in zip(old, rows[1 + selector][mask])]
    new = [sum(w * x for w, x in zip(WEIGHTS, col)) for col in zip(*rows[3])]
    assert not any(old[:150] + new[:150])
    closed = matrix(rows[0]).nullspace().to_Matrix()
    assert closed.rows == 77
    old_pairings = list(closed @ matrix([old[150:]]).to_Matrix().T)
    new_pairings = list(closed @ matrix([new[150:]]).to_Matrix().T)
    assert any(new_pairings) and all(
        x == CLASS_RATIO * y for x, y in zip(old_pairings, new_pairings)
    )
    print("six transverse lines, 60 polarizations, old/new =", CLASS_RATIO)

    # The five-line fixture uses opposite momenta on its repeated line.
    # A rational row combination must annihilate all incoming boundaries.
    fixed = rows[4]
    incoming = [r[:150] for r in fixed]
    pick = next(i for i, y in enumerate(new_pairings) if y)
    cocycle = list(map(int, closed.row(pick)))
    values = [sum(x * y for x, y in zip(r[150:], cocycle)) for r in fixed]
    columns = pivots(incoming)
    assert len(columns) == 33
    augmented = [row + [z] for row, z in zip(incoming, values)]
    assert pivots(augmented) == columns + [150]
    masks = pivots(zip(*([r[j] for j in columns + [150]] for r in augmented)))
    minor = matrix([[augmented[i][j] for i in masks] for j in columns + [150]])
    numer, den = minor.solve_den(matrix([[0]] * 33 + [[1]]))
    weights = list(map(int, numer.to_Matrix()))
    combined = [
        sum(w * fixed[i][j] for w, i in zip(weights, masks)) for j in range(255)
    ]
    assert not any(combined[:150])
    assert den and sum(x * y for x, y in zip(combined[150:], cocycle)) == den
    print(
        "five transverse lines, exact boundary cancellations", 150, "augmented rank", 34
    )


if __name__ == "__main__":
    verify(stdin)
