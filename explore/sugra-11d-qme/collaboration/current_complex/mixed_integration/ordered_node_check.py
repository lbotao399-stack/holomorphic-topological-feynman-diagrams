"""Independent rational integration of one full 94-row node on every shape."""

from collections import Counter
from pathlib import Path

from detector_node import single_fields
from ordered_shapes import determinant, expectation, numerator
from single_projection import cases, detector_rows
from theta_denominators import SHAPES

ORDERS = ("0146723589", "0146723589", "0146723589", "6901472358", "0691472358")


def exact_node(shape, seq, case=35):
    order = tuple(map(int, seq))
    _, br, sc, m, _ = next(v for v in cases() if v[0] == case)
    total = Counter()
    for weight, bits, eta in detector_rows():
        data, cols = single_fields(br, bits, 0, eta, sc, m)
        for key, v in numerator(data, shape, order, True).items():
            total[key] += weight * v
    d = determinant(cols, shape, order) * (-1) ** (9 - order.index(9))
    return d * expectation({k: v for k, v in total.items() if v}, shape, order, data[0])


if __name__ == "__main__":
    for (shape, _, _), seq in zip(SHAPES, ORDERS):
        q = exact_node(shape, seq)
        residue = q.numerator * pow(q.denominator, -1, 1009) % 1009
        file = (
            Path(__file__).resolve().parent
            / "ordered_detector"
            / ("".join(map(str, shape)) + "_" + seq + ".tsv")
        )
        if file.exists():
            rows = dict(
                tuple(map(int, r.split())) for r in file.read_text().splitlines()[1:]
            )
            if 35 in rows:
                assert rows[35] == residue
        print(shape, "case 35", q, "residue", residue, flush=True)
