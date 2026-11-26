"""Verify exact scalar restriction ranks on three and four ordinary momentum lines."""

from collections import defaultdict
from sys import stdin

import numpy as np

from mixed_probe_jets import matrix, pivots
from mixed_transverse_directions import solve_modulo


def verify(lines):
    data = defaultdict(list)
    for line in lines:
        tag, *row = map(int, line.split())
        data[tag].append(row)
    general = np.asarray([[v % 251 for v in r] for r in data[0]], dtype=np.int64)
    cochains, boundaries = general[:, :105], general[:, 105:]
    independent = pivots(cochains.T)
    assert len(independent) == 105
    differential = solve_modulo(cochains[independent], boundaries[independent])
    assert np.array_equal(cochains @ differential % 251, boundaries)
    for tag, expected in ((1, 8), (2, 50)):
        rows = list(set(map(tuple, data[tag])))
        assert matrix(rows).rank() == expected
        assert len(pivots(np.asarray(rows) @ differential % 251)) == expected
        print(
            "ordinary lines",
            tag + 2,
            "cochain and boundary restriction ranks",
            expected,
        )


if __name__ == "__main__":
    verify(stdin)
