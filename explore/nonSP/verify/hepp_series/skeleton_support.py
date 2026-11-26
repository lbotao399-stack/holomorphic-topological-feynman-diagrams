"""Strict edge-subset margins for the four-loop R x C^2 compressed numerator."""
import json
from collections import Counter
from fractions import Fraction
from itertools import product
from pathlib import Path

import numpy as np

from period import HeppSeries

ROOT = Path(__file__).resolve().parents[4]


def rows():
    graphs = json.loads((ROOT / 'deepen/o3/higher_orders/four-loop-graphs.json').read_text())
    exponents = np.array([a for a in product(range(3), repeat=9) if sum(a)==9], dtype=np.int16)
    subsets = np.array([[(mask >> j) & 1 for j in range(9)] for mask in range(1,511)], dtype=np.int16)
    original = exponents @ subsets.T
    assert len(exponents) == 3139
    core_cache = {}
    for i, graph in enumerate(graphs):
        chains = graph['chains']
        edges = tuple((c[0],c[-1]) for c in chains)
        if edges not in core_cache:
            core_cache[edges] = HeppSeries(edges, Fraction(11,2), {})
        engine = core_cache[edges]
        powers = np.array([len(c)-1 for c in chains], dtype=np.int16)
        margins = 2*(original + powers @ subsets.T) - 11*np.array(engine.cycles[1:-1])
        admissible = (margins>0).all(axis=1)
        yield i, engine.tree_count, int(admissible.sum())


if __name__ == '__main__':
    counts = Counter((trees, number) for i,trees,number in rows())
    for item, multiplicity in sorted(counts.items()):
        print(*item, multiplicity)
