"""Exact bipartite five-hard Ward period for a triple and pair on disjoint paths."""

from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations, islice
from math import factorial, prod
from pathlib import Path
from tempfile import TemporaryDirectory

from flint import fmpq
from five_hard_mixed_double import chamber_polynomials
from fourfold_bulk import symmetries
from fourfold_quartic_prism import project_slices
from mixed_bipartite_factor_contract import Periods, native, packed
from quartic_seven_plane_lift import cubic_cores
from quartic_seven_plane_prism import graph_data


MARKS = (0, 0, 0, 5, 5)
CHANNELS = (F(-21391605787354675307, 103526382821621760000000),
            F(-31533732279499, 57058191590400000), F(3396079, 4043239200))


def calculate(folder):
    core = cubic_cores()[4]
    data = graph_data(MARKS, core)
    roots = []
    for pair in combinations(range(5), 2):
        ds, voltage = chamber_polynomials(data, pair)
        assert all(d == ds[0] for d in ds)
        if not ds[0]:
            continue
        for row in roots:
            if ds[0] == row[0] or ds[0] == -row[0]:
                row[1] += voltage
                break
        else:
            roots.append([ds[0], voltage])
    adj = sum(map(sum, data[3]))
    periods = Periods(native(folder))
    weights = defaultdict(fmpq)
    bad = defaultdict(fmpq)
    for number, (root, voltage) in enumerate(roots):
        for poly in project_slices(data, [(root, adj)], voltage, 1, stars=()):
            terms = iter(poly.terms())
            while block := list(islice(terms, 8192)):
                raw = [packed(ex) for ex, _ in block]
                keys = periods.charts(raw)
                for key, r, (_, c) in zip(keys, raw, block, strict=True):
                    if key:
                        weights[int(key)] += c
                    else:
                        bad[r] += c
        print('bipartite root', number+1, 'of', len(roots), 'scalar powers', len(weights), flush=True)
    assert not any(bad.values())
    weights = {k: v for k, v in weights.items() if v}
    periods.read(weights)
    actions = symmetries(core)
    images = {tuple(sorted(a[h] for h in MARKS)) for a in actions}
    # 12 assignments per chosen pair; the forest voltage numerator is 4 U K.
    factor = F(3*len(images), len(actions)*factorial(12))
    total = [fmpq(0)]*3
    for k, c in weights.items():
        moment = prod(factorial(((k >> (4*i)) & 15)-1) for i in range(9))
        total = [v+c*moment*x for v, x in zip(total, periods.values[k], strict=True)]
    return tuple(factor*F(str(v)) for v in total)


if __name__ == '__main__':
    with TemporaryDirectory() as folder:
        value = calculate(Path(folder))
        assert value == CHANNELS
        print('bipartite', MARKS, *value, flush=True)
