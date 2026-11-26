"""Compare every current-density and displacement coefficient before integration."""
from pathlib import Path
import sys

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import q3lib
from vops import GRAPH_SPECS, Poly, bp_mul, density_and_X, det_lin


def verify(name, orientation):
    spec = GRAPH_SPECS[name]
    source = sp.symbols('p1:8')
    # The Kirchhoff solution in each spatial component is the same linear map.
    q3lib.LAM = [[p if v in spec['outputs'] else 0] for v, p in enumerate(source, 1)]
    chart = q3lib.Reduction(name)
    n = spec['nr']
    ratios = q3lib.R[:n]
    def polynomial(expr):
        return Poly(dict(sp.Poly(expr, *ratios).terms()) if expr != 0 else {}, n)
    def linear(current):
        expanded = sp.expand(current[0])
        return {v: polynomial(expanded.coeff(source[v-1])) for v in spec['outputs']}

    one = Poly.const(1, n)
    expected, shifts = density_and_X(name, {v: {v: one} for v in spec['outputs']}, n, 0)
    density = {((), ()): Poly.const(orientation, n)}
    for step in chart.steps[:-1]:
        rho = step['r'][0]
        weight = rho**2 if step['kind'] == 'chain' else rho*(1-rho)
        factor = det_lin(*(linear(step[key]) for key in ('Y', 'lv', 'lw')))
        density = bp_mul(density, {k: v*polynomial(weight) for k,v in factor.items()})
    density = bp_mul(density, det_lin(*(linear(y) for y in chart.steps[-1]['lam'])))
    assert set(density) == set(expected)
    assert all(density[key].t == expected[key].t for key in density)

    for j, (edge, sign) in enumerate(chart.spokes, 1):
        current = linear(chart.y[edge])
        for v in spec['outputs']:
            assert (current[v]*(-sign) - shifts[v][j]).iszero()
    print(name, len(density), 'bracket monomials;', len(chart.spokes)*len(shifts),
          'displacement entries')


if __name__ == '__main__':
    for name, sign in zip(('box', 'A', 'B', 'G1', 'G2', 'G3', 'G4', 'G5'),
                          (1,-1,1,1,-1,1,1,-1)):
        verify(name, sign)
