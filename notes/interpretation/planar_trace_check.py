"""Compare both ordered two-loop colour contractions with the trace and CE-boundary decomposition."""
from fractions import Fraction
from pathlib import Path
import sys

project = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(project / "code"))
import htk  # Pin the shared project engine before vgraph imports it.
sys.path.insert(0, str(project / "deepen/o3/verify"))
from vlib import LP, alt, bracket_poly, cubic, cycles_of, d, fpmul, from_raw, lin
from vgraph import colour, local_orders, theta_data, wick_sign

momenta = fpmul(bracket_poly(0, 1, 4), bracket_poly(2, 3, 4))


def cochain(traces, traceless):
    return alt(from_raw(4, [(c,word,momenta) for c,word in traces], traceless), 4, traceless)


def graph_cochain(shape, traceless):
    graph, external, _ = theta_data(shape)
    order = local_orders(graph, external, 'A')
    parity = wick_sign(graph, external, order)
    traces = [(parity*c,cycles_of(labels))
              for labels,c in colour(graph,external,order,{},traceless).items()]
    return cochain(traces,traceless)


for traceless in (False, True):
    beta = d(cubic('f', 2, traceless), 3, traceless)
    double = cochain([(1,((0,1),(2,3)))],traceless)
    triple = cochain([(1,((0,1),(2,),(3,)))],traceless)
    four = cochain([(1,((0,),(1,),(2,),(3,)))],traceless)
    assert not four
    assert double and bool(triple) != traceless
    assert not d(double,4,traceless) and not d(triple,4,traceless)
    f223 = graph_cochain((2,2,3),traceless)
    f133 = graph_cochain((1,3,3),traceless)
    assert not lin((1,f223),(24,beta),(6*LP.N(),double),(-12,triple))
    assert not lin((1,f133),(1,f223),(-2*LP.N(2),beta))
    bare = lin((Fraction(-3,240),f133),(Fraction(2,240),f223))
    assert not lin((1,bare),(Fraction(1,2)+LP.N(2)/40,beta),
                   (LP.N()/8,double),(Fraction(-1,4),triple))
    for polynomial in (f133,f223,bare):
        for (permutation,monomial),coefficient in polynomial.items():
            boundaries = len(cycles_of(permutation))
            for faces in coefficient.c:
                genus_twice = 3-boundaries-faces
                assert genus_twice in (0,2)
    print(('sl_N' if traceless else 'gl_N'), 'all formal-N coefficients agree:',
          len(f133), 'theta133 monomials,', len(f223), 'theta223 monomials;',
          'CE and ribbon-face checks passed')
