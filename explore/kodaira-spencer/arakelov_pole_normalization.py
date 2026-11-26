"""Check the planar graph height, equal-strength fourfold vertex and axial D3 pole."""
from fractions import Fraction as F
from itertools import combinations, product
from math import factorial
from pathlib import Path
from subprocess import run
import sys
from tempfile import TemporaryDirectory

import networkx as nx

from arakelov_graph_periods import graph_periods, metric_numerators
from bulk_three import evaluate as bulk_period, forest_data, integrate_core
from collinear_abel_moments import genus_three
from four_loop_arakelov import period_polynomial
from planar_bulk import evaluate as planar_period, transfer
from quartic_planar_cores import cubic_cores
from quartic_scalar_response import BIPARTITE, interval, total_channels
from square_three import tight_graphs
from worldsheet_fourfold_positivity import fourfold_coefficients
from zero_wave_descendant import IIB

sys.path.insert(0,str(IIB.parent/'twisted-holography-N4'/'code'))
from flux_three_tangent_sources import momentum_frame, reduced_traces, repeated_soft_fourfold


def check_metric_identity():
    cores = [((0,1),)*3, ((0,1),(0,1),(0,2),(1,3),(2,3),(2,3)),
             tuple(combinations(range(4),2)), *cubic_cores()]
    terms = []
    for core in cores:
        vertices = max(max(edge) for edge in core)+1
        genus = len(core)-vertices+1
        ctx,u,forest = forest_data(core,vertices)
        _,phi = metric_numerators(core,ctx,u,forest)
        planar = transfer(core,ctx,u,forest,genus)
        assert genus*planar == -u**2*phi
        for times in (list(range(1,len(core)+1)),[1]*(len(core)-1)+[1000]):
            tree = u(*times)
            probabilities = [t*u.derivative(i)(*times)/(genus*tree) for i,t in enumerate(times)]
            diagonal = sum(t*p*p for t,p in zip(times,probabilities))*genus*(genus-1)/6
            assert phi(*times)/tree >= diagonal and 54*genus*(genus-1)*diagonal >= sum(times)
        terms.append(len(planar))
    return terms


def axial_fivefold(genus, exe):
    # At A=B=0 after stripping A²B, take -2[x²y²] of the signed fivefold trace.
    grid = (-1,0,1)
    frames = [momentum_frame(genus,((0,0),(0,0),(x,y)),third=(0,0))
              for x,y in product(grid,repeat=2)]
    traces = reduced_traces(genus,frames,exe)
    second = (F(1,2),-1,F(1,2))
    return tuple(-2*sum(x*y*t for (x,y),t in zip(product(second,repeat=2),column))
                 for column in zip(*traces))


def four_loop_projection():
    # W=18 I5-3 I6 and the centered trace is 12 I6. Cancel W componentwise.
    for (k,j,w),centered in zip(((1,0,0),(0,1,0),(0,0,1)),fourfold_coefficients((1,1,-2))):
        six = centered/12
        five = (w+3*six)/18
        assert 3*six-22*five == F(11,576)*k+j
    coefficients = tuple(F(11,576)*k+j for k,j in total_channels(BIPARTITE))
    assert tuple(10**24*x for x in interval(coefficients)) == (-73152318930388,-73152318930387)
    return coefficients


def check_vertices(genus, exe):
    # Recorded one-descendant threefold simplex periods in the 1,zeta(3) channels.
    descendant = {2:(F(1,9720),F()),3:(F(411065203,94152554496000),-F(3617,1009008000))}[genus]
    rows = []
    for graph in tight_graphs(3,genus):
        if max(dict(nx.Graph(graph.edges).degree()).values()) != 3:
            continue
        core,poly = period_polynomial(graph.edges)
        auxiliary = integrate_core(core,poly,genus,2*genus+4)
        _,phi = graph_periods(graph.edges)
        planar = planar_period(graph.edges)
        assert phi == tuple(-2*j for j in planar)
        rows.append((auxiliary,bulk_period(graph.edges),phi,planar))
    assert len(rows) == {2:2,3:9}[genus]
    auxiliary,bulk,phi,planar = (tuple(map(sum,zip(*group))) for group in zip(*rows))
    equal = tuple(-F(factorial(genus-1)*(genus-1),2*factorial(genus)**2)*k
                  -F(factorial(genus-1),12)*(a-4*d)
                  for a,k,d in zip(auxiliary,bulk,descendant))
    five_mark = tuple((a+v-4*d)/8 - k/(2*factorial(genus)**2)
                      for a,v,d,k in zip(auxiliary,phi,descendant,bulk))
    axis = repeated_soft_fourfold(genus,(1,0,0,0))
    assert equal == (axis if genus==2 else genus_three(1,1))
    assert axis == tuple(factorial(genus-2)*v for v in five_mark)
    # At genus two the cube of the one-coordinate Betti form is zero.
    six_mark = (F(),F()) if genus==2 else genus_three(1,-1)
    combination = tuple(3*h-2*(3*genus-1)*s/factorial(genus-2) for h,s in zip(six_mark,axis))
    assert combination == tuple(F(3*genus-1,factorial(genus)**2)*k+j for k,j in zip(bulk,planar))
    pole = axial_fivefold(genus,exe)
    assert pole == tuple(4*factorial(genus-1)**2*v for v in five_mark[:len(pole)])
    print(f'Genus {genus}: equal-strength {equal}; axial D3 {pole}; five/six-mark {combination}',flush=True)


if __name__ == '__main__':
    print('Nonzero coefficients of the identical metric polynomials:',check_metric_identity(),flush=True)
    print('Genus-four five/six-mark difference:',four_loop_projection(),flush=True)
    with TemporaryDirectory() as tmp:
        exe = Path(tmp)/'nine'
        run(['c++','-std=c++17','-O3',str(IIB/'operator'/'nine_sextics.cpp'),'-o',str(exe)],check=True)
        for g in (2,3):
            check_vertices(g,exe)
