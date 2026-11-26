"""Exact original-kernel and Dirichlet checks for the H=2 odd triangle."""
from fractions import Fraction
from math import factorial, prod
from pathlib import Path
import sys

import networkx as nx
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "code"))
from enumerate_graphs import tight_ok
from htk import Graph, brute
from htk.simplex import dirichlet
from htk.linalg import det
from dimension_raise import split_times, subdivide


edges = [(0,3),(3,1),(0,4),(4,1),(1,5),(5,2),(1,6),(6,2),
         (2,7),(7,0),(2,8),(8,9),(9,0)]
graph = Graph(10, edges)


def sources(p):
    return [(-1,p-1),(-1,-p-1),(-1,-1),(1,0),(0,1),(1,0),(0,1),(1,0),(1,0),(-1,1)]


def check_density():
    coords = sp.symbols("x y w1 w2 w3 a1 b1 a2 b2 a3 b3 c3")
    x,y,w1,w2,w3,a1,b1,a2,b2,a3,b3,c3 = coords
    z = 1-x-y
    times = []
    for scale,share,a,second_path in [
        (x,w1,a1,(1-b1,b1)), (y,w2,a2,(1-b2,b2)),
        (z,w3,a3,(1-b3-c3,b3,c3))]:
        times.extend((scale*(1-a)/share, scale*a/share))
        times.extend(scale*b/(1-share) for b in second_path)
    times = sp.Matrix(times)
    jacobian = times.row_join(times.jacobian(coords))
    radicand = x*y*z/(w1*(1-w1)*w2*(1-w2)*w3*(1-w3))
    common = z*w1*(1-w1)*w2*(1-w2)*w3/4
    fixtures = [(2,3,4,6,8,3,5,6,2,4,1,2), (2,6,8,3,1,2,5,7,8,4,2,3)]
    points = [dict(zip(coords, (sp.Rational(i,10) for i in row))) for row in fixtures]
    sign_point = dict.fromkeys(coords,sp.Rational(1,4))
    sign_point.update({x:sp.Rational(1,5), y:sp.Rational(3,10),
                       w1:sp.Rational(1,2), w2:sp.Rational(1,2), w3:sp.Rational(1,2)})
    for point in [sign_point, *points]:
        time_values = [Fraction(value) for value in times.subs(point)]
        absolute_jacobian = abs(Fraction(jacobian.subs(point).det()))
        for p in (-3,0,1):
            expected = common*((1/x-1/z)*(b2-c3)-(1/y-1/z)*(p+b1-c3))
            rational = expected.subs(point)
            observed = brute.density(graph,time_values,sources(p),T=1).scaled(absolute_jacobian)
            assert observed.radicand == Fraction(radicand.subs(point))
            assert observed.sign() == sp.sign(rational)
            assert observed.square() == Fraction(rational**2 * radicand.subs(point))
            if point == sign_point and p == 1:
                assert observed.coeff == Fraction(-1,100)
                assert observed.radicand == Fraction(48,25)
    print("Original signed density agrees at 9 graph/source points; sign anchor = -sqrt(3)/240")


def check_period():
    half = Fraction(1,2)
    central = dirichlet((half,half,half))
    angular = dirichlet((half,half))
    third = dirichlet((half,-half))
    assert central == (Fraction(2,105),1)
    assert angular == (Fraction(1,8),1)
    assert third == (Fraction(1,2),1)
    # Each sin-squared angle contributes a Jacobian factor of 1/2.
    # The final path simplex contributes 1/2, and the curvature prefactor is 4.
    value = -4 * central[0] * (angular[0]/2)**2 * (third[0]/4)
    power = central[1] + 2*angular[1] + third[1]
    assert (value,power) == (Fraction(-1,26880),4)
    print(f"p=1: {value} pi^{power}")



def check_subdivisions():
    selected = (0,4,8,10)
    cycle = graph.cycle_matrix()
    assert det([[row[e] for e in selected] for row in cycle]) == -1
    # First chamber anchor from the original integral.
    times = [Fraction(v,20) for v in (6,2,6,2,9,3,9,3,15,5,10,5,5)]
    old_source = sources(1)[:-1]
    old = brute.density(graph,times,old_source,T=1)
    period = brute.symanzik_U(graph,times,cycle)
    raised = subdivide(graph,selected)
    lifted = [list(row)+[0] for row in old_source] + [[0,0,1] for _ in selected]
    new_times = split_times(times,selected,[Fraction(i,7) for i in range(1,5)])
    new = brute.density(raised,new_times,lifted,T=1)
    assert new.same_as(old.scaled(-1/period))
    assert prod(times[e] for e in selected)/period == Fraction(27,1024)
    half = Fraction(1,2)
    for k,known in enumerate((Fraction(-1,26880),Fraction(1,3870720),Fraction(-25,39239811072))):
        central,cp = dirichlet((half,half,k+half))
        angular,ap = dirichlet((half,k+half))
        third,tp = dirichlet((half,-half))
        value = ((-1)**(k+1)*4*(k+1)*central*(angular/2)**2*(third/2)
                 /((k+1)**4*(k+2)*factorial(k)**4))
        assert cp+2*ap+tp == 4 and value == known
        print(f"H={2+k}: explicit cotree family has K(1)/pi^4 = {value}")


if __name__ == "__main__":
    assert tight_ok(nx.Graph(edges), 3)
    check_density()
    check_period()
    check_subdivisions()
