"""Compute the five-dimensional sl_N order-four Taylor cochain in the heat-kernel gauge."""
import sys

from bipartite_cut import solve as bipartite_coefficients
from doubled_kite import doubled_coefficients
from odd_cut import HERE, replay
from reconstruct_check import series_parallel_coefficients
from reconstruct_prism.evaluate import prism_coefficients
import sympy as sp

sys.path.insert(0,str(HERE.parent))
from order_four import N, project


def solve(directory,executable):
    rows = [series_parallel_coefficients(),doubled_coefficients(),prism_coefficients(),
            bipartite_coefficients(directory,executable)]
    bare = sp.expand(sum(s*N**k for row in rows for k,s in row.items()))
    # With epsilon*U=1, the four-loop factor in units of (epsilon*pi*nu)^2
    # is 1/pi^2. The differential on the cyclic probe has value 16.
    coefficient = sp.expand(project()+bare/(16*sp.pi**2))
    expected = (sp.Rational(1,24)
                +N**2*(sp.Rational(1102559,5443200)-35701*sp.pi**2/1935360)
                +N**4*(sp.Rational(4216,212625)-130727*sp.pi**2/92897280))
    assert sp.expand(coefficient-expected)==0
    assert bare.coeff(N,0)==-172*sp.pi**2/945
    assert all(coefficient.coeff(N,power).subs(sp.pi**2,10)>0 for power in (0,2,4))
    print('bare',sp.collect(bare,N))
    print('primitive in units of (epsilon*pi*nu)^2',sp.collect(coefficient,N))


if __name__=='__main__':
    replay.run(solve)
