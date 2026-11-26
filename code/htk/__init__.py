"""htk: kernel engine for holomorphic-topological graph integrals on R^T x C^H (T in {0, 1}).

Layers (see CONVENTIONS.md):
    graph      Graph encoding, named graphs
    brute      exact brute-force layer at rational Schwinger times (no reductions)
    oracle     independent numerical integration over the Schwinger simplex
    spred      series-parallel reduction: ratio coordinates, polynomial currents, densities
    simplex    exact Dirichlet integrals over products of simplices
    series     constant-source-column contraction on arbitrary series cores
    dirichlet  exact rational sparse cubature on products of simplices
    fibres     series-fibre integration before core reconstruction
    support    graph-rank bounds on numerator monomials
    catalog    graph utilities (enumeration, series-parallel test, parallel-chain reducibility)
    conventions  conversion to the normalizations of the source notes
"""
from .graph import (Graph, polygon, theta, banana, fan, box, graph_A, graph_B, ladder,
                    glued_polygons, from_networkx)
from . import linalg, brute, simplex, spred, oracle, catalog, conventions
from .brute import Density
from .simplex import PiNumber
from .spred import sp_reduce, SPReduction, RatioDensity

__all__ = ["Graph", "polygon", "theta", "banana", "fan", "box", "graph_A", "graph_B", "ladder",
           "glued_polygons", "from_networkx", "linalg", "brute", "simplex", "spred", "oracle",
           "catalog", "conventions", "Density", "PiNumber", "sp_reduce", "SPReduction", "RatioDensity",
           "kernel_exact", "kernel_numeric"]


def kernel_exact(G, lam, T=0, powers=None, plan=None, check=2):
    """K_Gamma(lam; 0) (or the moment with insertion prod (y^i_e)^p, powers = {(e, i): p}) of a
    series-parallel graph by the reduction layer, self-checked against the brute-force layer at
    `check` random rational points.  Fraction (T = 0), PiNumber (T = 1), sympy (symbolic sources)."""
    red = sp_reduce(G, plan=plan, prefer="large" if T else "small")
    return red.integrate(lam, T, powers=powers, check=check)


def kernel_numeric(G, lam, T=0, powers=None, z=None, **kw):
    """(value, error estimate) of K_Gamma(lam; z) by the numerical oracle (no reductions)."""
    return oracle.kernel(G, lam, T, powers=powers, z=z, **kw)
