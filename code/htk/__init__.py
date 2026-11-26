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
    forest_model    marked product-forest recognition and canonical cubic trees
    forest_moments  exact local current messages for cubic product forests
    product_forests cubic messages or the existing all-source series-parallel method
    fan_moments     all-loop fan words and bounded-state displacement Taylor coefficients
    conventions  conversion to the normalizations of the source notes
"""
from .graph import (Graph, polygon, theta, banana, fan, box, graph_A, graph_B, ladder,
                    glued_polygons, from_networkx)
from . import linalg, brute, simplex, spred, oracle, catalog, conventions
from .brute import Density
from .simplex import PiNumber
from .spred import sp_reduce, SPReduction, RatioDensity
from .forest_model import (CubicForest, ProductForestStructure, NotProductForest,
                           NonCubicProductForest, UnsupportedForestSources,
                           recognize_product_forest, forest_graph, forest_counts)
from .fan_moments import (FanTaylorResult, fan_endpoint_moment, fan_current_moment,
                          fan_taylor)

__all__ = ["Graph", "polygon", "theta", "banana", "fan", "box", "graph_A", "graph_B", "ladder",
           "glued_polygons", "from_networkx", "linalg", "brute", "simplex", "spred", "oracle",
           "catalog", "conventions", "Density", "PiNumber", "sp_reduce", "SPReduction", "RatioDensity",
           "kernel_exact", "kernel_numeric", "CubicForest", "ProductForestStructure",
           "NotProductForest", "NonCubicProductForest", "UnsupportedForestSources",
           "recognize_product_forest", "forest_graph", "forest_counts", "ForestMoments",
           "forest_moment", "ForestBudgetExceeded", "FanTaylorResult", "fan_endpoint_moment",
           "fan_current_moment", "fan_taylor", "ProductForestMoments", "product_forest_moment",
           "UnsupportedProductForestReduction"]


def __getattr__(name):
    # Exact sparse forest messages use optional python-flint. Keep the
    # existing engine and the pure-Python fan transfer independent of it.
    if name in ("ForestMoments", "forest_moment", "ForestBudgetExceeded"):
        from . import forest_moments
        value = getattr(forest_moments, name)
        globals()[name] = value
        return value
    if name in ("ProductForestMoments", "product_forest_moment", "UnsupportedProductForestReduction"):
        from . import product_forests
        value = getattr(product_forests, name)
        globals()[name] = value
        return value
    raise AttributeError("module %r has no attribute %r" % (__name__, name))


def __dir__():
    return sorted(set(globals()) | set(__all__))


def kernel_exact(G, lam, T=0, powers=None, plan=None, check=2):
    """K_Gamma(lam; 0) (or the moment with insertion prod (y^i_e)^p, powers = {(e, i): p}) of a
    series-parallel graph by the reduction layer, self-checked against the brute-force layer at
    `check` random rational points.  Fraction (T = 0), PiNumber (T = 1), sympy (symbolic sources)."""
    red = sp_reduce(G, plan=plan, prefer="large" if T else "small")
    return red.integrate(lam, T, powers=powers, check=check)


def kernel_numeric(G, lam, T=0, powers=None, z=None, **kw):
    """(value, error estimate) of K_Gamma(lam; z) by the numerical oracle (no reductions)."""
    return oracle.kernel(G, lam, T, powers=powers, z=z, **kw)
