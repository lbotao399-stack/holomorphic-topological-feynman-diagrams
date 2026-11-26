"""Integrate every series-chain simplex exactly before any core integration.

This composes the source-preserving determinant identity in ``series`` with
the rational polynomial cubature in ``dirichlet``.  The result is an exact
function of the retained positive core resistances.  Only a one-chain graph
has no further projective core integral; in general this is not a completed
Feynman period.  T=0 and finite monomials in original edge currents only.
"""
from dataclasses import dataclass
from fractions import Fraction

from .dirichlet import product_simplex_rule
from .series import SeriesReduction

__all__ = ["FibreIntegral", "integrate_fibre"]


@dataclass(frozen=True)
class FibreIntegral:
    """Exact pushed-forward density and reproducible work counts."""
    value: Fraction
    nodes: int
    degree: int
    simplex_dimensions: tuple
    determinant_size: int
    zero_from_sources: bool


def integrate_fibre(reduction, resistances, powers=None, *,
                    max_nodes=100000, max_terms=1000000):
    """Return the integral over all chain ratios at fixed core resistances.

    ``reduction`` is a ``SeriesReduction(G, sources)``. ``powers`` maps
    (original edge ID, holomorphic component) to a nonnegative power.
    The ordinary simplex measures, all sources, orientation signs, current
    insertions and the factors R_chain**(length-1) are included exactly.
    The safe total polynomial degree is number_of_chains-1+sum(powers).

    A core resistance may be fixed to one for subsequent projective
    integration.  The caller must still perform that remaining integral.
    """
    if not isinstance(reduction, SeriesReduction):
        raise TypeError("reduction must be a SeriesReduction")
    resistances = tuple(Fraction(value) for value in resistances)
    dimensions = tuple(chain.length - 1 for chain in reduction.chains)
    # Validate the core point even when the source columns certify zero.
    centre = tuple(tuple(Fraction(1, chain.length) for _ in chain.edges)
                   for chain in reduction.chains)
    reduction.times(resistances, centre)
    degree = reduction.fibre_degree(powers)
    if reduction.zero:
        return FibreIntegral(Fraction(0), 0, degree, dimensions,
                             reduction.determinant_size, True)
    rule = product_simplex_rule(dimensions, degree,
                                max_nodes=max_nodes, max_terms=max_terms)
    value = rule.integrate(lambda point: reduction.fibre_density(resistances, point, powers))
    return FibreIntegral(value, len(rule), degree, dimensions,
                         reduction.determinant_size, False)
