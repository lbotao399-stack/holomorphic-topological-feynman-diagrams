"""Exact series-chain aggregation before scalar-core automorphisms.

This additive CoreCurrentPeriods subclass changes only the order of finite
rational additions. Every original Schwinger monomial receives the same
Dirichlet factorial as in the original implementation. Equal chain powers
are summed first, then graph automorphisms identify the remaining powers.
The source numerator, inserted currents, convergence checks, period tables
and complete-series fallback are inherited without modification.
"""
from fractions import Fraction
from math import factorial, prod

from flint import fmpq

from .forest_periods import CoreCurrentPeriods, ForestPeriodBudgetExceeded

__all__ = ["ChainFirstPeriods"]


class _SparseFactorials(dict):
    def __missing__(self, power):
        value = factorial(power)
        self[power] = value
        return value


class ChainFirstPeriods(CoreCurrentPeriods):
    """Aggregate complete chain coefficients before canonicalizing the core."""

    def _core_polynomials(self, polynomial):
        schwinger = self._schwinger(polynomial)
        self._chain_data()
        result, before_symmetry = {}, {}
        for q, numerator in schwinger.items():
            degree = self.G.loops * q - self.G.ne
            if degree < 0:
                raise ArithmeticError("a nonzero projective numerator has negative total degree")
            maximum = degree + self.G.ne
            # Avoid a large range allocation for a sparse high-power word.
            factorials = (tuple(factorial(k) for k in range(maximum + 1))
                          if maximum <= 4096 else _SparseFactorials())
            raw = {}
            for exponent, coefficient in numerator.terms():
                exponent = tuple(map(int, exponent))
                alpha = tuple(sum(exponent[edge] + 1 for edge in chain.edges)
                              for chain in self._chains)
                factor = coefficient * fmpq(
                    prod(factorials[k] for k in exponent),
                    prod(factorials[k - 1] for k in alpha))
                raw[alpha] = raw.get(alpha, fmpq(0)) + factor
            values = {}
            before_symmetry[q] = sum(bool(value) for value in raw.values())
            for alpha, coefficient in raw.items():
                if not coefficient:
                    continue
                orbit = min(tuple(alpha[edge] for edge in permutation)
                            for permutation in self._permutations)
                values[orbit] = values.get(orbit, fmpq(0)) + coefficient
            values = {alpha: Fraction(int(coefficient.p), int(coefficient.q))
                      for alpha, coefficient in values.items() if coefficient}
            if len(values) > self.max_period_terms:
                raise ForestPeriodBudgetExceeded("weighted core periods exceed max_period_terms")
            result[q] = values
        self.last_stats = {
            "backend": "chain_first_schwinger_numerator",
            "input_current_terms": len(polynomial),
            "bare_numerator_terms": len(self._bare),
            "inserted_numerator_terms": sum(len(value) for value in schwinger.values()),
            "core_edges": len(self._core),
            "used_automorphisms": len(self._permutations),
            "core_period_orbits": sum(len(value) for value in result.values()),
            "chain_power_terms_before_symmetry": before_symmetry,
            "denominator_powers": sorted(result),
            "independent_checks": self.check,
            "current_U_cancellations": dict(self._last_cancellations),
        }
        return result
