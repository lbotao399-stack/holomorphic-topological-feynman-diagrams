"""Exact full-current integration after the complete series-fibre integral.

For every original source and edge, y=y0-C^t M^{-1} C diag(t)y0.
The loop matrix depends only on total chain lengths.  Integrating the
polynomial P and every P*t_j over the chain simplices before multiplying
by adj(M) avoids a large intermediate polynomial on the original edges.
Every Dirichlet factor and the affine term y0 are retained.
"""
from fractions import Fraction
import hashlib
import json
from math import factorial, prod

from flint import fmpq, fmpq_mpoly_ctx

from bulk3_periods import _PeriodWorker, encoded, polynomial_record


class ChainCurrentPeriods:
    def __init__(self, engine):
        self.engine = engine
        self.context = fmpq_mpoly_ctx.get(tuple("R%d" % i for i in range(len(engine._chains))))
        self._ready = False
        self.integrations = []

    def prepare(self):
        if self._ready:
            return
        engine = self.engine
        chains = engine._chains
        count, loops = len(chains), engine.G.loops
        chain_index = {edge: i for i, chain in enumerate(chains) for edge in chain.edges}
        cycles = [tuple(chain.signs[0] * engine.C[r][chain.edges[0]] for r in range(loops))
                  for chain in chains]
        for chain, column in zip(chains, cycles):
            if any(tuple(sign * engine.C[r][edge] for r in range(loops)) != column
                   for edge, sign in zip(chain.edges, chain.signs)):
                raise ArithmeticError("The original cycle columns are not constant along a chain.")
        lengths = self.context.gens()
        matrix = [[sum((column[r] * column[s]) * length
                       for column, length in zip(cycles, lengths))
                   for s in range(loops)] for r in range(loops)]
        self.U = engine._determinant(matrix, self.context)
        self.adjugate = [[None] * loops for _ in range(loops)]
        for r in range(loops):
            for s in range(r, loops):
                minor = [[matrix[i][j] for j in range(loops) if j != r]
                         for i in range(loops) if i != s]
                value = (-1) ** (r + s) * engine._determinant(minor, self.context)
                self.adjugate[r][s] = self.adjugate[s][r] = value
        for r in range(loops):
            for s in range(loops):
                observed = sum(matrix[r][j] * self.adjugate[j][s] for j in range(loops))
                if observed != (self.U if r == s else self.context.constant(0)):
                    raise ArithmeticError("The complete polynomial loop adjugate identity failed.")

        factorials = tuple(factorial(i) for i in range(1 + 24 + engine.G.ne))
        base, moments = {}, [[{} for _ in range(3)] for _ in range(loops)]
        active = []
        for edge in range(engine.G.ne):
            terms = [(r, a, engine.C[r][edge] * engine.flow[edge][a])
                     for r in range(loops) for a in range(3)
                     if engine.C[r][edge] * engine.flow[edge][a]]
            if terms:
                active.append((edge, chain_index[edge], terms))
        for powers, coefficient in engine._bare.terms():
            powers = tuple(map(int, powers))
            alpha = tuple(sum(powers[e] + 1 for e in chain.edges) for chain in chains)
            exponent = tuple(a - 1 for a in alpha)
            weight = coefficient * fmpq(prod(factorials[p] for p in powers),
                                        prod(factorials[a - 1] for a in alpha))
            base[exponent] = base.get(exponent, fmpq(0)) + weight
            for edge, chain, terms in active:
                raised = list(exponent)
                raised[chain] += 1
                raised = tuple(raised)
                first = weight * fmpq(powers[edge] + 1, alpha[chain])
                for r, a, scalar in terms:
                    value = first * fmpq(scalar.numerator, scalar.denominator)
                    target = moments[r][a]
                    target[raised] = target.get(raised, fmpq(0)) + value
        self.base = self.context.from_dict(base)
        self.moments = [[self.context.from_dict(values) for values in row] for row in moments]
        self._ready = True

    def numerator(self, edge, component):
        if not 0 <= edge < self.engine.G.ne or component not in range(3):
            raise ValueError("Specify one original edge and one of the three complete source directions.")
        self.prepare()
        engine = self.engine
        value = engine.flow[edge][component]
        result = self.U * self.base * fmpq(value.numerator, value.denominator)
        for r in range(engine.G.loops):
            if engine.C[r][edge]:
                result -= engine.C[r][edge] * sum(
                    self.adjugate[r][s] * self.moments[s][component]
                    for s in range(engine.G.loops))
        return engine._bounded(result)

    def coefficients(self, numerator):
        combined = {}
        for powers, value in numerator.terms():
            alpha = tuple(int(p) + 1 for p in powers)
            canonical = min(tuple(alpha[e] for e in permutation)
                            for permutation in self.engine._permutations)
            combined[canonical] = combined.get(canonical, fmpq(0)) + value
        return {alpha: Fraction(int(value.p), int(value.q))
                for alpha, value in combined.items() if value}

    def integrate(self, edge, component):
        numerator = self.numerator(edge, component)
        q = 9
        while numerator:
            quotient, remainder = divmod(numerator, self.U)
            if remainder:
                break
            numerator = quotient
            q -= 1
        coefficients = self.coefficients(numerator)
        if any(sum(alpha) != 5 * q for alpha in coefficients):
            raise ArithmeticError("The complete pushed current has the wrong projective degree.")
        scalar_rows = [[list(alpha), str(value)] for alpha, value in sorted(coefficients.items())]
        entry = {"edge": edge, "source_direction": component,
                 "denominator_power": q,
                 "method": "complete_series_fibre_before_linear_current",
                 "complete_core_numerator": polynomial_record(numerator),
                 "weighted_core_periods": len(coefficients),
                 "complete_scalar_coefficients_sha256": hashlib.sha256(
                     json.dumps(scalar_rows, separators=(",", ":")).encode()).hexdigest()}
        if coefficients:
            if self.engine._worker is None:
                self.engine._worker = _PeriodWorker(self.engine.period_timeout)
            value, stats = self.engine._worker.integrate(self.engine._core, q, coefficients)
            entry["period_stats"] = stats
        else:
            value = {}
        entry["value"] = encoded(value)
        self.integrations.append(entry)
        return value
