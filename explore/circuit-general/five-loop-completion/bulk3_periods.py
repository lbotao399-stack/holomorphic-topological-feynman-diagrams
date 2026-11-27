"""Complete four-source H=3 five-loop kernels from an exact linear minor.

The canonical source frame represents the entire supported source polynomial:
K(lambda)=det(lambda_s0,lambda_s1,lambda_s2)^5 K(canonical), and every
first-current vector is the same determinant power times a full linear
combination of those three original sources. No source interpolation is used.
"""
from __future__ import annotations

from fractions import Fraction
import hashlib
from itertools import combinations
import json
from math import factorial, prod
from pathlib import Path
import sys
from types import SimpleNamespace

from flint import fmpq, fmpq_mpoly_ctx

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "code"))
sys.path.insert(0, str(ROOT / "explore/nonSP/verify/k4_periods"))
from htk import brute
from htk.graph import Graph
from htk.treeflow import tree_flow
from htk.forest_periods import CoreCurrentPeriods, _PeriodWorker
import current_minor


def encoded(value):
    return {key: str(coefficient) for key, coefficient in sorted(value.items()) if coefficient}


def add(first, second, factor=1):
    result = dict(first)
    for label, coefficient in second.items():
        result[label] = result.get(label, Fraction(0)) + factor * coefficient
    return {label: value for label, value in result.items() if value}


def polynomial_record(value):
    rows = [[list(map(int, powers)), str(coefficient)] for powers, coefficient in value.terms()]
    rows.sort()
    return {"terms": len(rows),
            "sha256": hashlib.sha256(json.dumps(rows, separators=(",", ":")).encode()).hexdigest()}


class Bulk3Periods:
    """Budgeted exact arithmetic; all original edges and three directions kept."""
    _bounded = CoreCurrentPeriods._bounded
    _determinant = CoreCurrentPeriods._determinant
    _chain_data = CoreCurrentPeriods._chain_data

    def __init__(self, row, *, tree=None, max_terms=2000000, period_timeout=180):
        self.G = Graph(row["nv"], row["edges"])
        if self.G.loops != 5 or self.G.ne != 16 or len(row["S0"]) != 4:
            raise ValueError("Expected an H=3 five-loop graph with four original source vertices.")
        self.source_vertices = tuple(row["S0"])
        self.sources = [[Fraction(0)] * 3 for _ in range(self.G.nv)]
        for component, vertex in enumerate(self.source_vertices[:3]):
            self.sources[vertex][component] = Fraction(1)
        self.sources[self.source_vertices[-1]] = [Fraction(-1)] * 3
        self.tree = tuple(self.G.spanning_tree() if tree is None else tree)
        self.C = self.G.cycle_matrix(self.tree)
        self.fundamental_cycles = self.G.fundamental_cycles(self.tree)
        self.flow = tree_flow(self.G, self.sources, self.tree)
        self.max_terms = max_terms
        self.period_timeout = period_timeout
        self._chains = self._core = self._permutations = None
        self._U = self._currents = None
        self._worker = None
        self.context = fmpq_mpoly_ctx.get(tuple("t%d" % e for e in range(self.G.ne)))
        proxy = SimpleNamespace(n=self.G.loops, N=self.G.ne, C=self.C)
        numerator = current_minor.numerator(
            proxy, self.sources, [[int(value) for value in row] for row in self.flow],
            determinant=self._determinant)
        self._bare = self.context.from_dict({
            tuple(map(int, powers)): fmpq(int(coefficient))
            for powers, coefficient in numerator.terms()})
        self._bounded(self._bare)
        self.density_checks = []
        for offset in (1, 3):
            point = [Fraction(edge + offset + 1, edge + 1) for edge in range(self.G.ne)]
            U = brute.symanzik_U(self.G, point)
            expected = brute.density(self.G, point, self.sources).rational() * U ** 8
            observed = self._evaluate(self._bare, point)
            if observed != expected:
                raise ArithmeticError("Full H=3 numerator differs from the original grounded-Laplacian density.")
            self.density_checks.append({"resistances": list(map(str, point)),
                                        "U": str(U), "complete_numerator": str(observed)})
        self._chain_data()
        self.integrations = []

    @staticmethod
    def _evaluate(polynomial, point):
        value = polynomial(*(fmpq(t.numerator, t.denominator) for t in point))
        return Fraction(int(value.p), int(value.q))

    def close(self):
        if self._worker is not None:
            self._worker.close()
            self._worker = None

    def tree_polynomials(self):
        if self._currents is not None:
            return
        scalar = {}
        values = [[{} for _ in range(3)] for _ in self.G.edges]
        for tree in combinations(range(self.G.ne), self.G.nv - 1):
            parent = list(range(self.G.nv))
            def find(vertex):
                while parent[vertex] != vertex:
                    parent[vertex] = parent[parent[vertex]]
                    vertex = parent[vertex]
                return vertex
            valid = True
            for edge in tree:
                left, right = map(find, self.G.edges[edge])
                if left == right:
                    valid = False
                    break
                parent[left] = right
            if not valid:
                continue
            powers = tuple(int(edge not in tree) for edge in range(self.G.ne))
            scalar[powers] = fmpq(1)
            flow = tree_flow(self.G, self.sources, tree)
            for edge in tree:
                for component in range(3):
                    value = flow[edge][component]
                    if value:
                        values[edge][component][powers] = fmpq(value.numerator, value.denominator)
        self._U = self.context.from_dict(scalar)
        self._currents = tuple(tuple(self.context.from_dict(direction) for direction in row)
                               for row in values)
        cotree = self.G.cotree(self.tree)
        for edge in range(self.G.ne):
            for component in range(3):
                reconstructed = self._U * fmpq(self.flow[edge][component].numerator,
                                                self.flow[edge][component].denominator)
                for r, chord in enumerate(cotree):
                    reconstructed += self.fundamental_cycles[r][edge] * self._currents[chord][component]
                if reconstructed != self._currents[edge][component]:
                    raise ArithmeticError("The complete affine original-current identity failed.")
        point = [Fraction(edge + 2, edge + 3) for edge in range(self.G.ne)]
        U = brute.symanzik_U(self.G, point)
        original = brute.currents(self.G, point, self.sources)
        if self._evaluate(self._U, point) != U:
            raise ArithmeticError("Spanning-tree Kirchhoff polynomial failed.")
        for edge in range(self.G.ne):
            for component in range(3):
                if self._evaluate(self._currents[edge][component], point) != U * original[edge][component]:
                    raise ArithmeticError("A complete original current numerator failed.")

    def scalar_coefficients(self, numerator):
        maximum = max((sum(powers) for powers, _ in numerator.terms()), default=0) + self.G.ne
        factorials = tuple(factorial(j) for j in range(maximum + 1))
        raw = {}
        for powers, coefficient in numerator.terms():
            powers = tuple(map(int, powers))
            alpha = tuple(sum(powers[edge] + 1 for edge in chain.edges) for chain in self._chains)
            factor = coefficient * fmpq(prod(factorials[k] for k in powers),
                                       prod(factorials[k - 1] for k in alpha))
            raw[alpha] = raw.get(alpha, fmpq(0)) + factor
        combined = {}
        for alpha, value in raw.items():
            if value:
                canonical = min(tuple(alpha[edge] for edge in permutation)
                                for permutation in self._permutations)
                combined[canonical] = combined.get(canonical, fmpq(0)) + value
        return {alpha: Fraction(int(value.p), int(value.q)) for alpha, value in combined.items() if value}

    def integrate(self, edge=None, component=None):
        q, numerator = 8, self._bare
        if edge is not None:
            if component not in range(3) or edge not in range(self.G.ne):
                raise ValueError("Every insertion must specify an original edge and a full source direction.")
            self.tree_polynomials()
            q += 1
            numerator = self._bounded(numerator * self._currents[edge][component])
            while numerator:
                quotient, remainder = divmod(numerator, self._U)
                if remainder:
                    break
                numerator = quotient
                q -= 1
        coefficients = self.scalar_coefficients(numerator)
        scalar_rows = [[list(alpha), str(value)] for alpha, value in sorted(coefficients.items())]
        if any(sum(alpha) != self.G.loops * q for alpha in coefficients):
            raise ArithmeticError("A complete scalar target has wrong projective degree.")
        entry = {"edge": edge, "source_direction": component, "denominator_power": q,
                 "complete_numerator": polynomial_record(numerator),
                 "weighted_core_periods": len(coefficients),
                 "complete_scalar_coefficients_sha256": hashlib.sha256(
                     json.dumps(scalar_rows, separators=(",", ":")).encode()).hexdigest()}
        if coefficients:
            if self._worker is None:
                self._worker = _PeriodWorker(self.period_timeout)
            value, stats = self._worker.integrate(self._core, q, coefficients)
            entry["period_stats"] = stats
        else:
            value = {}
        entry["value"] = encoded(value)
        self.integrations.append(entry)
        return value

    def complete_first_currents(self, bare, *, direct=False, progress=None):
        if direct:
            values = [[self.integrate(edge, component) for component in range(3)]
                      for edge in range(self.G.ne)]
        else:
            cotree = tuple(edge for edge in range(self.G.ne) if edge not in self.tree)
            chords = []
            for edge in cotree:
                chords.append([self.integrate(edge, component) for component in range(3)])
                if progress:
                    progress(edge)
            values = []
            for edge in range(self.G.ne):
                row = []
                for component in range(3):
                    value = {label: self.flow[edge][component] * coefficient
                             for label, coefficient in bare.items()
                             if self.flow[edge][component] * coefficient}
                    for r in range(self.G.loops):
                        value = add(value, chords[r][component], self.fundamental_cycles[r][edge])
                    row.append(value)
                values.append(row)
        for vertex in range(self.G.nv):
            for component in range(3):
                value = {}
                for edge, (left, right) in enumerate(self.G.edges):
                    value = add(value, values[edge][component], int(left == vertex) - int(right == vertex))
                expected = {label: self.sources[vertex][component] * coefficient
                            for label, coefficient in bare.items()
                            if self.sources[vertex][component] * coefficient}
                if value != expected:
                    raise ArithmeticError("A complete integrated original-source identity failed.")
        return values
