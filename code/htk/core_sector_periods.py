"""Boundary-preserving, certificate-producing scalar-period relations.

The scalar graph has a fixed Kirchhoff polynomial.  A Hepp sector turns
every convergent integer-weight monomial into a regular cube form.  On
that cube the complete relative relation schema is

    integral partial_i[x_i (1-x_i) x**a / U**q] = 0.

This module constructs those identities and searches a *specified finite*
set for a reduction to caller-supplied master forms.  A successful result
contains every rational coefficient and primitive and is checked by an
exact rational-function identity.  A failed finite search does not certify
master independence, stabilization, or a nonzero residual period.

No graph-period table, numerical integration, regulator, or arithmetic
independence assumption is used.  The scalar sector measure is positive;
any oriented graph/source coefficient remains an external exact scalar.
"""
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import comb
from numbers import Integral
from types import MappingProxyType

from .graph import Graph

__all__ = ["CoreSector", "CubeForm", "CubePeriodRelations", "CubeReduction",
           "SectorBudgetExceeded", "ReductionNotFound", "check_convergence"]


class SectorBudgetExceeded(RuntimeError):
    """A declared exact-computation budget was exceeded."""


class ReductionNotFound(RuntimeError):
    """The stated finite relation search did not eliminate the whole target."""

    def __init__(self, residual, statistics):
        super().__init__("finite relation search left a complete residual; "
                         "no period value or independence claim is returned")
        self.residual = residual
        self.statistics = statistics


def _integer(value, name, minimum=0):
    if isinstance(value, bool) or not isinstance(value, Integral) or value < minimum:
        raise ValueError(f"{name} must be an integer at least {minimum}")
    return int(value)


def _rational(value):
    if isinstance(value, (float, complex, bool)):
        raise TypeError("exact rational coefficients are required")
    if hasattr(value, "p") and hasattr(value, "q"):
        return Fraction(int(value.p), int(value.q))
    return Fraction(value)


def _clean(poly, variables):
    answer = {}
    for powers, coefficient in poly.items():
        powers = tuple(_integer(a, "monomial exponent") for a in powers)
        if len(powers) != variables:
            raise ValueError("one exponent is required per variable")
        coefficient = _rational(coefficient)
        if coefficient:
            answer[powers] = answer.get(powers, Fraction()) + coefficient
    return {a: c for a, c in answer.items() if c}


def _add(target, source, multiplier=Fraction(1), budget=None):
    for key, value in source.items():
        target[key] = target.get(key, Fraction()) + multiplier * value
        if not target[key]:
            del target[key]
        if budget is not None and len(target) > budget:
            raise SectorBudgetExceeded("polynomial addition exceeded max_terms")


def _shift(powers, coordinate, amount=1):
    result = list(powers)
    result[coordinate] += amount
    return tuple(result)


def _multiply(left, right, budget):
    result = {}
    for a, ca in left.items():
        for b, cb in right.items():
            key = tuple(x + y for x, y in zip(a, b))
            result[key] = result.get(key, Fraction()) + ca * cb
            if not result[key]:
                del result[key]
            if len(result) > budget:
                raise SectorBudgetExceeded("polynomial product exceeded max_terms")
    return result


def _monomials(variables, maximum):
    if not variables:
        yield ()
        return
    for first in range(maximum + 1):
        for rest in _monomials(variables - 1, maximum - first):
            yield (first,) + rest


def _poly_record(poly):
    return [[list(a), str(c)] for a, c in sorted(poly.items())]


def _read_poly(rows):
    result = {}
    for powers, coefficient in rows:
        powers = tuple(powers)
        if powers in result:
            raise ValueError("duplicate monomial in certificate")
        result[powers] = _rational(coefficient)
    return result


@dataclass(frozen=True)
class CubeForm:
    """An exact numerator over the containing cube's fixed U**power."""
    numerator: object
    power: int

    def record(self):
        return {"power": self.power, "numerator": _poly_record(self.numerator)}


@dataclass(frozen=True)
class CubeReduction:
    """The complete identity target = sum masters + sum derivatives."""
    denominator: object
    target: CubeForm
    masters: tuple
    coefficients: tuple
    # Each tuple is (rational coefficient, coordinate, exponent tuple, power).
    primitives: tuple
    statistics: object

    def record(self):
        return {"schema": "htk.cube-relative-reduction.v1",
                "denominator": _poly_record(self.denominator),
                "target": self.target.record(),
                "masters": [form.record() for form in self.masters],
                "coefficients": [str(c) for c in self.coefficients],
                "primitives": [[str(c), i, list(a), q]
                               for c, i, a, q in self.primitives],
                "statistics": dict(self.statistics)}


class CubePeriodRelations:
    """Exact relative identities for a fixed positive cube denominator.

    Nonnegative rational coefficients and a positive constant coefficient
    certify that U has no zero on the entire closed cube.  No positive
    real test-point heuristic substitutes for this condition.
    """

    def __init__(self, denominator, *, max_terms=100000):
        if not denominator:
            raise ValueError("the denominator must be a nonzero polynomial")
        self.variables = len(next(iter(denominator)))
        if not self.variables:
            raise ValueError("at least one cube coordinate is required")
        self.max_terms = _integer(max_terms, "max_terms", 1)
        self.zero = (0,) * self.variables
        denominator = _clean(denominator, self.variables)
        if len(denominator) > self.max_terms:
            raise SectorBudgetExceeded("denominator exceeded max_terms")
        if denominator.get(self.zero, 0) <= 0 or any(c < 0 for c in denominator.values()):
            raise ValueError("need nonnegative denominator coefficients and a positive constant")
        self.denominator = MappingProxyType(denominator)
        self._powers = {0: {self.zero: Fraction(1)}, 1: dict(denominator)}

    def form(self, numerator, power):
        power = _integer(power, "denominator power")
        numerator = _clean(numerator, self.variables)
        if len(numerator) > self.max_terms:
            raise SectorBudgetExceeded("numerator exceeded max_terms")
        return CubeForm(MappingProxyType(numerator), power)

    def monomial(self, powers, power=0, coefficient=1):
        return self.form({tuple(powers): coefficient}, power)

    def _power(self, power):
        if power not in self._powers:
            half = self._power(power // 2)
            value = _multiply(half, half, self.max_terms)
            if power % 2:
                value = _multiply(value, self.denominator, self.max_terms)
            self._powers[power] = value
        return self._powers[power]

    def lift(self, form, power):
        form = self.form(form.numerator, form.power)
        power = _integer(power, "common denominator power", form.power)
        return _multiply(form.numerator, self._power(power - form.power), self.max_terms)

    def boundary_derivative(self, powers, power, coordinate):
        """Derivative of x_i(1-x_i)x**powers/U**power, with all terms."""
        powers = tuple(_integer(a, "primitive exponent") for a in powers)
        if len(powers) != self.variables:
            raise ValueError("one primitive exponent is required per coordinate")
        power = _integer(power, "primitive denominator power")
        coordinate = _integer(coordinate, "coordinate")
        if coordinate >= self.variables:
            raise ValueError("coordinate outside cube")
        numerator = {}
        for b, coefficient in self.denominator.items():
            a = tuple(x + y for x, y in zip(powers, b))
            # Derivative of the polynomial primitive, and derivative of U**-q.
            _add(numerator, {a: coefficient * (powers[coordinate] + 1 - power * b[coordinate])},
                 budget=self.max_terms)
            _add(numerator, {_shift(a, coordinate): coefficient *
                             (power * b[coordinate] - powers[coordinate] - 2)},
                 budget=self.max_terms)
        return self.form(numerator, power + 1)

    def verify(self, certificate):
        """Check the complete rational-function identity, not reported metadata."""
        if dict(certificate.denominator) != dict(self.denominator):
            return False
        if len(certificate.masters) != len(certificate.coefficients):
            return False
        target = self.form(certificate.target.numerator, certificate.target.power)
        masters = tuple(self.form(x.numerator, x.power) for x in certificate.masters)
        derivatives = [(self.boundary_derivative(a, q, i), _rational(c))
                       for c, i, a, q in certificate.primitives]
        common = max([target.power] + [x.power for x in masters]
                     + [x.power for x, _ in derivatives])
        residual = self.lift(target, common)
        for coefficient, form in zip(certificate.coefficients, masters):
            _add(residual, self.lift(form, common), -_rational(coefficient), self.max_terms)
        for form, coefficient in derivatives:
            _add(residual, self.lift(form, common), -coefficient, self.max_terms)
        return not residual

    def read_certificate(self, record):
        if record.get("schema") != "htk.cube-relative-reduction.v1":
            raise ValueError("unknown reduction certificate schema")
        def form(row):
            return self.form(_read_poly(row["numerator"]), row["power"])
        result = CubeReduction(
            MappingProxyType(_clean(_read_poly(record["denominator"]), self.variables)),
            form(record["target"]), tuple(form(row) for row in record["masters"]),
            tuple(_rational(c) for c in record["coefficients"]),
            tuple((_rational(c), _integer(i, "coordinate"), tuple(a),
                   _integer(q, "primitive power")) for c, i, a, q in record["primitives"]),
            MappingProxyType(dict(record.get("statistics", {}))))
        if not self.verify(result):
            raise ArithmeticError("reduction certificate has a nonzero exact residual")
        return result

    def search(self, target, masters, *, max_primitive_power,
               max_primitive_degree, max_generators=20000,
               max_elimination_entries=1000000):
        """Find a fully certified reduction inside explicit finite bounds.

        This is not a master-basis discovery or a stabilization test.  The
        supplied masters need not be independent.  Exhausted search returns
        no answer; ReductionNotFound retains the entire rational residual.
        """
        target = self.form(target.numerator, target.power)
        masters = tuple(self.form(x.numerator, x.power) for x in masters)
        max_primitive_power = _integer(max_primitive_power, "max_primitive_power")
        max_primitive_degree = _integer(max_primitive_degree, "max_primitive_degree")
        max_generators = _integer(max_generators, "max_generators", 1)
        max_elimination_entries = _integer(max_elimination_entries, "max_elimination_entries", 1)
        count = (max_primitive_power + 1) * self.variables * comb(
            self.variables + max_primitive_degree, self.variables) + len(masters)
        if count > max_generators:
            raise SectorBudgetExceeded("finite relation family exceeded max_generators")
        common = max([target.power, max_primitive_power + 1] + [x.power for x in masters])
        basis, labels = {}, []
        stored = 0

        def eliminate(row, representation):
            if len(row) + len(representation) > max_elimination_entries:
                raise SectorBudgetExceeded("exact elimination input exceeded max_elimination_entries")
            for pivot in sorted(basis):
                coefficient = row.get(pivot, Fraction())
                if coefficient:
                    known, witness = basis[pivot]
                    _add(row, known, -coefficient, self.max_terms)
                    _add(representation, witness, -coefficient)
                    if len(row) + len(representation) > max_elimination_entries:
                        raise SectorBudgetExceeded("exact elimination exceeded max_elimination_entries")

        def add_generator(form, label):
            nonlocal stored
            index = len(labels)
            labels.append(label)
            row, representation = self.lift(form, common), {index: Fraction(1)}
            eliminate(row, representation)
            if not row:
                return
            pivot = min(row)
            scale = row[pivot]
            row = {a: c / scale for a, c in row.items()}
            representation = {a: c / scale for a, c in representation.items()}
            stored += len(row) + len(representation)
            if stored > max_elimination_entries:
                raise SectorBudgetExceeded("stored exact witnesses exceeded max_elimination_entries")
            basis[pivot] = (row, representation)

        for index, form in enumerate(masters):
            add_generator(form, ("master", index))
        for q in range(max_primitive_power + 1):
            for powers in _monomials(self.variables, max_primitive_degree):
                for coordinate in range(self.variables):
                    add_generator(self.boundary_derivative(powers, q, coordinate),
                                  ("primitive", coordinate, powers, q))
        row, representation = self.lift(target, common), {}
        eliminate(row, representation)
        statistics = {"generators": len(labels), "rank": len(basis),
                      "common_denominator_power": common,
                      "max_primitive_power": max_primitive_power,
                      "max_primitive_degree": max_primitive_degree,
                      "stored_elimination_entries": stored}
        if row:
            raise ReductionNotFound(self.form(row, common), statistics)
        coefficients = [Fraction()] * len(masters)
        primitives = []
        for index, coefficient in sorted(representation.items()):
            label = labels[index]
            coefficient = -coefficient
            if label[0] == "master":
                coefficients[label[1]] += coefficient
            else:
                primitives.append((coefficient, label[1], label[2], label[3]))
        answer = CubeReduction(self.denominator, target, masters, tuple(coefficients),
                               tuple(primitives), MappingProxyType(statistics))
        if not self.verify(answer):
            raise ArithmeticError("internal exact reduction witness failed verification")
        return answer


def _cycle_count(graph, selected):
    parent = list(range(graph.nv))
    def root(v):
        while parent[v] != v:
            parent[v] = parent[parent[v]]
            v = parent[v]
        return v
    loops = 0
    for edge in selected:
        a, b = map(root, graph.edges[edge])
        if a == b:
            loops += 1
        else:
            parent[a] = b
    return loops


def check_convergence(graph, alpha, power, *, max_subsets=1048576):
    """Check every strict edge-subset inequality and projective homogeneity."""
    if not isinstance(graph, Graph) or not graph.is_connected() or graph.ne < 1:
        raise ValueError("need a connected htk.Graph with at least one edge")
    alpha = tuple(_integer(a, "edge power", 1) for a in alpha)
    power = _integer(power, "denominator power", 1)
    max_subsets = _integer(max_subsets, "max_subsets", 1)
    if len(alpha) != graph.ne or sum(alpha) != graph.loops * power:
        raise ValueError("scalar period has incorrect projective homogeneity")
    total = (1 << graph.ne) - 2
    if total > max_subsets:
        raise SectorBudgetExceeded("convergence check exceeded max_subsets")
    minimum = None
    for mask in range(1, (1 << graph.ne) - 1):
        selected = tuple(i for i in range(graph.ne) if mask & (1 << i))
        margin = sum(alpha[i] for i in selected) - power * _cycle_count(graph, selected)
        if margin <= 0:
            raise ValueError(f"nonpositive convergence margin on edge subset {selected}")
        minimum = margin if minimum is None else min(minimum, margin)
    return {"proper_subsets": total, "minimum_integer_margin": minimum}


class CoreSector:
    """One complete ordered Hepp sector, retaining the original edge labels."""

    def __init__(self, graph, order, *, max_tree_candidates=1000000,
                 max_terms=100000):
        if not isinstance(graph, Graph) or not graph.is_connected() or graph.ne < 2:
            raise ValueError("need a connected htk.Graph with at least two edges")
        self.graph = graph
        self.order = tuple(_integer(i, "edge index") for i in order)
        if sorted(self.order) != list(range(graph.ne)):
            raise ValueError("sector order must contain every edge exactly once")
        max_tree_candidates = _integer(max_tree_candidates, "max_tree_candidates", 1)
        max_terms = _integer(max_terms, "max_terms", 1)
        if comb(graph.ne, graph.nv - 1) > max_tree_candidates:
            raise SectorBudgetExceeded("spanning-tree enumeration exceeded max_tree_candidates")
        self.prefix_loops = tuple(_cycle_count(graph, self.order[:j])
                                  for j in range(1, graph.ne))
        cotrees = []
        polynomial = {}
        for tree in combinations(range(graph.ne), graph.nv - 1):
            if _cycle_count(graph, tree):
                continue
            cotree = tuple(e for e in range(graph.ne) if e not in tree)
            cotrees.append(cotree)
            powers = tuple(sum(e in cotree for e in self.order[:j]) - self.prefix_loops[j - 1]
                           for j in range(1, graph.ne))
            if any(a < 0 for a in powers):
                raise ArithmeticError("cotree failed the exact prefix-rank bound")
            polynomial[powers] = polynomial.get(powers, 0) + 1
            if len(polynomial) > max_terms:
                raise SectorBudgetExceeded("sector denominator exceeded max_terms")
        self.cotrees = tuple(cotrees)
        self.relations = CubePeriodRelations(polynomial, max_terms=max_terms)
        if polynomial.get((0,) * (graph.ne - 1)) != 1:
            raise ArithmeticError("ordered sector has no unique greedy cotree")

    def pullback(self, alpha, power, *, check_all_subsets=True, max_subsets=1048576):
        """Transform a complete scalar monomial, with the exact Jacobian.

        The default checks global convergence.  Disabling it checks only
        this ordered sector and makes no assertion about other sectors.
        """
        alpha = tuple(_integer(a, "edge power", 1) for a in alpha)
        power = _integer(power, "denominator power", 1)
        if len(alpha) != self.graph.ne or sum(alpha) != self.graph.loops * power:
            raise ValueError("scalar period has incorrect projective homogeneity")
        if check_all_subsets:
            check_convergence(self.graph, alpha, power, max_subsets=max_subsets)
        exponents = tuple(sum(alpha[e] for e in self.order[:j])
                          - power * self.prefix_loops[j - 1] - 1
                          for j in range(1, self.graph.ne))
        if any(a < 0 for a in exponents):
            raise ValueError("the scalar monomial diverges in this sector")
        return self.relations.monomial(exponents, power)

    def shifted_weights(self, alpha, power):
        """All convergent q+1 summands of the exact U-multiplication identity."""
        alpha = tuple(_integer(a, "edge power", 1) for a in alpha)
        check_convergence(self.graph, alpha, power)
        return tuple(tuple(a + (e in cotree) for e, a in enumerate(alpha))
                     for cotree in self.cotrees)
