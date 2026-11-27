"""Exact sourced three-terminal Markov response, followed by full integration.

This implementation does not import BulkPolynomial, its cycle inverse,
the pointwise locality implementation, or the complete source answers.
Only the reduced grounded Laplacian is inverted.  The Dirichlet Gaussian
at the deleted vertex is retained before every polynomial division.
"""

from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations, permutations
from math import prod
from pathlib import Path
import sys

from flint import fmpq, fmpq_mpoly_ctx

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "code"))
from htk.graph import Graph


EDGES = {
    "prism": ((0, 1), (0, 2), (1, 2), (3, 4), (3, 5),
              (4, 5), (0, 3), (1, 4), (2, 5)),
    "K33": ((0, 3), (0, 4), (0, 5), (1, 3), (1, 4),
            (1, 5), (2, 3), (2, 4), (2, 5)),
}


def rational(value):
    if hasattr(value, "p"):
        return Fraction(int(value.p), int(value.q))
    return Fraction(value)


def scalar(value):
    value = rational(value)
    return fmpq(value.numerator, value.denominator)


def parity(order):
    return (-1) ** sum(a > b for i, a in enumerate(order) for b in order[i + 1:])


def integer_determinant(matrix):
    return sum(parity(p) * prod(matrix[i][j] for i, j in enumerate(p))
               for p in permutations(range(len(matrix))))


def evaluate(polynomial, values):
    return sum((rational(c) * prod(x**int(k) for x, k in zip(values, exponent))
                for exponent, c in polynomial.to_dict().items()), Fraction(0))


@dataclass
class ResponseDensity:
    numerator: object
    power: Fraction
    currents: tuple
    source: tuple
    exact_U_divisions: int


class StarResponse:
    """Delete one degree-three vertex with its full source and Gaussian law.

    The supported examples have six vertices, nine edges, and H=T=1.
    All five independent sources and all original edge currents survive.
    Polynomial budgets are checked on complete expressions, never by
    truncating their support.
    """

    def __init__(self, name="prism", *, centre=0, max_terms=500_000,
                 max_resident_terms=3_000_000):
        if name not in EDGES:
            raise ValueError("unknown exact four-loop response example")
        self.name = name
        self.graph = Graph(6, EDGES[name])
        self.edges = tuple(self.graph.edges)
        self.centre = int(centre)
        self.root = self.graph.root
        if self.centre == self.root or self.centre not in range(6):
            raise ValueError("the deleted vertex must be a non-root vertex")
        self.star = tuple(e for e, edge in enumerate(self.edges) if self.centre in edge)
        if len(self.star) != 3:
            raise ValueError("the subnetwork must be one three-leg star")
        self.max_terms = int(max_terms)
        self.max_resident_terms = int(max_resident_terms)
        if min(self.max_terms, self.max_resident_terms) < 1:
            raise ValueError("positive exact-polynomial budgets are required")
        self.context = fmpq_mpoly_ctx.get(tuple("x%d" % e for e in range(9)), ordering="lex")
        self.x = self.context.gens()
        self.zero = self.context.constant(0)
        self.one = self.context.constant(1)
        self.peak_terms = 0
        self.divisions = 0
        self.cycles = tuple(tuple(row) for row in self.graph.cycle_matrix())
        self._response()
        self._prepare_exterior_coefficients()

    def bounded(self, polynomial):
        self.peak_terms = max(self.peak_terms, len(polynomial))
        if len(polynomial) > self.max_terms:
            raise RuntimeError("complete response polynomial exceeds max_terms")
        return polynomial

    def resident(self, polynomials):
        if sum(len(p) for p in polynomials) > self.max_resident_terms:
            raise RuntimeError("complete response buffer exceeds max_resident_terms")

    def add(self, polynomials):
        result = self.zero
        for polynomial in polynomials:
            result = self.bounded(result + polynomial)
        return result

    def multiply(self, polynomials):
        result = self.one
        for polynomial in polynomials:
            result = self.bounded(result * polynomial)
        return result

    def divide(self, numerator, denominator):
        quotient, remainder = divmod(numerator, denominator)
        if remainder:
            raise ArithmeticError("a stated response-polynomial division is not exact")
        self.divisions += 1
        return self.bounded(quotient)

    def determinant(self, matrix):
        size = len(matrix)
        return self.add(parity(p) * self.multiply(matrix[i][j] for i, j in enumerate(p))
                        for p in permutations(range(size)))

    def reciprocal(self, polynomial):
        """Multiply p(1/t) by the product of all nine original times."""
        result = {}
        for exponent, coefficient in polynomial.to_dict().items():
            if any(k not in (0, 1) for k in exponent):
                raise ArithmeticError("the collected response is not multiaffine")
            result[tuple(1-int(k) for k in exponent)] = coefficient
        return self.bounded(self.context.from_dict(result))

    def _response(self):
        # x_e is first a conductance.  Reusing its polynomial context after
        # reciprocal() merely renames the variables as the original times.
        outer = tuple(v for v in range(6) if v not in (self.centre, self.root))
        self.outer = outer
        index = {v: i for i, v in enumerate(outer)}
        g = {v: self.zero for v in range(6)}
        for e in self.star:
            a, b = self.edges[e]
            other = b if a == self.centre else a
            g[other] = g[other] + self.x[e]
        G = self.add(self.x[e] for e in self.star)
        self.G = G
        outside = [[self.zero for _ in outer] for _ in outer]
        for e, (a, b) in enumerate(self.edges):
            if e in self.star:
                continue
            for u, v in ((a, b), (b, a)):
                if u in index:
                    outside[index[u]][index[u]] += self.x[e]
                    if v in index:
                        outside[index[u]][index[v]] -= self.x[e]
        # K = G Lambda_red.  No inverse of the original Laplacian or of
        # the original cycle matrix is formed anywhere in this module.
        matrix = [[self.bounded(G * outside[i][j]
                               + (G*g[u] if i == j else self.zero) - g[u]*g[v])
                   for j, v in enumerate(outer)] for i, u in enumerate(outer)]
        determinant = self.determinant(matrix)
        tree_polynomial = self.divide(determinant, G**3)
        inverse_numerator = [[self.divide((-1)**(i+j) * self.determinant(
            [[matrix[r][s] for s in range(4) if s != i]
             for r in range(4) if r != j]), G**2)
            for j in range(4)] for i in range(4)]
        self.resident(p for row in inverse_numerator for p in row)
        self.tree_polynomial = tree_polynomial
        self.U = self.reciprocal(tree_polynomial)
        self.reduced_inverse_numerator = inverse_numerator
        # Phi_outer = Lambda_red^-1 (lambda_outer + h lambda_centre).
        # Potential numerators below all have denominator K_G.
        potentials = [[self.zero for _ in range(5)] for _ in range(6)]
        for i, u in enumerate(outer):
            for source in range(5):
                if source == self.centre:
                    potentials[u][source] = self.divide(self.add(
                        inverse_numerator[i][j]*g[v] for j, v in enumerate(outer)), G)
                else:
                    potentials[u][source] = inverse_numerator[i][index[source]]
        for source in range(5):
            potentials[self.centre][source] = self.divide(self.add(
                g[v]*potentials[v][source] for v in outer)
                + (tree_polynomial if source == self.centre else self.zero), G)
        current_g = [[self.bounded(self.x[e]*(potentials[a][source]-potentials[b][source]))
                      for source in range(5)] for e, (a, b) in enumerate(self.edges)]
        self.current_coefficients = tuple(tuple(self.reciprocal(p) for p in row)
                                          for row in current_g)
        # A/G maps the reduced Gaussian field to all original currents;
        # b is the coefficient of the one independent Dirichlet field.
        A, b = [], []
        for e, (a, c) in enumerate(self.edges):
            if e in self.star:
                orientation = 1 if a == self.centre else -1
                other = c if orientation == 1 else a
                A.append([self.bounded(orientation*self.x[e]*(g[v]-(G if v == other else self.zero)))
                          for v in outer])
                b.append(orientation*self.x[e])
            else:
                A.append([G*self.x[e]*(int(a == v)-int(c == v)) for v in outer])
                b.append(self.zero)
        AH = [[self.add(A[e][i]*inverse_numerator[i][j] for i in range(4))
               for j in range(4)] for e in range(9)]
        covariance = {}
        contact_nonzero = 0
        for e, f in combinations(range(9), 2):
            harmonic = self.add(AH[e][i]*A[f][i] for i in range(4))
            contact = self.bounded(G*tree_polynomial*b[e]*b[f])
            contact_nonzero += bool(contact)
            combined = self.divide(self.bounded(harmonic + contact), G**2)
            # Independent source-response check of the combined Markov
            # covariance; no cycle-space covariance is used for this check.
            a, c = self.edges[e]
            u, v = self.edges[f]
            def potential_difference(source):
                return (potentials[a][source]-potentials[c][source]
                        if source != self.root else self.zero)
            expected = self.x[e]*self.x[f]*(potential_difference(u)-potential_difference(v))
            if combined != expected:
                raise ArithmeticError("reduced and Dirichlet covariance do not reconstruct the response")
            covariance[e, f] = self.reciprocal(combined)
        self.covariance_numerators = covariance  # twice covariance = R_ef/U
        self.contact_pairs = contact_nonzero
        self.resident([self.U] + [p for row in self.current_coefficients for p in row]
                      + list(covariance.values()))
        self.check_current_identities()
        self.check_tree_polynomial()

    def check_current_identities(self):
        for source in range(5):
            for v in range(6):
                divergence = self.add((int(a == v)-int(b == v))*self.current_coefficients[e][source]
                                      for e, (a, b) in enumerate(self.edges))
                expected = self.U*(int(v == source)-int(v == self.root))
                if divergence != expected:
                    raise ArithmeticError("an original source is missing from the response")
            for cycle in self.cycles:
                if self.add(cycle[e]*self.x[e]*self.current_coefficients[e][source]
                            for e in range(9)):
                    raise ArithmeticError("the recovered currents violate an original cycle")
        return True

    def check_tree_polynomial(self):
        result = self.zero
        for chosen in combinations(range(9), 5):
            parent = list(range(6))
            def find(v):
                while parent[v] != v:
                    v = parent[v]
                return v
            valid = True
            for e in chosen:
                a, b = map(find, self.edges[e])
                if a == b:
                    valid = False
                    break
                parent[a] = b
            if valid:
                result += self.multiply(self.x[e] for e in range(9) if e not in chosen)
        if result != self.U:
            raise ArithmeticError("the response determinant has the wrong original tree normalization")
        return True

    def _prepare_exterior_coefficients(self):
        self.exterior = []
        for top in combinations(range(8), 4):
            hol = tuple(e for e in range(8) if e not in top)
            top_minor = integer_determinant([[row[e] for e in top] for row in self.cycles])
            hol_minor = integer_determinant([[row[e] for e in hol] for row in self.cycles])
            if not top_minor or not hol_minor:
                continue
            a, b, c, d = top
            covariance = self.covariance_numerators
            pairs = ((a, b, c, d), (a, c, b, d), (a, d, b, c))
            wick = self.add(covariance[i, j]*covariance[k, l] for i, j, k, l in pairs)
            coefficient = scalar(Fraction(parity(top + hol)*top_minor*hol_minor, 4))
            self.exterior.append((hol, self.bounded(coefficient*wick), top, coefficient))
        self.resident([p for _, p, _, _ in self.exterior]
                      + [p for row in self.current_coefficients for p in row]
                      + list(self.covariance_numerators.values()))
        self.resident(self.stored_polynomials())

    def stored_polynomials(self):
        return ([self.U, self.G, self.tree_polynomial, self.zero, self.one] + list(self.x)
                + [p for row in self.reduced_inverse_numerator for p in row]
                + [p for row in self.current_coefficients for p in row]
                + list(self.covariance_numerators.values())
                + [p for _, p, _, _ in self.exterior])

    def density(self, source):
        source = tuple(rational(x) for x in source)
        if len(source) == 6:
            if sum(source):
                raise ValueError("the original six sources must sum to zero")
            source = source[:5]
        if len(source) != 5:
            raise ValueError("all five independent original sources are required")
        currents = tuple(self.add(scalar(source[v])*row[v] for v in range(5))
                         for row in self.current_coefficients)
        retained = self.stored_polynomials()+list(currents)
        raw = self.zero
        for hol, wick, _, _ in self.exterior:
            product = self.one
            for edge in hol:
                following = self.bounded(product*currents[edge])
                self.resident(retained+[raw, product, following])
                product = following
            term = self.bounded(wick*product)
            following = self.bounded(raw+term)
            self.resident(retained+[raw, product, term, following])
            raw = following
            del product, term, following
        raw = self.divide(raw, self.x[-1])
        power, divisions = Fraction(15, 2), 0
        while raw:
            quotient, remainder = divmod(raw, self.U)
            self.resident(retained+[raw, quotient, remainder])
            if remainder:
                break
            raw = self.bounded(quotient)
            divisions += 1
            power -= 1
        if raw and any(sum(map(int, exponent)) != 4*power-9 for exponent in raw.to_dict()):
            raise ArithmeticError("the response density has incorrect projective degree")
        return ResponseDensity(raw, power, currents, source, divisions)

    def contact_control(self, times, source):
        """Exact coefficient of U^(-1/2), with and without every contact.

        Only this control specializes times.  The density construction
        above keeps every time and every source coefficient exactly.
        """
        times, source = tuple(map(Fraction, times)), tuple(map(Fraction, source))
        if len(times) != 9 or min(times) <= 0 or len(source) != 5:
            raise ValueError("positive times and all five original sources are required")
        u = evaluate(self.U, times)
        currents = [sum((source[v]*evaluate(row[v], times) for v in range(5)), Fraction(0))/u
                    for row in self.current_coefficients]
        conductance_sum = sum((1/times[e] for e in self.star), Fraction(0))
        b = [Fraction(0)]*9
        for e in self.star:
            b[e] = (1 if self.edges[e][0] == self.centre else -1)/times[e]
        full, harmonic, once, twice = (Fraction(0) for _ in range(4))
        for hol, _, top, coefficient in self.exterior:
            a, c, d, e = top
            pairs = ((a, c, d, e), (a, d, c, e), (a, e, c, d))
            prefactor = rational(coefficient)*prod(currents[j] for j in hol)/(times[-1]*u)
            for i, j, k, l in pairs:
                q1 = evaluate(self.covariance_numerators[i, j], times)/u
                q2 = evaluate(self.covariance_numerators[k, l], times)/u
                c1, c2 = b[i]*b[j]/conductance_sum, b[k]*b[l]/conductance_sum
                full += prefactor*q1*q2
                harmonic += prefactor*(q1-c1)*(q2-c2)
                once += prefactor*((q1-c1)*c2+c1*(q2-c2))
                twice += prefactor*c1*c2
        if twice or full != harmonic+once:
            raise ArithmeticError("the three-terminal single-contact identity failed")
        return {"U": u, "complete_coefficient": full,
                "contact_omitted_coefficient": harmonic,
                "single_contact_coefficient": once,
                "double_contact_coefficient": twice,
                "original_currents": currents}

