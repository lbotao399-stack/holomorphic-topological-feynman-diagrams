"""Exact current moments of cubic marked-product forests.

Each rooted subtree is contracted on demand to a polynomial in three
boundary currents: root->A, root->B, and A->B.  Its six scalar coordinates
do not grow with the graph.  A separate integer records the total degree
of an outside marked-current word.  This degree shifts when two child
messages are composed; ignoring that shift would lose marked-edge
correlations.

Only the two or four ratios of the current local move are introduced.
They are integrated immediately by exact beta/Dirichlet moments.  Root
sources are then transferred by a finite, exact polynomial translation in
the two A-source coordinates.  No global ratio chart or precomputed tree
table is used.  This implements the source-free single/paired star and
sourced root-series identities in sources/LLC_main.tex.

The optional python-flint dependency supplies exact sparse rational
polynomials.  There is no fixed loop or current-degree cutoff.  Explicit
term and state budgets raise ForestBudgetExceeded without returning a
partial answer.
"""

from fractions import Fraction
from functools import cache
from math import factorial
from numbers import Integral

from .forest_model import CubicForest

try:
    from flint import fmpq, fmpq_mpoly_ctx
except ImportError as error:  # Keep unrelated htk functionality importable.
    _FLINT_ERROR = error
else:
    _FLINT_ERROR = None

__all__ = ["ForestMoments", "forest_moment", "ForestBudgetExceeded"]


class ForestBudgetExceeded(RuntimeError):
    """An explicit exact-computation budget was exceeded."""


def _integer(value, name, minimum=0):
    if isinstance(value, bool) or not isinstance(value, Integral) or value < minimum:
        raise ValueError("%s must be an integer at least %d" % (name, minimum))
    return int(value)


def _q(value):
    value = Fraction(value)
    return fmpq(value.numerator, value.denominator)


def _fraction(value):
    return Fraction(int(value.p), int(value.q))


def _add(left, right):
    return tuple(a + b for a, b in zip(left, right))


def _sub(left, right):
    return tuple(a - b for a, b in zip(left, right))


def _scale(vector, scalar):
    return tuple(value * scalar for value in vector)


def _bracket(left, right):
    return left[0] * right[1] - left[1] * right[0]


def _parity(sequence):
    return -1 if sum(a > b for i, a in enumerate(sequence)
                     for b in sequence[i + 1:]) % 2 else 1


def _edge(left, right):
    return (left, right) if left < right else (right, left)


@cache
def _factorial(value):
    return factorial(value)


@cache
def _paired_weight(first, second, third, fourth, marked_degree):
    # Integral r0^first r1^second (r0-r1) over the square, followed
    # by r2^third r3^fourth (1-r2-r3)^marked_degree over the triangle.
    if first == second:
        return fmpq(0)
    interval = fmpq(first - second,
                    (first + 1) * (first + 2) * (second + 1) * (second + 2))
    triangle = fmpq(_factorial(third) * _factorial(fourth) * _factorial(marked_degree),
                    _factorial(third + fourth + marked_degree + 2))
    return interval * triangle


def _plan_and_phase(model):
    """Use the same ordered local forms as the original ratio evaluator."""
    plan = []

    def descend(vertex):
        ports = model.canonical_nodes[vertex].ports
        children = [port for port in ports if not isinstance(port, str)]
        for child in children:
            descend(child)
        if len(children) == 1:
            boundary = next(port for port in ports if isinstance(port, str))
            near = model.A if boundary == "A" else model.B
            far = model.B if boundary == "A" else model.A
            plan.append(("single", (children[0], vertex, far, near)))
        elif len(children) == 2:
            plan.append(("paired", (children[0], children[1], vertex, model.A, model.B)))

    for root in model.canonical_roots:
        descend(root)
        plan.append(("series", (root, model.A, model.B)))
    edges = list(model.canonical_graph.edges)
    phase = (-1) ** (model.loops * (model.loops + 1) // 2)
    for kind, vertices in plan:
        if kind == "series":
            vertex, left, right = vertices
            before, after = ((vertex, left), (vertex, right)), ()
        elif kind == "single":
            vertex, parent, far, near = vertices
            before = ((vertex, parent), (vertex, far), (vertex, near),
                      (parent, near), (far, near))
            after = ((parent, far), (parent, near), (far, near))
        else:
            left, right, parent, A, B = vertices
            before = ((left, parent), (left, A), (left, B),
                      (right, parent), (right, A), (right, B), (A, B))
            after = ((parent, A), (parent, B), (A, B))
        before = [_edge(*item) for item in before]
        after = [_edge(*item) for item in after]
        rest = [item for item in edges if item not in before]
        phase *= _parity([edges.index(item) for item in before + rest])
        phase *= _parity(after + rest)
        edges = sorted(after + rest)
    if edges != [(model.A, model.B)]:
        raise ArithmeticError("forest reduction did not end at the marked edge")
    return tuple(plan), phase


class ForestMoments:
    """Source-complete finite current moments on a cubic marked forest.

    ``G`` is an htk Graph; ``marked=(A,B)`` orders the insertion vertices.
    ``lam`` follows htk.brute.sources: non-root rows, full balanced rows, or
    a vertex dictionary.  Every interaction source except each component's
    degree-two root must vanish, as required by the cubic product family.

    ``moment(powers)`` integrates the current monomial specified by
    ``{(original_edge_index, component): nonnegative_power}``.  Component
    indices are zero based.  The result uses the original htk graph order
    and orientations.  ``source_polynomial`` returns its exact polynomial
    in the A-source coordinates with the component-root vectors held fixed.
    """

    def __init__(self, G, lam, marked, *, max_terms=250000, max_states=100000):
        if _FLINT_ERROR is not None:
            raise ImportError("forest moments require the optional python-flint package") from _FLINT_ERROR
        self.model = CubicForest(G, marked=marked)
        self.sources = self.model.normalize_sources(lam)
        self.max_terms = _integer(max_terms, "max_terms", 1)
        self.max_states = _integer(max_states, "max_states", 1)
        self.plan, self.phase = _plan_and_phase(self.model)
        self.last_stats = {}
        self._results = {}
        self._message_cache = {}

    def _powers(self, powers):
        normalized = {}
        for key, power in ({} if powers is None else dict(powers)).items():
            if not isinstance(key, tuple) or len(key) != 2:
                raise ValueError("powers keys must be (edge index, component)")
            edge = _integer(key[0], "edge index")
            component = _integer(key[1], "component")
            power = _integer(power, "current power")
            if edge >= self.model.G.ne or component >= 2:
                raise ValueError("current edge/component index is out of range")
            if power:
                normalized[edge, component] = power
        return tuple(sorted(normalized.items()))

    def source_polynomial(self, powers=None):
        """Return {(A-source exponent pair): Fraction} in original conventions."""
        key = self._powers(powers)
        if key not in self._results:
            request = _Request(self, dict(key))
            polynomial = request.run()
            result = {tuple(map(int, exponent)): _fraction(value)
                      for exponent, value in polynomial.terms() if value}
            self._results[key] = result, dict(request.stats)
        result, stats = self._results[key]
        self.last_stats = dict(stats)
        return dict(result)

    def moment(self, powers=None):
        """Return the exact original-edge current moment as Fraction."""
        polynomial = self.source_polynomial(powers)
        source = self.sources[self.model.A]
        return sum((coefficient * source[0] ** exponents[0] * source[1] ** exponents[1]
                    for exponents, coefficient in polynomial.items()), Fraction(0))

    def moments(self, words):
        """Evaluate a finite list of power dictionaries, sharing subtree messages."""
        return tuple(self.moment(powers) for powers in words)


def forest_moment(G, lam, marked, powers=None, **budgets):
    """Convenience wrapper for ``ForestMoments(...).moment(powers)``."""
    return ForestMoments(G, lam, marked, **budgets).moment(powers)


class _Request:
    def __init__(self, engine, powers):
        self.engine, self.model = engine, engine.model
        self.current = fmpq_mpoly_ctx.get(("l0", "l1", "h0", "h1", "y0", "y1"))
        self.single = fmpq_mpoly_ctx.get(("l0", "l1", "h0", "h1", "y0", "y1", "r0", "r1"))
        self.paired = fmpq_mpoly_ctx.get(("l0", "l1", "h0", "h1", "y0", "y1", "r0", "r1", "r2", "r3"))
        self.root = fmpq_mpoly_ctx.get(("s", "r", "j0", "j1"))
        self.terminal = fmpq_mpoly_ctx.get(("s", "j0", "j1"))
        self.source = fmpq_mpoly_ctx.get(("a0", "a1"))
        self.shift = fmpq_mpoly_ctx.get(("a0", "a1", "s"))
        self.parents = {}
        for vertex, node in enumerate(self.model.canonical_nodes):
            for port in node.ports:
                if not isinstance(port, str):
                    self.parents[port] = vertex
        self.edge_powers = {}
        self.sign = engine.phase * self.model.kernel_sign
        for (edge, component), power in powers.items():
            mapped = self.model.original_to_canonical_edges[edge]
            endpoints = self.model.canonical_graph.edges[mapped]
            previous = self.edge_powers.setdefault(tuple(endpoints), [0, 0])
            previous[component] += power
            self.sign *= self.model.original_edge_signs[edge] ** power
        self.marked_powers = tuple(self.edge_powers.get((self.model.A, self.model.B), (0, 0)))
        self.signatures = {}
        # Preorder IDs put descendants after their parents.  A signature
        # retains the insertion pattern but forgets irrelevant vertex labels.
        for vertex in reversed(range(len(self.model.canonical_nodes))):
            node = self.model.canonical_nodes[vertex]
            incoming = (tuple(self.edge_powers.get(_edge(self.parents[vertex], vertex), (0, 0)))
                        if vertex in self.parents else (0, 0))
            outgoing = []
            for port in node.ports:
                if isinstance(port, str):
                    boundary = self.model.A if port == "A" else self.model.B
                    outgoing.append((port, tuple(self.edge_powers.get((vertex, boundary), (0, 0)))))
                else:
                    outgoing.append(self.signatures[port])
            self.signatures[vertex] = (node.encoding, incoming, tuple(outgoing))
        self.messages = {}
        self.grades = {}
        self.terminals = {}
        self.suffixes = {}
        self.translations = {}
        self.shifted = {}
        self.stats = dict(loops=self.model.loops, current_degree=sum(powers.values()),
                          tree_message_states=0, terminal_message_states=0,
                          suffix_moment_states=0, translation_states=0,
                          tree_message_cache_hits=0, source_shift_expansions=0,
                          largest_polynomial_terms=0, current_interface_coordinates=6,
                          largest_local_ratio_count=0, global_ratio_coordinates=0)

    def check(self, polynomial):
        terms = len(polynomial)
        self.stats["largest_polynomial_terms"] = max(self.stats["largest_polynomial_terms"], terms)
        if terms > self.engine.max_terms:
            raise ForestBudgetExceeded("local polynomial has %d terms, exceeding max_terms=%d"
                                       % (terms, self.engine.max_terms))
        return polynomial

    def state(self, name):
        self.stats[name] += 1
        total = sum(self.stats[key] for key in ("tree_message_states", "terminal_message_states",
                                                "suffix_moment_states", "translation_states"))
        if total > self.engine.max_states:
            raise ForestBudgetExceeded("forest contraction exceeds max_states=%d" % self.engine.max_states)

    def word(self, endpoints, vector, context):
        result = context.constant(1)
        for component, power in enumerate(self.edge_powers.get(_edge(*endpoints), (0, 0))):
            if power:
                result *= vector[component] ** power
        return self.check(result)

    def parts(self, vertex, degree):
        key = vertex, degree
        if key not in self.grades:
            grouped = {}
            for exponents, coefficient in self.message(vertex, degree).terms():
                marked = int(exponents[4] + exponents[5])
                grouped.setdefault(marked, {})[exponents] = coefficient
            self.grades[key] = {marked: self.current.from_dict(terms)
                                for marked, terms in grouped.items()}
        return self.grades[key]

    def _single_integral(self, polynomial, degree):
        coefficients = {}
        for exponents, coefficient in polynomial.terms():
            factor = fmpq(1, (int(exponents[6]) + 2) * (int(exponents[7]) + degree + 2))
            key = exponents[:6]
            coefficients[key] = coefficients.get(key, 0) + coefficient * factor
        return self.check(self.current.from_dict(coefficients))

    def _paired_integral(self, polynomial, degree):
        coefficients = {}
        for exponents, coefficient in polynomial.terms():
            factor = _paired_weight(*map(int, exponents[6:]), degree)
            if factor:
                key = exponents[:6]
                coefficients[key] = coefficients.get(key, 0) + coefficient * factor
        return self.check(self.current.from_dict(coefficients))

    def message(self, vertex, degree):
        key = vertex, degree
        if key in self.messages:
            return self.messages[key]
        shared_key = self.signatures[vertex], degree
        if shared_key in self.engine._message_cache:
            self.stats["tree_message_cache_hits"] += 1
            self.messages[key] = self.check(self.engine._message_cache[shared_key])
            return self.messages[key]
        self.state("tree_message_states")
        ports = self.model.canonical_nodes[vertex].ports
        children = [port for port in ports if not isinstance(port, str)]
        variables = self.current.gens()
        left, right, marked = variables[:2], variables[2:4], variables[4:6]
        result = self.current.constant(0)
        if not children:
            result = (self.word((vertex, self.model.A), left, self.current)
                      * self.word((vertex, self.model.B), right, self.current))
        elif len(children) == 1:
            self.stats["largest_local_ratio_count"] = max(self.stats["largest_local_ratio_count"], 2)
            direct = next(port for port in ports if isinstance(port, str))
            variables = self.single.gens()
            L, H, Y = variables[:2], variables[2:4], variables[4:6]
            first, second = variables[6:]
            if direct == "A":
                child_left = _sub(_scale(L, 1 - first), _scale(Y, 1 - second))
                child_right = _add(H, _scale(Y, 1 - second))
                direct_current, direct_vertex = _scale(L, first), self.model.A
                bracket = -_bracket(left, marked)
            else:
                child_left = _sub(L, _scale(Y, 1 - second))
                child_right = _add(_scale(H, 1 - first), _scale(Y, 1 - second))
                direct_current, direct_vertex = _scale(H, first), self.model.B
                bracket = _bracket(right, marked)
            original_spoke = self.word((vertex, direct_vertex), direct_current, self.single)
            for child_degree, polynomial in self.parts(children[0], degree).items():
                transformed = self.check(polynomial.compose(*(child_left + child_right + Y)))
                transformed = self.check(transformed * original_spoke)
                result += self._single_integral(transformed, degree + child_degree)
            result = self.check(result * bracket)
        else:
            self.stats["largest_local_ratio_count"] = 4
            variables = self.paired.gens()
            L, H, Y = variables[:2], variables[2:4], variables[4:6]
            first, second, third, fourth = variables[6:]
            left_currents = (_sub(_scale(L, first), _scale(Y, third))
                             + _add(_scale(H, second), _scale(Y, third)) + Y)
            right_currents = (_sub(_scale(L, 1 - first), _scale(Y, fourth))
                              + _add(_scale(H, 1 - second), _scale(Y, fourth)) + Y)
            for left_degree, first_polynomial in self.parts(children[0], degree).items():
                first_polynomial = self.check(first_polynomial.compose(*left_currents))
                for right_degree, second_polynomial in self.parts(children[1], degree + left_degree).items():
                    second_polynomial = self.check(second_polynomial.compose(*right_currents))
                    transformed = self.check(first_polynomial * second_polynomial)
                    result += self._paired_integral(transformed,
                                                    degree + left_degree + right_degree + 1)
            result = self.check(result * _bracket(left, marked) * _bracket(right, marked))
        if vertex in self.parents:
            result = self.check(result * self.word((self.parents[vertex], vertex),
                                                   _add(left, right), self.current))
        result = self.check(result)
        self.messages[key] = result
        self.engine._message_cache[shared_key] = result
        return result

    def terminal_message(self, index, degree):
        key = index, degree
        if key in self.terminals:
            return self.terminals[key]
        self.state("terminal_message_states")
        vertex = self.model.canonical_roots[index]
        source = tuple(_q(value) for value in self.engine.sources[vertex])
        split, parallel, J0, J1 = self.root.gens()
        J = (J0, J1)
        left = _sub(_scale(source, split), _scale(J, 1 - parallel))
        right = _add(_scale(source, 1 - split), _scale(J, 1 - parallel))
        coefficients = {}
        for marked_degree, polynomial in self.parts(vertex, degree).items():
            transformed = self.check(polynomial.compose(*(left + right + J)))
            for exponents, coefficient in transformed.terms():
                weight = fmpq(1, int(exponents[1]) + degree + marked_degree + 2)
                exponent = (int(exponents[0]), int(exponents[2]), int(exponents[3]))
                coefficients[exponent] = coefficients.get(exponent, 0) + coefficient * weight
        variables = self.terminal.gens()
        result = self.check(self.terminal.from_dict(coefficients)
                            * _bracket(source, variables[1:3]))
        self.terminals[key] = result
        return result

    def translated(self, index, powers, split_degree):
        key = index, powers, split_degree
        if key in self.translations:
            return self.translations[key]
        self.state("translation_states")
        shift_key = index, powers
        if shift_key not in self.shifted:
            polynomial = self.suffix(index + 1, powers)
            vertex = self.model.canonical_roots[index]
            source = tuple(_q(value) for value in self.engine.sources[vertex])
            A0, A1, split = self.shift.gens()
            self.shifted[shift_key] = self.check(polynomial.compose(A0 + split * source[0],
                                                                   A1 + split * source[1]))
            self.stats["source_shift_expansions"] += 1
        shifted = self.shifted[shift_key]
        coefficients = {}
        for exponents, coefficient in shifted.terms():
            weight = fmpq(1, int(exponents[2]) + split_degree + 1)
            exponent = exponents[:2]
            coefficients[exponent] = coefficients.get(exponent, 0) + coefficient * weight
        result = self.check(self.source.from_dict(coefficients))
        self.translations[key] = result
        return result

    def suffix(self, index, powers):
        key = index, powers
        if key in self.suffixes:
            return self.suffixes[key]
        self.state("suffix_moment_states")
        if index == len(self.model.canonical_roots):
            result = self.source.from_dict({powers: 1})
        else:
            result = self.source.constant(0)
            for exponents, coefficient in self.terminal_message(index, sum(powers)).terms():
                following = (powers[0] + int(exponents[1]), powers[1] + int(exponents[2]))
                result += coefficient * self.translated(index, following, int(exponents[0]))
            result = self.check(result)
        self.suffixes[key] = result
        return result

    def run(self):
        return self.check(self.suffix(0, self.marked_powers) * self.sign)
