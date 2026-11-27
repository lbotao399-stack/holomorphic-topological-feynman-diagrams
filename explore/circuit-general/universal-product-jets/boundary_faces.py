"""Every oriented time boundary in the K4 source-simplex Stokes identity.

All 62 nonempty proper edge subsets are retained. Forty-six have an
explicit source-independent null direction in the complete current map:
an independent component scale or an inert self-loop scale. The remaining
six single edges, four triangles, and six edge complements are evaluated
as complete six-current exterior forms in polynomial SP charts.

This does not replace a boundary by a scalar Laman kernel. Both affine
source-simplex variables, every shifted inner source, and all differentials
of all three unchanged chord currents are included in the six-by-six
Jacobian. Integration uses exact product-simplex factorial moments.
"""
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import factorial

from flint import fmpq, fmpq_mpoly_ctx

import jet_tensor as jets
from htk.graph import Graph
from htk.spred import SPReduction

EDGES = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
CHORDS = (3, 4, 5)


@lru_cache(maxsize=1)
def spanning_trees():
    result = tuple(tree for tree in combinations(range(6), 3)
                   if Graph(4, [EDGES[edge] for edge in tree]).is_connected())
    if len(result) != 16:
        raise ArithmeticError("the complete K4 spanning-tree census changed")
    return result


def _q(value):
    value = Fraction(value)
    return fmpq(value.numerator, value.denominator)


def _fraction(value):
    return Fraction(int(value.p), int(value.q))


def classify_faces():
    rows = []
    for size in range(1, 6):
        for selected in combinations(range(6), size):
            complement = tuple(edge for edge in range(6) if edge not in selected)
            vertices = sorted({vertex for edge in selected for vertex in EDGES[edge]})
            parent = {vertex: vertex for vertex in vertices}

            def find(vertex):
                while parent[vertex] != vertex:
                    vertex = parent[vertex]
                return vertex

            for edge in selected:
                a, b = EDGES[edge]
                parent[find(b)] = find(a)
            components = [tuple(vertex for vertex in vertices if find(vertex) == root)
                          for root in sorted({find(vertex) for vertex in vertices})]
            row = {"selected_edges": selected, "complement_edges": complement,
                   "edge_components": tuple(tuple(edge for edge in selected if EDGES[edge][0] in component)
                                            for component in components)}
            if len(components) > 1:
                row.update(kind="identically_zero", null_direction="relative scale of a complete selected-edge component",
                           null_edge_set=row["edge_components"][0], full_current_derivative=((0, 0),) * 6)
            else:
                contracted = set(vertices)
                loops = tuple(edge for edge in complement if set(EDGES[edge]) <= contracted)
                active = tuple(edge for edge in complement if edge not in loops)
                if loops and (active or len(loops) > 1):
                    row.update(kind="identically_zero", null_direction="relative scale of an inert quotient self-loop",
                               null_edge_set=(loops[0],), quotient_zero_current_edges=loops,
                               full_current_derivative=((0, 0),) * 6)
                else:
                    row.update(kind="polynomial_SP_chart", selected_vertices=tuple(vertices),
                               quotient_zero_current_edges=loops)
            tree_powers = tuple(len(set(selected).difference(tree)) for tree in spanning_trees())
            minimum = min(tree_powers)
            leading = tuple(index for index, power in enumerate(tree_powers) if power == minimum)
            row["all_leading_spanning_tree_indices"] = leading
            if row["kind"] == "identically_zero":
                null_edges = set(row["null_edge_set"])
                scale_block = set(selected if "component" in row["null_direction"] else complement)
                if not null_edges or not null_edges < scale_block:
                    raise ArithmeticError("the proposed null scale is only the removed projective scale")
                powers = tuple(len(null_edges.difference(spanning_trees()[index])) for index in leading)
                if len(set(powers)) != 1:
                    raise ArithmeticError("the complete leading current numerator is not invariant under the proposed scale")
                row["exact_source_independent_rank_proof"] = {
                    "all_leading_weight_scale_exponents": powers,
                    "common_weight_scale_exponent": powers[0],
                    "full_numerator_and_denominator_identity": "Every leading-tree weight has the same scale power; hence every full source-linear current numerator and its denominator acquire the identical factor.",
                    "complete_current_derivative": "(N' U-N U')/U^2=0 for every original current and every source coefficient",
                    "relative_projective_direction_is_nontrivial": True,
                }
            rows.append(row)
    if len(rows) != 62 or sum(row["kind"] == "polynomial_SP_chart" for row in rows) != 16:
        raise ArithmeticError("the complete K4 scale-face classification changed")
    return tuple(rows)


def _determinant(matrix, context, max_terms):
    if any(not any(row[column] for row in matrix) for column in range(6)):
        return context.constant(0)
    states = {0: context.constant(1)}
    for row in matrix:
        following = {}
        for used, value in states.items():
            for column, entry in enumerate(row):
                if used & (1 << column) or not entry:
                    continue
                sign = -1 if (used >> (column + 1)).bit_count() % 2 else 1
                target = used | (1 << column)
                following[target] = following.get(target, context.constant(0)) + sign * value * entry
                if len(following[target]) > max_terms:
                    raise RuntimeError("the complete boundary Jacobian exceeded max_terms")
        states = following
    return states.get(63, context.constant(0))


def _bounded_power(value, power, context, max_terms):
    result = context.constant(1)
    while power:
        if power & 1:
            result *= value
            if len(result) > max_terms:
                raise RuntimeError("the complete inserted boundary polynomial exceeded max_terms")
        power >>= 1
        if power:
            value *= value
            if len(value) > max_terms:
                raise RuntimeError("the complete inserted boundary polynomial exceeded max_terms")
    return result


def _response(red, context, offset):
    result = []
    for row in red.resp:
        converted = []
        for polynomial in row:
            terms = {}
            for monomial, coefficient in polynomial.items():
                powers = [0] * 6
                for position, power in enumerate(monomial):
                    powers[offset + position] = int(power)
                terms[tuple(powers)] = _q(Fraction(str(coefficient)))
            converted.append(context.from_dict(terms))
        result.append(tuple(converted))
    return tuple(result)


def _currents(response, sources, context):
    return tuple(tuple(sum((coefficient * sources[vertex][component]
                            for vertex, coefficient in enumerate(row)), context.constant(0))
                       for component in range(2)) for row in response)


def _blocks(red, offset):
    return tuple(tuple(offset + red._cidx[name] for name in block.indep) for block in red.blocks if block.indep)


@lru_cache(maxsize=16)
def _topology(selected):
    vertices = sorted({vertex for edge in selected for vertex in EDGES[edge]})
    inner_map = {vertex: i for i, vertex in enumerate(vertices)}
    inner = SPReduction(Graph(len(vertices), [(inner_map[EDGES[edge][0]], inner_map[EDGES[edge][1]])
                                              for edge in selected]))
    complement = tuple(edge for edge in range(6) if edge not in selected)
    groups = sorted([tuple(vertices)] + [(vertex,) for vertex in range(4) if vertex not in vertices], key=lambda group: group[0])
    outer_map = {vertex: position for position, group in enumerate(groups) for vertex in group}
    outer_edges = [(outer_map[EDGES[edge][0]], outer_map[EDGES[edge][1]]) for edge in complement]
    if any(a == b for a, b in outer_edges):
        if len(complement) != 1:
            raise ArithmeticError("an inert quotient scale was not retained in the zero-form classification")
        outer = None
    else:
        outer = SPReduction(Graph(len(groups), outer_edges))
    if len(inner.coords) + (len(outer.coords) if outer else 0) != 4:
        raise ArithmeticError("the complete four-dimensional time face was not retained")
    parity = (-1) ** sum(a > b for i, a in enumerate(selected + complement) for b in (selected + complement)[i + 1:])
    orientation = parity * (-1) ** len(selected) * inner.orientation * (outer.orientation if outer else 1)
    return vertices, groups, complement, inner, outer, orientation


class SourceTriangleTimeFaces:
    """Finite rational time-boundary functional for any complete chord polynomial."""

    def __init__(self, vectors, triangle="TWA", *, max_terms=1000000):
        if triangle not in ("TWA", "TBA"):
            raise ValueError("the source triangle must be TWA or TBA")
        if len(vectors) != 3 or any(len(vector) != 2 for vector in vectors):
            raise ValueError("all three complete original source spinors are required")
        if any(isinstance(value, (float, complex, bool)) for vector in vectors for value in vector):
            raise TypeError("all sources must be exact rational values")
        self.vectors = tuple(tuple(map(Fraction, vector)) for vector in vectors)
        self.triangle = triangle
        self.max_terms = jets.covariants.integer(max_terms, "max_terms", 1)
        self.context = fmpq_mpoly_ctx.get(tuple("ratio%d" % i for i in range(6)))
        zero = self.context.constant(0)
        variables = self.context.gens()
        source_first, source_second = variables[4:]
        P, Q, R = self.vectors
        sources = [[zero, zero], [self.context.constant(_q(value)) for value in P],
                   [self.context.constant(_q(value)) for value in R],
                   [self.context.constant(_q(-P[i] - Q[i] - R[i])) for i in range(2)]]
        triangle_vertices = (0, 1, 2) if triangle == "TWA" else (0, 3, 2)
        for vertex, weight in zip(triangle_vertices, (1 - source_first - source_second, source_first, source_second)):
            for component in range(2):
                sources[vertex][component] += weight * _q(Q[component])
        if any(sum((row[i] for row in sources), zero) for i in range(2)):
            raise ArithmeticError("the complete source simplex is not conserved")
        self.classification = classify_faces()
        self.faces = []
        for classification in self.classification:
            if classification["kind"] == "identically_zero":
                continue
            selected = classification["selected_edges"]
            vertices, groups, complement, inner, outer, orientation = _topology(selected)
            inner_dimension = len(inner.coords)
            if outer is None:
                outer_currents = ((zero, zero),)
            else:
                outer_sources = tuple(tuple(sum((sources[vertex][i] for vertex in group), zero) for i in range(2))
                                      for group in groups)
                outer_currents = _currents(_response(outer, self.context, inner_dimension), outer_sources, self.context)
            inner_sources = [[value for value in sources[vertex]] for vertex in vertices]
            for edge, current in zip(complement, outer_currents):
                a, b = EDGES[edge]
                for component in range(2):
                    if a in vertices:
                        inner_sources[vertices.index(a)][component] -= current[component]
                    if b in vertices:
                        inner_sources[vertices.index(b)][component] += current[component]
            if any(sum((row[i] for row in inner_sources), zero) for i in range(2)):
                raise ArithmeticError("an inner source lost its affine outer-current shift")
            inner_currents = _currents(_response(inner, self.context, 0), inner_sources, self.context)
            currents = [None] * 6
            for edge, current in zip(selected, inner_currents):
                currents[edge] = current
            for edge, current in zip(complement, outer_currents):
                currents[edge] = current
            chord_coordinates = tuple(currents[edge][component] for edge in CHORDS for component in range(2))
            jacobian = _determinant([[current.derivative(index) for index in range(6)]
                                     for current in chord_coordinates], self.context, self.max_terms)
            blocks = _blocks(inner, 0) + (_blocks(outer, inner_dimension) if outer else ()) + ((4, 5),)
            if sorted(index for block in blocks for index in block) != list(range(6)):
                raise ArithmeticError("a boundary integration direction was omitted")
            self.faces.append({"selected_edges": selected, "orientation": orientation,
                               "all_original_currents": tuple(currents), "chord_coordinates": chord_coordinates,
                               "complete_jacobian": jacobian, "simplex_blocks": blocks})
        self.last_stats = {"all_scale_subsets": 62, "source_independent_zero_forms": 46,
                           "complete_polynomial_face_charts": len(self.faces),
                           "nonzero_complete_jacobians": sum(bool(face["complete_jacobian"]) for face in self.faces)}

    def integrate(self, polynomial):
        terms = {}
        for powers, coefficient in polynomial.items():
            if len(powers) != 6:
                raise ValueError("every common-chord monomial has six exponents")
            powers = tuple(jets.covariants.integer(power, "current power") for power in powers)
            if isinstance(coefficient, (float, complex, bool)):
                raise TypeError("current polynomial coefficients must be exact")
            coefficient = Fraction(coefficient)
            if coefficient:
                terms[powers] = terms.get(powers, Fraction(0)) + coefficient
                if len(terms) > self.max_terms:
                    raise RuntimeError("the complete input boundary polynomial exceeded max_terms")
        result = Fraction(0)
        for face in self.faces:
            if not face["complete_jacobian"]:
                continue
            insertion = self.context.constant(0)
            for powers, coefficient in terms.items():
                term = self.context.constant(_q(coefficient))
                for current, power in zip(face["chord_coordinates"], powers):
                    if power:
                        term *= _bounded_power(current, power, self.context, self.max_terms)
                        if len(term) > self.max_terms:
                            raise RuntimeError("the complete inserted boundary polynomial exceeded max_terms")
                insertion += term
                if len(insertion) > self.max_terms:
                    raise RuntimeError("the complete inserted boundary polynomial exceeded max_terms")
            density = face["complete_jacobian"] * insertion
            if len(density) > self.max_terms:
                raise RuntimeError("the complete boundary density exceeded max_terms")
            for powers, coefficient in density.terms():
                weight = Fraction(face["orientation"])
                for block in face["simplex_blocks"]:
                    exponents = [int(powers[index]) for index in block]
                    denominator = factorial(sum(exponents) + len(block))
                    for exponent in exponents:
                        weight *= factorial(exponent)
                    weight /= denominator
                result += _fraction(coefficient) * weight
        return result
