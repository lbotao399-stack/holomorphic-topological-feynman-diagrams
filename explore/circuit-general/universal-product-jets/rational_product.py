"""Complete quartic product with arbitrary simultaneous finite input jets.

The two nonrational graph pairs are evaluated by their COMPLETE rational
source/time boundaries. The six remaining color words use source-free
contractions followed by polynomial SP integration. No individual K4
period or period table is requested.

All source spinors and all ten manuscript t-generator words are retained.
moment(powers) accepts (input_name,lower_spinor_entry)->nonnegative power.
tensor(orders) returns EVERY entry of the tensor separately symmetric in
each input's derivative indices, with no contraction between those indices.
Inputs are beta1,beta2,beta3,b,gamma5; gamma5 is the free neutral spectator.
The physical multiplier is +g^2/[128*(8*pi^2)^3] after source replacement
by its own surviving-field derivative, at every finite total jet order.
For gamma5, a jet means the physical spectator derivative applied AFTER
that source conversion. Its source factor here records that derivative;
it is not a propagator current and carries no propagator-head factor.
"""
from collections import defaultdict
from collections.abc import Mapping
from fractions import Fraction
from itertools import product
from math import prod
from pathlib import Path
import sys

from flint import fmpq_mpoly_ctx

import jet_tensor as jets
from boundary_faces import SourceTriangleTimeFaces, _q, _fraction
from routing import INPUTS, INPUT_NAMES

sys.path.insert(0, str(jets.HERE.parent / "forest-jet-identities"))
from outer_faces import SourceTriangleOuterFaces
from htk.forest_contractions import ForestContractions
from htk.forest_periods import CoreCurrentPeriods
from htk.spred import plan_generic

ALL_INPUTS = INPUT_NAMES + ("gamma5",)


def _powers(powers):
    if powers is None:
        return {}
    if not isinstance(powers, Mapping):
        raise TypeError("input current powers must be a mapping")
    result = {}
    for location, power in powers.items():
        if not isinstance(location, tuple) or len(location) != 2 or location[0] not in ALL_INPUTS:
            raise ValueError("current locations are (beta1|beta2|beta3|b|gamma5, lower-spinor entry)")
        component = jets.covariants.integer(location[1], "lower-spinor entry")
        if component not in (0, 1):
            raise ValueError("the lower spinor has exactly two entries")
        power = jets.covariants.integer(power, "input derivative order")
        if power:
            result[location[0], component] = power
    return result


def _bounded_power(value, power, context, max_terms):
    result = context.constant(1)
    while power:
        if power & 1:
            result *= value
            if len(result) > max_terms:
                raise RuntimeError("the complete transferred input word exceeded max_terms")
        power >>= 1
        if power:
            value *= value
            if len(value) > max_terms:
                raise RuntimeError("the complete transferred input word exceeded max_terms")
    return result


def _canonical_cross_word(record, powers):
    candidates = []
    for mapping in ((0, 1, 2, 3, 4), (0, 2, 1, 3, 4)):
        oriented = tuple((mapping[a], mapping[b]) for a, b in record["ordered_edges"])
        edges = tuple(sorted(tuple(sorted(edge)) for edge in oriented))
        permutation = tuple(edges.index(tuple(sorted(edge))) for edge in oriented)
        locations = tuple(mapping[next(vertex for species, flavor, vertex, slot in record["surviving_fields"]
                                       if species == name)] for name in ("gamma", "c", "gamma5"))
        sign = record["wick_sign"] * jets.base._parity(permutation)
        exponent = [0] * 14
        for name, field in zip(INPUT_NAMES, INPUTS):
            pair = next(pair for pair in record["pairs"] if field in pair)
            original_edge = record["ordered_edges"].index(tuple(sorted((pair[0][2], pair[1][2]))))
            edge_sign = 1 if oriented[original_edge][0] < oriented[original_edge][1] else -1
            for component in range(2):
                power = powers.get((name, component), 0)
                exponent[2 * permutation[original_edge] + component] += power
                sign *= edge_sign ** power
        # This is the same canonical graph/source choice as the immutable
        # Wick ledger, independently of the requested derivative orders.
        initial_current = record["current_edge"]
        reference_sign = record["wick_sign"] * jets.base._parity(permutation)
        reference_sign *= 1 if oriented[initial_current][0] < oriented[initial_current][1] else -1
        reference = ((edges, locations), permutation[initial_current], reference_sign)
        candidates.append((reference, (edges, locations), tuple(exponent), sign))
    unused, key, exponent, sign = min(candidates, key=lambda row: row[0])
    return key, exponent, sign


class RationalJetProduct:
    """Complete rational tensor functional, with explicit term and output budgets.

    max_terms bounds polynomial expansions. max_entries bounds the output
    of tensor() and the total number of cached color-word values. A single
    moment() always returns the fixed ten color words, including zeros;
    it can run with a cache budget smaller than ten, without caching.
    """

    def __init__(self, vectors, *, max_terms=1000000, max_entries=100000, check=1):
        if len(vectors) != 3 or any(len(vector) != 2 for vector in vectors):
            raise ValueError("all three original source spinors are required")
        if any(isinstance(value, (float, complex, bool)) for vector in vectors for value in vector):
            raise TypeError("all sources must be exact rational values")
        self.vectors = tuple(tuple(map(Fraction, vector)) for vector in vectors)
        self.max_terms = jets.covariants.integer(max_terms, "max_terms", 1)
        self.max_entries = jets.covariants.integer(max_entries, "max_entries", 1)
        self.check = jets.covariants.integer(check, "check count")
        self.outer = {triangle: SourceTriangleOuterFaces(self.vectors, triangle, max_terms=self.max_terms)
                      for triangle in ("TWA", "TBA")}
        self.time = {triangle: SourceTriangleTimeFaces(self.vectors, triangle, max_terms=self.max_terms)
                     for triangle in ("TWA", "TBA")}
        self.context = fmpq_mpoly_ctx.get(tuple("chord%d" % i for i in range(6)))
        self._sp = {}
        self._cache = {}
        self.records = jets.universal.grouped_matchings()[1]
        self.words = tuple(sorted({jets.universal.matching_action(record)[0] for record in self.records
                                   if record["cubic_types"] == ("M", "M")}))
        if len(self.words) != 10:
            raise ArithmeticError("the complete rational color-word census changed")
        self.last_stats = {}

    def close(self):
        for planner, engine in self._sp.values():
            engine.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, traceback):
        self.close()

    def _word(self, slot, powers):
        wa, wb, ab = (self.context.gens()[start:start + 2] for start in (0, 2, 4))
        P, Q, R = self.vectors
        other_A = tuple(ab[i] - wa[i] - _q(R[i]) for i in range(2))
        current_b = tuple(_q(P[i] + Q[i] + R[i]) - wb[i] - ab[i] for i in range(2))
        currents = {"beta1": other_A if slot == 1 else wa,
                    "beta2": wa if slot == 1 else other_A,
                    "beta3": wb, "b": current_b,
                    "gamma5": tuple(self.context.constant(_q(value)) for value in R)}
        value = self.context.constant(1)
        for (name, component), power in powers.items():
            value *= _bounded_power(currents[name][component], power, self.context, self.max_terms)
            if len(value) > self.max_terms:
                raise RuntimeError("the complete transferred input word exceeded max_terms")
        return {tuple(map(int, exponent)): _fraction(coefficient) for exponent, coefficient in value.terms()}

    def pair_moment(self, polynomial):
        """Both complete pair functionals; the three individual periods are not used."""
        return tuple(self.outer[triangle].integrate(polynomial) - self.time[triangle].integrate(polynomial)
                     for triangle in ("TWA", "TBA"))

    def _sp_engine(self, key):
        if key not in self._sp:
            edges, locations = key
            sources = jets.base._sources(locations, self.vectors)
            planner = ForestContractions(jets.universal.wick.Graph(5, edges), sources, marked=(3, 4),
                                         max_terms=self.max_terms, check=self.check)
            # An exact topological gate makes a period-worker fallback impossible.
            plan_generic(planner.reduced_graph)
            engine = CoreCurrentPeriods(planner.reduced_graph, planner.reduced_sources,
                                        planner.reduced_marked, check=self.check, max_terms=self.max_terms)
            self._sp[key] = planner, engine
        return self._sp[key]

    def moment(self, powers=None):
        """All ten color words in the derivative-normalized convention.

        gamma5 powers multiply by the physical spectator derivative after
        source conversion; they are not raw geometric current insertions.
        """
        powers = _powers(powers)
        cache_key = tuple(sorted(powers.items()))
        if cache_key in self._cache:
            return dict(self._cache[cache_key])
        result = {word: Fraction(0) for word in self.words}
        for slot in (1, 2):
            first_pair, second_pair = self.pair_moment(self._word(slot, powers))
            sign = 1 if slot == 1 else -1
            result["MM", ((slot, ("A", "B")),)] = 2 * sign * (first_pair - second_pair)
            result["MM", ((slot, ("B", "A")),)] = 2 * sign * second_pair
        grouped = defaultdict(Fraction)
        for record in self.records:
            word, coefficient = jets.universal.matching_action(record)
            if word[0] != "MM" or len(word[1]) != 2:
                continue
            key, exponent, sign = _canonical_cross_word(record, powers)
            grouped[key, exponent, word] += sign * coefficient
        spectator = prod(self.vectors[2][component] ** powers.get(("gamma5", component), 0)
                         for component in range(2))
        for (key, exponent, word), coefficient in grouped.items():
            planner, engine = self._sp_engine(key)
            edge_word = {(edge, component): exponent[2 * edge + component]
                         for edge in range(7) for component in range(2) if exponent[2 * edge + component]}
            contracted = planner.contract(edge_word)
            value = engine.integrate(contracted.polynomial, allow_series=False)
            if not isinstance(value, Fraction):
                raise ArithmeticError("a rational SP current chart returned a nonrational value")
            result[word] += coefficient * spectator * value
        if 10 * (len(self._cache) + 1) <= self.max_entries:
            self._cache[cache_key] = dict(result)
        self.last_stats = {"complete_color_words": 10, "original_Wick_matchings": 28,
                           "time_scale_faces_per_triangle": 62, "source_independent_zero_faces_per_triangle": 46,
                           "SP_reduction_topologies": len(self._sp), "individual_period_lookups": 0,
                           "total_input_jet_order": sum(powers.values())}
        return result

    def tensor(self, orders):
        if not isinstance(orders, Mapping) or any(name not in ALL_INPUTS for name in orders):
            raise ValueError("tensor orders must name beta1,beta2,beta3,b,gamma5")
        ordered = tuple(jets.covariants.integer(orders.get(name, 0), name + " derivative order") for name in ALL_INPUTS)
        count = prod(order + 1 for order in ordered)
        if 10 * count > self.max_entries:
            raise RuntimeError("the complete open-color tensor exceeds max_entries; no partial tensor is returned")
        result = {word: {} for word in self.words}
        for entry in product(*(range(order + 1) for order in ordered)):
            powers = {(name, component): power for name, order, second in zip(ALL_INPUTS, ordered, entry)
                      for component, power in ((0, order - second), (1, second)) if power}
            for word, value in self.moment(powers).items():
                result[word][entry] = value
        return {"source_spinors_gamma4_c_gamma5": self.vectors,
                "input_order": ALL_INPUTS, "derivative_orders": ordered,
                "complete_symmetric_entries_per_color_word": count,
                "complete_open_color_tensor": result,
                "physical_multiplier": "+g^2/[128*(8*pi^2)^3]; source conversion already included",
                "spectator_derivatives": "physical gamma5 derivatives applied after source conversion; no propagator-head factor"}
