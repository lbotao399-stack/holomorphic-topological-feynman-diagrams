"""Exact higher-loop lifts of the complete three-loop forest Ward identity.

Replace the two spokes at vertex 1 by nested source-free paired stars and
insert the deepest first A-spoke current.  Every source row, edge direction,
and edge position is retained.  Actual higher-loop quotient kernels are
evaluated by ForestContractions followed by CoreCurrentPeriods.  A separate
factorial transfer is checked as a full current polynomial, then integrated
on the rational diamond--triangle boundary by explicit beta moments.

Five- and seven-loop graphs have exhaustive tree/proper-face certificates.
Larger examples are exact inductive lifts of the base all-current-word Ward
identity; their certificates do not assert an unperformed full face census.
The object proved here is a graph identity, before Wick/color assembly.

Run all recorded examples:
    python3 explore/circuit-general/forest-leibniz/lifted_ward.py
Regenerate the deterministic certificate:
    python3 explore/circuit-general/forest-leibniz/lifted_ward.py --write
Run a selected source assignment and number of paired insertions:
    python3 explore/circuit-general/forest-leibniz/lifted_ward.py --pairs 2 --sources physical
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
from itertools import combinations
import json
from math import comb, factorial
from numbers import Integral
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "code"))

from flint import fmpq_mpoly_ctx
from htk.graph import Graph
from htk.forest_contractions import ForestContractions
from htk.forest_periods import CoreCurrentPeriods
from htk.linalg import perm_sign
from ward_three import (EDGES, MARKED, GENERIC_SOURCES, PHYSICAL_SOURCES,
                        _sources, _q, _fraction, coefficients,
                        boundary_first_spinor, diamond_first_spinor)


PATH_EDGES = ((0, 1), (0, 5), (2, 5))
SOURCE_CASES = {"generic": GENERIC_SOURCES, "physical": PHYSICAL_SOURCES}
CERTIFIED_CASES = (("generic", 1), ("generic", 2),
                   ("physical", 1), ("physical", 2),
                   ("physical", 3), ("physical", 5), ("physical", 8))


def _integer(value, name, minimum=0):
    if isinstance(value, bool) or not isinstance(value, Integral) or value < minimum:
        raise ValueError("%s must be an integer at least %d" % (name, minimum))
    return int(value)


def _bracket(left, right):
    return left[0] * right[1] - left[1] * right[0]


def _polynomial_dict(polynomial):
    return {tuple(map(int, exponent)): _fraction(coefficient)
            for exponent, coefficient in polynomial.terms() if coefficient}


def reconstruct_last(first, middle, rational_boundary, *, signs=(1, -1, 1)):
    """Recover the third complete period expression without integrating it.

    ``first`` and ``middle`` are raw kernel values (Fraction or exact period
    dictionaries). The three supplied signs are those of the actual ordered
    Ward relation. In the physical family the middle value is rational, so
    this computes the other period-valued forest from one such forest and
    two rational quantities. No period coefficient is discarded.
    """
    signs = tuple(signs)
    if len(signs) != 3 or any(isinstance(sign, bool) or not isinstance(sign, Integral)
                              or sign not in (-1, 1) for sign in signs):
        raise ValueError("three actual Ward signs, each +1 or -1, are required")
    boundary = coefficients(rational_boundary)
    if any(label != "1" for label in boundary):
        raise ValueError("this lifted family's boundary is rational")
    value = dict(boundary)
    for sign, expression in zip(signs[:2], (first, middle)):
        for label, scalar in coefficients(expression).items():
            value[label] = value.get(label, F(0)) - sign * scalar
    return {label: signs[2] * scalar for label, scalar in value.items() if scalar}


def transfer_coefficient(pairs):
    """Complete htk coefficient in the stated inherited edge order."""
    pairs = _integer(pairs, "paired insertions")
    return F(2 * (-1) ** pairs,
             factorial(pairs) * factorial(pairs + 1) ** 2
             * factorial(pairs + 2) * factorial(2 * pairs + 1))


def lifted_fixture(pairs, sources=GENERIC_SOURCES, component=0):
    pairs = _integer(pairs, "paired insertions")
    component = _integer(component, "spinor component")
    if component >= 2:
        raise ValueError("the spinor has two components")
    edges, rows, parent = list(EDGES), list(_sources(sources)), 1
    for _ in range(pairs):
        edges.remove((parent, MARKED[0]))
        edges.remove((parent, MARKED[1]))
        first, second = len(rows), len(rows) + 1
        rows.extend(((F(0), F(0)), (F(0), F(0))))
        edges.extend(((first, parent), (first, MARKED[0]), (first, MARKED[1]),
                      (second, parent), (second, MARKED[0]), (second, MARKED[1])))
        parent = first
    graph = Graph(len(rows), edges, "lifted-forest-Ward-%d" % (3 + 2 * pairs))
    return graph, tuple(rows), MARKED, {(edges.index((parent, MARKED[0])), component): 1}


def _quotient(graph, sources, position):
    left, right = graph.edges[position]
    vertices = tuple(v for v in range(graph.nv) if v != right)
    mapping = {v: vertices.index(left if v == right else v) for v in range(graph.nv)}
    edges = tuple((mapping[a], mapping[b]) for e, (a, b) in enumerate(graph.edges) if e != position)
    rows = [[F(0), F(0)] for _ in vertices]
    for vertex, row in enumerate(sources):
        for component in range(2):
            rows[mapping[vertex]][component] += row[component]
    return Graph(len(vertices), edges), tuple(map(tuple, rows)), mapping


def _overfull(nv, edges):
    """An explicit induced-subgraph witness for an identically zero kernel."""
    for size in range(2, nv):
        for vertices in combinations(range(nv), size):
            members = set(vertices)
            inside = tuple(e for e, edge in enumerate(edges) if set(edge) <= members)
            if len(inside) > 2 * size - 3:
                return {"vertices": vertices, "edge_indices": inside,
                        "edge_count": len(inside), "allowed_count": 2 * size - 3}
    return None


@lru_cache(maxsize=None)
def face_census(pairs):
    """Exhaustive census only for the small graphs explicitly advertised."""
    if pairs not in (1, 2):
        raise ValueError("the exhaustive face census is restricted to five and seven loops")
    graph, sources, _, _ = lifted_fixture(pairs)
    tree_faces = []
    for position in range(graph.ne):
        quotient, _, mapping = _quotient(graph, sources, position)
        tree_faces.append({"position": position, "edge": graph.edges[position],
                           "sign": (-1) ** position, "quotient_edges": tuple(quotient.edges),
                           "vertex_map": tuple(mapping[v] for v in range(graph.nv)),
                           "zero_witness": _overfull(quotient.nv, quotient.edges)})
    survivors = tuple(row["edge"] for row in tree_faces if row["zero_witness"] is None)
    if survivors != PATH_EDGES:
        raise ArithmeticError("unexpected higher-loop tree faces")

    proper = []
    for size in range(3, graph.nv):
        for cut in combinations(range(graph.nv), size):
            members = set(cut)
            inside = tuple(e for e, edge in enumerate(graph.edges) if set(edge) <= members)
            if len(inside) != 2 * size - 3:
                continue
            outside = tuple(e for e in range(graph.ne) if e not in inside)
            vertices = tuple(v for v in range(graph.nv) if v not in members) + (cut[0],)
            mapping = {v: vertices.index(cut[0] if v in members else v) for v in range(graph.nv)}
            inner = tuple((cut.index(graph.edges[e][0]), cut.index(graph.edges[e][1])) for e in inside)
            outer = tuple((mapping[a], mapping[b]) for e, (a, b) in enumerate(graph.edges) if e in outside)
            inner_loops, outer_loops = size - 2, graph.loops - size + 2
            proper.append({"cut": cut, "inner_original_edges": inside,
                           "outer_original_edges": outside, "inner_edges": inner, "outer_edges": outer,
                           "inner_zero_witness": _overfull(len(cut), inner),
                           "outer_zero_witness": _overfull(len(vertices), outer),
                           "sign": -(-1) ** (inner_loops * outer_loops) * perm_sign(inside + outside)})
    survivors = [row for row in proper if row["inner_zero_witness"] is None
                 and row["outer_zero_witness"] is None]
    expected = tuple(v for v in range(graph.nv) if v not in (0, 5))
    if len(survivors) != 1 or survivors[0]["cut"] != expected or survivors[0]["sign"] != 1:
        raise ArithmeticError("unexpected higher-loop proper faces")
    return {"tree_faces": tuple(tree_faces), "proper_faces": tuple(proper)}


def reference_polynomial(pairs, component, base_contracted_position):
    """Independent factorial transfer in the base quotient's inherited order.

    It is c_k y_13^a <y_13,y_34>^k <y_14,y_34>^k, with complete
    coefficient c_k. Its current degree is 4k+1, not 2k+1.
    """
    context = fmpq_mpoly_ctx.get(tuple("j%d" % i for i in range(14)))
    variables = context.gens()
    remaining = tuple(e for e in range(8) if e != base_contracted_position)
    left = variables[2 * remaining.index(2):2 * remaining.index(2) + 2]
    right = variables[2 * remaining.index(3):2 * remaining.index(3) + 2]
    marked = variables[2 * remaining.index(7):2 * remaining.index(7) + 2]
    polynomial = (_q(transfer_coefficient(pairs)) * left[component]
                  * _bracket(left, marked) ** pairs * _bracket(right, marked) ** pairs)
    return _polynomial_dict(polynomial)


def inner_boundary_spinor(pairs, P, Q, A):
    """Closed beta-integrated lifted diamond, in full spinor brackets.

    Only the following finite sums remain: integrating (v+t*u)^(2k+1)
    against 1 or t gives denominators j+1 or j+2. The square moments are
    beta(k+1,k+1), beta(k+2,k+1), beta(k+3,k+1). This expression does not
    call either forest engine or any period table.
    """
    pairs = _integer(pairs, "paired insertions")
    if any(isinstance(value, float) for vector in (P, Q, A) for value in vector):
        raise TypeError("the exact inner boundary requires rational sources")
    P, Q, A = (tuple(_q(value) if isinstance(value, (F, Integral)) else value for value in vector)
               for vector in (P, Q, A))
    u, v, w = _bracket(P, Q), _bracket(P, A), _bracket(Q, A)
    power, divisor = 2 * pairs + 1, 2 * pairs + 4
    beta0 = _q(F(factorial(pairs) ** 2, factorial(2 * pairs + 1)))
    beta1 = _q(F(factorial(pairs + 1) * factorial(pairs), factorial(2 * pairs + 2)))
    beta2 = _q(F(factorial(pairs + 2) * factorial(pairs), factorial(2 * pairs + 3)))
    integral0 = sum(_q(F(comb(power, j), j + 1)) * v ** (power - j) * u ** j
                    for j in range(power + 1))
    integral1 = sum(_q(F(comb(power, j), j + 2)) * v ** (power - j) * u ** j
                    for j in range(power + 1))
    scalar = _q(-transfer_coefficient(pairs) / ((2 * pairs + 2) * (2 * pairs + 3) * divisor))
    return tuple(scalar * ((divisor - 1) * P[a] * (w * beta1 - u * beta2) * integral0
                           - (w * beta0 - u * beta1) * (A[a] * integral0 + Q[a] * integral1))
                 for a in range(2))


def boundary_spinor(pairs, sources, *, max_terms=1000000):
    """Independent complete rational boundary, with no global ratio chart."""
    sources = _sources(sources)
    max_terms = _integer(max_terms, "max_terms", 1)
    context = fmpq_mpoly_ctx.get(("x", "z"))
    x, z = context.gens()
    first, second = tuple(map(_q, sources[0])), tuple(map(_q, sources[5]))
    outer_first = tuple((1 - x) * first[a] + z * second[a] for a in range(2))
    outer_last = tuple(-x * first[a] + (z - 1) * second[a] for a in range(2))
    P = tuple(_q(sources[1][a]) + outer_first[a] for a in range(2))
    Q = tuple(_q(sources[2][a]) - outer_last[a] for a in range(2))
    A = tuple(context.constant(_q(sources[3][a])) for a in range(2))
    polynomials = inner_boundary_spinor(pairs, P, Q, A)
    values = []
    for polynomial in polynomials:
        if len(polynomial) > max_terms:
            raise RuntimeError("the exact outer triangle polynomial exceeds max_terms")
        value = F(0)
        for exponent, coefficient in polynomial.terms():
            a, b = map(int, exponent)
            value += _fraction(coefficient) * F(factorial(a) * factorial(b), factorial(a + b + 2))
        values.append(_bracket(sources[0], sources[5]) * value)
    return tuple(values)


def _compare_transfers(actual, reference, higher_map, base_map):
    """Compare complete polynomials after only vertex/edge conventions change."""
    vertex_map = {}
    for vertex in range(6):
        a = actual.original_to_reduced_vertices[higher_map[vertex]]
        b = reference.original_to_reduced_vertices[base_map[vertex]]
        if (a is None) != (b is None):
            raise ArithmeticError("the analytic and actual lifts retain different base vertices")
        if a is not None:
            vertex_map[a] = b
    if set(vertex_map) != set(range(actual.graph.nv)) or set(vertex_map.values()) != set(range(reference.graph.nv)):
        raise ArithmeticError("a passive inserted vertex survived the planned lift")
    if any(actual.sources[a] != reference.sources[b] for a, b in vertex_map.items()):
        raise ArithmeticError("the transfer comparison changed a source row")
    if tuple(vertex_map[v] for v in actual.marked) != reference.marked:
        raise ArithmeticError("the transfer comparison changed a marked vertex")
    context = fmpq_mpoly_ctx.get(tuple("j%d" % i for i in range(2 * reference.graph.ne)))
    variables, order, replacements = context.gens(), [], []
    for a, b in actual.graph.edges:
        edge = (vertex_map[a], vertex_map[b])
        if edge in reference.graph.edges:
            location, sign = reference.graph.edges.index(edge), 1
        elif edge[::-1] in reference.graph.edges:
            location, sign = reference.graph.edges.index(edge[::-1]), -1
        else:
            raise ArithmeticError("the actual reduction changed a core edge")
        order.append(location)
        replacements.extend(sign * variables[2 * location + component] for component in range(2))
    original = fmpq_mpoly_ctx.get(tuple("j%d" % i for i in range(2 * actual.graph.ne)))
    actual_value = original.from_dict({a: _q(c) for a, c in actual.polynomial.items()}).compose(*replacements)
    expected = context.from_dict({a: _q(c) for a, c in reference.polynomial.items()})
    if perm_sign(order) * actual_value != expected:
        raise ArithmeticError("the exact paired-star transfer disagrees with the factorial formula")
    return {"vertex_map": tuple(vertex_map[v] for v in range(actual.graph.nv)),
            "edge_permutation": tuple(order), "kernel_permutation_sign": perm_sign(order),
            "complete_current_polynomial_equal": True}


class LiftedWard:
    """Cache all three actual quotient engines and both spinor insertions."""
    def __init__(self, pairs, sources, *, check=1, max_terms=1000000):
        self.pairs = _integer(pairs, "paired insertions", 1)
        self.sources = _sources(sources)
        self.check = _integer(check, "check count")
        self.max_terms = _integer(max_terms, "max_terms", 1)
        self.graph, self.full_sources, self.marked, first_word = lifted_fixture(self.pairs, self.sources)
        self.inserted_edge = next(iter(first_word))[0]
        self.census = face_census(self.pairs) if self.pairs <= 2 else None
        self.boundary = boundary_spinor(self.pairs, self.sources, max_terms=self.max_terms)
        self.quotients, self._results = [], {}
        base = Graph(6, EDGES)
        for edge in PATH_EDGES:
            position, base_position = self.graph.edges.index(edge), EDGES.index(edge)
            quotient, rows, mapping = _quotient(self.graph, self.full_sources, position)
            base_quotient, base_rows, base_mapping = _quotient(base, self.sources, base_position)
            marked = tuple(mapping[v] for v in MARKED)
            reduction = ForestContractions(quotient, rows, marked, check=self.check, max_terms=self.max_terms)
            reference = ForestContractions(base_quotient, base_rows, tuple(base_mapping[v] for v in MARKED),
                                           check=self.check, max_terms=self.max_terms)
            engine = CoreCurrentPeriods(reduction.reduced_graph, reduction.reduced_sources,
                                        reduction.reduced_marked, check=self.check, max_terms=self.max_terms)
            self.quotients.append({"position": position, "base_position": base_position,
                                   "contracted_edge": edge, "sign": (-1) ** position,
                                   "graph": quotient, "sources": rows, "mapping": mapping,
                                   "base_mapping": base_mapping, "reduction": reduction,
                                   "reference": reference, "engine": engine})

    def close(self):
        for quotient in self.quotients:
            quotient["engine"].close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, traceback):
        self.close()

    def reconstruct_last(self, component, first, middle):
        """Use this graph's actual signs and complete independent boundary."""
        component = _integer(component, "spinor component")
        if component >= 2:
            raise ValueError("the spinor has two components")
        return reconstruct_last(first, middle, self.boundary[component],
                                signs=tuple(quotient["sign"] for quotient in self.quotients))

    def verify(self, component):
        component = _integer(component, "spinor component")
        if component >= 2:
            raise ValueError("the spinor has two components")
        if component in self._results:
            return self._results[component]
        terms, total = [], {}
        for quotient in self.quotients:
            position = quotient["position"]
            current_edge = self.inserted_edge - (self.inserted_edge > position)
            actual = quotient["reduction"].contract({(current_edge, component): 1})
            formula = reference_polynomial(self.pairs, component, quotient["base_position"])
            reference = quotient["reference"].contract(polynomial=formula)
            comparison = _compare_transfers(actual, reference, quotient["mapping"], quotient["base_mapping"])
            value = coefficients(quotient["engine"].integrate(actual.polynomial, allow_series=False))
            signed = {key: quotient["sign"] * scalar for key, scalar in value.items()}
            for key, scalar in signed.items():
                total[key] = total.get(key, F(0)) + scalar
            terms.append({"position": position, "base_position": quotient["base_position"],
                          "contracted_edge": quotient["contracted_edge"], "sign": quotient["sign"],
                          "edges": tuple(quotient["graph"].edges), "sources": quotient["sources"],
                          "vertex_map": tuple(quotient["mapping"][v] for v in range(self.graph.nv)),
                          "reduced_edges": tuple(actual.graph.edges), "reduced_sources": actual.sources,
                          "reduced_marked": actual.marked,
                          "current_polynomial": tuple((a, c) for a, c in sorted(actual.polynomial.items())),
                          "contraction_steps": actual.steps, "contraction_stats": actual.stats,
                          "analytic_transfer_comparison": comparison,
                          "integration_stats": dict(quotient["engine"].last_stats),
                          "coefficients": value, "signed_coefficients": signed})
        total = {key: scalar for key, scalar in total.items() if scalar}
        if total != coefficients(self.boundary[component]):
            raise ArithmeticError("the complete lifted period expression disagrees with its rational boundary")
        reconstructed = self.reconstruct_last(component, terms[0]["coefficients"], terms[1]["coefficients"])
        if reconstructed != terms[2]["coefficients"]:
            raise ArithmeticError("the last-quotient reconstruction changed a period coefficient")
        result = {"original_word": ((self.inserted_edge, component, 1),),
                  "terms": tuple(terms), "signed_sum": total,
                  "rational_boundary": self.boundary[component]}
        self._results[component] = result
        return result


def _json(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _json(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json(item) for item in value]
    return value


def certificate(cases=CERTIFIED_CASES):
    records = []
    for name, pairs in cases:
        with LiftedWard(pairs, SOURCE_CASES[name]) as ward:
            checks = []
            for component in range(2):
                result = ward.verify(component)
                checks.append(result)
                print("%s %d loops, deepest A-spoke component %d: %s" %
                      (name, ward.graph.loops, component, result["rational_boundary"]), flush=True)
            records.append({"sources_name": name, "paired_insertions": pairs, "loops": ward.graph.loops,
                            "edges": tuple(ward.graph.edges), "sources": ward.full_sources,
                            "marked": ward.marked, "transfer_coefficient": transfer_coefficient(pairs),
                            "raw_transferred_current_degree": 4 * pairs + 1,
                            "face_analysis": ("exhaustive graph-level Stokes faces" if ward.census is not None
                                              else "exact inductive lift of the base all-current-word Ward identity"),
                            "face_census": ward.census, "complete_spinor_boundary": ward.boundary,
                            "checks": tuple(checks)})
    paths = (Path(__file__), HERE / "ward_three.py", ROOT / "code/htk/forest_contractions.py",
             ROOT / "code/htk/forest_periods.py", ROOT / "code/htk/_forest_period_worker.py")
    return _json({"scope": "Exact graph identities for fully retained sources; Wick/color/operator assembly is not asserted.",
                  "max_terms": 1000000, "records": records,
                  "source_hashes": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                                    for path in paths}})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--pairs", type=int, nargs="+")
    parser.add_argument("--sources", choices=tuple(SOURCE_CASES), default="physical")
    arguments = parser.parse_args()
    cases = CERTIFIED_CASES if arguments.pairs is None else tuple((arguments.sources, k) for k in arguments.pairs)
    actual = certificate(cases)
    path = HERE / "lifted_ward_certificate.json"
    if arguments.write:
        if arguments.pairs is not None:
            parser.error("--write regenerates the complete certificate; do not combine it with --pairs")
        path.write_text(json.dumps(actual, indent=2, sort_keys=True) + "\n")
        print("Wrote the complete lifted Ward certificate.", flush=True)
    elif arguments.pairs is None:
        if actual != json.loads(path.read_text()):
            raise ArithmeticError("the lifted Ward certificate changed")
        print("Every complete lifted Ward certificate field replays.", flush=True)


if __name__ == "__main__":
    main()
