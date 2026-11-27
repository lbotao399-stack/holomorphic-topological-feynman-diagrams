"""Complete supported-source first-current spinors on a five-loop quartic forest.

Vertices 1, 2 and 5 carry arbitrary rational spinors, vertex 6 balances
their sum, and vertices 0, 3 and 4 have zero sources.  Every original edge
is retained.  Finite coefficients use the existing exact prism tables.
"""
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "code"))
sys.path.insert(0, str(HERE.parent / "symmetry"))

from htk import Graph
from htk.forest_chain_periods import ChainFirstPeriods
from htk.forest_contractions import normalize_forest_sources
from htk.forest_periods import ExactPeriod
from htk.treeflow import tree_flow
from tree_minor import prepare_source_numerator

EDGES = ((0, 1), (1, 3), (0, 2), (2, 4), (0, 5), (0, 6),
         (3, 5), (3, 6), (4, 5), (4, 6), (5, 6))
MARKED = (5, 6)
SOURCE_SUPPORT = (1, 2, 5, 6)
ZERO_SOURCE_VERTICES = (0, 3, 4)
DEGREE = 6
BASIS = tuple((i, degree - i) for degree in range(DEGREE + 1)
              for i in range(degree + 1))
POINTS = BASIS
DIRECT_POINTS = tuple((a, b) for a, b in POINTS if a <= b)
SWAP_VERTICES = (0, 2, 1, 4, 3, 5, 6)
SWAP_EDGES = (2, 3, 0, 1, 4, 5, 8, 9, 6, 7, 10)
LABELS = ("1", "z3", "z5")
OPTIONS = {"check": 1, "max_terms": 4000000,
           "max_period_terms": 1000000, "period_timeout": 300}


def graph():
    return Graph(7, EDGES, "five-loop-supported-source-current-family")


def canonical(value):
    return json.loads(json.dumps(value, sort_keys=True))


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True,
                                    separators=(",", ":")).encode()).hexdigest()


def coefficients(value):
    if isinstance(value, (int, Fraction)):
        return {"1": Fraction(value)} if value else {}
    if not isinstance(value, ExactPeriod):
        raise ArithmeticError("a complete finite period was required")
    return dict(value.coefficients)


def record(value):
    return {label: str(coefficient) for label, coefficient in sorted(value.items())
            if coefficient}


def add(terms):
    result = {}
    for weight, values in terms:
        weight = Fraction(weight)
        for label, coefficient in values.items():
            result[label] = result.get(label, Fraction()) + weight * Fraction(coefficient)
    return {label: coefficient for label, coefficient in result.items() if coefficient}


def normalized_sources(sources):
    sources = normalize_forest_sources(graph(), sources)
    if any(any(sources[vertex]) for vertex in ZERO_SOURCE_VERTICES):
        raise ValueError("this source family requires zero sources at vertices 0, 3 and 4")
    return sources


def sources_at(a, b):
    a, b = Fraction(a), Fraction(b)
    return ((0, 0), (1, 0), (0, 1), (0, 0), (0, 0),
            (a, b), (-1 - a, -1 - b))


def polynomial_record(polynomial):
    rows = [[list(map(int, exponent)), str(coefficient)]
            for exponent, coefficient in sorted(polynomial.terms())]
    return {"terms": len(rows), "sha256": digest(rows),
            "degrees": sorted({sum(exponent) for exponent, _ in rows})}


def reconstruct_currents(sources, bare, chord_values):
    G = graph()
    tree = G.spanning_tree()
    chords = G.cotree(tree)
    if set(chord_values) != set(chords):
        raise ArithmeticError("all five complete original chord spinors are required")
    cycles = G.fundamental_cycles(tree)
    flow = tree_flow(G, sources, tree)
    return [[add([(flow[edge][entry], bare)] +
                 [(cycles[k][edge], chord_values[chord][entry])
                  for k, chord in enumerate(chords)])
             for entry in range(2)] for edge in range(G.ne)]


def check_conservation(sources, bare, currents):
    for vertex in range(7):
        for entry in range(2):
            observed = add((int(left == vertex) - int(right == vertex), currents[edge][entry])
                           for edge, (left, right) in enumerate(EDGES))
            expected = add([(sources[vertex][entry], bare)])
            if observed != expected:
                raise ArithmeticError("an integrated original-source identity failed")


def swap_frame_record(source):
    """The exact graph automorphism and determinant-minus-one source map."""
    if any(tuple(SWAP_VERTICES[v] for v in edge) != EDGES[SWAP_EDGES[k]]
           for k, edge in enumerate(EDGES)):
        raise ArithmeticError("the oriented source-exchange automorphism changed")
    inversions = sum(a > b for i, a in enumerate(SWAP_EDGES) for b in SWAP_EDGES[i + 1:])
    if inversions % 2:
        raise ArithmeticError("the source-exchange edge permutation is not even")
    a, b = map(Fraction, source["full_balanced_sources"][5])
    bare = add([(-1, source["bare"])])
    old = source["all_original_first_current_spinors"]
    currents = [[add([(-1, old[SWAP_EDGES[edge]][1 - entry])])
                 for entry in range(2)] for edge in range(len(EDGES))]
    check_conservation(sources_at(b, a), bare, currents)
    return {"full_balanced_sources": [[str(value) for value in row] for row in sources_at(b, a)],
            "bare": record(bare), "all_original_first_current_spinors":
            [[record(value) for value in spinor] for spinor in currents],
            "exact_source_exchange_transport": {"from_point": [str(a), str(b)],
                "full_source_record_sha256": digest(source),
                "vertex_permutation": SWAP_VERTICES, "edge_permutation": SWAP_EDGES,
                "source_matrix": [[0, 1], [1, 0]], "source_matrix_determinant": -1},
            "integrated_original_source_identities": 14,
            "all_values_finite": True, "residual_periods": []}


def calculate_sources(sources, *, direct=False, tree=None, original=False,
                      scalar_archive=False, progress=None):
    """Integrate a whole sample; never use a partial spinor or a residual value."""
    from htk.forest_periods import CoreCurrentPeriods
    from flint import fmpq
    G = graph()
    sources = normalized_sources(sources)
    cls = CoreCurrentPeriods if original else ChainFirstPeriods
    core = cls(G, sources, MARKED, **OPTIONS)
    try:
        chosen = tuple(G.spanning_tree()) if tree is None else tuple(tree)
        plan = prepare_source_numerator(core, tree_candidates=[chosen],
                                        max_tree_subsets=1, max_plans=1)
        core._tree_polynomials()
        for vertex in range(G.nv):
            for entry in range(2):
                actual = sum(((int(left == vertex) - int(right == vertex)) *
                              core._currents[2 * edge + entry]
                              for edge, (left, right) in enumerate(EDGES)),
                             core._context.constant(0))
                source = sources[vertex][entry]
                if actual != fmpq(source.numerator, source.denominator) * core._U:
                    raise ArithmeticError("a complete original current polynomial changed source")
        bare = coefficients(core.integrate({(0,) * (2 * G.ne): Fraction(1)}, allow_series=False))
        stats = [dict(core.last_stats)]
        selected = tuple(range(G.ne)) if direct else tuple(G.cotree())
        integrated, archives = {}, []
        for edge in selected:
            integrated[edge] = []
            for entry in range(2):
                power = [0] * (2 * G.ne)
                power[2 * edge + entry] = 1
                polynomial = {tuple(power): Fraction(1)}
                if scalar_archive:
                    scalar = core._core_polynomials(polynomial)
                    archives.append({"edge": edge, "entry": entry,
                        "components": [[q, [[list(alpha), str(value)]
                                             for alpha, value in sorted(values.items())]]
                                       for q, values in sorted(scalar.items())]})
                value = coefficients(core.integrate(polynomial, allow_series=False))
                integrated[edge].append(value)
                stats.append(dict(core.last_stats))
            if progress:
                progress(edge, integrated[edge])
        currents = [integrated[edge] for edge in range(G.ne)] if direct else reconstruct_currents(
            sources, bare, integrated)
        check_conservation(sources, bare, currents)
        result = {"full_balanced_sources": [[str(value) for value in row] for row in sources],
                  "bare": record(bare),
                  "all_original_first_current_spinors": [[record(value) for value in spinor]
                                                          for spinor in currents],
                  "source_minor_plan": plan,
                  "original_numerator": polynomial_record(core._bare),
                  "original_current_numerators": [polynomial_record(value) for value in core._currents],
                  "all_original_density_and_current_point_checks": True,
                  "coefficientwise_original_source_identities": 14,
                  "integrated_original_source_identities": 14,
                  "independently_integrated_original_edges": selected,
                  "used_cotree_reconstruction": not direct,
                  "used_chain_first_aggregation": not original,
                  "statistics": stats, "all_values_finite": True, "residual_periods": []}
        if scalar_archive:
            result["complete_scalar_numerators"] = archives
        return canonical(result)
    finally:
        core.close()


def interpolate(records):
    from flint import fmpq_mat
    if len(records) != len(POINTS):
        raise ArithmeticError("all 28 exact interpolation samples are required")
    matrix = fmpq_mat([[a ** i * b ** j for i, j in BASIS] for a, b in POINTS])
    inverse = matrix.inv()
    weights = [[Fraction(int(inverse[k, n].p), int(inverse[k, n].q))
                for n in range(len(POINTS))] for k in range(len(BASIS))]
    def recover(values):
        return [record(add(zip(row, values))) for row in weights]
    bare = recover([row["bare"] for row in records])
    currents = [[recover([row["all_original_first_current_spinors"][edge][entry]
                          for row in records]) for entry in range(2)]
                for edge in range(len(EDGES))]
    if any(bare[k] for k, (i, j) in enumerate(BASIS) if i + j == 6):
        raise ArithmeticError("the scalar source polynomial exceeded bracket degree five")
    covariants = []
    for first, second in currents:
        high = {}
        for k, (i, j) in enumerate(BASIS):
            if i + j != 6:
                continue
            if not i and first[k]:
                raise ArithmeticError("the degree-six first entry is not divisible by a")
            if not j and second[k]:
                raise ArithmeticError("the degree-six second entry is not divisible by b")
            if i:
                high[i - 1, j] = first[k]
        for k, (i, j) in enumerate(BASIS):
            if i + j == 6 and j and second[k] != high.get((i, j - 1), {}):
                raise ArithmeticError("degree-six spinor covariance b*f6=a*g6 failed")
        covariants.append([
            [[list(power), values] for power, values in zip(BASIS, entry)
             if sum(power) <= 5 and values] for entry in (first, second)] +
            [[[list(power), values] for power, values in sorted(high.items()) if values]])
    return {"frame_monomials": BASIS, "bare_frame_coefficients": bare,
            "full_current_frame_coefficients": currents,
            "denominator_free_original_source_covariants": covariants,
            "degree_six_covariance_verified_on_every_original_edge": True}


def bracket(first, second):
    return first[0] * second[1] - first[1] * second[0]


def evaluate(reconstruction, sources):
    """Evaluate all 11 full spinors, including every degenerate source frame."""
    sources = normalized_sources(sources)
    first, second, marked = sources[1], sources[2], sources[5]
    determinant = bracket(first, second)
    first_weight, second_weight = bracket(marked, second), bracket(first, marked)
    def monomial(i, j):
        return determinant ** (5 - i - j) * first_weight ** i * second_weight ** j
    bare = add((monomial(i, j), values) for (i, j), values in
               zip(reconstruction["frame_monomials"], reconstruction["bare_frame_coefficients"])
               if i + j <= 5)
    currents = []
    for covariant in reconstruction["denominator_free_original_source_covariants"]:
        direction_values = [add((monomial(i, j), values) for (i, j), values in rows)
                            for rows in covariant]
        currents.append([add((source[entry], values) for source, values in
                             zip((first, second, marked), direction_values))
                         for entry in range(2)])
    check_conservation(sources, bare, currents)
    return {"bare": record(bare), "all_original_first_current_spinors":
            [[record(value) for value in spinor] for spinor in currents]}


def five_loop_source_currents(sources, *, certificate=None):
    """Return the certified bare value and every complete original spinor."""
    document = json.loads((HERE / "certificate.json").read_text()) if certificate is None else certificate
    if (document.get("schema_version") != 1 or document.get("edges") != [list(edge) for edge in EDGES]
            or document.get("source_support") != list(SOURCE_SUPPORT)
            or not document.get("complete")):
        raise ValueError("the certificate does not describe the complete supported source family")
    return evaluate(document["reconstruction"], sources)
