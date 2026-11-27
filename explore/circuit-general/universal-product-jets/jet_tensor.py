"""Complete higher-input-jet tensor with all sources and color indices.

This extends the immutable universal-quartic-product Wick ledger. The q
derivatives all act on the original beta2 field at marked A. The complete
rank-q symmetric spinor tensor is retained, as are all twelve color words
before and all ten words after the representation relation. The physical
manuscript t-word multiplier remains +g^2/[128*(8*pi^2)^3] for every q.

The full normalized source degree is at most q+3, so (q+4)*(q+5)/2
exact triangular nodes determine the entire covariant tensor. The default
certificate calculation reintegrates all geometric data; no old source or
certificate is modified and no period table is built.
"""
from collections import defaultdict
from fractions import Fraction
import argparse
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE.parent / "universal-quartic-product"))
import quartic_tensor as universal
import covariants

base = universal.base


def _clean(values):
    return {word: tuple({label: coefficient for label, coefficient in entry.items() if coefficient}
                        for entry in tensor) for word, tensor in sorted(values.items())}


def evaluate_sources(vectors, orders, *, check=1, max_terms=1000000):
    orders = tuple(sorted({covariants.integer(order, "jet order") for order in orders}))
    if not orders:
        raise ValueError("at least one complete jet order is required")
    if len(vectors) != 3 or any(len(vector) != 2 for vector in vectors):
        raise ValueError("all three complete source spinors are required")
    if any(isinstance(value, float) for vector in vectors for value in vector):
        raise TypeError("source spinors must be exact")
    vectors = tuple(tuple(map(Fraction, vector)) for vector in vectors)
    groups, records, rejected = universal.grouped_matchings()
    if any(record["ordered_edges"][record["current_edge"]][1] != 3 for record in records):
        raise ArithmeticError("the original differentiated head was not retained")
    values = {order: defaultdict(lambda order=order: [defaultdict(Fraction) for k in range(order + 1)])
              for order in orders}
    graphs = []
    for (edges, locations), edge_words in groups.items():
        sources = base._sources(locations, vectors)
        engine = universal.wick.CoreCurrentPeriods(universal.wick.Graph(5, edges), sources,
                                                  marked=(3, 4), check=check, max_terms=max_terms)
        moments = {}
        try:
            for order in orders:
                moments[order] = {}
                for edge, words in sorted(edge_words.items()):
                    entries = []
                    for k in range(order + 1):
                        exponent = [0] * 14
                        exponent[2 * edge], exponent[2 * edge + 1] = order - k, k
                        complete = base._periods(engine.integrate({tuple(exponent): Fraction(1)}, allow_series=False))
                        entries.append(complete)
                        for word, coefficient in words.items():
                            values[order][word]
                            for label, value in complete.items():
                                values[order][word][k][label] += coefficient * value
                    moments[order][edge] = tuple(entries)
        finally:
            engine.close()
        graphs.append({"ordered_edges": edges, "source_vertices_gamma4_c_gamma5": locations,
                       "full_balanced_sources": sources, "all_original_current_tensors": moments})
    return {"source_spinors_gamma4_c_gamma5": vectors, "all_geometric_tensors": graphs,
            "complete_formal_tensors": {order: _clean(value) for order, value in values.items()}}


def assemble_stored_sample(sample, order):
    """Reassemble every tensor word from all stored original-edge moments."""
    groups = universal.grouped_matchings()[0]
    graph_rows = {(tuple(map(tuple, row["ordered_edges"])), tuple(row["source_vertices_gamma4_c_gamma5"])): row
                  for row in sample["all_geometric_tensors"]}
    values = defaultdict(lambda: [defaultdict(Fraction) for k in range(order + 1)])
    for key, edges in groups.items():
        row = graph_rows[key]
        moments = row["all_original_current_tensors"].get(order, row["all_original_current_tensors"].get(str(order)))
        if moments is None:
            raise ValueError("the complete requested jet order is absent from this sample")
        for edge, words in edges.items():
            tensor = moments.get(edge, moments.get(str(edge)))
            if tensor is None or len(tensor) != order + 1:
                raise ValueError("an original current tensor entry was omitted")
            for word, coefficient in words.items():
                values[word]
                for k, entry in enumerate(tensor):
                    for label, value in entry.items():
                        values[word][k][label] += coefficient * Fraction(value)
    return _clean(values)


def restore_tensor(samples, order, *, max_terms=1000000):
    order = covariants.integer(order, "jet order")
    grid = {(a, b) for a in range(order + 4) for b in range(order + 4 - a)}
    values = {tuple(map(Fraction, sample["source_spinors_gamma4_c_gamma5"][2])): assemble_stored_sample(sample, order)
              for sample in samples}
    if set(values) != grid:
        raise ValueError("the entire proven-degree source grid is required")
    words = sorted(next(iter(values.values())))
    if len(words) != 12 or any(sorted(tensor) != words for tensor in values.values()):
        raise ArithmeticError("the complete open-color tensor was not retained")
    labels = {label for tensor in values.values() for entries in tensor.values() for entry in entries for label in entry}
    result = {}
    for word in words:
        channels = {}
        for label in sorted(labels):
            polynomials = tuple(base.interpolate({node: tensor[word][k].get(label, 0) for node, tensor in values.items()},
                                                 degree=order + 3) for k in range(order + 1))
            restored = covariants.restore(polynomials, order, max_terms=max_terms)
            if covariants.normalized_polynomials(restored, order) != polynomials:
                raise ArithmeticError("covariant restoration lost an entry of the original complete tensor")
            if restored:
                channels[label] = restored
        result[word] = channels
    return result


def representation_relation(tensor):
    result = defaultdict(lambda: defaultdict(lambda: defaultdict(Fraction)))
    for word, channels in tensor.items():
        replacements = ((word, 1),) if word[0] == "MM" else (
            (("MM", ((word[1], ("A", "B")),)), 1), (("MM", ((word[1], ("B", "A")),)), -1))
        for new_word, sign in replacements:
            for label, terms in channels.items():
                for counts, powers, coefficient in terms:
                    result[new_word][label][counts, powers] += sign * coefficient
    return {word: {label: tuple((counts, powers, coefficient) for (counts, powers), coefficient in sorted(terms.items())
                                if coefficient) for label, terms in sorted(channels.items()) if any(terms.values())}
            for word, channels in sorted(result.items())}


def _serial_tensor(tensor):
    return [{"generator_word": word, "complete_covariant_tensor": channels} for word, channels in sorted(tensor.items())]


def certificate_from_samples(order, samples, replays):
    tensor = restore_tensor(samples, order)
    reduced = representation_relation(tensor)
    defects = []
    for slot in (1, 2):
        ab = ("MM", ((slot, ("A", "B")),))
        ba = ("MM", ((slot, ("B", "A")),))
        fr = ("fR", slot)
        for label, terms in tensor[ab].items():
            if label == "1":
                continue
            opposite = tuple((counts, powers, -coefficient) for counts, powers, coefficient in terms)
            if tensor[ba].get(label) == opposite and tensor[fr].get(label) == opposite:
                defects.append({"W_slot": slot, "period_label": label,
                                "complete_covariant_tensor_coefficient": terms,
                                "exact_color_factor": "[t_A,t_B]-i*C^C_AB*t_C"})
    identities = []
    for sample in replays:
        direct = assemble_stored_sample(sample, order)
        vectors = sample["source_spinors_gamma4_c_gamma5"]
        for word, channels in tensor.items():
            predicted = [dict() for k in range(order + 1)]
            for label, terms in channels.items():
                for k, value in enumerate(covariants.evaluate(terms, vectors, order)):
                    if value:
                        predicted[k][label] = value
            if tuple(predicted) != direct[word]:
                raise ArithmeticError("an independent complete tensor/source replay disagrees")
        identities.append({"source_spinors_gamma4_c_gamma5": vectors,
                           "every_original_color_word_and_symmetric_entry_matches": True})
    paths = (HERE / "jet_tensor.py", HERE / "covariants.py", universal.HERE / "quartic_tensor.py",
             universal.PREVIOUS / "nonabelian_wick.py", universal.PREVIOUS / "quartic_wick.py",
             ROOT / "code/htk/forest_periods.py", ROOT / "code/htk/_forest_period_worker.py",
             ROOT / "explore/circuit-general/fiveloop/per5.py", ROOT / "explore/circuit-general/fiveloop/msys5.py",
             ROOT / "explore/circuit-general/cores/wheel.py")
    return {"schema": 1, "jet_order": order,
            "scope": "complete pure-quartic W coefficient with every input/output color index and every lower spinor index retained",
            "input_derivative_routing": "all ordered holomorphic derivatives act on original beta2 at marked A; no derivative is reassigned",
            "source_total_degree": 6 + order, "spectator_source_degree_bound": 3 + order,
            "complete_grid_count": (order + 4) * (order + 5) // 2,
            "canonical_covariant_coefficient_bound_per_word": 2 * (order + 1) * (order + 5),
            "covariant_source_order": ["lambda_gamma4", "lambda_c", "lambda_gamma5"],
            "covariant_bracket_order": ["<lambda_gamma4,lambda_c>", "<lambda_gamma5,lambda_c>", "<lambda_gamma4,lambda_gamma5>"],
            "covariant_term_semantics": "(three source multiplicities, three bracket powers, exact coefficient); tensor products are averaged over all permutations of all free lower +dot a indices",
            "normalization": "+g^2/[128*(8*pi^2)^3] in manuscript t words after each source is replaced by its own surviving-field derivative; source conversion already included",
            "implementation_sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                                      for path in paths},
            "complete_source_samples": samples, "independent_full_source_replays": replays,
            "full_replay_identities": identities,
            "complete_unreduced_tensor": _serial_tensor(tensor), "complete_reduced_tensor": _serial_tensor(reduced),
            "complete_representation_defect_coefficients": defects,
            "complete_representation_cancellation": all(set(channels) <= {"1"} for channels in reduced.values()),
            "all_color_words_retained": len(tensor) == 12 and len(reduced) == 10,
            "period_independence_assumed": False}


def calculate(orders, *, progress=False, max_samples=1000):
    orders = tuple(sorted({covariants.integer(order, "jet order") for order in orders}))
    if not orders:
        raise ValueError("at least one jet order is required")
    maximum = max(orders)
    count = (maximum + 4) * (maximum + 5) // 2
    if count > covariants.integer(max_samples, "max_samples", 1):
        raise RuntimeError("the complete source grid exceeds max_samples; no partial result is returned")
    samples = {order: [] for order in orders}
    nodes = [(a, b) for a in range(maximum + 4) for b in range(maximum + 4 - a)]
    for index, node in enumerate(nodes, 1):
        active = [order for order in orders if sum(node) <= order + 3]
        sample = evaluate_sources(((1, 0), (0, 1), node), active)
        for order in active:
            samples[order].append(sample)
        if progress:
            print("Complete source node %d/%d, orders %s" % (index, len(nodes), active), flush=True)
    replays = []
    for vectors in (((1, 2), (3, -1), (-2, 4)), ((2, -1), (1, 3), (4, -2)),
                    ((1, 2), (2, 4), (3, -1)), ((1, 2), (0, 0), (-2, 4)),
                    ((1, 2), (3, -1), (0, 0))):
        replays.append(evaluate_sources(vectors, orders))
        if progress:
            print("Replayed every full tensor at sources " + str(vectors), flush=True)
    return {order: certificate_from_samples(order, samples[order], replays) for order in orders}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orders", nargs="+", type=int, default=[2, 3])
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--max-samples", type=int, default=1000)
    args = parser.parse_args()
    for order, result in calculate(args.orders, progress=True, max_samples=args.max_samples).items():
        result = base._jsonable(result)
        path = HERE / ("jet_order_%d_certificate.json" % order)
        if args.write:
            path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
            print("Wrote complete jet-order %d certificate; representation cancellation: %s" %
                  (order, result["complete_representation_cancellation"]), flush=True)
        elif not path.is_file() or json.loads(path.read_text(encoding="utf-8")) != result:
            raise ArithmeticError("the complete regenerated jet-order %d certificate differs" % order)
        else:
            print("Verified the complete jet-order %d certificate." % order, flush=True)


if __name__ == "__main__":
    main()
