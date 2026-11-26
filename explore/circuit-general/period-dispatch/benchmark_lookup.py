"""Measure exact K5 orbit lookup, with table loading reported separately."""
import argparse
from collections import Counter
import json
from math import factorial, prod
import random
from time import perf_counter

import k5_bridge as bridge


def benchmark(count):
    start = perf_counter()
    lookup = bridge.K5Lookup(bridge.E5, 8)
    load_seconds = perf_counter() - start
    rows = json.loads((bridge.TABLE_DIRECTORY / "K5_q8.json").read_text())["values"]
    rng = random.Random(277)
    # Distinct relabelings test the canonicalizer, not repeated cached keys.
    queries = []
    for index in range(count):
        powers = rows[index % len(rows)][0]
        permutation = rng.choice(bridge._PERMUTATIONS)
        queries.append(tuple(powers[e] for e in permutation))
    bridge.signature_key.cache_clear()
    start = perf_counter()
    grouped = Counter()
    for powers in queries:
        grouped[lookup.key(powers)] += 1
    lookup_seconds = perf_counter() - start
    cache = bridge.signature_key.cache_info()._asdict()
    workload = Counter()
    for powers, _, _ in rows:
        signatures = [tuple(sorted(powers[e] for e in incident))
                      for incident in bridge._INCIDENT]
        workload[prod(factorial(size) for size in Counter(signatures).values())] += 1
    # A cold reference timing exposes the cost of the original 120 permutations.
    reference_count = min(10000, count)
    bridge.canonical_powers.cache_clear()
    start = perf_counter()
    for powers in queries[:reference_count]:
        bridge.canonical_powers(powers)
    reference_seconds = perf_counter() - start
    return {
        "q": 8,
        "table_orbits": len(lookup.values),
        "queries": count,
        "distinct_input_vectors": len(set(queries)),
        "distinct_result_orbits": len(grouped),
        "load_and_index_seconds": load_seconds,
        "validated_lookup_and_integer_grouping_seconds": lookup_seconds,
        "lookups_per_second": count / lookup_seconds,
        "signature_cache": cache,
        "candidate_permutations_to_orbits": dict(sorted(workload.items())),
        "full_120_reference_queries": reference_count,
        "full_120_reference_seconds": reference_seconds,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--count", type=int, default=100000)
    args = parser.parse_args()
    if args.count < 1:
        parser.error("--count must be positive")
    print(json.dumps(benchmark(args.count), indent=2))
