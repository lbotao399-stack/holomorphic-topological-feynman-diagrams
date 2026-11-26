#!/usr/bin/env python3
"""Exact all-order color-topology dictionary for a connected boundary fan.

Three independent constructions are compared:
  (1) fundamental trace splitting/joining;
  (2) oriented ribbon half-edge face cycles;
  (3) a mod-two intersection form after contracting a spanning tree.

No existing repository evaluator is imported. All arithmetic is integral.
The theorem in REPORT.md is all-order; these finite cases audit its signs,
genus, trace counts, and composition conventions.
"""
from collections import Counter
from itertools import permutations, product
from pathlib import Path
import json


HERE = Path(__file__).resolve().parent


def cyclic(word):
    word = tuple(word)
    return min(word[i:] + word[:i] for i in range(len(word))) if word else ()


def canonical(words):
    return tuple(sorted(cyclic(word) for word in words))


def set_partitions(values):
    if not values:
        yield ()
        return
    first, *rest = values
    for blocks in set_partitions(rest):
        yield ((first,), *blocks)
        for i, block in enumerate(blocks):
            yield (*blocks[:i], (first, *block), *blocks[i + 1:])


def trace_surgery(blocks, path, choices):
    """Expand one nested commutator term, then use exact completeness.

    The generators have trace norm one here. The omitted common factor is
    tau^(n+1), before converting the structure-constant chain to the trace.
    """
    inner = ("D",)
    sign = 1
    for letter, choice in reversed(tuple(zip(path, choices))):
        if choice:
            inner = inner + (letter,)
            sign = -sign
        else:
            inner = (letter,) + inner
    words = [(*block, f"S{j}") for j, block in enumerate(blocks)]
    words.append((*inner, "E"))
    empty_faces = 0
    for letter in path:
        positions = [(i, j) for i, word in enumerate(words)
                     for j, value in enumerate(word) if value == letter]
        assert len(positions) == 2
        (i, u), (j, v) = positions
        if i == j:
            word = words.pop(i)
            if u > v:
                u, v = v, u
            pieces = (word[u + 1:v], word[v + 1:] + word[:u])
        else:
            first, second = words[i], words[j]
            after_first = first[u + 1:] + first[:u]
            after_second = second[v + 1:] + second[:v]
            words = [word for k, word in enumerate(words) if k not in (i, j)]
            pieces = (after_first + after_second,)
        for word in pieces:
            if word:
                words.append(word)
            else:
                empty_faces += 1
    return sign, empty_faces, canonical(words)


def ribbon_data(blocks, path, choices):
    """Original trivalent chain and input polygons, before contractions."""
    vertices, pairs, external = [], [], {}
    active_half = {}
    for j, block in enumerate(blocks):
        cycle = []
        for letter in block:
            half = f"in:{letter}"
            active_half[letter] = half
            cycle.append(half)
        marker = f"source:{j}"
        external[marker] = f"S{j}"
        cycle.append(marker)
        vertices.append(cycle)
    for i, (letter, choice) in enumerate(zip(path, choices)):
        spoke, inner, outer = f"v{i}:spoke", f"v{i}:inner", f"v{i}:outer"
        vertices.append([inner, spoke, outer] if choice else [spoke, inner, outer])
        pairs.append((active_half[letter], spoke))
        if i:
            pairs.append((f"v{i - 1}:inner", outer))
        else:
            external[outer] = "E"
        if i == len(path) - 1:
            external[inner] = "D"
    return vertices, pairs, external


def ribbon_faces(vertices, pairs, external):
    next_half = {}
    for vertex in vertices:
        next_half.update({a: b for a, b in zip(vertex, vertex[1:] + vertex[:1])})
    across = {half: half for half in next_half}
    for first, second in pairs:
        across[first], across[second] = second, first
    permutation = {half: next_half[across[half]] for half in next_half}
    unseen, words, empty = set(permutation), [], 0
    while unseen:
        start = min(unseen)
        half, word = start, []
        while half in unseen:
            unseen.remove(half)
            if half in external:
                word.append(external[half])
            half = permutation[half]
        assert half == start
        if word:
            words.append(tuple(word))
        else:
            empty += 1
    return empty, canonical(words)


def rank_mod_two(rows):
    rows = list(rows)
    pivot = 0
    for column in range(len(rows)):
        chosen = next((i for i in range(pivot, len(rows))
                       if (rows[i] >> column) & 1), None)
        if chosen is None:
            continue
        rows[pivot], rows[chosen] = rows[chosen], rows[pivot]
        for i in range(len(rows)):
            if i != pivot and (rows[i] >> column) & 1:
                rows[i] ^= rows[pivot]
        pivot += 1
    return pivot


def intersection_genus(vertices, pairs, external, blocks, path):
    """Find genus from crossing rank, without counting any face cycles."""
    cycles = [list(vertex) for vertex in vertices]
    tree = [(f"v{i - 1}:inner", f"v{i}:outer") for i in range(1, len(path))]
    position = {letter: i for i, letter in enumerate(path)}
    tree += [(f"in:{block[0]}", f"v{position[block[0]]}:spoke") for block in blocks]
    tree_sets = {frozenset(edge) for edge in tree}
    for first, second in tree:
        located = [(i, vertex.index(half))
                   for half in (first, second)
                   for i, vertex in enumerate(cycles) if half in vertex]
        assert len(located) == 2
        (i, a), (j, b) = located
        assert i != j
        left, right = cycles[i], cycles[j]
        merged = left[a + 1:] + left[:a] + right[b + 1:] + right[:b]
        cycles = [vertex for k, vertex in enumerate(cycles) if k not in (i, j)]
        cycles.append(merged)
    assert len(cycles) == 1
    circle = [half for half in cycles[0] if half not in external]
    location = {half: i for i, half in enumerate(circle)}
    chords = [tuple(sorted((location[a], location[b])))
              for a, b in pairs if frozenset((a, b)) not in tree_sets]
    assert len(chords) == len(path) - len(blocks)
    rows = []
    for a, b in chords:
        row = 0
        for j, (c, d) in enumerate(chords):
            if (a < c < b) != (a < d < b):
                row |= 1 << j
        rows.append(row)
    rank = rank_mod_two(rows)
    assert rank % 2 == 0
    return rank // 2, len(chords)


def check_case(blocks, path):
    census = Counter()
    coefficient = Counter()
    for choices in product((0, 1), repeat=len(path)):
        sign, empty, words = trace_surgery(blocks, path, choices)
        vertices, pairs, external = ribbon_data(blocks, path, choices)
        assert ribbon_faces(vertices, pairs, external) == (empty, words)
        genus, cycle_rank = intersection_genus(vertices, pairs, external, blocks, path)
        n, r, b = len(path) - 1, len(blocks), len(words)
        assert len(vertices) == len(path) + r
        assert len(pairs) == 2 * len(path) - 1
        assert empty == n + 2 - r - b - 2 * genus
        assert 0 <= genus <= (n + 1 - r) // 2
        assert cycle_rank == n + 1 - r
        census[n, r, b, genus, empty] += 1
        # D and E are traceless physical output matrices.
        if ("D",) not in words and ("E",) not in words:
            coefficient[empty, words] += sign
    return census, {key: value for key, value in coefficient.items() if value}


def serialize_polynomial(coefficients):
    return [{"N_power": power, "traces": [list(word) for word in words],
             "integer_coefficient": value}
            for (power, words), value in sorted(coefficients.items())]


def matrix_product(left, right):
    size = len(left)
    return tuple(tuple(sum(left[i][k]*right[k][j] for k in range(size))
                       for j in range(size)) for i in range(size))


def matrix_trace(matrix):
    return sum(matrix[i][i] for i in range(len(matrix)))


def matrix_word(matrices, size):
    value = tuple(tuple(int(i == j) for j in range(size)) for i in range(size))
    for matrix in matrices:
        value = matrix_product(value, matrix)
    return value


def commutator(left, right):
    first, second = matrix_product(left, right), matrix_product(right, left)
    return tuple(tuple(a-b for a,b in zip(row1,row2)) for row1,row2 in zip(first,second))


def matrix_controls():
    """Direct finite-rank sums in the dual matrix-unit basis E_ij,E_ji."""
    records = []
    for size in (2,3):
        def specimen(seed, traceless=False):
            value = [[((i+1)*(j+2)+seed*(i-j+2)) % 7 - 3
                      for j in range(size)] for i in range(size)]
            if traceless:
                value[-1][-1] = -sum(value[i][i] for i in range(size-1))
            return tuple(map(tuple,value))
        end_d, end_e = specimen(1,True), specimen(3,True)
        units = {
            (i,j): tuple(tuple(int(a==i and b==j) for b in range(size))
                         for a in range(size))
            for i in range(size) for j in range(size)
        }
        for m in range(1,4):
            path = tuple(f"a{i}" for i in range(m))
            for blocks in set_partitions(path):
                _, coefficient = check_case(blocks,path)
                external = {"D":end_d,"E":end_e}
                external.update({f"S{j}":specimen(5+j) for j in range(len(blocks))})
                expected = 0
                for (faces,words),value in coefficient.items():
                    term = value * size**faces
                    for word in words:
                        term *= matrix_trace(matrix_word([external[x] for x in word],size))
                    expected += term
                direct = 0
                for colors in product(units,repeat=m):
                    assignments = dict(zip(path,colors))
                    scalar = 1
                    for j,block in enumerate(blocks):
                        factors = [units[assignments[letter]] for letter in block]
                        factors.append(external[f"S{j}"])
                        scalar *= matrix_trace(matrix_word(factors,size))
                        if scalar == 0:
                            break
                    if scalar == 0:
                        continue
                    nested = end_d
                    for letter in reversed(path):
                        i,j = assignments[letter]
                        nested = commutator(units[j,i],nested)
                    direct += scalar * matrix_trace(matrix_product(nested,end_e))
                assert direct == expected
                records.append({"rank":size,"boundary_order":m-1,
                                "input_blocks":[list(block) for block in blocks],
                                "exact_value":direct})
    assert any(row["exact_value"] for row in records)
    return records


def composition_checks():
    count = 0
    examples = []
    for first_genus, second_genus in product(range(4), repeat=2):
        for r1, b1, r2, b2 in product(range(1, 5), repeat=4):
            for glued in range(1, min(b1, r2) + 1):
                r, b = r1 + r2 - glued, b1 + b2 - glued
                genus = first_genus + second_genus + glued - 1
                left = (2-r1-b1-2*first_genus)+(2-r2-b2-2*second_genus)
                right = 2-r-b-2*genus
                assert left == right
                count += 1
                if (r1, b1, r2, b2, first_genus, second_genus) == (2, 2, 2, 1, 0, 0):
                    examples.append({"glued_traces": glued, "input_traces": r,
                                     "output_traces": b, "genus": genus,
                                     "N_power_after_boundary_loop_reweighting": right})
    return {"exact_gluing_checks": count, "examples": examples}


def main():
    census, term_count, case_count = Counter(), 0, 0
    cases_by_order = Counter()
    examples = {}
    for m in range(1, 6):
        labels = tuple(f"a{i}" for i in range(m))
        for blocks in set_partitions(labels):
            for path in permutations(labels):
                data, coefficients = check_case(blocks, path)
                census.update(data)
                term_count += 2 ** m
                case_count += 1
                cases_by_order[m-1] += 1
                if path == labels and blocks in ((labels,), tuple((a,) for a in labels)):
                    key = f"n={m-1},input_traces={len(blocks)}"
                    examples[key] = serialize_polynomial(coefficients)
    # Larger-order controls exercise every input trace count, both path
    # orientations, and adjacent/nonadjacent spokes without factorial growth.
    for m in (6, 7, 8):
        labels = tuple(f"a{i}" for i in range(m))
        for r in range(1, m + 1):
            blocks = tuple(tuple(labels[j::r]) for j in range(r))
            paths = (labels, labels[::-1], labels[::2]+labels[1::2])
            for path in dict.fromkeys(paths):
                data, coefficients = check_case(blocks, path)
                census.update(data)
                term_count += 2 ** m
                case_count += 1
                cases_by_order[m-1] += 1
    result = {
        "status": "PASS",
        "arithmetic": "integers and finite-field ranks; standard library only",
        "source_revision": "29ef4be3c92b8055ac51067357524c1048e6cda0",
        "independent_methods": [
            "fundamental trace splitting and joining",
            "ribbon half-edge face cycles",
            "spanning-tree contraction and intersection-form rank over F2",
        ],
        "complete_color_cases": case_count,
        "individual_commutator_terms": term_count,
        "cases_by_boundary_order": dict(sorted(cases_by_order.items())),
        "all_topologies": [
            {"boundary_order": n, "input_traces": r, "output_traces": b,
             "ribbon_genus": genus, "internal_color_faces": faces,
             "commutator_terms": count}
            for (n,r,b,genus,faces),count in sorted(census.items())
        ],
        "selected_complete_color_polynomials": examples,
        "direct_matrix_unit_controls": matrix_controls(),
        "gluing": composition_checks(),
        "exact_scope": "Direct connected fan color topology. Cochain transfer and bulk amplitudes are separate.",
    }
    destination = HERE / "ribbon_dictionary.json"
    destination.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")
    print(json.dumps({key: result[key] for key in (
        "status", "complete_color_cases", "individual_commutator_terms", "cases_by_boundary_order")}))


if __name__ == "__main__":
    main()
