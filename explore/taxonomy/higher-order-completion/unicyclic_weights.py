"""Exact six-dimensional order-two weight determination with all Wick histories.

Eight complete plane-wave functionals are retained.  Four flavor counts
certify the coefficientwise interpolation to sixteen fundamental multiplets;
an explicit matter-line component count bounds its degree by two.  Color
contractions use Python integers and all simplex integrations use FLINT Q.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction
from hashlib import sha256
import itertools
import json
from pathlib import Path
import random
import sys
import time

import numpy as np
from flint import fmpq, fmpq_mat

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
ENGINE = ROOT / "explore/taxonomy/qn-holomorphic"
sys.path[:0] = [str(ENGINE), str(ROOT / "code")]
import theories
from wick import Algebra
from explicit import graph_ops
from fmode import Context, FOp, apply_types
import exact_colour
import exact_holomorphic_kernel as exact_kernel

OPERATOR_TYPES = ("gv", "gv", "u")
FUNCTIONAL_TYPES = ("g", "g", "g", "g", "u")
JETS = {("O", 0): (0, 1), ("O", 1): (1, 2), ("O", 2): (2,)}
SEED, FUNCTIONAL_COUNT, POOL_SIZE, MOMENTUM_CLASSES = 117, 8, 3, 5


def atomic_json(path, data):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, sort_keys=True, indent=2) + "\n")
    temporary.replace(path)


def blocks_for(sl):
    cover = {}
    for block in sorted(sl.blocks, key=lambda b: (len(sl.blocks[b]), b)):
        indices = sl.blocks[block]
        if len({int(sl.pi[i]) for i in indices}) != 1 or any(i in cover for i in indices):
            continue
        cover.update({i: block for i in indices})
    out = {}
    for i in range(sl.N):
        out.setdefault(cover.get(i, "par%d" % int(sl.pi[i])), []).append(i)
    return out


class Functionals:
    def __init__(self, ctx, rng, allowed):
        self.types = FUNCTIONAL_TYPES
        self.N, self.F, self.R = 5, FUNCTIONAL_COUNT, POOL_SIZE
        self.pool = {}
        for kind in sorted(set(self.types)):
            values = np.zeros((self.R, ctx.alg.N), dtype=object)
            for r in range(self.R):
                for a in allowed(kind):
                    values[r, a] = rng.randint(-3, 3)
            self.pool[kind] = values
        self.ridx = np.zeros((self.F, self.N), dtype=np.int64)
        self.u = np.zeros((self.F, self.N, 3), dtype=np.int64)
        for f in range(self.F):
            momenta = [[rng.randint(-5, 5) for _ in range(3)] for _ in range(MOMENTUM_CLASSES)]
            for kind in sorted(set(self.types)):
                slots = [i for i, label in enumerate(self.types) if label == kind]
                selected = rng.sample([(r, m) for r in range(self.R) for m in range(MOMENTUM_CLASSES)], len(slots))
                for i, (r, m) in zip(slots, selected):
                    self.ridx[f, i], self.u[f, i] = r, momenta[m]
        self.perms = np.array(list(itertools.permutations(range(self.N))), dtype=np.int64)
        patterns = defaultdict(list)
        signs = []
        for row, permutation in enumerate(self.perms):
            patterns[tuple(self.types[i] for i in permutation)].append(row)
            odd = [int(i) for i in permutation if ctx.bpar[self.types[i]]]
            signs.append((-1) ** sum(a > b for j, a in enumerate(odd) for b in odd[j + 1:]))
        self.pat_of = {key: np.array(value, dtype=np.int64) for key, value in patterns.items()}
        self.sg = np.array(signs, dtype=object)
        self.prow, rows = {}, []
        for kind in sorted(set(self.types)):
            self.prow[kind] = np.arange(len(rows), len(rows) + self.R, dtype=np.int64)
            rows.extend(self.pool[kind])
        self.pseudo = type("ExactFunctionalPool", (), {})()
        self.pseudo.F = 1
        self.pseudo.Af = np.array(rows, dtype=object)[None, :, :]


def inputs(flavour):
    H, sl = theories.build("6d10", "sl2", "fund*%d" % flavour)
    alg = Algebra(sl, H)
    ctx = Context(alg, blocks_for(sl))
    ctx.certify_matter_flavour_cycles = True

    def allowed(kind):
        return [ctx.blocks[kind][i] for i in (0, 1, 2 * flavour, 2 * flavour + 1)] if kind == "u" else ctx.blocks[kind]

    rng = random.Random(SEED)
    operator = tuple((("O", i), kind) for i, kind in enumerate(OPERATOR_TYPES))
    colors = {}
    for tag, kind in operator:
        vector = np.zeros(alg.N, dtype=object)
        for a in allowed(kind):
            vector[a] = rng.randint(-2, 2)
        colors[tag] = vector
    functionals = Functionals(ctx, rng, allowed)
    graphs = {graph.name: graph for n in (1, 2)
              for graph in graph_ops(alg, n, unicyclic=True, drop_zero=False)}
    metadata = {"flavours": flavour, "basis_names": list(sl.names), "basis_parities": list(map(int, sl.pi)),
                "blocks": ctx.blocks, "operator_colors": [list(map(int, colors[("O", i)])) for i in range(3)],
                "functional_color_pools": {key: [[int(x) for x in row] for row in value] for key, value in functionals.pool.items()},
                "functional_color_indices": functionals.ridx.tolist(), "functional_momenta": functionals.u.tolist()}
    return ctx, operator, colors, functionals, graphs, metadata


def color_coefficients(bundle, outer_name, inner_name):
    ctx, operator, colors, functionals, graphs, _ = bundle
    op1 = FOp(ctx.alg, None if inner_name == "Q0" else graphs[inner_name], "i")
    op2 = FOp(ctx.alg, None if outer_name == "Q0" else graphs[outer_name], "o")
    original_tags = [tag for tag, _ in operator]
    aggregate, information = {}, {}
    histories = 0
    for binding1 in itertools.permutations(original_tags, op1.nin):
        for sign1, word1, blocks1, q1 in apply_types(op1, ctx, operator, binding1):
            for binding2 in itertools.permutations([tag for tag, _ in word1], op2.nin):
                for sign2, word2, blocks2, q2 in apply_types(op2, ctx, word1, binding2):
                    slots = [tag for tag, _ in word2]
                    pattern = tuple(block for _, block in word2)
                    if pattern not in functionals.pat_of:
                        continue
                    histories += 1
                    types = dict(operator)
                    types.update(zip(op1.outtags(), blocks1))
                    types.update(zip(op2.outtags(), blocks2))
                    locations = []
                    for tag in original_tags:
                        if tag in binding1:
                            locations.append((tag, ("i", binding1.index(tag))))
                        elif tag in binding2:
                            locations.append((tag, ("o", binding2.index(tag))))
                        else:
                            locations.append((tag, ("s", slots.index(tag))))
                    locations = tuple(locations)
                    smap = tuple(tag[1] if tag[0] == op1.prefix else None for tag in binding2)
                    key = (outer_name, inner_name, smap, locations)
                    information[key] = {"ops": (op1, op2), "smap": smap, "jets": JETS, "jetloc": locations}
                    value = exact_colour.colour_values(
                        ctx, (op1, op2), (binding1, binding2), types, colors, slots,
                        [functionals.prow[pattern[s]] for s in range(5)], functionals.pseudo, (q1, q2))[0]
                    if not value.any():
                        continue
                    rows = functionals.pat_of[pattern]
                    permutations = functionals.perms[rows]
                    indices = tuple(functionals.ridx[:, permutations[:, s]] for s in range(5))
                    coefficients = value[indices] * sign1 * sign2 * functionals.sg[rows][None, :]
                    accumulation = aggregate.setdefault(key, {}).setdefault(pattern, np.zeros(coefficients.shape, dtype=object))
                    accumulation += coefficients
    grouped = {}
    for key, patterns in aggregate.items():
        result = defaultdict(int)
        for pattern, coefficients in patterns.items():
            rows = functionals.pat_of[pattern]
            for f in range(FUNCTIONAL_COUNT):
                for column, row in enumerate(rows):
                    coefficient = int(coefficients[f, column])
                    if coefficient:
                        momenta = tuple(tuple(map(int, functionals.u[f, functionals.perms[row, slot]])) for slot in range(5))
                        result[(f, momenta)] += coefficient
        grouped[key] = {key: value for key, value in result.items() if value}
    norm = Fraction(1, ctx.alg.Gden ** (op1.nG + op2.nG))
    if op1.gop is None:
        norm /= 2
    if op2.gop is None:
        norm /= 2
    return grouped, information, norm, histories


def combine_flavours(data):
    keys = set().union(*(set(entry) for entry in data))
    combined = {}
    checked = 0
    for key in keys:
        rows = set().union(*(set(entry.get(key, {})) for entry in data))
        result = {}
        for row in rows:
            values = [entry.get(key, {}).get(row, 0) for entry in data]
            assert values[3] == values[0] - 3 * values[1] + 3 * values[2], (key, row, values)
            checked += 1
            coefficient = 91 * values[0] - 195 * values[1] + 105 * values[2]
            if coefficient:
                result[row] = coefficient
        combined[key] = result
    return combined, checked


def dependencies():
    paths = [Path(__file__), HERE / "exact_colour.py", HERE / "exact_holomorphic_kernel.py"]
    # Bind the actual import closure.  Unrelated experimental bulk workers
    # can evolve without invalidating this purely holomorphic calculation.
    paths += [ENGINE / name for name in ("theories.py", "superalg.py", "wick.py", "polyops.py",
                                        "explicit.py", "fmode.py", "kern_h.py", "rgraphs.py")]
    paths += [ROOT / "code/htk" / name for name in
              ("__init__.py", "graph.py", "linalg.py", "brute.py", "simplex.py", "spred.py",
               "oracle.py", "catalog.py", "conventions.py", "forest_model.py", "fan_moments.py")]
    return {str(path.relative_to(ROOT)): sha256(path.read_bytes()).hexdigest() for path in paths}


def compute(output, only=None, benchmark=False):
    started = time.monotonic()
    bundles = [inputs(flavour) for flavour in (1, 2, 3, 4)]
    print("four exact color algebras and complete inputs prepared", round(time.monotonic() - started, 2), flush=True)
    graph_names = ["G2.1", "G2.2", "G2.3u", "G2.4u", "G2.5u"]
    pairs = [("Q0", graph) for graph in graph_names] + [("G1.1", "G1.1")] + [(graph, "Q0") for graph in graph_names]
    if only:
        pairs = [pair for pair in pairs if "|".join(pair) in only]
    certificate = {"schema": "six-dimensional-complete-order-two-weight-equations-v1",
                   "theory": "6d (1,0), sl2 with sixteen fundamental cotangent hypermultiplets",
                   "operator_types": OPERATOR_TYPES, "functional_types": FUNCTIONAL_TYPES,
                   "operator_jets": [list(JETS[("O", i)]) for i in range(3)],
                   "input_data": [bundle[-1] for bundle in bundles], "seed": SEED,
                   "graphs": {name: {"edges": graph.edges, "bulk_vertices": graph.nb,
                                      "root_spokes": graph.spokes, "output_slots": graph.outv,
                                      "automorphism_order": graph.aut, "kernel_orientation": graph.gk.sign_n}
                              for name, graph in bundles[0][4].items()},
                   "implementation_sha256": dependencies(), "compositions": {}, "complete": False}
    atomic_json(output, certificate)
    for outer, inner in pairs:
        stage = time.monotonic()
        color_data, information, norms, histories = [], [], [], []
        for bundle in bundles:
            colors, info, norm, count = color_coefficients(bundle, outer, inner)
            color_data.append(colors)
            information.append(info)
            norms.append(norm)
            histories.append(count)
        assert len(set(norms)) == len(set(histories)) == 1
        combined, checks = combine_flavours(color_data)
        print(outer, inner, "integer colors and flavour checks", checks,
              "seconds", round(time.monotonic() - stage, 2), flush=True)
        results = [fmpq(0) for _ in range(FUNCTIONAL_COUNT)]
        scalar_count, nonzero_count = 0, 0
        for key in sorted(combined, key=repr):
            info = information[0][key]
            for (f, momenta), coefficient in sorted(combined[key].items()):
                value = exact_kernel.key_kernel(info, momenta)
                results[f] += coefficient * value
                scalar_count += 1
                nonzero_count += int(bool(value))
                if scalar_count % 500 == 0:
                    print(outer, inner, "exact integrals", scalar_count,
                          "seconds", round(time.monotonic() - stage, 2), flush=True)
                if benchmark and scalar_count == 500:
                    print("BENCHMARK_ONLY", outer, inner, "500 exact integrals",
                          round(time.monotonic() - stage, 2), flush=True)
                    return
        norm = fmpq(norms[0].numerator, norms[0].denominator)
        results = [value * norm for value in results]
        certificate["compositions"][outer + "|" + inner] = {
            "complete_rational_values": list(map(str, results)), "all_Wick_histories": histories[0],
            "coefficientwise_fourth_flavour_checks": checks, "complete_kernel_integrals": scalar_count,
            "nonzero_kernel_integrals": nonzero_count, "normalization": str(norm)}
        certificate["exact_color_statistics"] = dict(exact_colour.STATISTICS)
        atomic_json(output, certificate)
        print(outer, inner, "values", list(map(str, results)),
              "seconds", round(time.monotonic() - stage, 2), flush=True)
    certificate["complete"] = len(certificate["compositions"]) == 11
    certificate["seconds"] = round(time.monotonic() - started, 2)
    atomic_json(output, certificate)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE / "unicyclic_weights.json")
    parser.add_argument("--only", default="")
    parser.add_argument("--benchmark", action="store_true")
    arguments = parser.parse_args()
    compute(arguments.output, set(arguments.only.split(",")) if arguments.only else None, arguments.benchmark)
