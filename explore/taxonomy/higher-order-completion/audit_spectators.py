"""Concrete regression for an untouched input color in a real Wick history."""
from hashlib import sha256
import itertools
from pathlib import Path

import numpy as np
import unicyclic_weights as weights
import exact_colour
from fmode import FOp, apply_types, colour_values as inherited_colour

HERE = Path(__file__).resolve().parent


def old_call(arguments, parities):
    arguments = list(arguments)
    arguments[4] = {key: np.array(value, dtype=np.float64) for key, value in arguments[4].items()}
    functionals = type("InheritedFloatPool", (), {})()
    functionals.F = arguments[7].F
    functionals.Af = np.array(arguments[7].Af, dtype=np.float64)
    arguments[7] = functionals
    return inherited_colour(*arguments, qs=parities)


def history(bundle, graph_name, require_spectator):
    ctx, operator, colors, functionals, graphs, _ = bundle
    inner = FOp(ctx.alg, graphs[graph_name], "i")
    outer = FOp(ctx.alg, None, "o")
    tags = [tag for tag, _ in operator]
    for binding1 in itertools.permutations(tags, inner.nin):
        for sign1, word1, blocks1, q1 in apply_types(inner, ctx, operator, binding1):
            for binding2 in itertools.permutations([tag for tag, _ in word1], outer.nin):
                for sign2, word2, blocks2, q2 in apply_types(outer, ctx, word1, binding2):
                    pattern = tuple(block for _, block in word2)
                    if pattern not in functionals.pat_of:
                        continue
                    spectators = set(tags) - set(binding1) - set(binding2)
                    if bool(spectators) != require_spectator:
                        continue
                    slots = [tag for tag, _ in word2]
                    types = dict(operator)
                    types.update(zip(inner.outtags(), blocks1))
                    types.update(zip(outer.outtags(), blocks2))
                    slot_functions = [functionals.prow[kind] for kind in pattern]
                    args = (ctx, (inner, outer), (binding1, binding2), types, colors,
                            slots, slot_functions, functionals.pseudo)
                    old = old_call(args, (q1, q2))
                    fresh = exact_colour.colour_values(*args, qs=(q1, q2))
                    if not fresh.any():
                        continue
                    if require_spectator and np.array_equal(old, fresh):
                        continue
                    result = {"inner_graph": graph_name, "outer_graph": "Q0",
                              "bindings": [binding1, binding2], "remaining_slots": slots,
                              "output_pattern": pattern, "spectator_tags": sorted(spectators),
                              "loop_parities": [q1, q2], "Wick_sign": sign1 * sign2,
                              "shape": list(fresh.shape),
                              "inherited_complete_color_tensor": list(map(int, old.flat)),
                              "corrected_complete_color_tensor": list(map(int, fresh.flat))}
                    if require_spectator:
                        tag = sorted(spectators)[0]
                        zero_colors = dict(colors)
                        zero_colors[tag] = np.zeros(ctx.alg.N, dtype=object)
                        altered = list(args)
                        altered[4] = zero_colors
                        wrong_zero = old_call(altered, (q1, q2))
                        correct_zero = exact_colour.colour_values(*altered, qs=(q1, q2))
                        assert np.array_equal(wrong_zero, old) and wrong_zero.any()
                        assert not correct_zero.any()
                        twice_colors = dict(colors)
                        twice_colors[tag] = 2 * colors[tag]
                        altered[4] = twice_colors
                        assert np.array_equal(exact_colour.colour_values(*altered, qs=(q1, q2)), 2 * fresh)
                        result["zero_original_spectator_gives_zero_corrected_tensor"] = True
                        result["inherited_tensor_wrongly_independent_of_original_spectator"] = True
                        result["exact_linearity_under_doubling_spectator"] = True
                    else:
                        assert np.array_equal(old, fresh)
                        result["full_tensor_unchanged_without_original_spectator"] = True
                    return result
    raise AssertionError("required nonzero regression history not found")


def main():
    bundle = weights.inputs(1)
    report = {"schema": "original-spectator-color-regression-v1", "result": "passed",
              "complete_input": bundle[-1], "with_spectator": history(bundle, "G2.3u", True),
              "without_spectator": history(bundle, "G2.1", False),
              "implementation_sha256": weights.dependencies(),
              "verifier_sha256": sha256(Path(__file__).read_bytes()).hexdigest()}
    weights.atomic_json(HERE / "spectator_audit.json", report)
    print("PASS: original spectator must remain multilinear; complete non-spectator tensor unchanged")


if __name__ == "__main__":
    main()
