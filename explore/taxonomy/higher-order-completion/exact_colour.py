"""Integer color contractions, including every untouched original letter.

The combinatorial block assignments and ordered Wick signs are supplied by
the original operator engine.  Every tensor product and summation here uses
Python integers.  Original input color vectors are inserted after all final
slots have been named, so spectator letters retain their original colors.
"""
import itertools
import numpy as np
from fmode import leg_assignments

PATHS = {}
STATISTICS = {"histories_with_spectators": 0, "integer_contractions": 0,
              "maximum_closed_matter_flavour_cycles": 0}


def integers(array):
    array = np.asarray(array)
    out = np.array([int(x) for x in array.flat], dtype=object).reshape(array.shape)
    assert (out == array).all(), "color input was not an exact integer"
    return out


def colour_values(ctx, ops, bindings, wtypes, O_cols, slots, slot_funcs, func,
                  qs=(None, None)):
    bd = ctx.bd
    assignments = []
    for k, (op, binding) in enumerate(zip(ops, bindings)):
        external = {("out", i): wtypes[(op.prefix, i)]
                    for i in range(len(op.outv))}
        external.update({("spoke", j): wtypes[binding[j]] for j in range(op.nin)})
        assignments.append(leg_assignments(bd, op, k, external, qs[k]))
    names = {}

    def name(key):
        return names.setdefault(key, len(names))

    structure = []
    for k, (op, binding) in enumerate(zip(ops, bindings)):
        spokes = {e: j for j, e in enumerate(op.spokes)}
        for v in range(op.nb):
            indices, blocks = [], []
            for leg in op.ga.legs[v]:
                if leg[0] == "out":
                    tag = (op.prefix, op.outv.index((v, leg[2])))
                    indices.append(name(("col", tag)))
                    blocks.append(("fixed", wtypes[tag]))
                else:
                    key = ("leg", k, leg[1], v)
                    indices.append(name(key))
                    blocks.append(("var", key))
            structure.append(("V", blocks, indices))
        for e in op.internal:
            tail, head = op.edges[e]
            keys = [("leg", k, e, v) for v in (tail, head)]
            structure.append(("G", [("var", key) for key in keys],
                              [name(key) for key in keys]))
        for e in op.spokes:
            tail, _ = op.edges[e]
            tag = binding[spokes[e]]
            key = ("leg", k, e, tail)
            structure.append(("G", [("var", key), ("fixed", wtypes[tag])],
                              [name(key), name(("col", tag))]))
    bound = set(bindings[0]) | set(bindings[1])
    if any(tag not in bound for tag in O_cols):
        STATISTICS["histories_with_spectators"] += 1
    final_slots, output = [], []
    for slot, tag in enumerate(slots):
        functional = name(("fun", slot))
        final_slots.append((slot, wtypes[tag], [functional, name(("col", tag))]))
        output.append(functional)
    # At this point *every* original letter, including spectators, is named.
    original = []
    for tag, vector in O_cols.items():
        assert ("col", tag) in names
        original.append((integers(vector[bd.ids[wtypes[tag]]]),
                         [name(("col", tag))]))
    assert len(names) <= 52
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
    subscripts = [entry[2] for entry in structure]
    subscripts += [indices for _, indices in original]
    subscripts += [indices for _, _, indices in final_slots]
    specification = ",".join("".join(alphabet[i] for i in axes) for axes in subscripts)
    specification += "->" + "".join(alphabet[i] for i in output)
    total = None
    for combination in itertools.product(*assignments):
        assignment = {key: value for part in combination for key, value in part.items()}
        arrays = []
        matter_edges = []
        for kind, blocks, _ in structure:
            labels = tuple(value if mode == "fixed" else assignment[value]
                           for mode, value in blocks)
            arrays.append(integers(bd.Vs[labels] if kind == "V" else bd.Gs[labels]))
        if getattr(ctx, "certify_matter_flavour_cycles", False):
            for (_, blocks, indices) in structure:
                labels = tuple(value if mode == "fixed" else assignment[value]
                               for mode, value in blocks)
                matter = [axis for axis, label in zip(indices, labels) if label == "u"]
                assert len(matter) in (0, 2)
                if matter:
                    matter_edges.append(tuple(matter))
            boundary = {name(("col", tag)) for tag in O_cols if wtypes[tag] == "u"}
            boundary |= {name(("col", tag)) for tag in slots if wtypes[tag] == "u"}
            adjacency = {}
            for a, b in matter_edges:
                adjacency.setdefault(a, set()).add(b)
                adjacency.setdefault(b, set()).add(a)
            seen, closed = set(), 0
            for start in adjacency:
                if start in seen:
                    continue
                component, stack = set(), [start]
                while stack:
                    vertex = stack.pop()
                    if vertex in component:
                        continue
                    component.add(vertex)
                    stack.extend(adjacency[vertex] - component)
                seen.update(component)
                if not component & boundary:
                    closed += 1
            assert closed <= 2, "flavour interpolation degree bound exceeded"
            STATISTICS["maximum_closed_matter_flavour_cycles"] = max(
                STATISTICS["maximum_closed_matter_flavour_cycles"], closed)
        arrays.extend(vector for vector, _ in original)
        values = []
        for f in range(func.F):
            external = [integers(func.Af[f][slot_funcs[slot], :][:, bd.ids[block]])
                        for slot, block, _ in final_slots]
            operands = arrays + external
            key = (specification, tuple(value.shape for value in operands))
            if key not in PATHS:
                PATHS[key] = np.einsum_path(specification, *operands, optimize="greedy")[0]
            values.append(np.einsum(specification, *operands, optimize=PATHS[key]))
            STATISTICS["integer_contractions"] += 1
        result = np.array(values, dtype=object)
        total = result if total is None else total + result
    if total is None:
        return np.zeros([func.F] + [len(x) for x in slot_funcs], dtype=object)
    return total
