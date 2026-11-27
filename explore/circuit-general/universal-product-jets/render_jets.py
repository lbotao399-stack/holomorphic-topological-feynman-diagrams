"""Render and exactly verify every full-index higher-jet tensor coefficient."""
from collections import defaultdict
from fractions import Fraction
from functools import reduce
from math import lcm
import argparse
import json

import jet_tensor as jets
import render_tensor as original_renderer

SOURCES = (r"\gamma_4", "c", r"\gamma_5")
INDICES = (r"\dot a", r"\dot b", r"\dot c")


def nested_tuple(value):
    return tuple(nested_tuple(item) for item in value) if isinstance(value, list) else value


def free_tensor(counts):
    labels = [source for source, count in zip(SOURCES, counts) for unused in range(count)]
    if len(labels) not in (2, 3):
        raise ValueError("the printed complete tables have two or three free indices")
    factors = []
    for position, (source, index) in enumerate(zip(labels, INDICES)):
        if position == 0:
            index = "(" + index
        if position == len(labels) - 1:
            index += ")"
        factors.append(r"\lambda_{%s,+%s}" % (source, index))
    return "".join(factors)


def coefficient_tex(terms):
    grouped = defaultdict(lambda: defaultdict(Fraction))
    for counts, powers, coefficient in terms:
        grouped[tuple(counts)][tuple(powers)] += Fraction(coefficient)
    output = []
    for counts, coefficients in sorted(grouped.items(), reverse=True):
        coefficients = {powers: value for powers, value in coefficients.items() if value}
        if not coefficients:
            continue
        denominator = reduce(lcm, (value.denominator for value in coefficients.values()), 1)
        integers = {powers: int(value * denominator) for powers, value in coefficients.items()}
        if any(Fraction(integers[powers], denominator) != value for powers, value in coefficients.items()):
            raise ArithmeticError("rendering changed an exact tensor coefficient")
        factor = free_tensor(counts)
        if denominator != 1:
            factor = r"\frac{%s}{%d}" % (factor, denominator)
        output.append(factor + r"\left[" + original_renderer.integer_polynomial(integers) + r"\right]")
    return "&" + "\\\\\n&+".join(output) if output else "&0"


def render(order):
    path = jets.HERE / ("jet_order_%d_certificate.json" % order)
    data = json.loads(path.read_text())
    if order not in (2, 3) or data["jet_order"] != order:
        raise ValueError("a complete second- or third-jet certificate is required")
    derivative = "".join(r"\partial_{+%s}" % index for index in INDICES[:order])
    lines = [
        "# Complete order-%d input-jet tensor" % order, "",
        "Generated and checked coefficient by coefficient by [render_jets.py](render_jets.py).", "",
        r"The polarized superpotential is \(\bar W(\gamma)=\bar W^{i_1i_2i_3i_4}\gamma_{1i_1}\gamma_{2i_2}\gamma_{3i_3}\gamma_{4i_4}\). The inputs are", "",
        r"\[\beta^{1i_1}(%s\beta^{2i_2})\gamma_5,\qquad b_A\beta^{3i_3}.\]" % derivative, "",
        r"The output fields are \(\gamma_{4i_4},c^B,\gamma_5\), with every input/output color index and every displayed spinor index retained. The spectator \(\gamma_5\) is free, gauge neutral, and absent from \(\bar W\).", "",
        r"Sources carry the manuscript's lower indices. Before using angle brackets, fix their convention:", "",
        r"\[\langle\lambda_{\gamma_4},\lambda_c\rangle=\epsilon^{\dot b\dot c}\lambda_{\gamma_4,+\dot b}\lambda_{c,+\dot c},\qquad \epsilon^{\dot+\dot-}=1=-\epsilon^{\dot-\dot+}.\]", "",
        "The same antisymmetric pairing defines the other brackets. Parentheses on the free spinor indices mean the average over every permutation of those indices. These are the actual symmetric derivative indices; all their values are retained. For example,", "",
        r"\[\lambda_{\gamma_4,+(\dot a}\lambda_{c,+\dot b)}=\frac{\lambda_{\gamma_4,+\dot a}\lambda_{c,+\dot b}+\lambda_{\gamma_4,+\dot b}\lambda_{c,+\dot a}}{2}.\]", "",
        "For three indices the average contains all six index permutations. Each color tensor below multiplies its entire displayed spinor tensor. Replace each source spinor by the derivative acting on its own surviving field, and multiply their complete sum by", "",
        r"\[\frac{g^2}{128(8\pi^2)^3}.\]", "",
        r"This is the complete \(g^2\bar W\) coefficient at \(\hbar^3\). The factors from contracted input derivatives and from \(\lambda=-\partial/2\) have already been included in this common multiplier. See [DERIVATION.md](DERIVATION.md).", "",
        "## Complete rational tensor after the representation identity", "",
    ]
    if len(data["complete_reduced_tensor"]) != 10:
        raise ArithmeticError("a complete open-color word is missing")
    for number, row in enumerate(data["complete_reduced_tensor"], 1):
        word = nested_tuple(row["generator_word"])
        channels = row["complete_covariant_tensor"]
        if set(channels) - {"1"}:
            raise ArithmeticError("the complete reduced tensor is not rational")
        lines += ["### Color tensor %d" % number, "", r"\[" + original_renderer.tensor_tex(word) + r"\]", "",
                  "Full coefficient:", "", "\\[\n\\begin{aligned}\n" + coefficient_tex(channels.get("1", ())) +
                  "\n\\end{aligned}\n\\]", ""]
    lines += ["## Complete coefficients of the representation defects", "",
              r"The complete unreduced twelve-word tensor is retained in the certificate. Its entire \(\zeta(3)\) contribution consists of the following two full spinor tensors multiplying the displayed exact representation defects. No independence of period values is used.", ""]
    defects = data["complete_representation_defect_coefficients"]
    if len(defects) != 2 or any(row["period_label"] != "z3" for row in defects):
        raise ArithmeticError("the complete period-defect census changed")
    for row in defects:
        slot = row["W_slot"]
        indices = ["i_1", "i_2", "i_3", "i_4"]
        indices[slot - 1] = "j_%d" % slot
        factor = (r"\zeta(3)\bigl([t_A,t_B]-iC^C{}_{AB}t_C\bigr)^{i_%d}{}_{j_%d}\bar W^{%s}"
                  % (slot, slot, " ".join(indices)))
        lines += ["### Superpotential slot %d" % slot, "", r"\[" + factor + r"\]", "", "Full coefficient:", "",
                  "\\[\n\\begin{aligned}\n" + coefficient_tex(row["complete_covariant_tensor_coefficient"]) +
                  "\n\\end{aligned}\n\\]", ""]
    lines += ["## Exact completeness", "",
              "Every rational and period coefficient before and after the representation relation, every original geometric tensor, every source interpolation value, and all further full-source replays are retained in [%s](%s)." % (path.name, path.name), "",
              "The all-finite-jet rational boundary evaluator independently reproduces every entry of this tensor at every proven-degree source interpolation node. Its proof and the strict scope are in [DERIVATION.md](DERIVATION.md) and [the source-simplex derivation](../forest-jet-identities/DERIVATION.md).", ""]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    for order, filename in ((2, "FULL_SECOND_JET.md"), (3, "FULL_THIRD_JET.md")):
        result = render(order)
        path = jets.HERE / filename
        if args.write:
            path.write_text(result, encoding="utf-8")
        elif not path.is_file() or path.read_text(encoding="utf-8") != result:
            raise ArithmeticError("the complete printed order-%d tensor differs from its exact certificate" % order)
        print("Verified every full-index order-%d coefficient: ten rational tensors and both complete period defects." % order)


if __name__ == "__main__":
    main()
