"""Render and verify the complete covariant tensor table without auxiliary spinors."""
from collections import defaultdict
from fractions import Fraction
from functools import reduce
from math import lcm
import argparse

import quartic_tensor as universal

VECTORS = (r"\lambda_{\gamma_4}", r"\lambda_c", r"\lambda_{\gamma_5}")
FREE_INDEX_VECTORS = (r"\lambda_{\gamma_4,+\dot a}", r"\lambda_{c,+\dot a}",
                      r"\lambda_{\gamma_5,+\dot a}")
BRACKETS = (r"\langle\lambda_{\gamma_4},\lambda_c\rangle",
            r"\langle\lambda_{\gamma_5},\lambda_c\rangle",
            r"\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle")


def tensor_tex(word):
    indices = ["i_1", "i_2", "i_3", "i_4"]
    factors = []
    actions = word[1] if word[0] == "MM" else ((word[1], ("C",)),)
    if word[0] == "fR":
        factors.append(r"iC^C{}_{AB}")
    for flavor, letters in actions:
        chain = ["i_%d" % flavor, "j_%d" % flavor, "k_%d" % flavor]
        for step, letter in enumerate(letters):
            factors.append(r"(t_%s)^{%s}{}_{%s}" % (letter, chain[step], chain[step + 1]))
        indices[flavor - 1] = chain[len(letters)]
    return "".join(factors) + r"\bar W^{%s}" % " ".join(indices)


def monomial_tex(powers):
    return "".join(bracket + ("^{%d}" % power if power != 1 else "")
                   for bracket, power in zip(BRACKETS, powers) if power) or "1"


def integer_polynomial(terms):
    output = []
    for powers, coefficient in sorted(terms.items(), reverse=True):
        if not coefficient:
            continue
        magnitude = abs(coefficient)
        scalar = str(magnitude) if magnitude != 1 else ""
        monomial = monomial_tex(powers)
        if monomial == "1":
            scalar, monomial = str(magnitude), ""
        output.append(("-" if coefficient < 0 else "+") + scalar + monomial)
    if not output:
        return "0"
    output[0] = output[0].lstrip("+")
    lines = ["".join(output[start:start + 3]) for start in range(0, len(output), 3)]
    return (r"\begin{aligned}" + "\n&" + ("\\\\\n&".join(lines))
            + "\n" + r"\end{aligned}")


def spinor_tex(terms):
    grouped = defaultdict(lambda: defaultdict(Fraction))
    for vector, powers, coefficient in terms:
        grouped[vector][powers] += coefficient
    result = []
    for vector, coefficients in sorted(grouped.items()):
        denominator = reduce(lcm, (coefficient.denominator for coefficient in coefficients.values()), 1)
        integers = {powers: int(coefficient * denominator) for powers, coefficient in coefficients.items()}
        numerator = integer_polynomial(integers)
        factor = FREE_INDEX_VECTORS[vector]
        if denominator != 1:
            factor = r"\frac{%s}{%d}" % (factor, denominator)
        result.append(factor + r"\left[" + numerator + r"\right]")
    return ("&" + "\\\\\n&+".join(result)) if result else "&0"


def compact_spinor(coefficients):
    output = []
    for vector, coefficient in enumerate(coefficients):
        if not coefficient:
            continue
        if coefficient.denominator != 1:
            raise ArithmeticError("the advertised compact common denominator is insufficient")
        value = coefficient.numerator
        output.append(("-" if value < 0 else "+") + (str(abs(value)) if abs(value) != 1 else "")
                      + FREE_INDEX_VECTORS[vector])
    return "".join(output).lstrip("+") or "0"


def render():
    tensor = universal.covariant_tensor(universal.read_geometry())
    reduced = universal.apply_representation_relation(tensor)
    lines = [
        "# Complete open-color, full-source three-loop tensor", "",
        "Generated exactly by `render_tensor.py` from the full tensor certificate's computation.", "",
        "The coefficient tensor of the polarized superpotential is defined by", "",
        r"\[\bar W(\gamma)=\bar W^{i_1i_2i_3i_4}\gamma_{1i_1}\gamma_{2i_2}\gamma_{3i_3}\gamma_{4i_4}.\]", "",
        r"The inputs are \(\beta^{1i_1}(\partial_{+\dot a}\beta^{2i_2})\gamma_5\) and \(b_A\beta^{3i_3}\). All six input/output color indices and the free spinor index remain open. The field \(\gamma_5\) is a separate free gauge singlet, absent from \(\bar W\).", "",
        r"Sources and currents have the manuscript's lower index \(+\dot a\). The angle bracket is its antisymmetric spinor pairing, with", "",
        r"\[\langle\lambda_{\gamma_4},\lambda_c\rangle=\epsilon^{\dot b\dot c}\lambda_{\gamma_4,+\dot b}\lambda_{c,+\dot c},\qquad\epsilon^{\dot+\dot-}=1=-\epsilon^{\dot-\dot+}.\]", "",
        "The same pairing defines the brackets involving the other original source spinor.", "",
        "Each of the ten color tensors below multiplies its displayed complete source spinor polynomial.",
        r"Their sum acts on \(\gamma_{4i_4}\), \(c^B\), and \(\gamma_5\). Replacing each source by the derivative on its own surviving field and multiplying by", "",
        r"\[\frac{g^2}{128(8\pi^2)^3}\]", "",
        r"gives the complete \(g^2\bar W\) coefficient of \(\hbar^3\). The factor from \(\lambda=-\partial/2\) is already included in this multiplier; it must not be included a second time. See [DERIVATION.md](DERIVATION.md) for all signs, normalization, and the exact source-degree proof.", "",
        "## The complete rational tensor", "",
    ]
    for number, (word, channels) in enumerate(sorted(reduced.items()), 1):
        if set(channels) != {"1"}:
            raise ArithmeticError("the complete tensor is not rational")
        lines += ["### %d" % number, "", r"\[" + tensor_tex(word) + r"\]", "", "Full coefficient:", "",
                  "\\[\n\\begin{aligned}\n" + spinor_tex(channels["1"])
                  + "\n\\end{aligned}\n\\]", ""]
    lines += ["## Complete compact specialization with the free spectator set to one", "",
              r"Set \(\gamma_5=1\), hence \(\lambda_{\gamma_5}=0\). The table retains all ten color words, including the two that become zero. Factor out", "",
              r"\[\frac{\langle\lambda_{\gamma_4},\lambda_c\rangle^3}{4320}.\]", "",
              "| Color tensor | Remaining full spinor coefficient |", "| --- | --- |"]
    for word, channels in sorted(reduced.items()):
        coefficients = [4320 * sum((coefficient for vector, powers, coefficient in channels["1"]
                                     if vector == index and powers == (3, 0, 0)), Fraction(0))
                        for index in range(2)]
        lines.append(r"| \(%s\) | \(%s\) |" % (tensor_tex(word), compact_spinor(coefficients)))
    lines += ["", r"In the physical derivative formula the corresponding common multiplier, including this compact denominator, is \(g^2/[552960(8\pi^2)^3]\).", "",
              "## Full coefficients of the two representation defects", "",
              r"Before applying \([t_A,t_B]=iC^C{}_{AB}t_C\), the complete \(\zeta(3)\) summand consists of the following two displayed spinors multiplying their displayed tensor defects. These spinors retain all three original source spinors. The tensor defects vanish exactly.", ""]
    for flavor in (1, 2):
        indices = ["i_1", "i_2", "i_3", "i_4"]
        indices[flavor - 1] = "j_%d" % flavor
        defect = (r"\bigl([t_A,t_B]-iC^C{}_{AB}t_C\bigr)^{i_%d}{}_{j_%d}\bar W^{%s}"
                  % (flavor, flavor, " ".join(indices)))
        lines += ["### Defect on superpotential slot %d" % flavor, "", r"\[" + defect + r"\]", "",
                  "Complete spinor coefficient:", "",
                  "\\[\n\\begin{aligned}\n"
                  + spinor_tex(tensor["MM", ((flavor, ("A", "B")),)]["z3"])
                  + "\n\\end{aligned}\n\\]", ""]
    lines += ["The complete unreduced twelve-word tensor, including all rational terms before using the representation relation, is retained in [quartic_tensor_certificate.json](quartic_tensor_certificate.json). No color contraction or special source assignment is required for either cancellation.", ""]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    target = universal.HERE / "FULL_SOURCE_TENSOR.md"
    result = render()
    if args.write:
        target.write_text(result, encoding="utf-8")
        print("Wrote all ten complete rational spinors and both complete representation-defect spinors.")
    elif not target.is_file() or target.read_text(encoding="utf-8") != result:
        raise ArithmeticError("the complete printed tensor table differs from the exact coefficients")
    else:
        print("Verified every printed full-source tensor coefficient.")


if __name__ == "__main__":
    main()
