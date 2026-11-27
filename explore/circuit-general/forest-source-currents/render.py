"""Print every nonzero covariant coefficient of every original free spinor."""
import argparse
from fractions import Fraction
import json

from source_currents import HERE


def rational(value):
    value = Fraction(value)
    if value.denominator == 1:
        return str(value.numerator)
    sign = "-" if value < 0 else ""
    return sign + r"\frac{%d}{%d}" % (abs(value.numerator), value.denominator)


def power(base, exponent):
    if not exponent:
        return ""
    return base if exponent == 1 else base + "^{%d}" % exponent


def monomial(i, j):
    return "".join((power(r"\langle\lambda_1,\lambda_2\rangle", 5 - i - j),
                    power(r"\langle\lambda_5,\lambda_2\rangle", i),
                    power(r"\langle\lambda_1,\lambda_5\rangle", j))) or "1"


def row(prefix, values):
    return "| " + " | ".join(prefix + ["$" + rational(values.get(label, 0)) + "$"
                                      for label in ("1", "z3", "z5")]) + " |"


def render():
    document = json.loads((HERE / "certificate.json").read_text())
    reconstruction = document["reconstruction"]
    lines = ["# Complete supported-source bare kernel and eleven current spinors", "",
        "The independent original spinors are $\\lambda_1$, $\\lambda_2$, and $\\lambda_5$. "
        "The other original sources are $\\lambda_6=-\\lambda_1-\\lambda_2-\\lambda_5$ "
        "and $\\lambda_0=\\lambda_3=\\lambda_4=0$.", "",
        "Every bracket is the determinant of its two displayed source spinors. "
        "In each row, multiply the displayed monomial and full original-source spinor "
        "by the rational column plus the second column times $\\zeta(3)$ plus the third "
        "column times $\\zeta(5)$. Add every row of the relevant original edge. "
        "All unlisted monomial coefficients are exactly zero. Every displayed bracket "
        "power is nonnegative, so these formulas include degenerate source frames.", "",
        "## Complete bare kernel", "",
        "| Source monomial | Rational | $\\zeta(3)$ | $\\zeta(5)$ |",
        "|---|---:|---:|---:|"]
    for (i, j), values in zip(reconstruction["frame_monomials"], reconstruction["bare_frame_coefficients"]):
        if values:
            lines.append(row(["$" + monomial(i, j) + "$"], values))
    for edge, ((left, right), covariants) in enumerate(zip(document["edges"],
            reconstruction["denominator_free_original_source_covariants"])):
        lines += ["", "## Original edge %d: %d to %d" % (edge, left, right), "",
            "| Complete free-spinor monomial | Rational | $\\zeta(3)$ | $\\zeta(5)$ |",
            "|---|---:|---:|---:|"]
        count = 0
        for vertex, coefficients in zip((1, 2, 5), covariants):
            for (i, j), values in coefficients:
                lines.append(row(["$" + monomial(i, j) + r"\lambda_{%d,+\dot a}" % vertex + "$"], values))
                count += 1
        if not count:
            lines.append(row(["$0$"], {}))
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    value = render()
    path = HERE / "FULL_RESULTS.md"
    if args.write:
        path.write_text(value)
    elif not path.exists() or path.read_text() != value:
        raise ArithmeticError("the complete printed free-spinor coefficients changed")
    print("Every bare and original free-spinor coefficient is printed exactly.")


if __name__ == "__main__":
    main()
