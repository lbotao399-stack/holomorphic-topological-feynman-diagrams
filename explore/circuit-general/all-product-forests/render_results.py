"""Print every certified first-current spinor in an original-source basis."""
from fractions import Fraction
import json
from pathlib import Path

from certify_integrals import add_records

HERE = Path(__file__).resolve().parent


def latex(value):
    terms = []
    for label, coefficient in sorted(value.items(), key=lambda item: (item[0] != "1", item[0])):
        coefficient = Fraction(coefficient)
        if not coefficient:
            continue
        magnitude = abs(coefficient)
        rational = str(magnitude.numerator) if magnitude.denominator == 1 else (
            r"\frac{%d}{%d}" % (magnitude.numerator, magnitude.denominator))
        if label != "1":
            rational = ("" if magnitude == 1 else rational) + "".join(
                r"\zeta(%s)" % factor[1:] for factor in label.split("*"))
        terms.append(("-" if coefficient < 0 else "+", rational))
    if not terms:
        return "0"
    return ("-" if terms[0][0] == "-" else "") + terms[0][1] + "".join(
        sign + value for sign, value in terms[1:])


def render():
    certificate = json.loads((HERE / "integral_certificate.json").read_text())
    out = ["# Every bare kernel and first-current spinor through four loops", "",
           "The ordered insertions and complete source assignments below are fixed inputs. "
           "The tables contain every original edge. Each row is the coefficient of the "
           "two indicated original source spinors in the full integral of "
           r"$y_{e,+\dot a}$; the free lower spinor index is retained. "
           "The two columns together give the entire spinor, including its sign. "
           "These are geometric kernels in the repository convention, before physical "
           "Wick factors and displacement signs. They are exact values at the displayed "
           "source assignments, not an interpolation claiming arbitrary-source coefficients.", "",
           "Both the compositional forest solver and direct original-graph integration "
           "verify every entry. Machine-readable complete results are in "
           "[integral_certificate.json](integral_certificate.json).", ""]
    for record in certificate["records"]:
        A, B = record["marked"]
        sources = [tuple(map(Fraction, row)) for row in record["full_balanced_sources"]]
        source = next(vertex for vertex in range(A)
                      if sources[vertex][0] * sources[A][1] - sources[vertex][1] * sources[A][0])
        denominator = sources[source][0] * sources[A][1] - sources[source][1] * sources[A][0]
        out.extend(["## %d loops, forest %d" % (record["loops"], record["index"]), "",
                    "Marked vertices: `%d,%d`. Ordered edges: `%s`." %
                    (A, B, ", ".join("%d→%d" % tuple(edge) for edge in record["edges"])), "",
                    "Full source rows in vertex order: `%s`." %
                    ", ".join("(%s,%s)" % tuple(row) for row in record["full_balanced_sources"]), "",
                    "Bare kernel:", "", "$$", latex(record["bare"]), "$$", "",
                    r"| Original edge | Coefficient of $\lambda_{%d,+\dot a}$ | Coefficient of $\lambda_{%d,+\dot a}$ |" % (source, A),
                    "|---|---|---|"])
        for edge, spinor in zip(record["edges"], record["all_original_first_current_spinors"]):
            first = add_records(((sources[A][1] / denominator, spinor[0]),
                                 (-sources[A][0] / denominator, spinor[1])))
            second = add_records(((-sources[source][1] / denominator, spinor[0]),
                                  (sources[source][0] / denominator, spinor[1])))
            for entry in range(2):
                if add_records(((sources[source][entry], first), (sources[A][entry], second))) != spinor[entry]:
                    raise ArithmeticError("the complete source-spinor form changed a certified entry")
            out.append("| %d→%d | $%s$ | $%s$ |" % (edge[0], edge[1], latex(first), latex(second)))
        out.append("")
    return "\n".join(out)


if __name__ == "__main__":
    path = HERE / "FULL_RESULTS.md"
    path.write_text(render())
    print("Wrote every bare kernel and all first-current spinors in their original source bases.")
