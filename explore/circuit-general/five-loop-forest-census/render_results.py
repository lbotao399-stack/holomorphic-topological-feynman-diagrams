"""Write every finite five-loop kernel with the full current spinor retained."""
from fractions import Fraction
import json

from certify import HERE, PREVIOUS, add_records
import importlib.util

spec = importlib.util.spec_from_file_location("_previous_forest_renderer", PREVIOUS / "render_results.py")
previous = importlib.util.module_from_spec(spec)
spec.loader.exec_module(previous)


def render():
    certificate = json.loads((HERE / "certificate.json").read_text())
    out = ["# Every five-loop bare kernel and first-current spinor", "",
           "These are complete geometric kernels at the recorded source assignments. "
           "All 200 topologies and all original edges are present. Each row gives the "
           "two coefficients multiplying the indicated original source spinors; the "
           r"free lower index $+\dot a$ is retained. The source assignments are fixed "
           "inputs, not an assertion of arbitrary-source interpolation. Wick factors, "
           "couplings and endpoint displacement signs are not included.", "",
           "[Complete exact coefficient archive](certificate.json).", ""]
    for record in certificate["records"]:
        A, B = record["marked"]
        sources = [tuple(map(Fraction, row)) for row in record["full_balanced_sources"]]
        source = next(vertex for vertex in range(A)
                      if sources[vertex][0] * sources[A][1] - sources[vertex][1] * sources[A][0])
        denominator = sources[source][0] * sources[A][1] - sources[source][1] * sources[A][0]
        out.extend(["## Forest %d" % record["index"], "",
                    "Ordered marked vertices: `%d,%d`. Ordered edges: `%s`." %
                    (A, B, ", ".join("%d→%d" % tuple(edge) for edge in record["edges"])), "",
                    "Full balanced source rows in vertex order: `%s`." %
                    ", ".join("(%s,%s)" % tuple(row) for row in record["full_balanced_sources"]), "",
                    "Bare kernel:", "", "$$", previous.latex(record["bare"]), "$$", "",
                    r"| Original edge | Coefficient of $\lambda_{%d,+\dot a}$ | Coefficient of $\lambda_{%d,+\dot a}$ |" % (source, A),
                    "|---|---|---|"])
        for edge, spinor in zip(record["edges"], record["all_original_first_current_spinors"]):
            first = add_records(((sources[A][1] / denominator, spinor[0]),
                                 (-sources[A][0] / denominator, spinor[1])))
            second = add_records(((-sources[source][1] / denominator, spinor[0]),
                                  (sources[source][0] / denominator, spinor[1])))
            for entry in range(2):
                if add_records(((sources[source][entry], first), (sources[A][entry], second))) != spinor[entry]:
                    raise ArithmeticError("a complete original spinor changed during rendering")
            out.append("| %d→%d | $%s$ | $%s$ |" % (edge[0], edge[1], previous.latex(first), previous.latex(second)))
        out.append("")
    return "\n".join(out)


if __name__ == "__main__":
    (HERE / "FULL_RESULTS.md").write_text(render())
    print("Printed all 200 bare kernels and all 2200 original current spinors.")
