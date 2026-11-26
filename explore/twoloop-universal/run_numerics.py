"""Direct numerical integrals of the brute-force density, compared with the closed values.
  own:    numint.integrate (chart Gauss-Legendre, brute-force Kirchhoff integrand), m nodes per direction;
  oracle: htk.oracle.kernel (log-coordinate lattice in Schwinger space, engine layer c).
Usage: python3 run_numerics.py own|oracle   -> numerics_<mode>.json"""
import json
import math
import sys
import time
from fractions import Fraction as Fr

sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/twoloop-universal")
from twoloop import theta_graph, Theta, full_sources, kernel_C3_closed, kernel_htk_closed, sgn_detBC  # noqa
import numint

# sources at the non-root vertices (N = 0 and interior), fixed rationals
SRC = {
    1: [[Fr(3, 2)], [Fr(-2)], [Fr(1, 3)], [Fr(5, 4)], [Fr(-1, 2)]],
    2: [[Fr(1), Fr(2)], [Fr(-1), Fr(1, 2)], [Fr(2), Fr(-3)], [Fr(1, 3), Fr(1)], [Fr(-2), Fr(1, 2)]],
    3: [[Fr(1), Fr(0), Fr(1, 2)], [Fr(-1), Fr(1), Fr(2)], [Fr(2), Fr(-1), Fr(1)], [Fr(0), Fr(1), Fr(-1)],
        [Fr(1, 2), Fr(1), Fr(1)]],
}
CASES = [(1, 1, (1, 2, 2)), (2, 1, (1, 3, 3)), (2, 1, (2, 2, 3)), (2, 0, (1, 2, 2)), (3, 0, (1, 3, 3)),
         (3, 0, (2, 2, 3))]


def digits(a, b):
    if a == b:
        return 16.0
    return round(-math.log10(abs(a - b) / abs(b)), 1)


def main(mode):
    out = []
    for H, T, pqr in CASES:
        nv, edges = theta_graph(*pqr)
        th = Theta(nv, edges)
        lam = SRC[H][:nv - 1]
        lf = full_sources(nv, nv - 1, {v: lam[v] for v in range(nv - 1)})
        KC3 = kernel_C3_closed(th, lf, H, T)
        Kh = kernel_htk_closed(th, lf, H, T)
        exC3 = float(KC3) * (math.pi if T else 1.0)
        exh = float(Kh) * (math.pi if T else 1.0)
        lamf = [[float(c) for c in r] for r in lam]
        rec = {"H": H, "T": T, "pqr": pqr, "sources_nonroot": [[str(c) for c in r] for r in lam],
               "K_C3_exact": str(KC3) + ("*pi" if T else ""), "K_C3_float": exC3,
               "K_htk_exact": str(Kh) + ("*pi" if T else ""), "K_htk_float": exh}
        if mode == "own":
            dim = len(edges) - 1
            ms = (8, 12, 16, 24) if dim == 4 else (6, 8, 10, 12)
            runs = []
            for m in ms:
                t0 = time.time()
                val, npts = numint.integrate(pqr, lamf, H, T, m)
                runs.append({"m": m, "points": npts, "value": val, "digits_vs_exact": digits(val, exC3),
                             "time_s": round(time.time() - t0, 1)})
                print(H, T, pqr, runs[-1], flush=True)
            rec["own_chart_GL"] = runs
        else:
            import htk
            G = htk.Graph(nv, edges)
            t0 = time.time()
            val, err = htk.kernel_numeric(G, lamf, T)
            rec["oracle"] = {"value": val, "error_estimate": err, "digits_vs_exact": digits(val, exh),
                             "time_s": round(time.time() - t0, 1)}
            print(H, T, pqr, rec["oracle"], "exact", exh, flush=True)
        out.append(rec)
    with open((__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/twoloop-universal/numerics_%s.json" % mode, "w") as f:
        json.dump(out, f, indent=1)


if __name__ == "__main__":
    main(sys.argv[1])
