"""Assemble results.json: (d) tables of K(lambda; 0) for all theta graphs with D <= 6, T in {0, 1}, as sums
of products of two H x H determinants, plus the status of the checks (outputs of run_*.py).
The tables are re-verified here against the closed determinant formula at random rational sources."""
import json
import random
import sys

sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/twoloop-universal")
from twoloop import *  # noqa

DIR = (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/twoloop-universal/"
rng = random.Random(4242)


def rf():
    return Fr(rng.randint(-9, 9), rng.randint(1, 4))


def sign_formula_T1(p, q, H, Bp):
    """epsilon(B') in D_theta = sum eps [Lam_A Lam_B'][Lam_B'' Lam_C], |B'| = H - p + 1 (B' 1-based in B)"""
    return sgn_power(H * (H + 1) // 2 + p * (p - 1) // 2 + (H - p + 1) * (p - 1) + sum(Bp))


def sign_formula_T0(p, q, H, Bp):
    """epsilon_0(B') in Q(v) = sum eps [v Lam_A Lam_B'][v Lam_B'' Lam_C], |B'| = H - p"""
    return sgn_power(H * (H + 1) // 2 + (p + 1) * (p + 2) // 2 + (H - p) * (p + 1) + sum(Bp))


def check_sign_formulas(th, H, T):
    """laplace_expansion signs (programmatic) against the closed sign formulas"""
    p, q, r = th.pqr
    nm = vertex_names(th)
    if T == 1:
        spec = [(nm[u], "top") for u in th.interior["A"]] + [(nm[u], "both") for u in th.interior["B"]] + \
               [(nm[u], "bot") for u in th.interior["C"]]
    else:
        spec = [("@L", "top"), ("@L", "bot")] + [(nm[u], "top") for u in th.interior["A"]] + \
               [(nm[u], "both") for u in th.interior["B"]] + [(nm[u], "bot") for u in th.interior["C"]]
    ok = True
    for sgn, top, bot in laplace_expansion(spec, H):
        Bp = [int(x[1:]) for x in top if x.startswith("B")]
        f = sign_formula_T1(p, q, H, Bp) if T == 1 else sign_formula_T0(p, q, H, Bp)
        ok &= (f == sgn)
    return ok


tables = []
allok = True
for T in (0, 1):
    for D in range(1, 7):
        H = D - T
        if H < 0 or (T == 0 and H < 1):
            continue
        for pqr in theta_list(D):
            p, q, r = pqr
            nv, edges = theta_graph(*pqr)
            th = Theta(nv, edges)
            nm = vertex_names(th)
            mu = mu_weights(th)
            if T == 0:
                c = c0(pqr, H)
                sC3 = sgn_power(H * (H + 1) // 2 + 1)
                terms = expand_T0(th, H)
            else:
                c = gamma_half_ratio(pqr, H)
                sC3 = sgn_power(H * (H + 1) // 2)
                terms = expand_T1(th, H)
            shtk = sC3 * sgn_power((q + 1) * D)
            assert sgn_detBC(nv, edges, nv - 1, th.C3()) ** D == sgn_power((q + 1) * D)
            # verification of the table against the closed formula
            ok = check_sign_formulas(th, H, T)
            for _ in range(3):
                lam = full_sources(nv, nv - 1, {v: [rf() for _ in range(H)] for v in range(nv - 1)})
                val = shtk * c * eval_terms(terms, th, lam)
                ok &= (val == kernel_htk_closed(th, lam, H, T))
            allok &= ok
            entry = {
                "T": T, "H": H, "D": D, "pqr": [p, q, r], "edges": len(edges), "vertices": nv,
                "vertex_labels_htk_theta": {nm[v]: v for v in sorted(nm)},
                "c": str(c), "c_meaning": ("prod_P (H-n_P)! / ((H+1)! prod_P (n_P-1)!)" if T == 0 else
                                           "prod_P Gamma(H+3/2-n_P) / (2 Gamma(H+3/2) prod_P (n_P-1)!) / pi"),
                "sign_C3_basis": sC3, "sign_htk": shtk,
                "formula": ("K_htk(lam;0) = sign_htk * c * Q(Lbar)" if T == 0 else
                            "K_htk(lam;0) = sign_htk * c * pi * D_theta"),
                "Lbar_weights": {nm[v]: str(mu[v]) for v in mu if v != th.Np} if T == 0 else None,
                "polynomial": ("Q(Lbar)" if T == 0 else "D_theta") + " = sum coef [top][bottom]",
                "terms": [[str(cf), a, b] for cf, a, b in terms],
                "K_htk_terms": [[str(shtk * c * cf) + ("" if T == 0 else "*pi"), a, b] for cf, a, b in terms],
                "verified_against_closed_formula": ok,
            }
            if (H, T) in ((2, 0), (3, 0), (5, 1)):
                entry["source_kernel"] = {(2, 0): "I^S1 = -K_htk", (3, 0): "I^S2 = -K_htk",
                                          (5, 1): "K^S3 = -K_htk"}[(H, T)]
            tables.append(entry)
            print("T=%d H=%d %s c=%s sign_C3=%+d sign_htk=%+d terms=%d verified=%s" % (
                T, H, pqr, c, sC3, shtk, len(terms), ok), flush=True)

res = {"topic": "twoloop-universal",
       "conventions": {
           "graph": "htk.graph.theta(p,q,r): N = 0, N' = root = p+q+r-2, interior vertices along A, B, C; edge order A, B, C; edges N -> N'",
           "names": "N, Np = trivalent vertices (Np = root), Pk = k-th interior vertex of path P from N; [v1..vH] = det of the H x H matrix with columns lam_v1..lam_vH",
           "K_htk": "reduced kernel of P/code/htk/CONVENTIONS.md (positively oriented cycle basis)",
           "C3": "cycle matrix (A-B, B-C); K[C3] = sign_C3_basis * c * (...); K_htk = (-1)^{(q+1)D} K[C3]",
           "root_at_interior_vertex": "K as a function of the full source vector (sum zero) is the same for C3; substitute lam_root = -sum_{v != root} lam_v; htk sign: K_htk[G relabelled by the transposition rho <-> N'] = (-1)^D K_htk[G]",
           "S_normalizations": "I^S1 = (-1)^{n(n+1)/2} K_htk (H=2,T=0); I^S2 = same (H=3,T=0); K^S3 = same (H=5,T=1); n = 2 gives -1",
       },
       "tables_verified": allok,
       "tables": tables}

checks = {}
for name in ("exact_results", "sources_results", "moments_results", "numerics_own", "numerics_oracle",
             "formscheck_results"):
    try:
        with open(DIR + name + ".json") as f:
            checks[name] = json.load(f)
    except FileNotFoundError:
        checks[name] = None
res["checks"] = checks
with open(DIR + "results.json", "w") as f:
    json.dump(res, f, indent=1, default=str)
print("tables verified:", allok)
