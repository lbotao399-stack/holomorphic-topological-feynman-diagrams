"""General-source kernel K_n(lam; 0) at (H,T) = (1,1) as an explicit finite sum (independent
implementation of the configuration formula; compared with integrate_general.py and the engine).

  K_n(lam; 0) = 2^{-n/2} sum_{O} N(O) sum_{w, chi} prod_{i in O} lam_{w_i} prod_{j in C} lam_{chi_j}
               * prod_{l=1}^{n-1} B(e_l + 1, f_l + 1) * prod_{m=1}^{n-1} B(3/2, m/2 - c_m) * B(e_n + 1, f_n + p + 1) / p

  O: openers (smaller elements of the pairs of a perfect matching of 1..n), C its complement, N(O) the
  number of matchings with opener set O;  w_i in {1..i};  chi_j in {(w,k): j < w <= k <= n} u {inf},
  lam_{(w,k)} = lam_w, lam_inf = lam_{n+1};  p = #{chi_j = inf};  f_l = #{chi_j = (w,l)};
  c_m = #{chi_j = (w,k): k <= m};  e_l = (n-1-l)/2 - #{i: w_i > l} + #{j: chi_j = (w,k), w <= l < k}.
"""
import sys
import itertools
from fractions import Fraction as Fr

sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/fans-ht/code")
from integrate_general import beta2


def opener_sets(n):
    for O in itertools.combinations(range(1, n + 1), n // 2):
        Os = set(O)
        avail, cnt, ok = 0, 1, True
        for k in range(1, n + 1):
            if k in Os:
                avail += 1
            else:
                if avail == 0:
                    ok = False
                    break
                cnt *= avail
                avail -= 1
        if ok:
            yield list(O), [k for k in range(1, n + 1) if k not in Os], cnt


def kernel_configurations(n):
    """{lam exponent tuple (length n+1): {k: q}}, value = sum q sqrt(pi)^k."""
    res = {}
    for O, C, N in opener_sets(n):
        wch = [list(range(1, i + 1)) for i in O]
        cch = [[(w, k) for k in range(j + 1, n + 1) for w in range(j + 1, k + 1)] + ["inf"] for j in C]
        for ws in itertools.product(*wch):
            for chis in itertools.product(*cch):
                p = sum(1 for c in chis if c == "inf")
                lamexp = [0] * (n + 1)
                for w in ws:
                    lamexp[w - 1] += 1
                for c in chis:
                    lamexp[(n if c == "inf" else c[0] - 1)] += 1
                q, kpi = Fr(N, 2 ** (n // 2)), 0
                for l in range(1, n + 1):
                    e2 = (n - 1 - l) - 2 * sum(1 for w in ws if w > l) + \
                        2 * sum(1 for c in chis if c != "inf" and c[0] <= l < c[1])
                    f = sum(1 for c in chis if c != "inf" and c[1] == l)
                    if l < n:
                        dq, dk = beta2(e2 + 2, 2 * f + 2)
                    else:
                        dq, dk = beta2(e2 + 2, 2 * (f + p) + 2)
                        dq = dq / p
                    q *= dq
                    kpi += dk
                for m in range(1, n):
                    cm = sum(1 for c in chis if c != "inf" and c[1] <= m)
                    dq, dk = beta2(3, m - 2 * cm)
                    q *= dq
                    kpi += dk
                key = tuple(lamexp)
                slot = res.setdefault(key, {})
                slot[kpi] = slot.get(kpi, 0) + q
    return res


if __name__ == "__main__":
    import sympy as sp
    from integrate_general import kernel_general
    for n in [int(x) for x in sys.argv[1:]] or [2, 4, 6]:
        a = kernel_configurations(n)
        b, _, _ = kernel_general(n)
        lam = sp.symbols("l1:%d" % (n + 2))
        def expr(res):
            out = 0
            for lk, d in res.items():
                for k, q in d.items():
                    if q:
                        out += sp.Rational(q.numerator, q.denominator) * sp.sqrt(sp.pi) ** k * sp.Mul(*[x ** pp for x, pp in zip(lam, lk)])
            return sp.expand(out)
        ea, eb = expr(a), expr(b)
        assert sp.expand(ea - eb) == 0
        print("n=%d: configuration formula = expanded Beta integration (%d lam-monomials)" % (n, len(sp.Add.make_args(ea))))
