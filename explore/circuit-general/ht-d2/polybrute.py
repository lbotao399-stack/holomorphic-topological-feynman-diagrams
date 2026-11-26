"""Fast symbolic density with polynomial numerators.  D(t) = Num(t) / (t_{e0} U^k) * U^{-T/2}.
Forms are dicts {tuple(sorted edges): poly}; all polys in ring QQ[t_0..t_{N-1}, extra symbols]."""
import itertools
import sympy as sp
from sympy import QQ
from sympy.polys.rings import ring as mkring


def _wedge(a, b):
    out = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            if set(ka) & set(kb):
                continue
            perm = list(ka + kb)
            inv = sum(1 for i in range(len(perm)) for j in range(i + 1, len(perm)) if perm[i] > perm[j])
            key = tuple(sorted(perm))
            v = va * vb
            out[key] = out.get(key, 0) + (-v if inv % 2 else v)
    return {k: v for k, v in out.items() if v != 0}


def _add(a, b, c=1):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0) + c * v
    return {k: v for k, v in out.items() if v != 0}


def density_poly(G, lam_syms, H, T, C=None, extra=()):
    """lam_syms: list over non-root vertices of H-lists of symbol names or ints.
    Returns (R, gens, num, upow) with  D = num / (t_{N-1} * U^upow) * U^{-T/2},  U = det M (poly)."""
    C = G.cycle_matrix() if C is None else C
    n, N = len(C), G.ne
    names = ["t%d" % e for e in range(N)] + list(extra)
    R, *gens = mkring(",".join(names), QQ)
    t = gens[:N]
    sym = dict(zip(names, gens))
    g = G.nv - 1
    tree = G.spanning_tree()
    B = G.incidence()
    # tree solution y0 of the current law: B_T^t y0_T = lam
    BT = sp.Matrix([B[e] for e in tree]).T          # g x g
    BTinv = BT.inv()
    def conv(x):
        if isinstance(x, str):
            return sym[x]
        return R(x)
    lamR = [[conv(x) for x in vec] for vec in lam_syms]
    y0 = [[R(0)] * H for _ in range(N)]
    for k, e in enumerate(tree):
        for i in range(H):
            s = R(0)
            for v in range(g):
                c = BTinv[k, v]
                if c != 0:
                    s += R(QQ(int(sp.fraction(c)[0]), int(sp.fraction(c)[1]))) * lamR[v][i]
            y0[e][i] = s
    M = [[sum((t[e] * C[r][e] * C[s][e] for e in range(N)), R(0)) for s in range(n)] for r in range(n)]
    Ms = sp.Matrix(n, n, lambda r, s: M[r][s].as_expr())
    U = R(sp.expand(Ms.det()))
    adj = Ms.adjugate()
    adjR = [[R(sp.expand(adj[r, s])) for s in range(n)] for r in range(n)]
    # loop currents numerators: L = -adj C t y0  (per component)
    Ynum = [[None] * H for _ in range(N)]
    for i in range(H):
        Cty0 = [sum((C[r][e] * t[e] * y0[e][i] for e in range(N)), R(0)) for r in range(n)]
        L = [-sum((adjR[r][s] * Cty0[s] for s in range(n)), R(0)) for r in range(n)]
        for e in range(N):
            Ynum[e][i] = U * y0[e][i] + sum((C[r][e] * L[r] for r in range(n)), R(0))
    # psi_r^i numerators (divide by U)
    psi = [{(e,): C[r][e] * Ynum[e][i] for e in range(N) if C[r][e] != 0} for i in range(H) for r in range(n)]
    hol = {(): R(1)}
    for f in psi:
        hol = _wedge(hol, f)
    upow = H * n  # from psi denominators
    if T == 1:
        dM = [[{(e,): R(C[r][e] * C[a][e]) for e in range(N) if C[r][e] * C[a][e] != 0} for a in range(n)] for r in range(n)]
        A = [[{} for _ in range(n)] for _ in range(n)]
        for r in range(n):
            for s in range(r + 1, n):
                acc = {}
                for a in range(n):
                    for b in range(n):
                        if adjR[a][b] != 0:
                            acc = _add(acc, _wedge(dM[r][a], dM[s][b]), adjR[a][b] * QQ(-1, 2))
                A[r][s] = acc
                A[s][r] = {k: -v for k, v in acc.items()}
        def pf(ix):
            if not ix:
                return {(): R(1)}
            if len(ix) % 2:
                return {}
            out = {}
            for j in range(1, len(ix)):
                term = _wedge(A[ix[0]][ix[j]], pf(ix[1:j] + ix[j + 1:]))
                out = _add(out, term, (-1) ** (j - 1))
            return out
        top = pf(list(range(n)))
        upow += n // 2   # A has 1/U
        form = _wedge(top, hol)
    else:
        form = hol
    upow += H  # U^{-H}
    key = tuple(range(N - 1))
    num = form.get(key, R(0))
    if (N - 1) % 2:
        num = -num
    return R, gens, num, upow, U
