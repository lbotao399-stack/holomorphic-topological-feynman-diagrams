"""Two-loop (theta-graph) kernels on R^T x C^H, T in {0, 1}: own exact brute force and closed forms.

Conventions (P/code/htk/CONVENTIONS.md, Sections 1-2):
  graph (nv, edges, root): vertices 0..nv-1, edges (e(0), e(1)) oriented, edge order = list order;
  Ohm t_e y_e = phi_{e(0)} - phi_{e(1)}, phi_root = 0; KCL sum_{e(0)=v} y_e - sum_{e(1)=v} y_e = lam_v.
  Omega = U^{-H} rho_2(M)^T ^ prod_{i=1..H} (psi^i_1 ^ psi^i_2) = D(t) iota_E dt_0 ^ ... ^ dt_{N-1},
  psi^i_r = sum_e C_re y^i_e dt_e,  M = C t C^T,  rho_2 = U^{-1/2} Pf(-1/2 dM M^{-1} dM),
  K = int_{RP^{N-1}_{>0}} D(t) e^{-sum y.z} |iota_E vol|.

Theta structure: trivalent vertices N, N'; paths P in (A, B, C) from N to N' with n_P edges
P_1..P_{n_P}; s_e = +1 if e is oriented along N -> N'; interior vertices P_1..P_{n_P - 1}
(vertex P_k sits between the edges P_k and P_{k+1}).  Cycle matrix C3 (S3 Section 7):
row 1 = A - B, row 2 = B - C (path vectors with the signs s_e).

Closed forms derived in REPORT.md (Section 1), for C3 and the standard edge order
(A_1..A_p, B_1..B_q, C_1..C_r):
  R_P = sum_j t_{P_j},  U = R_A R_B + R_B R_C + R_C R_A,  sigma^P_k = sum_{j>k} t_{P_j} / R_P,
  L = lam_N + sum_P sum_k sigma^P_k lam^P_k,
  Q(v) = det[(v;0), (0;v), (Lam_A;0), (Lam_B;Lam_B), (0;Lam_C)],
  D_theta = det[(Lam_A;0), (Lam_B;Lam_B), (0;Lam_C)],
  T = 0:  D(t) = (-1)^{H(H+1)/2+1} R_A R_B R_C Q(L) / U^{H+2},
  T = 1:  D(t) = (-1)^{H(H+1)/2} D_theta / (2 U^{H+3/2}).
"""
from fractions import Fraction as Fr
import itertools
import math

# ---------------------------------------------------------------------------- exact linear algebra


def det(M):
    n = len(M)
    if n == 0:
        return Fr(1)
    A = [[Fr(x) for x in row] for row in M]
    s = Fr(1)
    for c in range(n):
        piv = None
        for r in range(c, n):
            if A[r][c] != 0:
                piv = r
                break
        if piv is None:
            return Fr(0)
        if piv != c:
            A[c], A[piv] = A[piv], A[c]
            s = -s
        pv = A[c][c]
        s *= pv
        Ac = A[c]
        for r in range(c + 1, n):
            f = A[r][c] / pv
            if f:
                Ar = A[r]
                for k in range(c, n):
                    Ar[k] -= f * Ac[k]
    return s


def solve(A, B):
    """A X = B, A square (Fractions); B: list of columns -> list of solution columns."""
    n = len(A)
    M = [[Fr(x) for x in A[i]] + [Fr(b[i]) for b in B] for i in range(n)]
    m = len(B)
    for c in range(n):
        piv = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[piv] = M[piv], M[c]
        pv = M[c][c]
        M[c] = [x / pv for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return [[M[i][n + j] for i in range(n)] for j in range(m)]


def inv(A):
    n = len(A)
    cols = solve(A, [[Fr(int(i == j)) for i in range(n)] for j in range(n)])
    return [[cols[j][i] for j in range(n)] for i in range(n)]


def perm_sign(seq):
    seq = list(seq)
    s = 1
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                s = -s
    return s


def sgn_power(k):
    return -1 if k % 2 else 1


# ---------------------------------------------------------------------------- graphs


def theta_graph(p, q, r):
    """Same labelling as htk.graph.theta: N = 0, N' = root = nv-1, interior vertices numbered
    along A, then B, then C; edge order A, B, C; all edges oriented N -> N'."""
    nv = p + q + r - 1
    edges, nxt = [], 1
    for n_ in (p, q, r):
        prev = 0
        for j in range(n_):
            cur = nv - 1 if j == n_ - 1 else nxt
            if j < n_ - 1:
                nxt += 1
            edges.append((prev, cur))
            prev = cur
    return nv, edges


class Theta:
    """Theta structure of a labelled graph.  paths: {P: [(edge, s_e)]} in order from N to N';
    interior: {P: [vertices]} in order from N."""

    def __init__(self, nv, edges, root=None, N=None, Np=None, order=None):
        self.nv, self.edges = nv, list(edges)
        self.root = nv - 1 if root is None else root
        deg = [0] * nv
        for a, b in edges:
            deg[a] += 1
            deg[b] += 1
        tri = [v for v in range(nv) if deg[v] == 3]
        if nv == 2:
            tri = [0, 1]
        assert len(tri) == 2 and all(d in (2, 3) or nv == 2 for d in deg), (deg,)
        if N is None:
            if self.root in tri:
                Np = self.root
                N = tri[0] if tri[1] == Np else tri[1]
            else:
                N, Np = tri
        self.N, self.Np = N, Np
        inc = {v: [] for v in range(nv)}
        for e, (a, b) in enumerate(edges):
            inc[a].append(e)
            inc[b].append(e)
        paths = []
        for e0 in inc[N]:
            path, cur, e, ints = [], N, e0, []
            while True:
                a, b = edges[e]
                s = 1 if a == cur else -1
                nxt = b if a == cur else a
                path.append((e, s))
                if nxt == Np:
                    break
                ints.append(nxt)
                cur = nxt
                e = next(f for f in inc[cur] if f != e)
            paths.append((path, ints))
        if order is None:
            paths.sort(key=lambda pi: (len(pi[0]), min(e for e, _ in pi[0])))
        else:
            paths = [paths[i] for i in order]
        self.paths = {P: paths[i][0] for i, P in enumerate("ABC")}
        self.interior = {P: paths[i][1] for i, P in enumerate("ABC")}
        self.n = {P: len(self.paths[P]) for P in "ABC"}
        self.pqr = tuple(self.n[P] for P in "ABC")
        self.std_order = [e for P in "ABC" for e, _ in self.paths[P]]
        # iota_E vol in the graph's edge order = perm_sign * iota_E vol in the standard order
        self.order_sign = perm_sign(self.std_order)

    def C3(self):
        Ne = len(self.edges)
        C = [[0] * Ne, [0] * Ne]
        for e, s in self.paths["A"]:
            C[0][e] = s
        for e, s in self.paths["B"]:
            C[0][e] = -s
            C[1][e] = s
        for e, s in self.paths["C"]:
            C[1][e] = -s
        return C


def full_sources(nv, root, lam):
    """lam: dict v -> H-vector for some vertices (root optional) -> full list summing to zero."""
    H = len(next(iter(lam.values()))) if lam else 0
    out = [[Fr(0)] * H for _ in range(nv)]
    for v, vec in lam.items():
        if v != root:
            out[v] = [Fr(x) for x in vec]
    out[root] = [-sum(out[v][i] for v in range(nv) if v != root) for i in range(H)]
    return out


# ---------------------------------------------------------------------------- brute force (own)


def incidence(nv, edges, root):
    g = [v for v in range(nv) if v != root]
    idx = {v: k for k, v in enumerate(g)}
    B = [[0] * len(g) for _ in edges]
    for e, (a, b) in enumerate(edges):
        if a != root:
            B[e][idx[a]] += 1
        if b != root:
            B[e][idx[b]] -= 1
    return B, g


def currents(nv, edges, root, t, lamfull):
    """Kirchhoff currents y[e][i] from the grounded Laplacian (exact)."""
    B, g = incidence(nv, edges, root)
    H = len(lamfull[0])
    m = len(g)
    Lap = [[Fr(0)] * m for _ in range(m)]
    for e in range(len(edges)):
        ce = 1 / Fr(t[e])
        for i in range(m):
            if B[e][i]:
                for j in range(m):
                    if B[e][j]:
                        Lap[i][j] += B[e][i] * B[e][j] * ce
    if H == 0:
        return [[] for _ in edges]
    phis = solve(Lap, [[lamfull[v][i] for v in g] for i in range(H)])
    y = []
    for e in range(len(edges)):
        y.append([sum((B[e][k] * phis[i][k] for k in range(m)), Fr(0)) / Fr(t[e]) for i in range(H)])
    return y


def check_kirchhoff(nv, edges, root, t, lamfull, y, C):
    H = len(lamfull[0])
    for v in range(nv):
        for i in range(H):
            s = sum((y[e][i] for e, (a, b) in enumerate(edges) if a == v), Fr(0)) - \
                sum((y[e][i] for e, (a, b) in enumerate(edges) if b == v), Fr(0))
            if s != lamfull[v][i]:
                return False
    for row in C:
        for i in range(H):
            if sum(row[e] * Fr(t[e]) * y[e][i] for e in range(len(edges))) != 0:
                return False
    return True


def period(C, t):
    return [[sum(Fr(ra[e] * rb[e]) * Fr(t[e]) for e in range(len(t))) for rb in C] for ra in C]


def brute_density(nv, edges, root, t, lamfull, H, T, C):
    """D(t) = coeff * radicand^{-1/2} (radicand = U for T = 1, 1 for T = 0), from the forms."""
    t = [Fr(x) for x in t]
    Ne = len(edges)
    y = currents(nv, edges, root, t, lamfull) if H else [[] for _ in edges]
    n = len(C)
    assert n == 2 and n * (H + T) == Ne - 1
    Psi = [[C[r][e] * y[e][i] for e in range(Ne)] for i in range(H) for r in range(n)]
    M = period(C, t)
    U = det(M)
    e0 = Ne - 1
    cols = [e for e in range(Ne) if e != e0]
    if T == 0:
        val = det([[row[e] for e in cols] for row in Psi])
        return sgn_power(e0) * val / (t[e0] * U ** H), Fr(1)
    Mi = inv(M)
    # A12 = -1/2 sum_ab Minv_ab dM_1a ^ dM_b2 ; coefficient of dt_e ^ dt_f (e < f)
    def a12(e, f):
        s = Fr(0)
        for a in range(n):
            for b in range(n):
                if Mi[a][b]:
                    s += Mi[a][b] * (C[0][e] * C[a][e] * C[b][f] * C[1][f] - C[0][f] * C[a][f] * C[b][e] * C[1][e])
        return -s / 2
    top = Fr(0)
    for i, e in enumerate(cols):
        for j in range(i + 1, len(cols)):
            f = cols[j]
            c = a12(e, f)
            if c == 0:
                continue
            rest = [x for x in cols if x not in (e, f)]
            d = det([[row[x] for x in rest] for row in Psi]) if Psi else Fr(1)
            top += sgn_power(i + j - 1) * c * d
    return sgn_power(e0) * top / (t[e0] * U ** H), U


# ---------------------------------------------------------------------------- closed forms


def Q_form(th, lamfull, v):
    """Q(v) = det[(v;0), (0;v), (Lam_A;0), (Lam_B;Lam_B), (0;Lam_C)] (2H x 2H)."""
    H = len(v)
    I = th.interior
    cols = [[x for x in v] + [Fr(0)] * H, [Fr(0)] * H + [x for x in v]]
    for u in I["A"]:
        cols.append(list(lamfull[u]) + [Fr(0)] * H)
    for u in I["B"]:
        cols.append(list(lamfull[u]) + list(lamfull[u]))
    for u in I["C"]:
        cols.append([Fr(0)] * H + list(lamfull[u]))
    assert len(cols) == 2 * H
    return det([[cols[c][r] for c in range(2 * H)] for r in range(2 * H)])


def D_theta(th, lamfull):
    H = len(lamfull[0])
    if H == 0:
        return Fr(1)
    I = th.interior
    cols = []
    for u in I["A"]:
        cols.append(list(lamfull[u]) + [Fr(0)] * H)
    for u in I["B"]:
        cols.append(list(lamfull[u]) + list(lamfull[u]))
    for u in I["C"]:
        cols.append([Fr(0)] * H + list(lamfull[u]))
    assert len(cols) == 2 * H
    return det([[cols[c][r] for c in range(2 * H)] for r in range(2 * H)])


def mu_weights(th):
    """mu_v = (edges from v to N' along its path) / n_P; mu_N = 1, mu_N' = 0."""
    mu = {th.N: Fr(1), th.Np: Fr(0)}
    for P in "ABC":
        nP = th.n[P]
        for k, u in enumerate(th.interior[P], start=1):
            mu[u] = Fr(nP - k, nP)
    return mu


def Lbar(th, lamfull):
    """bar L = sum_v mu_v lam_v."""
    H = len(lamfull[0])
    mu = mu_weights(th)
    return [sum((mu[v] * lamfull[v][i] for v in mu), Fr(0)) for i in range(H)]


def L_of_t(th, t, lamfull):
    H = len(lamfull[0])
    L = [Fr(x) for x in lamfull[th.N]]
    for P in "ABC":
        es = [e for e, _ in th.paths[P]]
        R = sum(Fr(t[e]) for e in es)
        for k, u in enumerate(th.interior[P], start=1):
            sig = sum(Fr(t[e]) for e in es[k:]) / R
            for i in range(H):
                L[i] += sig * lamfull[u][i]
    return L


def closed_density_t(th, t, lamfull, H, T):
    """D(t) for the cycle matrix th.C3() in the graph's edge order: (coeff, radicand)."""
    t = [Fr(x) for x in t]
    R = {P: sum(t[e] for e, _ in th.paths[P]) for P in "ABC"}
    U = R["A"] * R["B"] + R["B"] * R["C"] + R["C"] * R["A"]
    s = th.order_sign
    if T == 1:
        return s * sgn_power(H * (H + 1) // 2) * D_theta(th, lamfull) / (2 * U ** (H + 1)), U
    L = L_of_t(th, t, lamfull)
    return (s * sgn_power(H * (H + 1) // 2 + 1) * R["A"] * R["B"] * R["C"] * Q_form(th, lamfull, L)
            / U ** (H + 2)), Fr(1)


def c0(pqr, H):
    """T = 0: prod_P (H - n_P)! / ((H+1)! prod_P (n_P - 1)!)."""
    num = 1
    den = math.factorial(H + 1)
    for n_ in pqr:
        num *= math.factorial(H - n_)
        den *= math.factorial(n_ - 1)
    return Fr(num, den)


def gamma_half_ratio(pqr, H):
    """T = 1: prod_P Gamma(H + 3/2 - n_P) / (2 Gamma(H + 3/2) prod_P (n_P - 1)!) = c1 * pi."""
    # Gamma(k + 1/2) = (2k)! / (4^k k!) sqrt(pi), k >= 0;  Gamma(1/2 - k) = (-4)^k k! / (2k)! sqrt(pi)
    def g(k):  # Gamma(k + 1/2) / sqrt(pi)
        if k >= 0:
            return Fr(math.factorial(2 * k), 4 ** k * math.factorial(k))
        m = -k
        return Fr((-4) ** m * math.factorial(m), math.factorial(2 * m))
    num = Fr(1)
    for n_ in pqr:
        num *= g(H + 1 - n_)
    den = 2 * g(H + 1)
    for n_ in pqr:
        den *= math.factorial(n_ - 1)
    return num / den


def kernel_C3_closed(th, lamfull, H, T):
    """K[C3] at z = 0 in the graph's edge order: Fraction (T = 0) or coefficient of pi (T = 1)."""
    if T == 0:
        return th.order_sign * sgn_power(H * (H + 1) // 2 + 1) * c0(th.pqr, H) * Q_form(th, lamfull, Lbar(th, lamfull))
    return th.order_sign * sgn_power(H * (H + 1) // 2) * gamma_half_ratio(th.pqr, H) * D_theta(th, lamfull)


def sgn_detBC(nv, edges, root, C):
    B, g = incidence(nv, edges, root)
    rows = [[B[e][k] for e in range(len(edges))] for k in range(len(g))] + [list(r) for r in C]
    d = det(rows)
    assert d != 0
    return 1 if d > 0 else -1


def kernel_htk_closed(th, lamfull, H, T):
    """K of htk (positively oriented cycle basis for the grounding at th.root)."""
    D = H + T
    s = sgn_detBC(th.nv, th.edges, th.root, th.C3()) ** D
    return s * kernel_C3_closed(th, lamfull, H, T)


# ---------------------------------------------------------------------------- chart (w, x)


def chart_t(th, w, x):
    """w = (w_A, w_B, w_C) (sum 1), x = {P: [x_1..x_n]} (sum 1): t_{P_j} = x^P_j / w_P."""
    t = [None] * len(th.edges)
    for i, P in enumerate("ABC"):
        for (e, _), xj in zip(th.paths[P], x[P]):
            t[e] = Fr(xj) / Fr(w[i])
    return t


def chart_jacobian(th, w, x):
    """det[t; dt/dxi_1; ...] for xi = (w_A, w_B, x^A_2..x^A_p, x^B_2.., x^C_2..), w_C = 1 - w_A - w_B,
    x^P_1 = 1 - sum_{j>=2} x^P_j."""
    Ne = len(th.edges)
    t = chart_t(th, w, x)
    rows = [t]
    # w_A, w_B
    for a in range(2):
        v = [Fr(0)] * Ne
        for i, P in enumerate("ABC"):
            dw = Fr(1) if i == a else (Fr(-1) if i == 2 else Fr(0))
            for (e, _), xj in zip(th.paths[P], x[P]):
                v[e] = -Fr(xj) / Fr(w[i]) ** 2 * dw
        rows.append(v)
    for i, P in enumerate("ABC"):
        es = [e for e, _ in th.paths[P]]
        for j in range(1, len(es)):
            v = [Fr(0)] * Ne
            v[es[j]] = 1 / Fr(w[i])
            v[es[0]] = -1 / Fr(w[i])
            rows.append(v)
    return det(rows)


def closed_chart_density(th, w, x, lamfull, H, T):
    """g(w, x) with K[C3] = int g |d^2w dx| (graph edge order): (coeff, radicand) with
    g = coeff * radicand^{1/2} for T = 1 (radicand = w_A w_B w_C) and radicand 1 for T = 0."""
    s = th.order_sign
    if T == 1:
        mono = Fr(1)
        for i, P in enumerate("ABC"):
            mono *= Fr(w[i]) ** (H - th.n[P])
        return s * sgn_power(H * (H + 1) // 2) * D_theta(th, lamfull) / 2 * mono, Fr(w[0]) * Fr(w[1]) * Fr(w[2])
    L = [Fr(xx) for xx in lamfull[th.N]]
    for P in "ABC":
        xs = x[P]
        for k, u in enumerate(th.interior[P], start=1):
            sig = sum(Fr(xx) for xx in xs[k:])
            for i in range(H):
                L[i] += sig * lamfull[u][i]
    mono = Fr(1)
    for i, P in enumerate("ABC"):
        mono *= Fr(w[i]) ** (H - th.n[P])
    return s * sgn_power(H * (H + 1) // 2 + 1) * mono * Q_form(th, lamfull, L), Fr(1)


# ---------------------------------------------------------------------------- displacement


def X_map(th, w, x, z):
    """X_v (v != N') with -sum_e y_e . z_e = sum_v lam_v . X_v (root N'); z: list over edges of
    H-vectors in the edge orientation.  Z_P = sum_j s_e z_e along P, W = sum_P w_P Z_P,
    X_N = -W, X_{P_k} = sigma^P_k (Z_P - W) - sum_{j>k} s z_{P_j}."""
    H = len(z[0])
    zp = {P: [[s * Fr(c) for c in z[e]] for e, s in th.paths[P]] for P in "ABC"}
    Z = {P: [sum(zz[i] for zz in zp[P]) for i in range(H)] for P in "ABC"}
    W = [sum(Fr(w[k]) * Z[P][i] for k, P in enumerate("ABC")) for i in range(H)]
    X = {th.N: [-c for c in W]}
    for P in "ABC":
        xs = x[P]
        for k, u in enumerate(th.interior[P], start=1):
            sig = sum(Fr(xx) for xx in xs[k:])
            tail = [sum(zz[i] for zz in zp[P][k:]) for i in range(H)]
            X[u] = [sig * (Z[P][i] - W[i]) - tail[i] for i in range(H)]
    X[th.Np] = [Fr(0)] * H
    return X


# ---------------------------------------------------------------------------- symbolic expansion


def laplace_expansion(colspec, H):
    """colspec: list of (name, kind) with kind in {'top', 'bot', 'both'}; the 2H x 2H determinant
    det[col_1..col_2H], col = (v;0), (0;v), (v;v).  Returns [(sign, top names, bottom names)]."""
    n2 = len(colspec)
    assert n2 == 2 * H
    forced_top = [j for j, (_, k) in enumerate(colspec) if k == "top"]
    forced_bot = [j for j, (_, k) in enumerate(colspec) if k == "bot"]
    both = [j for j, (_, k) in enumerate(colspec) if k == "both"]
    need = H - len(forced_top)
    out = []
    if need < 0 or need > len(both):
        return out
    for Bp in itertools.combinations(both, need):
        S = sorted(forced_top + list(Bp))
        Sc = sorted(j for j in range(n2) if j not in S)
        sign = sgn_power(sum(range(1, H + 1)) + sum(j + 1 for j in S))
        out.append((sign, [colspec[j][0] for j in S], [colspec[j][0] for j in Sc]))
    return out


def canon_bracket(names):
    """Sort a bracket of vector names; return (sign, tuple) or (0, None) for a repeated name."""
    if len(set(names)) < len(names):
        return 0, None
    order = sorted(range(len(names)), key=lambda i: names[i])
    return perm_sign(order), tuple(names[i] for i in order)


def vertex_names(th):
    nm = {th.N: "N", th.Np: "Np"}
    for P in "ABC":
        for k, u in enumerate(th.interior[P], start=1):
            nm[u] = "%s%d" % (P, k)
    return nm


def expand_T1(th, H):
    """D_theta = sum coef [top][bottom], names of interior vertices: list of (Fraction, top, bot)."""
    nm = vertex_names(th)
    spec = [(nm[u], "top") for u in th.interior["A"]] + [(nm[u], "both") for u in th.interior["B"]] + \
           [(nm[u], "bot") for u in th.interior["C"]]
    terms = {}
    for sgn, top, bot in laplace_expansion(spec, H):
        s1, a = canon_bracket(top)
        s2, b = canon_bracket(bot)
        if s1 == 0 or s2 == 0:
            continue
        key = (a, b)
        terms[key] = terms.get(key, Fr(0)) + sgn * s1 * s2
    return [(c, list(a), list(b)) for (a, b), c in sorted(terms.items()) if c != 0]


def expand_T0(th, H):
    """Q(bar L) = sum coef [top][bottom] in the names of the non-root sources (root N'):
    bar L = lam_N + sum mu_u lam_u."""
    nm = vertex_names(th)
    mu = mu_weights(th)
    spec = [("@L", "top"), ("@L", "bot")] + [(nm[u], "top") for u in th.interior["A"]] + \
           [(nm[u], "both") for u in th.interior["B"]] + [(nm[u], "bot") for u in th.interior["C"]]
    Lterms = [(Fr(1), "N")] + [(mu[u], nm[u]) for P in "ABC" for u in th.interior[P]]
    terms = {}
    for sgn, top, bot in laplace_expansion(spec, H):
        for ca, na in Lterms:
            ta = [na if x == "@L" else x for x in top]
            s1, a = canon_bracket(ta)
            if s1 == 0:
                continue
            for cb, nb in Lterms:
                tb = [nb if x == "@L" else x for x in bot]
                s2, b = canon_bracket(tb)
                if s2 == 0:
                    continue
                key = (a, b)
                terms[key] = terms.get(key, Fr(0)) + sgn * s1 * s2 * ca * cb
    return [(c, list(a), list(b)) for (a, b), c in sorted(terms.items()) if c != 0]


def eval_terms(terms, th, lamfull):
    nm = vertex_names(th)
    vec = {nm[v]: lamfull[v] for v in nm}
    tot = Fr(0)
    for c, a, b in terms:
        tot += c * det([[vec[x][i] for x in a] for i in range(len(a))]) * det([[vec[x][i] for x in b] for i in range(len(b))])
    return tot


def theta_list(D):
    return [(a, b, 2 * D + 1 - a - b) for a in range(1, D + 1) for b in range(a, D + 1)
            if b <= 2 * D + 1 - a - b <= D]
