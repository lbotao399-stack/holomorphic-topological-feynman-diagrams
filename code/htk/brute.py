"""Layer b): exact brute-force evaluation at rational Schwinger times.  No reductions.

Circuit:  Ohm  t_e y_e = phi_{e(0)} - phi_{e(1)},  phi_root = 0;
          Kirchhoff  sum_{e(0)=v} y_e - sum_{e(1)=v} y_e = lambda_v  (v != root).
Integrand of the reduced kernel (CONVENTIONS.md, Section 2):

    exp(-sum_e y_e.z_e)  U^{-H}  rho_n(M)^T  ^  prod_{i=1..H} prod_{r=1..n} psi^i_r  =  D(t) iota_E prod_e dt_e,

    psi^i_r = sum_e C_re y^i_e dt_e,   M = C t C^t,   U = det M,
    rho_n(M) = det(M)^{-1/2} Pf(-1/2 dM M^{-1} dM).
"""
import math
from fractions import Fraction as Fr

from .linalg import (det, solve, inverse, one_form, wedge, wedge_all, form_add, form_scale,
                     pfaffian_forms, basic_coefficient)

__all__ = [
    "Density", "sources", "laplacian", "response", "potentials", "currents", "period_matrix",
    "symanzik_U", "psi_forms", "topological_form", "density", "density_T0_det", "jacobian_density",
    "exponent", "monomial", "integrand", "check_kirchhoff", "theta_D",
]


class Density:
    """The number  coeff * radicand**(-1/2) * exp(expo)  with coeff, radicand > 0, expo rational.
    T = 0: radicand = 1, expo = 0.  T = 1: radicand = U = det M; expo = sum_e t_e wbar_e^2 is
    nonzero only with topological sources kappa."""
    __slots__ = ("coeff", "radicand", "expo")

    def __init__(self, coeff, radicand=1, expo=0):
        self.coeff = Fr(coeff)
        self.radicand = Fr(radicand)
        self.expo = Fr(expo)
        if self.radicand <= 0:
            raise ValueError("radicand must be positive")

    def is_rational(self):
        return self.radicand == 1 and self.expo == 0

    def rational(self):
        if not self.is_rational():
            raise ValueError("not rational: %r" % self)
        return self.coeff

    def sign(self):
        return (self.coeff > 0) - (self.coeff < 0)

    def square(self):
        """(value * exp(-expo))^2 = coeff^2 / radicand, a rational number."""
        return self.coeff ** 2 / self.radicand

    def scaled(self, c):
        return Density(self.coeff * Fr(c), self.radicand, self.expo)

    def same_as(self, other):
        return (self.sign() == other.sign() and self.square() == other.square()
                and self.expo == other.expo)

    def __float__(self):
        return float(self.coeff) / math.sqrt(float(self.radicand)) * math.exp(float(self.expo))

    def __repr__(self):
        s = repr(self.coeff)
        if self.radicand != 1:
            s += " * (%r)**(-1/2)" % self.radicand
        if self.expo != 0:
            s += " * exp(%r)" % self.expo
        return "Density(%s)" % s


# ---------------------------------------------------------------------- circuit

def sources(G, lam, H=None):
    """Normalize sources to a list over the non-root vertices of H-vectors of Fractions.
    lam: list of nv-1 vectors (non-root vertices), or of nv vectors (must sum to zero), or a dict
    {vertex: vector} (missing non-root vertices carry 0; a root entry is checked).  H = 0: lam = None."""
    g = G.nv - 1
    if lam is None:
        return [[] for _ in range(g)]
    if isinstance(lam, dict):
        if H is None:
            H = len(next(iter(lam.values()))) if lam else 0
        out = [[Fr(x) for x in lam.get(v, [0] * H)] for v in range(g)]
        if G.root in lam:
            if any(sum(out[v][i] for v in range(g)) + Fr(lam[G.root][i]) != 0 for i in range(H)):
                raise ValueError("sources do not sum to zero")
        return out
    lam = [list(v) for v in lam]
    if len(lam) == G.nv:
        Hh = len(lam[0])
        if any(sum(Fr(lam[v][i]) for v in range(G.nv)) != 0 for i in range(Hh)):
            raise ValueError("sources do not sum to zero")
        lam = lam[:g]
    if len(lam) != g:
        raise ValueError("need sources at the %d non-root vertices" % g)
    return [[Fr(x) for x in v] for v in lam]


def laplacian(G, t):
    """Lambda = B^t t^{-1} B  ((nv-1) x (nv-1));  det Lambda = U / prod t_e."""
    g = G.nv - 1
    L = [[Fr(0)] * g for _ in range(g)]
    for (a, b), te in zip(G.edges, t):
        c = 1 / Fr(te)
        if a < g:
            L[a][a] += c
        if b < g:
            L[b][b] += c
        if a < g and b < g:
            L[a][b] -= c
            L[b][a] -= c
    return L


def response(G, t):
    """c[e][v] = d y_e / d lambda_v (scalar response), v over the non-root vertices:
    y_e = sum_v c[e][v] lambda_v."""
    g = G.nv - 1
    Linv = inverse(laplacian(G, t))
    out = []
    for (a, b), te in zip(G.edges, t):
        ra = Linv[a] if a < g else [Fr(0)] * g
        rb = Linv[b] if b < g else [Fr(0)] * g
        out.append([(x - y) / Fr(te) for x, y in zip(ra, rb)])
    return out


def potentials(G, t, lam):
    """phi[v][i], with phi[root] = 0."""
    lam = sources(G, lam)
    H = len(lam[0]) if lam else 0
    if H == 0:
        return [[] for _ in range(G.nv)]
    phi = solve(laplacian(G, t), lam)
    return phi + [[Fr(0)] * H]


def currents(G, t, lam):
    """Edge currents y[e][i] (H-vectors of Fractions)."""
    lam = sources(G, lam)
    c = response(G, t)
    H = len(lam[0]) if lam else 0
    return [[sum((c[e][v] * lam[v][i] for v in range(G.nv - 1)), Fr(0)) for i in range(H)]
            for e in range(G.ne)]


def check_kirchhoff(G, t, lam, y, C=None):
    """True iff y satisfies the current law with sources lam and the voltage law on the cycles C."""
    lam = sources(G, lam)
    H = len(lam[0]) if lam else 0
    B = G.incidence()
    for v in range(G.nv - 1):
        for i in range(H):
            if sum(B[e][v] * y[e][i] for e in range(G.ne)) != lam[v][i]:
                return False
    C = G.cycle_matrix() if C is None else C
    for row in C:
        for i in range(H):
            if sum(row[e] * Fr(t[e]) * y[e][i] for e in range(G.ne)) != 0:
                return False
    return True


def period_matrix(G, t, C=None):
    """M = C t C^t."""
    C = G.cycle_matrix() if C is None else C
    return [[sum((Fr(ra[e] * rb[e]) * Fr(t[e]) for e in range(G.ne)), Fr(0)) for rb in C] for ra in C]


def symanzik_U(G, t, C=None):
    """U = det M (for a fundamental or any unimodular cycle basis: the spanning-tree polynomial
    sum_T prod_{e not in T} t_e)."""
    return det(period_matrix(G, t, C))


def psi_forms(G, t, y, C):
    """The one-forms psi^i_r = sum_e C_re y^i_e dt_e in the order i = 1..H (outer), r = 1..n."""
    H = len(y[0]) if y else 0
    return [one_form([C[r][e] * y[e][i] for e in range(G.ne)]) for i in range(H) for r in range(len(C))]


def topological_form(G, t, C=None, kappa=None):
    """< prod_{r=1..n} psi^tau_r >  as an n-form (dict), psi^tau_r = sum_e C_re w_e dt_e, with
    < w_e w_e' >_conn = (1/2)[delta_ee'/t_e - (C^t M^{-1} C)_ee'] and mean wbar = Kirchhoff
    solution with sources kappa (S3 eq-ww, Section 6).  For kappa = None:
    Pf(-1/2 dM M^{-1} dM) = det(M)^{1/2} rho_n(M).  Returns (form, sum_e t_e wbar_e^2)."""
    C = G.cycle_matrix() if C is None else C
    n = len(C)
    t = [Fr(x) for x in t]
    Minv = inverse(period_matrix(G, t, C)) if n else []
    dM = [[one_form([C[r][e] * C[a][e] for e in range(G.ne)]) for a in range(n)] for r in range(n)]
    A = [[{} for _ in range(n)] for _ in range(n)]
    for r in range(n):
        for s in range(r + 1, n):
            acc = {}
            for a in range(n):
                for b in range(n):
                    if Minv[a][b] != 0:
                        acc = form_add(acc, form_scale(wedge(dM[r][a], dM[s][b]), -Minv[a][b] / 2))
            A[r][s] = acc
            A[s][r] = form_scale(acc, -1)
    if kappa is None or all(k == 0 for k in kappa):
        return pfaffian_forms(A), Fr(0)
    wbar = [row[0] for row in currents(G, t, [[k] for k in kappa])]
    mean = [one_form([C[r][e] * wbar[e] for e in range(G.ne)]) for r in range(n)]

    def wick(idx):
        if not idx:
            return {0: Fr(1)}
        i0, rest = idx[0], idx[1:]
        out = wedge(mean[i0], wick(rest))
        for j in range(len(rest)):
            term = wedge(A[i0][rest[j]], wick(rest[:j] + rest[j + 1:]))
            out = form_add(out, form_scale(term, -1 if j % 2 else 1))
        return out

    return wick(list(range(n))), sum(te * w * w for te, w in zip(t, wbar))


def density_T0_det(G, t, y, C, e0=None):
    """T = 0:  D = (-1)^{e0} det[C_re y^i_e]_{(i,r); e != e0} / (t_{e0} U^H)   (S3 eq-calD)."""
    n, N = len(C), G.ne
    H = len(y[0]) if y else 0
    if n * H != N - 1:
        raise ValueError("degree mismatch: n H = %d, |I| - 1 = %d" % (n * H, N - 1))
    e0 = N - 1 if e0 is None else e0
    cols = [e for e in range(N) if e != e0]
    Psi = [[C[r][e] * y[e][i] for e in cols] for i in range(H) for r in range(n)]
    U = symanzik_U(G, t, C)
    d = det(Psi)
    return (-d if e0 & 1 else d) / (Fr(t[e0]) * U ** H)


def density(G, t, lam, T=0, kappa=None, C=None, tree=None, check_basic=False):
    """Coefficient D(t) of the integrand at z = 0,
        U^{-H} rho_n(M)^T ^ prod_{i,r} psi^i_r = D(t) iota_E prod_e dt_e,
    as a Density.  lam: sources (see `sources`), T in {0, 1}.
    C: cycle matrix (default G.cycle_matrix(tree), positively oriented).
    kappa (T = 1 only): topological sources at the non-root vertices; then
        D = U^{-H-1/2} exp(sum_e t_e wbar_e^2) <prod_r psi^tau_r> ^ prod psi^i_r / iota_E vol
    (S3 Section 6); this form has weight 1/2 under t -> c t and is used on a gauge slice."""
    t = [Fr(x) for x in t]
    if any(x <= 0 for x in t):
        raise ValueError("Schwinger times must be positive")
    C = G.cycle_matrix(tree) if C is None else C
    lam = sources(G, lam)
    y = currents(G, t, lam)
    H = len(lam[0]) if lam else 0
    n, N = len(C), G.ne
    if n * (H + T) != N - 1:
        raise ValueError("degree mismatch: n (H + T) = %d, |I| - 1 = %d" % (n * (H + T), N - 1))
    if T == 0:
        D = density_T0_det(G, t, y, C)
        if check_basic:
            for e0 in range(N):
                assert density_T0_det(G, t, y, C, e0) == D
        return Density(D)
    if T != 1:
        raise ValueError("T must be 0 or 1")
    top, expo = topological_form(G, t, C, kappa)
    if not top:
        return Density(0, 1, 0)
    form = wedge(top, wedge_all(psi_forms(G, t, y, C)))
    U = symanzik_U(G, t, C)
    coeff = basic_coefficient(form, t, None if check_basic else N - 1)
    return Density(coeff / U ** H, U, expo)


def jacobian_density(G, t, lam, cotree=None):
    """D_dy(t) with  ^_{r=1..n} d^H y_{c_r} = D_dy iota_E prod dt_e,  c_r the cotree edges in the
    given order:  D_dy = (-1)^{N-1} det[ d y^i_{c_r} / d t_j ]_{(r,i); j < N-1} / t_{N-1}.
    Exact derivative:  d y_e / d t_j = b_e^t Lambda^{-1} b_j y_j / (t_e t_j) - delta_ej y_e / t_e.
    (Convention of S4 check_cd.py and of the kernels of S1, S2, up to their prefactors.)"""
    t = [Fr(x) for x in t]
    lam = sources(G, lam)
    H = len(lam[0]) if lam else 0
    cotree = G.cotree() if cotree is None else list(cotree)
    N, g = G.ne, G.nv - 1
    if len(cotree) * H != N - 1:
        raise ValueError("degree mismatch")
    B = G.incidence()
    Linv = inverse(laplacian(G, t))
    y = currents(G, t, lam)
    LB = [[sum((Linv[v][u] * B[j][u] for u in range(g)), Fr(0)) for v in range(g)] for j in range(N)]
    rows = []
    for c in cotree:
        for i in range(H):
            row = []
            for j in range(N - 1):
                val = sum((B[c][v] * LB[j][v] for v in range(g)), Fr(0)) * y[j][i] / (t[c] * t[j])
                if j == c:
                    val -= y[c][i] / t[c]
                row.append(val)
            rows.append(row)
    d = det(rows)
    return (d if (N - 1) % 2 == 0 else -d) / t[N - 1]


# ---------------------------------------------------------------------- insertions

def exponent(y, z):
    """- sum_e y_e . z_e  for currents y[e][i] and displacements z[e][i] (z: list or dict {e: vector})."""
    if isinstance(z, dict):
        return -sum((sum((yy * Fr(zz) for yy, zz in zip(y[e], ze)), Fr(0)) for e, ze in z.items()), Fr(0))
    return -sum((sum((yy * Fr(zz) for yy, zz in zip(ye, ze)), Fr(0)) for ye, ze in zip(y, z)), Fr(0))


def monomial(y, powers):
    """prod_{(e,i)} (y^i_e)^{powers[(e,i)]}  (i is 0-based)."""
    out = Fr(1)
    for (e, i), p in powers.items():
        out *= y[e][i] ** p
    return out


def integrand(G, t, lam, T=0, powers=None, z=None, **kw):
    """D(t) prod (y^i_e)^{powers[(e,i)]} exp(-sum_e y_e.z_e) as a Density (the rational exponent
    -sum y.z is added to .expo).  kw: kappa, C, tree as in `density`."""
    D = density(G, t, lam, T, **kw)
    y = currents(G, [Fr(x) for x in t], lam)
    c = D.coeff * (monomial(y, powers) if powers else 1)
    e = D.expo + (exponent(y, z) if z is not None else 0)
    return Density(c, D.radicand, e)


def theta_D(lam, interior):
    """D_theta of S3 eq. (R-Ktheta): det [[Lam_A, Lam_B, 0], [0, Lam_B, Lam_C]], Lam_P the
    H x (n_P - 1) matrix of the sources at the interior vertices of the path P (path order)."""
    H = len(lam[0])
    A, Bv, Cv = interior
    rows = []
    for i in range(H):
        rows.append([lam[v][i] for v in A] + [lam[v][i] for v in Bv] + [0] * len(Cv))
    for i in range(H):
        rows.append([0] * len(A) + [lam[v][i] for v in Bv] + [lam[v][i] for v in Cv])
    return det(rows) if rows else Fr(1)
