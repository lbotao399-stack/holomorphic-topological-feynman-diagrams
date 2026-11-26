"""Relations between the reduced kernel K of htk and the kernels of the source notes
(CONVENTIONS.md, Sections 3 and 4), and exact checks of the sign rules.

    S1 (Li-Liu-Chang, C^2):   I^{S1}_Gamma = (-1)^{n(n+1)/2} K_Gamma          (eq-operatope)
    S2 (6d note, C^3):        I^{S2}_Gamma = (-1)^{n(n+1)/2} K_Gamma          (eq-kernel-def)
    S3 (11d note, R x C^5):   K^{S3}_Gamma = varsigma_Gamma K_Gamma,
                              varsigma_Gamma = (-1)^{5n(5n+1)/2} = (-1)^{n(n+1)/2}   (eq-master)
    S4 (check_cd):            I^{S4}_{Gamma,T} = int ^_{c notin T} d^d y_c
                                               = (-1)^{nd + C(n,2)C(d,2)} sigma_T^d K_Gamma
"""
import random
from fractions import Fraction as Fr

from .linalg import det, wedge_all, one_form, basic_coefficient

__all__ = ["sign_n", "sign_dy", "S1_kernel", "S2_kernel", "S3_kernel", "S4_integral",
           "K_from_dy", "raw_factor", "extraction_sign", "extraction_sign_formula", "varsigma_S3",
           "edge_reversal_sign", "edge_permutation_sign", "vertex_relabelling_sign"]


def sign_n(n):
    """(-1)^{n(n+1)/2}."""
    return -1 if (n * (n + 1) // 2) % 2 else 1


def sign_dy(n, H):
    """prod_{i=1..H} prod_{r=1..n} psi^i_r = sign_dy(n, H) U^H wedge_{r=1..n} d^H y_{c_r} for the
    fundamental cycle matrix of a tree with cotree (c_1 < ... < c_n):
    (-1)^{nH + C(n,2) C(H,2)}."""
    return -1 if (n * H + (n * (n - 1) // 2) * (H * (H - 1) // 2)) % 2 else 1


def S1_kernel(K, n):
    """Kernel I_Gamma of S1 (H = 2, T = 0), eq. (eq-operatope), from the reduced kernel K."""
    return sign_n(n) * K


def S2_kernel(K, n):
    """Kernel I_Gamma of S2 (H = 3, T = 0), eq. (eq-kernel-def), from the reduced kernel K."""
    return sign_n(n) * K


def S3_kernel(K, n):
    """Kernel K_Gamma of S3 (H = 5, T = 1), eq. (eq-master), from the reduced kernel K."""
    return sign_n(n) * K


def S4_integral(K, G, cotree, H):
    """int_{omega_m} wedge_{c in cotree} d^H y_c (S4 check_cd.py; cotree in increasing order)."""
    tree = [e for e in range(G.ne) if e not in set(cotree)]
    if list(cotree) != sorted(cotree):
        raise ValueError("cotree must be in increasing order")
    return sign_dy(G.loops, H) * G.sigma(tree) ** H * K


def K_from_dy(I_dy, G, cotree, H, T=0):
    """Reduced kernel from int rho^T ^ wedge_{c in cotree} d^H y_c."""
    tree = [e for e in range(G.ne) if e not in set(cotree)]
    return sign_dy(G.loops, H) * G.sigma(tree) ** (H + T) * I_dy


def raw_factor(source, n):
    """Factor f with (position-space graph integral of the source) = f * (kernel of the source):
    S1 eq-1-7: (-8 pi^2)^{-n};  S2 eq-B-norm: (-1/(8 pi^3))^n;  S3 eq-master: (-i/(64 pi^{11/2}))^n."""
    import sympy as sp
    if source == "S1":
        return (-8 * sp.pi ** 2) ** (-n)
    if source == "S2":
        return (-1 / (8 * sp.pi ** 3)) ** n
    if source == "S3":
        return (-sp.I / (64 * sp.pi ** sp.Rational(11, 2))) ** n
    raise ValueError("source must be 'S1', 'S2' or 'S3'")


# ---------------------------------------------------------------------- transformation rules

def edge_reversal_sign(D, k=1):
    """K -> (-1)^{D k} K when k edges are reversed (and z_e -> -z_e on them)."""
    return -1 if (D * k) % 2 else 1


def edge_permutation_sign(D, perm_sign_):
    """K -> sgn^{D+1} K under a permutation of the edge order."""
    return perm_sign_ if (D + 1) % 2 else 1


def vertex_relabelling_sign(D, perm_sign_):
    """K -> sgn^D K under a relabelling of all vertices (the root is the vertex labelled last)."""
    return perm_sign_ if D % 2 else 1


# ---------------------------------------------------------------------- extraction identity (S3 eq-extraction)

def extraction_sign_formula(G, D):
    """Predicted sign in
        wedge_e (theta^0_e ^ ... ^ theta^{D-1}_e)
          = sign * prod_a (c_a)^n [prod_v dX^0_v ... prod_v dX^{D-1}_v] ^ prod_a prod_r psi^a_r + ...,
    theta^a_e = dX^a_e - c_a J^a_e dt_e,  psi^a_r = sum_e C_re J^a_e dt_e, for a positively oriented
    cycle matrix:  (-1)^{C(|I|,2) C(D,2) + n D + (|V|-1) n C(D,2)}."""
    N, n, V = G.ne, G.loops, G.nv
    cD = D * (D - 1) // 2
    return -1 if ((N * (N - 1) // 2) * cD + n * D + (V - 1) * n * cD) % 2 else 1


def extraction_sign(G, D, C=None, seed=0, c=None, k=None):
    """Exact evaluation (random rational positions and Schwinger times) of the ratio
        F / (prod_a c_a^n D_coef),
        dg ^ wedge_e wedge_a theta^a_e = F (-1)^{dim X} g vol_X ^ vol_t,   g = sum_e t_e,
        prod_a prod_r psi^a_r = D_coef iota_E vol_t,
    with J^a_e = X^a_e / (k_a t_e), theta^a_e = dX^a_e - c_a J^a_e dt_e, psi^a_r = sum_e C_re J^a_e dt_e
    (X^a_e = X^a_{e(0)} - X^a_{e(1)}, X_root = 0; default c_a = k_a = 1).  The result is +-1 and
    equals extraction_sign_formula(G, D) for the positively oriented cycle matrix.
    S3 eq-extraction (H = 5, T = 1): c = (1, 4, 4, 4, 4, 4), k = (2, 4, 4, 4, 4, 4); the result is
    varsigma^ext (S3 Appendix A, item 3: -1 for C_7, +1 for the five theta graphs)."""
    rng = random.Random(seed)
    C = G.cycle_matrix() if C is None else C
    N, g, n = G.ne, G.nv - 1, len(C)
    c = [Fr(1)] * D if c is None else [Fr(x) for x in c]
    k = [Fr(1)] * D if k is None else [Fr(x) for x in k]
    B = G.incidence()
    while True:
        X = [[Fr(rng.randint(-9, 9), rng.randint(1, 5)) for _ in range(g)] for _ in range(D)]
        t = [Fr(rng.randint(1, 20), rng.randint(1, 7)) for _ in range(N)]
        Xe = [[sum(B[e][v] * X[a][v] for v in range(g)) for e in range(N)] for a in range(D)]
        J = [[Xe[a][e] / (k[a] * t[e]) for e in range(N)] for a in range(D)]
        nvar = D * g + N
        rows = [[Fr(0)] * (D * g) + [Fr(1)] * N]
        for e in range(N):
            for a in range(D):
                row = [Fr(0)] * nvar
                for v in range(g):
                    row[a * g + v] = Fr(B[e][v])
                row[D * g + e] = -c[a] * J[a][e]
                rows.append(row)
        F = det(rows) / sum(t)
        if (D * g) % 2:
            F = -F
        psis = [one_form([C[r][e] * J[a][e] for e in range(N)]) for a in range(D) for r in range(n)]
        Dc = basic_coefficient(wedge_all(psis), t, 0)
        if Dc != 0:
            break
    pref = Fr(1)
    for a in range(D):
        pref *= c[a] ** n
    ratio = F / (pref * Dc)
    if ratio not in (1, -1):
        raise ArithmeticError("extraction ratio %s is not a sign" % ratio)
    return int(ratio)


def varsigma_S3(n, ext=None):
    """varsigma_Gamma = (-1)^n (-1)^{(5n+1) 5n / 2} varsigma^ext of S3 Theorem (Master formula);
    with varsigma^ext = (-1)^n (extraction_sign_formula at D = 6, |I| = 6n+1, |V| = 5n+2) this is
    (-1)^{5n(5n+1)/2} = (-1)^{n(n+1)/2}."""
    if ext is None:
        ext = -1 if n % 2 else 1
    return (-1 if n % 2 else 1) * (-1 if ((5 * n + 1) * 5 * n // 2) % 2 else 1) * ext
