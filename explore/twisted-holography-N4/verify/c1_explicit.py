"""C1 literal checks: fan colour tensor from structure constants, no ribbon reduction.

L(D,E,F) = sum_sigma sum tr(t_{A_1}...t_{A_m} F) C_{A_s(1) D R_1} C_{A_s(2) R_1 R_2} ... C_{A_s(m) R_n E},
m = n+1, C_{XYZ} = -i tr([t_X,t_Y] t_Z)/tau.  Claim: L = (-i)^(n+1) tau^n P_n(N) [tr(DEF)+(-1)^(n+1) tr(DFE)].
"""
import sys
from fractions import Fraction
from itertools import permutations, product
import numpy as np


def perm_sum(T, M):
    """sum_sigma sum_a T[a] M[a o sigma] for two m-index tensors."""
    m = T.ndim
    total = 0
    for p in permutations(range(m)):
        total += int(np.sum(T * M.transpose(p)))
    return total


# ---------- SU(2), t_i = sigma_i/2, literal Levi-Civita ----------
def su2_component():
    s = [np.array([[0, 1], [1, 0]], dtype=complex), np.array([[0, -1j], [1j, 0]]),
         np.array([[1, 0], [0, -1]], dtype=complex)]
    # tr(t_A t_B t_C t_3) * 16 is a Gaussian integer; check structure constants = epsilon
    def eps(i, j, k):
        return 0 if len({i, j, k}) < 3 else (1 if (i, j, k) in ((0, 1, 2), (1, 2, 0), (2, 0, 1)) else -1)
    for i, j, k in product(range(3), repeat=3):
        c = -1j * np.trace((s[i] @ s[j] - s[j] @ s[i]) @ s[k]) / 8 / 0.5   # t=s/2, tau=1/2
        assert abs(c - eps(i, j, k)) < 1e-12
    total = Fraction(0)
    for A, B, C in product(range(3), repeat=3):
        tr16 = np.trace(s[A] @ s[B] @ s[C] @ s[2])
        assert abs(tr16.imag) < 1e-12 and abs(tr16.real - round(tr16.real)) < 1e-12
        tr = Fraction(int(round(tr16.real)), 16)
        for x, y, z in permutations((A, B, C)):
            ch = sum(eps(x, 0, R) * eps(y, R, S) * eps(z, S, 1) for R, S in product(range(3), repeat=2))
            total += tr * ch
    return total


# ---------- SU(3) Gell-Mann, n = 3, D = E = t_1, F = diag(1,1,-2) ----------
def gellmann_n3():
    l = np.zeros((8, 3, 3), dtype=complex)
    l[0][0, 1] = l[0][1, 0] = 1
    l[1][0, 1], l[1][1, 0] = -1j, 1j
    l[2][0, 0], l[2][1, 1] = 1, -1
    l[3][0, 2] = l[3][2, 0] = 1
    l[4][0, 2], l[4][2, 0] = -1j, 1j
    l[5][1, 2] = l[5][2, 1] = 1
    l[6][1, 2], l[6][2, 1] = -1j, 1j
    l[7] = np.diag([1, 1, -2])            # sqrt(3) * lambda_8
    w = np.array([3] * 7 + [1])           # 3 * (orthonormal weight): 1 for A<8, 1/3 for lambda_8'
    # C'_{XYZ} = -i tr([t,t]t)/tau with t = l/2, tau = 1/2  ->  -i tr([l,l]l)/4
    C4 = np.zeros((8, 8, 8))              # 4*C'
    for x, y, z in product(range(8), repeat=3):
        v = -1j * np.trace((l[x] @ l[y] - l[y] @ l[x]) @ l[z])
        assert abs(v.imag) < 1e-9
        C4[x, y, z] = round(v.real)
    D = l[0]                              # 2*t_1
    F = np.diag([1, 1, -2]).astype(complex)
    # 4*C'(x, D, R) and 4*C'(x, R, E) with D = E = 2 t_1 -> index 0 of C4 (times 1)
    K1 = C4[:, 0, :]                      # C(x, D, R)
    K3 = C4[:, :, 0]                      # C(x, R, E)
    # chain: sum_R C(x1,D,R1) w C(x2,R1,R2) w C(x3,R2,R3) w C(x4,R3,E)
    Y = np.einsum('ar,r->ar', K1, w)
    Y = np.einsum('ar,brs,s->abs', Y, C4, w)
    Y = np.einsum('abr,crs,s->abcs', Y, C4, w)
    M = np.einsum('abcr,dr->abcd', Y, K3)
    # spoke weights and trace: tr(t t t t F) = tr(l l l l F)/16
    T = np.einsum('aij,bjk,ckl,dlm,mi->abcd', l, l, l, l, F)
    assert np.max(np.abs(T.imag - np.round(T.imag))) < 1e-9
    Wt = np.einsum('a,b,c,d->abcd', w, w, w, w)
    Tr = np.round(T.real).astype(np.int64) * Wt
    Ti = np.round(T.imag).astype(np.int64) * Wt
    Mi = np.round(M).astype(np.int64)
    re, im = perm_sum(Tr, Mi), perm_sum(Ti, Mi)
    # scale: T/16, each C' = C4/4 (four of them; index 0 of C4 is the generator t_1 itself),
    # weights: 3^7 for 7 contractions
    scale = Fraction(1, 16) * Fraction(1, 4**4) * Fraction(1, 3**7)
    return Fraction(re) * scale, Fraction(im) * scale


# ---------- general SU(N): integral basis, scaled dual basis, literal R sums ----------
def sl_basis(N):
    e = []
    for i in range(N):
        for j in range(N):
            if i != j:
                x = np.zeros((N, N), dtype=np.int64); x[i, j] = 1; e.append(x)
    for i in range(N - 1):
        x = np.zeros((N, N), dtype=np.int64); x[i, i] = 1; x[i + 1, i + 1] = -1; e.append(x)
    e = np.array(e)
    G = np.einsum('aij,bji->ab', e, e)
    # N * G^{-1} is integral for the Cartan block; off-diagonal block G = permutation
    Ginv = np.array([[Fraction(int(v)) for v in row] for row in G])
    import sympy
    Gi = sympy.Matrix(G.tolist()).inv() * N
    assert all(v.is_integer for v in Gi)
    Gi = np.array(Gi.tolist(), dtype=np.int64)
    dual = np.einsum('ab,bij->aij', Gi, e)            # N * dual basis
    assert (np.einsum('aij,bji->ab', e, dual) == N * np.eye(len(e), dtype=np.int64)).all()
    return e, dual


def comm(x, y):
    return x @ y - y @ x


def literal_SN(n, N, D, E, F):
    """S' = L / ((-i)^(n+1) tau^n), computed with explicit R sums. Returns Fraction."""
    m = n + 1
    e, du = sl_basis(N)
    d = len(e)
    # K1[x,r] = tr([du_x, D] e_r); K2[x,r,s] = tr([du_x, du_r] e_s); K3[x,s] = tr([du_x, du_s] E)
    cD = np.einsum('xij,jk->xik', du, D) - np.einsum('ij,xjk->xik', D, du)
    K1 = np.einsum('xik,rki->xr', cD, e)
    cc = np.einsum('xij,rjk->xrik', du, du) - np.einsum('rij,xjk->xrik', du, du)
    K2 = np.einsum('xrik,ski->xrs', cc, e)
    K3 = np.einsum('xsik,ki->xs', cc, E)
    Y = K1
    for _ in range(m - 2):
        Y = np.tensordot(Y, K2, axes=([Y.ndim - 1], [1]))      # (..., r) x (x, r, s) -> (..., x, s)
    M = np.tensordot(Y, K3, axes=([Y.ndim - 1], [1]))          # (..., x)
    # T[a1..am] = tr(e_a1 ... e_am F)
    P = e
    for _ in range(m - 2):
        P = np.einsum('...ij,ajk->...aik', P, e)
    eF = np.einsum('ajk,ki->aji', e, F)                        # (e_a F)_{ji}
    T = np.tensordot(P, eF, axes=([P.ndim - 2, P.ndim - 1], [2, 1]))
    val = perm_sum(T, M)
    return Fraction(val, N**(2 * n + 1))


def P_claim(n, N):
    return {2: N**2 + 6, 3: N**3 + 12 * N, 4: N**4 + 36 * N**2 + 120, 5: N**5 + 63 * N**3 + 360 * N}[n]


def random_traceless(N, rng):
    x = rng.integers(-2, 3, size=(N, N)).astype(np.int64)
    x[N - 1, N - 1] = -np.trace(x[:N - 1, :N - 1])
    return x


if __name__ == '__main__':
    v = su2_component()
    print('SU(2) (D,E,F)=(1,2,3) literal Levi-Civita sum:', v, '  claim -5/4:', v == Fraction(-5, 4))
    re, im = gellmann_n3()
    print('SU(3) Gell-Mann n=3, D=E=t1, F=diag(1,1,-2):', re, '+ i*', im, '  claim 63/8:', re == Fraction(63, 8) and im == 0)
    rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 7)
    cases = [(2, 2), (2, 3), (2, 4), (2, 5), (3, 3), (3, 4), (3, 5), (4, 2), (4, 3), (4, 4), (5, 3), (5, 4)]
    if len(sys.argv) > 2:
        cases = [tuple(map(int, c.split(','))) for c in sys.argv[2:]]
    for n, N in cases:
      for trial in range(int(__import__('os').environ.get('TRIALS', '1'))):
          while True:
              D, E, F = (random_traceless(N, rng) for _ in range(3))
              tdef = int(np.trace(D @ E @ F)); tdfe = int(np.trace(D @ F @ E))
              basis = tdef + (-1)**(n + 1) * tdfe
              if basis:
                  break
          S = literal_SN(n, N, D, E, F)
          ok = (S == P_claim(n, N) * basis)
          print(f'n={n} N={N}: S={S}, tr(DEF)+(-1)^(n+1)tr(DFE)={basis}, ratio={S / basis if basis else None}, '
                f'claim P={P_claim(n, N)}, match={ok}', flush=True)
