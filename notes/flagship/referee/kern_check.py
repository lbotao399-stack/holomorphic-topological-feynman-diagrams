"""Referee check of K_133, K_223 at (H,T)=(2,1) from the definition (flagship eq-K), own code.
Exact pointwise density (Fractions) at random rational t and random rational sources;
then the Schwinger integral of U^{-7/2} by a separate mpmath quadrature.
No import of htk or any project code."""
from fractions import Fraction as Fr
import random, itertools, sys
import mpmath as mp

random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 20261121)

# ---------- exterior algebra on dt_0..dt_{N-1}: dict {sorted tuple: coeff}
def wedge(a, b):
    out = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            if set(ka) & set(kb):
                continue
            lst = list(ka) + list(kb)
            # sign of sorting permutation
            sgn = 1
            arr = lst[:]
            for i in range(len(arr)):
                for j in range(len(arr) - 1 - i):
                    if arr[j] > arr[j + 1]:
                        arr[j], arr[j + 1] = arr[j + 1], arr[j]
                        sgn = -sgn
            k = tuple(arr)
            out[k] = out.get(k, 0) + sgn * va * vb
    return {k: v for k, v in out.items() if v != 0}

def add(a, b, cb=1):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0) + cb * v
    return {k: v for k, v in out.items() if v != 0}

def scal(a, c):
    return {k: c * v for k, v in a.items() if c * v != 0}

def solve(A, b):
    """exact Gaussian elimination, A square list of lists of Fractions"""
    n = len(A)
    M = [list(map(Fr, A[i])) + [Fr(b[i])] for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[p] = M[p], M[c]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c] / M[c][c]
                M[r] = [M[r][k] - f * M[c][k] for k in range(n + 1)]
    return [M[i][n] / M[i][i] for i in range(n)]

def det(A):
    n = len(A)
    M = [list(map(Fr, r)) for r in A]
    d = Fr(1)
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != 0), None)
        if p is None:
            return Fr(0)
        if p != c:
            M[c], M[p] = M[p], M[c]
            d = -d
        d *= M[c][c]
        for r in range(c + 1, n):
            f = M[r][c] / M[c][c]
            M[r] = [M[r][k] - f * M[c][k] for k in range(n)]
    return d

H = 2
graphs = {
    '133': dict(edges=[(0, 5), (0, 1), (1, 2), (2, 5), (0, 3), (3, 4), (4, 5)],
                C=[[-1, 1, 1, 1, 0, 0, 0], [-1, 0, 0, 0, 1, 1, 1]],
                paths={'A': [0], 'B': [1, 2, 3], 'C': [4, 5, 6]}),
    '223': dict(edges=[(0, 1), (1, 5), (0, 2), (2, 5), (0, 3), (3, 4), (4, 5)],
                C=[[1, 1, -1, -1, 0, 0, 0], [-1, -1, 0, 0, 1, 1, 1]],
                paths={'A': [0, 1], 'B': [2, 3], 'C': [4, 5, 6]}),
}
root, nV = 5, 6

def br(a, b):
    return a[0] * b[1] - a[1] * b[0]

def density(g, t, lam):
    """returns (D*U^{7/2}, U, list of D*U^{7/2} from each hat-j component)"""
    E = g['edges']; C = g['C']; N = len(E); n = 2
    nonroot = [v for v in range(nV) if v != root]
    B = [[(1 if v == e[0] else 0) - (1 if v == e[1] else 0) for v in nonroot] for e in E]  # N x 5
    # Kirchhoff: B^T y = lam (5 eqs), C t y = 0 (2 eqs)
    A = [[B[e][iv] for e in range(N)] for iv in range(len(nonroot))] + \
        [[C[r][e] * t[e] for e in range(N)] for r in range(n)]
    y = []
    for i in range(H):
        rhs = [lam[v][i] for v in nonroot] + [0] * n
        y.append(solve(A, rhs))           # y[i][e]
    # checks: currents
    for i in range(H):
        for iv, v in enumerate(nonroot):
            assert sum(B[e][iv] * y[i][e] for e in range(N)) == lam[v][i]
        for r in range(n):
            assert sum(C[r][e] * t[e] * y[i][e] for e in range(N)) == 0
    M = [[sum(C[r][e] * t[e] * C[s][e] for e in range(N)) for s in range(n)] for r in range(n)]
    U = det(M)
    Minv = [[M[1][1] / U, -M[0][1] / U], [-M[1][0] / U, M[0][0] / U]]
    dM = [[{(e,): Fr(C[r][e] * C[s][e]) for e in range(N) if C[r][e] * C[s][e] != 0} for s in range(n)] for r in range(n)]
    # A_12 = -1/2 sum_{k,l} dM_1k Minv_kl dM_l2 ; Pf for n=2 is A_12 (flagship eq-rho)
    A12 = {}
    for k in range(n):
        for l in range(n):
            A12 = add(A12, scal(wedge(dM[0][k], dM[l][1]), Fr(-1, 2) * Minv[k][l]))
    psi = [[{(e,): C[r][e] * y[i][e] for e in range(N) if C[r][e] * y[i][e] != 0} for r in range(n)] for i in range(H)]
    form = A12          # rho_2 * U^{1/2}
    for i in range(H):
        for r in range(n):
            form = wedge(form, psi[i][r])
    # omega = U^{-H} rho_2 ^ prod psi = U^{-2} U^{-1/2} form ; D = omega_{hat j} / ((-1)^j t_j)
    vals = []
    for j in range(N):
        key = tuple(k for k in range(N) if k != j)
        cj = form.get(key, Fr(0))
        vals.append(cj / ((-1) ** j * t[j]))   # = D * U^{5/2} ; times U gives D*U^{7/2}
    vals = [v * U for v in vals]
    return vals[0], U, vals, M

out = []
def P(*a):
    s = ' '.join(str(x) for x in a)
    print(s); out.append(s)

for name, g in graphs.items():
    E = g['edges']; N = len(E)
    nonroot = [v for v in range(nV) if v != root]
    B = [[(1 if v == e[0] else 0) - (1 if v == e[1] else 0) for v in nonroot] for e in E]
    BC = [[B[e][iv] for e in range(N)] for iv in range(5)] + g['C']
    P('graph', name, 'det[B^T;C] =', det(BC))
    ratios = []
    for trial in range(4):
        t = [Fr(random.randint(1, 40), random.randint(1, 40)) for _ in range(N)]
        lam = {v: (Fr(random.randint(-30, 30), random.randint(1, 9)), Fr(random.randint(-30, 30), random.randint(1, 9)))
               for v in range(5)}
        lam[root] = (-sum(lam[v][0] for v in range(5)), -sum(lam[v][1] for v in range(5)))
        DU, U, vals, M = density(g, t, lam)
        R = {P_: sum(t[e] for e in es) for P_, es in g['paths'].items()}
        Uform = R['A'] * R['B'] + R['B'] * R['C'] + R['C'] * R['A']
        dth = br(lam[1], lam[2]) * br(lam[3], lam[4])
        ratios.append(DU / dth)
        P('  trial', trial, 'all hat-j equal:', len(set(vals)) == 1,
          '| U == R_AR_B+R_BR_C+R_CR_A:', U == Uform,
          '| D*U^(7/2)/([p1p2][p3p4]) =', DU / dth, '| lambda_0 =', lam[0])
    P('  ratio constant over trials:', len(set(ratios)) == 1, 'value', ratios[0])
    g['c'] = ratios[0]

# ---------- Schwinger integral of U^{-7/2} over the 6-simplex, path reduction:
# int_{Delta^6} f(R) dt = int_{Delta^2} prod R_P^{n_P-1}/(n_P-1)! f(R) dR_A dR_B
mp.mp.dps = 30
def schwinger(nA, nB, nC):
    def f(u, v):
        RA = u; RB = (1 - u) * v; RC = (1 - u) * (1 - v)
        Uv = RA * RB + RB * RC + RC * RA
        return (1 - u) * RA ** (nA - 1) * RB ** (nB - 1) * RC ** (nC - 1) / (
            mp.factorial(nA - 1) * mp.factorial(nB - 1) * mp.factorial(nC - 1)) * Uv ** mp.mpf(-3.5)
    return mp.quad(f, [0, mp.mpf(1) / 2, 1], [0, mp.mpf(1) / 2, 1])

for name, nn, target in (('133', (1, 3, 3), -mp.pi / 20), ('223', (2, 2, 3), mp.pi / 30)):
    I = schwinger(*nn)
    K = graphs[name]['c'] * I       # K/([p1p2][p3p4])
    P('graph', name, ' int U^{-7/2} =', mp.nstr(I, 22), '| c*int =', mp.nstr(mp.mpf(K.numerator) / K.denominator * I if isinstance(K, Fr) else K, 22),
      '| flagship value', mp.nstr(target, 22), '| rel diff', mp.nstr(abs((graphs[name]['c'].numerator * I / graphs[name]['c'].denominator) - target) / abs(target), 3))

open((__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/notes/flagship/referee/kern_check.log', 'w').write('\n'.join(out) + '\n')
