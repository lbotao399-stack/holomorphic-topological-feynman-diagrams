"""Theories for hteng (field content: ESW 2002.10517 Part II).
  3d  (H=1): N=2 HT twist data (c, B, chirals phi/psi, W)  [ESW Thm 11.1]; 3d N=4 = N=2 vector + adjoint
             chiral Phi + hypers (Q, Qt), W = Qt Phi Q  [ESW Thm 11.3: T^*[-1]Sect(C x R, (U x K^1/2)//g)];
             3d N=8 = three adjoints, W = Phi[X,Y]  [ESW Thm 11.8, g/g on C_Dol x R].
  5d  (H=2): BF R/g  [ESW Thm 9.3]; U//g with symplectic U (u-u propagator)  [ESW Thm 9.2];
             5d N=2 = BF g/g  [ESW 9.2.1].
  7d  (H=3): BF Bg  [ESW Thm 7.2].
Gauge bases: so(3) real (f = eps); sl(2) (e,f,h): [h,e]=2e, [h,f]=-2f, [e,f]=h (rational, fundamentals)."""
from fractions import Fraction as Fr
from hteng import Theory


def eps(a, b, c):
    if len({a, b, c}) < 3:
        return 0
    return 1 if (a, b, c) in [(0, 1, 2), (1, 2, 0), (2, 0, 1)] else -1


F_SO3 = {(a, b, c): eps(a, b, c) for a in range(3) for b in range(3) for c in range(3) if eps(a, b, c)}

# sl2: index 0 = e, 1 = f, 2 = h
F_SL2 = {(2, 0, 0): 2, (0, 2, 0): -2, (2, 1, 1): -2, (1, 2, 1): 2, (0, 1, 2): 1, (1, 0, 2): -1}
FUND_SL2 = [{(0, 1): 1}, {(1, 0): 1}, {(0, 0): 1, (1, 1): -1}]  # rho(t_a)^i_j


def adj_rho(F, G):
    # rho(t_a)^b_c = f_ac^b
    return [{(b, c): v for (a2, c, b), v in F.items() if a2 == a} for a in range(G)]


def dual_rho(rho):
    return [{(j, i): -v for (i, j), v in r.items()} for r in rho]


def direct_sum(rhos, dims):
    out = [dict() for _ in range(len(rhos[0]))]
    off = 0
    for r, d in zip(rhos, dims):
        for a in range(len(r)):
            for (i, j), v in r[a].items():
                out[a][(i + off, j + off)] = v
        off += d
    return out


def check_rep(F, rho, dim):
    """[rho_a, rho_b] = f_ab^c rho_c"""
    G = len(rho)
    def mat(r):
        M = [[Fr(0)] * dim for _ in range(dim)]
        for (i, j), v in r.items():
            M[i][j] = Fr(v)
        return M
    Ms = [mat(r) for r in rho]
    for a in range(G):
        for b in range(G):
            C = [[sum(Ms[a][i][k] * Ms[b][k][j] - Ms[b][i][k] * Ms[a][k][j] for k in range(dim)) for j in range(dim)] for i in range(dim)]
            R = [[sum(Fr(F.get((a, b, c), 0)) * Ms[c][i][j] for c in range(G)) for j in range(dim)] for i in range(dim)]
            if C != R:
                return False
    return True


# ------------------------------------------------------------------ 3d (H = 1)
def u1_3d(charges, W=None, name=None):
    bf = [(1, [{(0, 0): Fr(q)}]) for q in charges]
    return Theory(1, 1, {}, bf=bf, W=W, name=name or "3d U(1) charges %s" % (charges,))


def sqed_3d(Nf, q=1):
    """3d N=4 U(1) with Nf hypers of charge q: multiplets 0 = Phi (q=0), 1..Nf = Q_i (+q),
    Nf+1..2Nf = Qt_i (-q); W = sum_i Qt_i Phi Q_i."""
    charges = [0] + [q] * Nf + [-q] * Nf
    W = {}
    for i in range(Nf):
        key = tuple(sorted([(1 + Nf + i, 0), (0, 0), (1 + i, 0)]))
        W[key] = W.get(key, 0) + 1
    th = u1_3d(charges, W=W, name="3d N=4 SQED Nf=%d q=%d" % (Nf, q))
    th.rcharges = [Fr(1)] + [Fr(1, 2)] * (2 * Nf)
    th.flavour = [(0,) * (Nf + 1)] + [tuple(1 if k == i else 0 for k in range(Nf)) + (1,) for i in range(Nf)] \
        + [tuple(-1 if k == i else 0 for k in range(Nf)) + (1,) for i in range(Nf)]
    # flavour: U(1)^Nf (Q_i +1, Qt_i -1) is mostly broken by W to the diagonal; keep axial (last entry)
    # sectors are labelled by (Q_i - Qt_i charges) and axial (N_Q + N_Qt - 2 N_Phi); define below
    th.flavour = None
    return th


def su2_N4_3d(Nf):
    """3d N=4 SU(2) (sl2 basis) with Nf fundamental hypers: multiplet 0 = Phi (adjoint),
    1..Nf = Q_i (fund), Nf+1..2Nf = Qt_i (antifund); W = sum_i Phi^a Qt_i rho_a Q_i."""
    adj = adj_rho(F_SL2, 3)
    bf = [(3, adj)] + [(2, FUND_SL2) for _ in range(Nf)] + [(2, dual_rho(FUND_SL2)) for _ in range(Nf)]
    W = {}
    for i in range(Nf):
        for a in range(3):
            for (k, l), v in FUND_SL2[a].items():
                key = tuple(sorted([(0, a), (1 + Nf + i, k), (1 + i, l)]))
                W[key] = W.get(key, 0) + v
    W = {k: v for k, v in W.items() if v}
    th = Theory(1, 3, F_SL2, bf=bf, W=W, name="3d N=4 SU(2) Nf=%d" % Nf)
    th.rcharges = [Fr(1)] + [Fr(1, 2)] * (2 * Nf)
    return th


def su2_pure_N4_3d():
    th = Theory(1, 3, F_SO3, bf=[(3, adj_rho(F_SO3, 3))], name="3d N=4 SU(2) pure")
    th.rcharges = [Fr(1)]
    return th


def su2_N8_3d():
    adj = adj_rho(F_SO3, 3)
    W = {}
    for a in range(3):
        for b in range(3):
            for c in range(3):
                if eps(a, b, c):
                    key = ((0, a), (1, b), (2, c))
                    W[key] = W.get(key, 0) + eps(a, b, c)
    th = Theory(1, 3, F_SO3, bf=[(3, adj)] * 3, W=W, name="3d N=8 SU(2)")
    th.rcharges = [Fr(2, 3)] * 3
    return th


# ------------------------------------------------------------------ 5d (H = 2), 7d (H = 3)
def bf_theory(H, gauge, reps, name):
    """BF R/g: gauge in {'so3','sl2','u1'}; reps: list of (dim, rho)"""
    if gauge == 'so3':
        G, F = 3, F_SO3
    elif gauge == 'sl2':
        G, F = 3, F_SL2
    elif gauge == 'u1':
        G, F = 1, {}
    return Theory(H, G, F, bf=list(reps), name=name)


def u1_reps(charges):
    return [(1, [{(0, 0): Fr(q)}]) for q in charges]


def sym_theory(H, gauge, sym, name, signs=None):
    """U//g: sym = list of (dimU, Om upper {(al,be): v}, rho); om lower computed"""
    if gauge == 'sl2':
        G, F = 3, F_SL2
    elif gauge == 'so3':
        G, F = 3, F_SO3
    elif gauge == 'u1':
        G, F = 1, {}
    oml = []
    for (dimU, Om, rho) in sym:
        # invert Om
        M = [[Fr(Om.get((i, j), 0)) for j in range(dimU)] for i in range(dimU)]
        inv = mat_inv(M)
        # om_{al be} Om^{be ga} = delta  -> om = inv(Om) (as matrices: om . Om = 1)
        oml.append({(i, j): inv[i][j] for i in range(dimU) for j in range(dimU) if inv[i][j]})
    return Theory(H, G, F, sym=list(sym), name=name, omlow=oml, signs=signs)


def mat_inv(M):
    n = len(M)
    A = [row[:] + [Fr(int(i == j)) for j in range(n)] for i, row in enumerate(M)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[p] = A[p], A[c]
        pv = A[c][c]
        A[c] = [x / pv for x in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]
                A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    return [row[n:] for row in A]


def half_hyper_sl2():
    """U = C^2 fundamental of sl2 = sp(2): Om^{01} = 1, Om^{10} = -1"""
    return (2, {(0, 1): 1, (1, 0): -1}, FUND_SL2)


def full_hyper_sym(rho, dim):
    """U = R + R^*: u = (phi^i, psi_i), Om^{i, dim+i} = 1, Om^{dim+i, i} = -1"""
    Om = {}
    for i in range(dim):
        Om[(i, dim + i)] = 1
        Om[(dim + i, i)] = -1
    rr = direct_sum([rho, dual_rho(rho)], [dim, dim])
    return (2 * dim, Om, rr)
