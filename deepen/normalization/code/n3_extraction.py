"""(b) Position space -> Schwinger space: extraction sign, reordering signs, constant.

theta^a_e = dX^a_e - c_a J^a_e dt_e  (a = 0: tau, c = 1, J = w = tau_e/(2 t_e);  a = 1,2: zb^i, c = 4, J = y^i = zb^i_e/(4 t_e)),
so that deta_e = theta^0_e/(2 sqrt t_e), dy^i_e = theta^i_e/(4 t_e).
Claim (htk CONVENTIONS Sec. 4, S3 eq-extraction):
  wedge_e (theta^0_e theta^1_e theta^2_e) = vs_ext prod_a c_a^n [prod_v dX^0_v prod_v dX^1_v prod_v dX^2_v] ^ prod_a prod_r psi^a_r + lower,
  psi^a_r = sum_e C_re J^a_e dt_e, C = positive cycle basis of htk.
Computed here: the coefficient of [all dX] ^ dt_{0..N-1 without e0} as a determinant (rows = one-forms in wedge order,
columns = generators in monomial order), exact rationals, random rational positions and Schwinger times.
"""
import sys
import random
from fractions import Fraction as Fr
import sympy as sp
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
import htk
from htk import conventions


def perm_sign(seq):
    seq = list(seq)
    s = 1
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                s = -s
    return s


def det(M):
    return sp.Matrix(M).det()


def extraction_ratio(G, cs, seed, e0=None):
    """returns ratio = coefficient(LHS) / coefficient(prod psi) for the monomial [all dX] dt_{hat e0}."""
    R = random.Random(seed)
    D = len(cs)
    nv, N = G.nv, G.ne
    m = nv - 1
    B = G.incidence()                 # N x m
    C = G.cycle_matrix()              # n x N
    n = len(C)
    if e0 is None:
        e0 = R.randrange(N)
    t = [Fr(R.randint(1, 40), R.randint(1, 40)) for _ in range(N)]
    X = [[Fr(R.choice([-1, 1]) * R.randint(1, 30), R.randint(1, 9)) for _ in range(m)] for a in range(D)]
    Xe = [[sum(B[e][v] * X[a][v] for v in range(m)) for e in range(N)] for a in range(D)]
    J = [[Xe[a][e] / ((2 if a == 0 and cs[0] == 1 else 4) * t[e]) for e in range(N)] for a in range(D)]
    if cs == (4,):                     # holomorphic only (T = 0, H = 1)
        J = [[Xe[0][e] / (4 * t[e]) for e in range(N)]]
    ngen = D * m + N
    rows = []
    for e in range(N):
        for a in range(D):
            row = [Fr(0)] * ngen
            for v in range(m):
                row[a * m + v] = Fr(B[e][v])
            row[D * m + e] = -cs[a] * J[a][e]
            rows.append(row)
    cols = list(range(D * m)) + [D * m + e for e in range(N) if e != e0]
    lhs = det([[row[c] for c in cols] for row in rows])
    psirows = []
    for a in range(D):
        for r in range(n):
            psirows.append([C[r][e] * J[a][e] for e in range(N) if e != e0])
    rhs = det(psirows)
    return sp.Rational(lhs) / sp.Rational(rhs), n, m, N


if __name__ == "__main__":
    cases = [("theta(1,3,3)", htk.theta(1, 3, 3)[0]), ("theta(2,2,3)", htk.theta(2, 2, 3)[0]),
             ("graph_A (S2 labelling)", htk.graph_A()), ("graph_B (S2 labelling)", htk.graph_B())]
    for p, q, r in [(1, 3, 3), (2, 2, 3)]:
        nv = p + q + r - 1
        # theta_rootN of msym-odd: interior first, N' = nv-2, root N = nv-1
        N_, Np = nv - 1, nv - 2
        edges, nxt = [], 0
        for n_ in (p, q, r):
            prev = N_
            for j in range(n_):
                cur = Np if j == n_ - 1 else nxt
                if j < n_ - 1:
                    nxt += 1
                edges.append((prev, cur))
                prev = cur
        cases.append(("theta_rootN(%d,%d,%d)" % (p, q, r), htk.Graph(nv, edges, "rootN")))
    print("(H,T) = (2,1), c = (1,4,4):")
    for name, G in cases:
        vals = set()
        for seed in range(3):
            ratio, n, m, N = extraction_ratio(G, (1, 4, 4), seed)
            vals.add(ratio)
        vs = [v / 16 ** n for v in vals]
        s1 = (-1) ** ((3 * n) * (3 * m))            # [dX top] ^ Psi -> Psi ^ [dX top]
        order = [a * m + v for v in range(m) for a in range(3)]  # prod_v (dX^0_v dX^1_v dX^2_v) in canonical labels
        s2 = perm_sign(order)
        print("  %-24s n=%d |V|-1=%d |I|=%d: LHS/(16^n psi) at 3 random points = %s ; htk formula %d ; "
              "s1 = (-1)^{9n(|V|-1)} = %d ; s2 (reorder to prod_v) = %d ; vs_ext*s1*s2 = %s"
              % (name, n, m, N, sorted(vs), conventions.extraction_sign_formula(G, 3), s1, s2,
                 [v * s1 * s2 for v in sorted(vs)]))
    # D = 1 cases used by the end-to-end checks (n4)
    print("(H,T) = (0,1), c = (1,):")
    for k in (3, 5):
        G = htk.banana(k)
        vals = {extraction_ratio(G, (1,), s)[0] for s in range(3)}
        print("  banana(%d): vs_ext = %s ; htk formula %d" % (k, vals, conventions.extraction_sign_formula(G, 1)))
    print("(H,T) = (1,0), c = (4,):")
    for k in (2, 3):
        G = htk.banana(k)
        n = k - 1
        vals = {extraction_ratio(G, (4,), s)[0] / 4 ** n for s in range(3)}
        print("  banana(%d): vs_ext = %s ; htk formula %d" % (k, vals, conventions.extraction_sign_formula(G, 1)))
    # constant
    n = sp.Symbol('n')
    pi = sp.pi
    te = sp.Symbol('t_e', positive=True)
    for nn in (1, 2, 3):
        N = 3 * nn + 1
        m = 2 * nn + 1
        const = (1 / (4 * pi ** sp.Rational(5, 2))) ** N * (4 * (4 * pi) ** sp.Rational(5, 2)) ** m * 4 ** (2 * nn) \
            * (te ** sp.Rational(5, 2) / (2 * sp.sqrt(te) * (4 * te) ** 2)) ** N
        print("constant n=%d: (4pi^{5/2})^{-|I|} (4 (4pi)^{5/2})^{|V|-1} 4^{2n} prod_e t^{5/2}/(2 sqrt t (4t)^2) = %s ; "
              "(8 pi^{5/2})^{-n} = %s" % (nn, sp.simplify(const), sp.simplify((8 * pi ** sp.Rational(5, 2)) ** (-nn))))
