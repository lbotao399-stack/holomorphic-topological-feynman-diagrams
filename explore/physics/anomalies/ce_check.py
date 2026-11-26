"""Direct exact check of the closedness lemma in the translation-reduced CE complex of g[[z_1,z_2,z_3]]
(plane-wave evaluation, cochains as polynomials in momenta restricted to sum p = 0).
omega(x_1 e^{p_1 z},...,x_4 e^{p_4 z}) = Delta(p_1,p_2,p_3)^n T(x_1,...,x_4),
d omega(y_0..y_4) = sum_{i<j} (-1)^{i+j} omega([y_i,y_j], y_0..^i..^j..y_4),  [x e^{pz}, y e^{qz}] = [x,y] e^{(p+q)z}.
Algebras: matrix Lie algebras over Q; T built from traces. Random integer samples, exact Fractions."""
import random, itertools
from fractions import Fraction as Fr
rng = random.Random(5)

def mm(a, b):
    n = len(a)
    return [[sum(a[i][k] * b[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
def br(a, b):
    ab, ba = mm(a, b), mm(b, a)
    return [[ab[i][j] - ba[i][j] for j in range(len(a))] for i in range(len(a))]
def tr(a):
    return sum(a[i][i] for i in range(len(a)))
def det3(a, b, c):
    return (a[0]*(b[1]*c[2]-b[2]*c[1]) - a[1]*(b[0]*c[2]-b[2]*c[0]) + a[2]*(b[0]*c[1]-b[1]*c[0]))
def vadd(a, b):
    return [x + y for x, y in zip(a, b)]
def psign(p):
    s = 1
    for i in range(len(p)):
        for j in range(i+1, len(p)):
            if p[i] > p[j]: s = -s
    return s
PERM4 = [(s, psign(s)) for s in itertools.permutations(range(4))]

def rand_mat(N, traceless):
    a = [[rng.randint(-3, 3) for _ in range(N)] for _ in range(N)]
    if traceless:
        a[N-1][N-1] -= tr(a)
    return a

# invariant 4-tensors
def str4(x):            # symmetrized tr(x1x2x3x4)
    return sum(tr(mm(mm(x[s[0]], x[s[1]]), mm(x[s[2]], x[s[3]]))) for s, _ in PERM4)
def sq22(x):            # symmetrized tr(x1x2)tr(x3x4)
    return sum(tr(mm(x[s[0]], x[s[1]])) * tr(mm(x[s[2]], x[s[3]])) for s, _ in PERM4)
def lam4(x):            # product of traces (factors through g/[g,g] for gl_N)
    r = 1
    for a in x: r *= tr(a)
    return r
def alt_lam_h3(x):      # alternation of tr(x1) tr(x2 [x3,x4])  (z^* (x) H^3 type, gl_N)
    return sum(sg * tr(x[s[0]]) * tr(mm(x[s[1]], br(x[s[2]], x[s[3]]))) for s, sg in PERM4)

def omega(T, n, args):
    xs = [a[0] for a in args]; ps = [a[1] for a in args]
    return Fr(det3(ps[0], ps[1], ps[2]))**n * T(xs)

def d_omega(T, n, args):
    tot = Fr(0)
    for i, j in itertools.combinations(range(5), 2):
        rest = [args[k] for k in range(5) if k not in (i, j)]
        y = (br(args[i][0], args[j][0]), vadd(args[i][1], args[j][1]))
        tot += (-1)**(i + j) * omega(T, n, [y] + rest)
    return tot

def sample(N, traceless):
    ps = [[rng.randint(-4, 4) for _ in range(3)] for _ in range(4)]
    ps.append([-sum(p[a] for p in ps) for a in range(3)])
    return [(rand_mat(N, traceless), ps[k]) for k in range(5)]

def antisym_check(T, n, N, traceless):
    """omega is S4-antisymmetric modulo sum p = 0 (needed for it to be a cochain)."""
    ps = [[rng.randint(-4, 4) for _ in range(3)] for _ in range(3)]
    ps.append([-sum(p[a] for p in ps) for a in range(3)])
    args = [(rand_mat(N, traceless), ps[k]) for k in range(4)]
    ref = omega(T, n, args)
    return all(omega(T, n, [args[s[i]] for i in range(4)]) == sg * ref for s, sg in PERM4)

cases = [('sl2', 2, True), ('sl3', 3, True), ('gl2', 2, False), ('gl3', 3, False)]
tens = [('Str4', str4, 'sym'), ('(tr xx)^2 sym', sq22, 'sym'), ('lambda^4', lam4, 'sym'), ('alt lambda.H3', alt_lam_h3, 'alt')]
for name, N, tl in cases:
    for tn, T, kind in tens:
        if tl and tn in ('lambda^4', 'alt lambda.H3'):
            continue
        for n in range(0, 6):
            if (kind == 'sym') != (n % 2 == 1):
                continue
            ok_cochain = antisym_check(T, n, N, tl)
            vals = [d_omega(T, n, sample(N, tl)) for _ in range(3)]
            nz = [omega(T, n, sample(N, tl)[:4]) != 0 for _ in range(2)]
            print('%-4s %-14s n=%d: cochain(S4-antisym) %s; omega nonzero %s; d omega == 0 on 3 samples: %s' %
                  (name, tn, n, ok_cochain, any(nz), all(v == 0 for v in vals)), flush=True)
