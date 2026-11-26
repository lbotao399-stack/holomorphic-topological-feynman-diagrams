"""C4: one-flux insertion in the flat genus-one vertex (1/24) Tr a1 det(d_i d_j a0); genus-two substitution."""
from itertools import permutations
import sympy as sp

z1, z2, w1, w2, w3, N = sp.symbols('z1 z2 w1 w2 w3 N')
X = [z1, z2, w1, w2, w3]
f = 1 / (w1 * w2 * w3)                       # normal coefficient of F = f dz1^dz2 (Laurent representative)
t = sp.symbols('t0:5')


def hess(g):
    return sp.Matrix(5, 5, lambda i, j: sp.diff(g, X[i], X[j]))


def flux_linear(waves):
    """coefficient of t0 t1 t2 t3 t4 in det(Hess(t0 f + sum t_i wave_i))."""
    H = t[0] * hess(f)
    for ti, wv in zip(t[1:], waves):
        H += ti * hess(wv)
    d = sp.expand(H.det(method='berkowitz'))
    c = sp.Poly(d, *t).coeff_monomial(t[0] * t[1] * t[2] * t[3] * t[4])
    return sp.simplify(c)


def residue(expr):
    """coefficient of (w1 w2 w3)^-1 in a Laurent monomial sum."""
    e = sp.expand(expr * w1 * w2 * w3)
    return sp.expand(e).subs({w1: 0, w2: 0, w3: 0}) if not e.free_symbols & {w1, w2, w3} else \
        sp.Poly(sp.expand(e * (w1 * w2 * w3)**10), w1, w2, w3).coeff_monomial((w1 * w2 * w3)**10)


def dictionary(expr):
    """Laurent sum -> {(a1,a2,a3): coefficient of (dtheta)^a}, using w^{-a-1} <-> (-1)^{|a|} (dtheta)^a / a!."""
    e = sp.expand(expr)
    out = {}
    for term in sp.Add.make_args(e):
        c, mon = term.as_coeff_Mul()
        # degree of a reciprocal: use as_powers_dict
        pd = mon.as_powers_dict()
        a = tuple(-int(pd.get(v, 0)) - 1 for v in (w1, w2, w3))
        rest = mon / sp.Mul(*[v**pd.get(v, 0) for v in (w1, w2, w3)])
        assert rest == 1, term
        coeff = c * rest * (-1)**sum(a) / sp.Mul(*[sp.factorial(x) for x in a])
        out[a] = out.get(a, 0) + coeff
    return {k: sp.nsimplify(v) for k, v in out.items() if v != 0}


if __name__ == '__main__':
    # 1. quadratic jets
    c1 = flux_linear([z1**2 / 2, z2**2 / 2, w2**2 / 2, w3**2 / 2])
    out1 = sp.simplify(N * c1 / 24)
    print('quadratic jets: multilinear coefficient =', c1, ';  (N/24)*it =', out1,
          '; equals N/(12 w1^3 w2 w3):', sp.simplify(out1 - N / (12 * w1**3 * w2 * w3)) == 0)
    print('   equals (N/24) d_w1^2 F coefficient:', sp.simplify(out1 - N * sp.diff(f, w1, 2) / 24) == 0)
    print('   dictionary:', dictionary(out1 / N), '(times N)')
    r1 = sp.expand(out1 * w1**2 * w1 * w2 * w3)
    print('   residue against w1^2 dw1^dw2^dw3: coefficient of (w1w2w3)^-1 =', r1)
    # 2. mixed jets
    c2 = flux_linear([z1 * w2, z1 * w3, z2 * w2, z2 * w3])
    out2 = sp.simplify(N * c2 / 24)
    r2 = sp.expand(out2 * w1**2 * w1 * w2 * w3)
    print('mixed jets z1w2,z1w3,z2w2,z2w3: multilinear coefficient =', c2, '; (N/24)*it =', out2, '; residue =', r2)
    s1, s2, s3, s4 = sp.symbols('s1:5')
    K = sp.Matrix([[0, 0, 0, s1, s2], [0, 0, 0, s3, s4], [0, 0, 1, 0, 0], [s1, s3, 0, 0, 0], [s2, s4, 0, 0, 0]])
    minor = sp.factor(K.extract([0, 1, 3, 4], [0, 1, 3, 4]).det())
    print('   complementary 4-row minor for a0 = s1 z1w2 + s2 z1w3 + s3 z2w2 + s4 z2w3:', minor,
          '; s1s2s3s4 coefficient:', sp.Poly(sp.expand(minor), s1, s2, s3, s4).coeff_monomial(s1 * s2 * s3 * s4))
    # 3. six directions
    probes = [((0, 1, 0), (0, 0, 1)), ((0, 0, 1), (1, 0, 0)), ((1, 0, 0), (0, 1, 0)),
              ((0, 0, 1), (1, -1, 0)), ((0, 1, 0), (-1, 0, 1)), ((1, 0, 0), (0, 1, -1))]
    basis = [(2, 0, 0), (0, 2, 0), (0, 0, 2), (1, 1, 0), (1, 0, 1), (0, 1, 1)]   # (dth_i)(dth_j), i<=j
    rows = []
    for p, q in probes:
        pw = sum(pi * wi for pi, wi in zip(p, (w1, w2, w3)))
        qw = sum(qi * wi for qi, wi in zip(q, (w1, w2, w3)))
        c = flux_linear([z1**2 / 2, z2**2 / 2, pw**2 / 2, qw**2 / 2])
        v = sp.Matrix(p).cross(sp.Matrix(q))
        quad = sum(v[i] * v[j] * sp.diff(f, (w1, w2, w3)[i], (w1, w2, w3)[j]) for i in range(3) for j in range(3))
        dct = dictionary(c)
        rows.append([dct.get(b, 0) for b in basis])
        print(f'   probe p={p} q={q}: v=p x q={list(v)}, response = v^T Hess(f) v: {sp.simplify(c - quad) == 0}, '
              f'coords {rows[-1]}')
    M = sp.Matrix(rows)
    print('six-direction response determinant (without N/24):', M.det(), ' with (N/24) per row:', sp.factor(M.det() * (N / 24)**6))

    # 4. genus-two seven-slot trace symbol on the stated configuration
    a, b, c, d, e, s, tt, q1, q2, q3 = sp.symbols('a b c d e s t q1 q2 q3')

    def P_symbol(mom):
        """mom: list of 7 momenta (5-vectors) summing to zero; formula of iib-bcov README."""
        cs = []
        for i in range(6):
            rows_ = [mom[j] for j in range(6) if j != i]
            cs.append((-1)**i * sp.Matrix(rows_).det())
        cs.append(0)
        tot = sum(cs)
        A = [7 * ci - tot for ci in cs]
        return sp.expand(-sp.Rational(31, 2987228160) * sum(x**2 for x in A)**2 + sp.Rational(1, 71124480) * sum(x**4 for x in A))

    inputs = [[0, 0, a, 0, 0], [0, 0, b, 0, 0], [0, 0, 0, c, 0], [0, d, 0, 0, s], [e, 0, 0, 0, tt], [0, 0, q1, q2, q3]]
    bal = [-sum(v[k] for v in inputs) for k in range(5)]
    P = P_symbol(inputs + [bal])
    claim = -a**2 * b**2 * c**4 * d**4 * e**4 * q3**4 / 31104
    print('genus-two P on configuration:', sp.factor(P), '; equals claim:', sp.expand(P - claim) == 0)
    # symmetry of the symbol under permutations of the seven slots (random rational momenta)
    import random
    random.seed(3)
    mom = [[sp.Rational(random.randint(-5, 5), random.randint(1, 3)) for _ in range(5)] for _ in range(6)]
    mom.append([-sum(v[k] for v in mom) for k in range(5)])
    vals = set()
    for perm in random.sample(list(permutations(range(7))), 30):
        vals.add(P_symbol([mom[i] for i in perm]))
    print('   P invariant under 30 random slot permutations of a generic conserved 7-tuple:', len(vals) == 1)
    # output coefficient: P-coefficient * N * d_w3^4 f * w3^2 along e d_z2 (coefficient of a^2 b^2 c^4 d^4 e^5 s t)
    Pc = sp.Poly(P, a, b, c, d, e, q3).coeff_monomial(a**2 * b**2 * c**4 * d**4 * e**4 * q3**4)
    vec = sp.expand(Pc * N * sp.diff(f, w3, 4) * w3**2)
    print('   P coefficient =', Pc, '; P-coeff * N * d_w3^4(f) * w3^2 =', sp.factor(vec),
          '; equals -N/(1296 w1 w2 w3^3):', sp.simplify(vec + N / (1296 * w1 * w2 * w3**3)) == 0)
    print('   P has no s,t dependence:', not P.has(s) and not P.has(tt), '; P depends on q only via q3:', not P.has(q1) and not P.has(q2))
