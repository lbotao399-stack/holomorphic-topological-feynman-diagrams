"""(b) Classical bracket as the tree graph with one heat edge, Wick rule of normalization b5 (= b11 of its REPORT):
  vertices v0 (cochain a, integrated), v1 (cochain b, root), edge e = (v0, v1), externals x1, x2 at v0, x3, x4 at v1;
  sigma = sgn(source -> target), source = local order at v0, then at v1; target = (e,0), (e,1), x1, x2, x3, x4;
  vertex factor = cochain in local order; momenta into the vertex: x_i -> p_i, (e,0) -> -y_e, (e,1) -> +y_e,
  y_e = Kirchhoff current = lambda_{v0} = p1 + p2;  raw tree integral G = lim int_X Omega ^ e^{lam z} K_eps = 1 (n = 0:
  vs_ext = s1 = s2 = 1, K_tree = 1);  |Aut| = 4 (legs at v0, legs at v1);  colour: sum_e X_e (x) X^e (dual basis, tr).
  tree cochain = (U_tree/4) alt( sigma a(L0) b(L1) ).
Compared (exact, explicit matrices, explicit dual basis, no ce bracket code) with the shuffle bracket {a,b}
(evalx.bracket_explicit and ce.bracket), for all 36 pairs of local orders.  (1/2){I,I} -> {I_a, I_b} = (1/4) w_tree:
  {I_a, I_b} = U_tree F_{(tree/U_tree)} , i.e. eps = U_tree x (tree/U_tree : {a,b}).
Also: Feynman rule of the J_1 vertex from (1/3) int dz d_abc d_{z1}A^a d_{z2}A^b A^c (dictionary F <-> alt(PW))."""
import itertools, random
from common import *
import evalx

rng = random.Random(7)

# J_1 vertex: (1/3) d_abc d_{z1}A^a d_{z2}A^b A^c  ->  w = alt(P W), P = (1/3) p1^(1) p2^(2), W = d
raw = ce.raw_from({(1, 0, 0, 1, 0, 0): Fr(1, 3)}, {((0, 1, 2),): {0: Fr(1)}, ((0, 2, 1),): {0: Fr(1)}})
w = ce.alt_raw(raw, 3)
print("alt((1/3) p1^1 p2^2 d(X1,X2,X3)) == D = [p1p2] tr(1{23}):", is_zero(ce.lin((1, w), (-1, J_cochain('d')))))
print("(1/6) alt([p1p2] d) == D:", is_zero(ce.lin((Fr(1, 6), ce.alt_raw(ce.raw_from(ce.bracket_poly(0, 1, 3),
      {((0, 1, 2),): {0: Fr(1)}, ((0, 2, 1),): {0: Fr(1)}}), 3)), (-1, J_cochain('d')))))
# 3! leg assignments: sum_s sgn(s) q^1_{h_s1} q^2_{h_s2} = [q_h1 q_h2] + [q_h2 q_h3] + [q_h3 q_h1] = 3 [q_h1 q_h2] on sum q = 0
q = [[sp.Symbol('q%d_%d' % (i, c)) for c in range(2)] for i in range(3)]
s_ = 0
for perm in itertools.permutations(range(3)):
    s_ += ce.perm_sign(perm) * q[perm[0]][0] * q[perm[1]][1]
s_ = sp.expand(s_.subs({q[2][0]: -q[0][0] - q[1][0], q[2][1]: -q[0][1] - q[1][1]}))
br01 = q[0][0] * q[1][1] - q[0][1] * q[1][0]
print("sum_s sgn(s) q^1_{s1} q^2_{s2} == 3 [q1 q2] on sum q = 0:", sp.expand(s_ - 3 * br01) == 0)


def tree_value(a, b, mats, moms, N, L0, L1):
    """sigma * sum_e a(L0) b(L1) at fixed slots x1..x4 (no alternation)."""
    target = [('e', 0), ('e', 1), 'x1', 'x2', 'x3', 'x4']
    src = list(L0) + list(L1)
    sig = ce.perm_sign([target.index(h) for h in src])
    Bs, Ds = evalx.basis_and_dual(N, False)
    y = [moms[0][c] + moms[1][c] for c in range(2)]
    tot = Fr(0)
    for Xe, Xd in zip(Bs, Ds):
        def arg(h):
            if h == ('e', 0):
                return Xe, [-y[0], -y[1]]
            if h == ('e', 1):
                return Xd, [y[0], y[1]]
            i = int(h[1]) - 1
            return mats[i], moms[i]
        A = [arg(h) for h in L0]
        B = [arg(h) for h in L1]
        va = evalx.eval_cochain(a, [m for m, _ in A], [p for _, p in A], N)
        if va == 0:
            continue
        vb = evalx.eval_cochain(b, [m for m, _ in B], [p for _, p in B], N)
        tot += va * vb
    return sig * tot


def tree_cochain_value(a, b, mats, moms, N, L0, L1):
    tot = Fr(0)
    for perm in itertools.permutations(range(4)):
        tot += ce.perm_sign(perm) * tree_value(a, b, [mats[i] for i in perm], [moms[i] for i in perm], N, L0, L1)
    return tot / 4


N = 3
pairs = [("D", J_cochain('d'), "[p1p2]^2 f", f_bracket_sq()),
         ("D", J_cochain('d'), "D", J_cochain('d')),
         ("s_0", ce.s3_cochain(), "[p1p2]^2 f", f_bracket_sq()),
         ("S1", J_cochain('s1'), "[p1p2]^2 f", f_bracket_sq())]
L0s = list(itertools.permutations(['x1', 'x2', ('e', 0)]))
L1s = list(itertools.permutations([('e', 1), 'x3', 'x4']))
for na, a, nb, b in pairs:
    ratios = set()
    for trial in range(3):
        mats = [evalx.rand_mat(N, 3, rng=rng) for _ in range(4)]
        moms = evalx.rand_moms(4, 4, rng=rng)
        ref = evalx.bracket_explicit(a, 3, b, 3, mats, moms, N)
        ref2 = evalx.eval_cochain(ce.bracket(a, 3, b, 3), mats, moms, N)
        vals = set()
        for L0 in L0s:
            for L1 in L1s:
                vals.add(tree_cochain_value(a, b, mats, moms, N, L0, L1))
        assert len(vals) == 1, vals
        v = vals.pop()
        ratios.add((v / ref) if ref != 0 else ('ref=0', v))
        print("  {%s, %s} at gl_%d, sample %d: tree/U_tree = %s (all 36 local orders equal), shuffle (explicit) = %s, ce.bracket = %s"
              % (na, nb, N, trial, v, ref, ref2), flush=True)
    print("{%s, %s}: tree/U_tree : shuffle bracket = %s" % (na, nb, ratios))
print("=> {F_a, F_b} = U_tree F_{{a,b}}:  eps = U_tree")
