"""Own operator engine for the nilpotency relations of hCS on C^{3|2} (S2 conventions), independent
of opcheck.py / salg.py.

Conventions (S2):
  superfield  C^A(theta) = c^A + theta_1 g1^A + theta_2 g2^A + theta_1 theta_2 b^A  (c = i c_text,
              S2 eq-res-letters), c, b odd, g1, g2 even;
  delta/delta C^A(theta) = theta_1 theta_2 d_c + theta_2 d_g1 - theta_1 d_g2 + d_b  (left derivatives,
              S2 eq-fd), int d^2theta theta_1 theta_2 = 1 (S2 eq-measure), thetas moved to the left;
  hat Q_0 = -i Q_0, odd derivation, hat Q_0 C^A = (1/2) f_{ABC} C^B C^C, Q_0 theta = -theta Q_0
              (S2 eq-Q0), commuting with d_a:
              hat Q_0 c = (1/2) f c c,  hat Q_0 g_i = f c g_i,  hat Q_0 b = f c b + f g1 g2;
  hat Q^I_Gamma = int d^2theta int dz I_Gamma(d; z) T_Gamma(theta) delta_{A_1}(x+z_1) ... delta_{A_m}(x+z_m),
              T_Gamma = prod_v f_{legs(v)} prod_{outputs, vertex order} C^{D_v}(x + X_v, theta),
              legs(v) = (in-edges in edge order, output, out-edges in edge order) (S2 (1.4),(1.6),(1.7));
              delta_m acts first; input jets -> z-derivatives of the kernel at z = 0 (S2 eq-moment).
  Q_0 = i hatQ_0, Q_1 = (i/16pi^3) hatQ^I_box, Q_2 = (i/128pi^6)(hatQ^I_A + hatQ^I_B)  (S2 (1.4),(1.5),(7.15)).

Evaluation: a final polynomial P in letters is paired with test functionals Phi_i = (tau_i, a_i, u_i),
i = 1..N: the algebra homomorphism letter(type tau, colour D, jet m) -> sum_{i: tau_i = tau} eps_i a_i^D u_i^m,
eps_i of the parity of tau, coefficient of eps_1 ... eps_N.  Letters are carried symbolically with tags;
colours by einsum contraction of the network of structure constants; jets by momenta (a letter consumed
by input j of a graph operator carries p_j = sum_v (dX_v/dz_j) lambda_v; by hat Q_0: the sum of the two
output momenta); kernels as bracket polynomials in the final momenta, integrated exactly over r.
"""
import itertools
import math
from fractions import Fraction as Fr

import numpy as np

PAR = {'c': 1, 'b': 1, 'g1': 0, 'g2': 0}
THETA_OUT = {'c': (), 'g1': (1,), 'g2': (2,), 'b': (1, 2)}          # C = c + th1 g1 + th2 g2 + th1 th2 b
THETA_DEL = {'c': ((1, 2), 1), 'g1': ((2,), 1), 'g2': ((1,), -1), 'b': ((), 1)}


# =============================================================================== Lie algebra
def so_n(N):
    idx = {}
    for i in range(N):
        for j in range(i + 1, N):
            idx[(i, j)] = len(idx)
    dim = len(idx)
    gens = []
    for (i, j) in idx:
        M = np.zeros((N, N), dtype=int)
        M[i, j], M[j, i] = 1, -1
        gens.append(M)
    f = np.zeros((dim, dim, dim), dtype=np.int64)
    for B in range(dim):
        for C in range(dim):
            Z = gens[B] @ gens[C] - gens[C] @ gens[B]
            for (i, j), A in idx.items():
                f[A, B, C] = Z[i, j]
    assert np.all(f == -f.transpose(1, 0, 2)) and np.all(f == -f.transpose(0, 2, 1))
    # Jacobi
    jac = (np.einsum('abe,ecd->abcd', f, f) + np.einsum('ace,edb->abcd', f, f) + np.einsum('ade,ebc->abcd', f, f))
    assert np.all(jac == 0)
    return dim, f


# =============================================================================== polynomials in r
class Poly:
    __slots__ = ('t', 'n')

    def __init__(self, t, n):
        self.t, self.n = t, n

    @staticmethod
    def const(c, n):
        return Poly({(0,) * n: Fr(c)} if c else {}, n)

    @staticmethod
    def var(i, n):
        e = [0] * n
        e[i] = 1
        return Poly({tuple(e): Fr(1)}, n)

    def __add__(self, o):
        if not isinstance(o, Poly):
            o = Poly.const(o, self.n)
        t = dict(self.t)
        for k, v in o.t.items():
            w = t.get(k, 0) + v
            if w:
                t[k] = w
            else:
                t.pop(k, None)
        return Poly(t, self.n)
    __radd__ = __add__

    def __neg__(self):
        return Poly({k: -v for k, v in self.t.items()}, self.n)

    def __sub__(self, o):
        return self + (-o if isinstance(o, Poly) else Poly.const(-Fr(o), self.n))

    def __rsub__(self, o):
        return (-self) + o

    def __mul__(self, o):
        if not isinstance(o, Poly):
            o = Fr(o)
            return Poly({k: v * o for k, v in self.t.items()} if o else {}, self.n)
        t = {}
        for k1, v1 in self.t.items():
            for k2, v2 in o.t.items():
                k = tuple(a + b for a, b in zip(k1, k2))
                t[k] = t.get(k, 0) + v1 * v2
        return Poly({k: v for k, v in t.items() if v}, self.n)
    __rmul__ = __mul__

    def __pow__(self, k):
        r = Poly.const(1, self.n)
        for _ in range(k):
            r = r * self
        return r

    def iszero(self):
        return not self.t


def integrate_poly(P, blocks):
    """blocks: list of (kind, [var indices]); kinds: 'u' single var in (0,1); 'o2' 1>b>c>0; 'o3' 1>a>b>c>0."""
    tot = Fr(0)
    for e, c in P.t.items():
        w = Fr(1)
        for kind, vs in blocks:
            ks = [e[i] for i in vs]
            if kind == 'u':
                w /= (ks[0] + 1)
            elif kind == 'o2':
                b, cc = ks
                w /= (cc + 1) * (b + cc + 2)
            else:
                a, b, cc = ks
                w /= (cc + 1) * (b + cc + 2) * (a + b + cc + 3)
        tot += c * w
    return tot


# =============================================================================== bracket polynomials
def lin_add(*terms):
    """terms: (Poly or number, LinMom) -> LinMom (dict label -> Poly)."""
    out = {}
    n = None
    for c, L in terms:
        for k, v in L.items():
            n = v.n
            out[k] = (out[k] + v * c) if k in out else v * c
    return {k: v for k, v in out.items() if not v.iszero()}


def tri_sort(t):
    s = sorted(t)
    perm = [s.index(x) for x in t]
    sg = 1
    for i in range(3):
        for j in range(i + 1, 3):
            if perm[i] > perm[j]:
                sg = -sg
    return tuple(s), sg


def det_lin(L1, L2, L3):
    out = {}
    for a, ca in L1.items():
        for b, cb in L2.items():
            if b == a:
                continue
            cab = ca * cb
            for c, cc in L3.items():
                if c == a or c == b:
                    continue
                key, sg = tri_sort((a, b, c))
                v = cab * cc * sg
                out[key] = out[key] + v if key in out else v
    return {((k,), ()): v for k, v in out.items() if not v.iszero()}


def bp_mul(P, Q):
    out = {}
    for (t1, l1), c1 in P.items():
        for (t2, l2), c2 in Q.items():
            key = (tuple(sorted(t1 + t2)), tuple(sorted(l1 + l2)))
            v = c1 * c2
            out[key] = out[key] + v if key in out else v
    return {k: v for k, v in out.items() if not (isinstance(v, Poly) and v.iszero()) and v != 0}


def bp_scale(P, c):
    return {k: v * c for k, v in P.items()}


def lin_factor(L, comp):
    """linear factor (L)_comp as a bracket polynomial."""
    return {((), ((k, comp),)): v for k, v in L.items()}


def bp_integrate(P, blocks):
    out = {}
    for k, v in P.items():
        w = integrate_poly(v, blocks)
        if w:
            out[k] = w
    return out


def bp_substitute(K, sub, n):
    """K: bracket poly (Fraction coefficients) in labels; sub: label -> LinMom (constant Polys)."""
    out = {}
    for (tris, lins), c in K.items():
        acc = {((), ()): Poly.const(c, n)}
        for t in tris:
            acc = bp_mul(acc, det_lin(*(sub[x] for x in t)))
        for (lab, comp) in lins:
            acc = bp_mul(acc, lin_factor(sub[lab], comp))
        for k, v in acc.items():
            out[k] = out[k] + v if k in out else v
    return {k: integrate_poly(v, []) for k, v in out.items() if not v.iszero()}


# =============================================================================== graph operators
GRAPH_SPECS = {
    'box': dict(edges=[('O', 1), (1, 2), (2, 3), (3, 'O')], doms=['o3']),
    'A': dict(edges=[('O', 1), (1, 2), (2, 3), ('O', 3), (3, 4), (4, 5), (5, 'O')], doms=['chain', 'o3']),
    'B': dict(edges=[('O', 1), (1, 3), ('O', 2), (2, 3), (3, 4), (4, 5), (5, 'O')], doms=['bubble', 'o3']),
    'G1': dict(edges=[('O', 1), (1, 2), (2, 3), ('O', 3), (3, 4), (4, 5), ('O', 5), (5, 6), (6, 7), (7, 'O')],
               doms=['chain', 'chain', 'o3']),
    'G2': dict(edges=[('O', 1), (1, 2), (2, 3), ('O', 3), (3, 4), (4, 5), (5, 6), (6, 'O'), (5, 7), (7, 'O')],
               doms=['chain', 'bubble', 'o3']),
    'G3': dict(edges=[('O', 1), (1, 2), (2, 3), ('O', 3), (3, 4), (4, 5), (5, 'O'), (4, 6), (6, 7), (7, 'O')],
               doms=['chain', 'bubble', 'o3']),
    'G4': dict(edges=[('O', 1), (1, 3), ('O', 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 'O'), (5, 7), (7, 'O')],
               doms=['bubble', 'bubble', 'o3']),
    'G5': dict(edges=[('O', 1), (1, 3), ('O', 2), (2, 3), (3, 4), (4, 5), (5, 'O'), (4, 6), (6, 7), (7, 'O')],
               doms=['bubble', 'bubble', 'o3']),
}
def build_spec(_g):
    E = _g['edges']
    verts = sorted({x for e in E for x in e if x != 'O'})
    _g['verts'] = verts
    _g['spokes'] = [i for i, e in enumerate(E) if 'O' in e]
    deg = {v: sum(v in e for e in E) for v in verts}
    _g['outputs'] = [v for v in verts if deg[v] == 2]
    legs = {}
    for v in verts:
        ins = [i for i, e in enumerate(E) if e[1] == v]
        outs = [i for i, e in enumerate(E) if e[0] == v]
        L = [('e', i) for i in ins] + ([('out', v)] if deg[v] == 2 else []) + [('e', i) for i in outs]
        assert len(L) == 3
        legs[v] = L
    _g['legs'] = legs
    _g['nr'] = 3 * len(_g['doms'])
    return _g


for _n, _g in GRAPH_SPECS.items():
    build_spec(_g)


def blocks_for(name, off):
    out = []
    for m, kind in enumerate(GRAPH_SPECS[name]['doms']):
        i = off + 3 * m
        if kind == 'chain':
            out += [('u', [i]), ('o2', [i + 1, i + 2])]
        elif kind == 'bubble':
            out += [('u', [i]), ('u', [i + 1]), ('u', [i + 2])]
        else:
            out += [('o3', [i, i + 1, i + 2])]
    return out


def density_and_X(name, lam, n, off):
    """J_Gamma with the I-sign (I_Gamma = int J) as bracket poly with Poly coefficients, and
    c[v][j] = dX_v/dz_j as Polys.  lam: dict output vertex -> LinMom."""
    spec = GRAPH_SPECS[name]
    if 'base' in spec:
        # relabelled graph: I' = eps * I_base, vertex map old -> new, spoke map new j -> old j
        lam_old = {v: lam[spec['vmap'][v]] for v in GRAPH_SPECS[spec['base']]['outputs']}
        J, c = density_and_X(spec['base'], lam_old, n, off)
        J = bp_scale(J, spec['eps'])
        cn = {spec['vmap'][v]: {jn + 1: cv[spec['smap'][jn] + 1] for jn in range(len(spec['smap']))}
              for v, cv in c.items()}
        return J, cn
    R = [Poly.var(off + i, n) for i in range(GRAPH_SPECS[name]['nr'])]
    one = Poly.const(1, n)
    if name == 'box':
        J = det_lin(lam[1], lam[2], lam[3])
        c = {w: {1: R[w - 1], 2: one - R[w - 1]} for w in (1, 2, 3)}
        return J, c
    if name in ('A', 'B'):
        r1, r2, r3, r4, r5, r6 = R
        d1 = det_lin(lam[1], lam[2], lin_add((r5, lam[4]), (r6, lam[5])))
        d2 = det_lin(lin_add((one - r2, lam[1]), (one - r3, lam[2])), lam[4], lam[5])
        pref = r1 * r1 if name == 'A' else r1 * (one - r1) * (-1)       # I_B = - int J_B
        J = bp_mul(bp_mul({((), ()): pref}, d1), d2)
        # S2 eq-res-X: W = (1-r1) z1 + r1 z2 - z3, Y = z3 + r4 W
        W = {1: one - r1, 2: r1, 3: -one}
        Y = {1: r4 * W[1], 2: r4 * W[2], 3: one + r4 * W[3]}
        c = {1: {j: r2 * (j == 1) + (one - r2) * Y[j] for j in (1, 2, 3)},
             4: {j: Poly.const(int(j == 3), n) + r5 * W[j] for j in (1, 2, 3)},
             5: {j: Poly.const(int(j == 3), n) + r6 * W[j] for j in (1, 2, 3)}}
        if name == 'A':
            c[2] = {j: r3 * (j == 1) + (one - r3) * Y[j] for j in (1, 2, 3)}
        else:
            c[2] = {j: r3 * (j == 2) + (one - r3) * Y[j] for j in (1, 2, 3)}
        return J, c
    r1, r2, r3, r4, r5, r6, r7, r8, r9 = R
    l = lam
    mu = lin_add((one - r2, l[1]), (one - r3, l[2]))
    if name == 'G1':
        pref = r1 * r1 * r4 * r4
        dets = [det_lin(l[1], l[2], lin_add((r6 + (one - r4) * (one - r6) * r7, l[4]), ((one - r4) * r8, l[6]), ((one - r4) * r9, l[7]))),
                det_lin(mu, l[4], lin_add((r8, l[6]), (r9, l[7]))),
                det_lin(lin_add((one - r5, mu), (one - r6, l[4])), l[6], l[7])]
    elif name in ('G2', 'G4'):
        pref = r1 * r1 * r4 * (one - r4) if name == 'G2' else r1 * (one - r1) * r4 * (one - r4) * (-1)
        dets = [det_lin(l[1], l[2], lin_add((r8, l[4]), (r9 * r5, l[6]), (r9 * r6, l[7]))),
                det_lin(lin_add((one - r7, mu), (one - r8, l[4])), l[6], l[7]),
                det_lin(mu, l[4], lin_add((r5, l[6]), (r6, l[7])))]
    else:
        pref = r1 * r1 * r4 * (one - r4) * (one - r4) if name == 'G3' else r1 * (one - r1) * r4 * (one - r4) * (one - r4) * (-1)
        dets = [det_lin(l[1], l[2], lin_add((r7 * (one - r6), l[5]), (r8, l[6]), (r9, l[7]))),
                det_lin(mu, l[5], lin_add((r8, l[6]), (r9, l[7]))),
                det_lin(lin_add((one - r5, mu), (one - r6, l[5])), l[6], l[7])]
    J = {((), ()): pref}
    for d in dets:
        J = bp_mul(J, d)
    # C4 displacement maps: coefficient of z_j in X_v
    Zc = {j: {jj: Poly.const(int(jj == j), n) for jj in (1, 2, 3, 4)} for j in (1, 2, 3, 4)}

    def comb(*terms):
        return {j: sum((c * v[j] for c, v in terms), Poly.const(0, n)) for j in (1, 2, 3, 4)}
    z1, z2, z3, z4 = Zc[1], Zc[2], Zc[3], Zc[4]
    W1 = comb((one - r1, z1), (r1, z2))
    if name == 'G1':
        W2 = comb((one - r4, W1), (r4, z3))
        P5 = comb((one - r7, z4), (r7, W2))
        P3 = comb((r5, W1), (one - r5, P5))
        X = {1: comb((r2, z1), (one - r2, P3)), 2: comb((r3, z1), (one - r3, P3)), 4: comb((r6, W1), (one - r6, P5)),
             6: comb((one - r8, z4), (r8, W2)), 7: comb((one - r9, z4), (r9, W2))}
    elif name in ('G2', 'G4'):
        W2 = comb((one - r4, z3), (r4, z4))
        P3 = comb((r7, W1), (one - r7, W2))
        P5 = comb((r9, W1), (one - r9, W2))
        X = {1: comb((r2, z1), (one - r2, P3)),
             2: comb((r3, z1 if name == 'G2' else z2), (one - r3, P3)),
             4: comb((r8, W1), (one - r8, W2)), 6: comb((r5, P5), (one - r5, z3)), 7: comb((r6, P5), (one - r6, z4))}
    else:
        W2 = comb((one - r4, W1), (r4, z3))
        P4 = comb((one - r7, z4), (r7, W2))
        P3 = comb((r5, W1), (one - r5, P4))
        X = {1: comb((r2, z1), (one - r2, P3)),
             2: comb((r3, z1 if name == 'G3' else z2), (one - r3, P3)),
             5: comb((r6, z3), (one - r6, P4)), 6: comb((one - r8, z4), (r8, W2)), 7: comb((one - r9, z4), (r9, W2))}
    return J, X


# =============================================================================== Grassmann layer
def canon(factors):
    """factors: list of ('th', i) / ('L', tag, type).  Moves thetas to the left (Koszul sign);
    returns (sign, letters) for theta-part theta_1 theta_2, else (0, None)."""
    sign = 1
    nodd_letters = 0
    ths = []
    letters = []
    for f in factors:
        if f[0] == 'th':
            if nodd_letters & 1:
                sign = -sign
            ths.append(f[1])
        else:
            letters.append((f[1], f[2]))
            nodd_letters += PAR[f[2]]
    if ths == [1, 2]:
        return sign, letters
    if ths == [2, 1]:
        return -sign, letters
    return 0, None


def apply_graph(name, prefix, mono, binding):
    """hat Q^I_name (Grassmann part) on the ordered monomial mono = [(tag, type)], spoke j bound to
    the letter with tag binding[j].  Returns list of (sign, letters) with output tags (prefix, v)."""
    g = GRAPH_SPECS[name]
    cur = [('L', t, ty) for t, ty in mono]
    sign = 1
    for tag in reversed(binding):
        pos = next(i for i, f in enumerate(cur) if f[0] == 'L' and f[1] == tag)
        ty = cur[pos][2]
        if PAR[ty]:
            before = sum(1 for f in cur[:pos] if f[0] == 'th') + sum(PAR[f[2]] for f in cur[:pos] if f[0] == 'L')
            if before & 1:
                sign = -sign
        del cur[pos]
        ths, s = THETA_DEL[ty]
        sign *= s
        cur = [('th', i) for i in ths] + cur
    nth_in = sum(1 for f in cur if f[0] == 'th')
    out = []
    outs = g['outputs']
    for types in itertools.product(('c', 'g1', 'g2', 'b'), repeat=len(outs)):
        if nth_in + sum(len(THETA_OUT[t]) for t in types) != 2:
            continue
        fac = []
        for v, t in zip(outs, types):
            fac += [('th', i) for i in THETA_OUT[t]] + [('L', (prefix, v), t)]
        s, letters = canon(fac + cur)
        if s:
            out.append((sign * s, letters))
    return out


Q0_RULE = {'c': [(Fr(1, 2), 'c', 'c')], 'g1': [(Fr(1), 'c', 'g1')], 'g2': [(Fr(1), 'c', 'g2')],
           'b': [(Fr(1), 'c', 'b'), (Fr(1), 'g1', 'g2')]}


def apply_q0(prefix, mono, tag):
    pos = next(i for i, (t, _) in enumerate(mono) if t == tag)
    s = -1 if sum(PAR[ty] for _, ty in mono[:pos]) & 1 else 1
    ty = mono[pos][1]
    out = []
    for c, t1, t2 in Q0_RULE[ty]:
        new = mono[:pos] + [((prefix, 1), t1), ((prefix, 2), t2)] + mono[pos + 1:]
        out.append((s * c, new))
    return out


# =============================================================================== compositions
class Op:
    """An operator instance in a composition: kind 'q0' or a graph name; prefix for output tags."""

    def __init__(self, kind, prefix):
        self.kind, self.prefix = kind, prefix

    def nin(self):
        return 1 if self.kind == 'q0' else len(GRAPH_SPECS[self.kind]['spokes'])

    def outtags(self):
        if self.kind == 'q0':
            return [(self.prefix, 1), (self.prefix, 2)]
        return [(self.prefix, v) for v in GRAPH_SPECS[self.kind]['outputs']]

    def apply(self, mono, binding):
        if self.kind == 'q0':
            return apply_q0(self.prefix, mono, binding[0])
        res = apply_graph(self.kind, self.prefix, mono, binding)
        return [(Fr(s), l) for s, l in res]


def histories(ops, O):
    """ops: list applied in order (first acts first on O).  O: list of (tag, type).
    Yields (bindings, terms) with terms = list of (coef, final monomial)."""
    def rec(k, terms, bindings):
        if k == len(ops):
            yield bindings, terms
            return
        op = ops[k]
        tags = [t for t, _ in terms[0][1]]
        for bnd in itertools.permutations(tags, op.nin()):
            new = []
            for coef, mono in terms:
                for c, m in op.apply(mono, list(bnd)):
                    new.append((coef * c, m))
            # combine identical monomials
            acc = {}
            for c, m in new:
                key = tuple(m)
                acc[key] = acc.get(key, 0) + c
            new = [(c, list(k2)) for k2, c in acc.items() if c != 0]
            if new:
                yield from rec(k + 1, new, bindings + [bnd])
    # all monomials at each stage share the same tag set; group by tag set
    yield from rec(0, [(Fr(1), list(O))], [])


def colour_network(ops, bindings, O_colours, final_tags, dim, f, A):
    """Returns G[i_1..i_N] = contraction with A[i, D] on the final letters (order final_tags)."""
    operands, subs = [], []
    names = {}

    def nm(x):
        if x not in names:
            names[x] = len(names)
        return names[x]
    for k, (op, bnd) in enumerate(zip(ops, bindings)):
        if op.kind == 'q0':
            operands.append(f)
            subs.append([nm(('col', bnd[0])), nm(('col', (op.prefix, 1))), nm(('col', (op.prefix, 2)))])
            continue
        g = GRAPH_SPECS[op.kind]
        spoke_of = {e: j for j, e in enumerate(g['spokes'])}
        for v in g['verts']:
            idx = []
            for leg in g['legs'][v]:
                if leg[0] == 'out':
                    idx.append(nm(('col', (op.prefix, leg[1]))))
                elif leg[1] in spoke_of:
                    idx.append(nm(('col', bnd[spoke_of[leg[1]]])))
                else:
                    idx.append(nm(('int', k, leg[1])))
            operands.append(f)
            subs.append(idx)
    # For a fixed tuple of output probes, bound the absolute colour sum.
    # Each new vertex sums its unassigned indices at fixed earlier indices.
    # Taking the largest such slice gives an upper bound for every prefix.
    pending = list(subs)
    known = set()
    bound = 1
    absolute_f = np.abs(f)
    while pending:
        indices = max(pending, key=lambda group: sum(i in known for i in group))
        pending.remove(indices)
        free = tuple(j for j, i in enumerate(indices) if i not in known)
        bound *= int(absolute_f.sum(axis=free).max())
        known.update(indices)
    outsub = []
    for tag in final_tags:
        operands.append(A)
        subs.append([nm(('fun', tag)), nm(('col', tag))])
        outsub.append(nm(('fun', tag)))
    for tag, col in O_colours.items():
        assert ('col', tag) in names
        if isinstance(col, (int, np.integer)):
            vec = np.zeros(dim, dtype=np.int64)
            vec[col] = 1
        else:
            vec = np.asarray(col)
        operands.append(vec)
        subs.append([nm(('col', tag))])
    bound *= dim ** len(set(names.values()) - known - set(outsub))
    for operand in operands:
        if operand is not f:
            bound *= int(np.abs(operand).max())
    # int64 contraction is arithmetic modulo 2**64. The bound makes the
    # signed representative unique, even if an intermediate wrapped.
    dtype = np.int64 if bound < 2 ** 63 else object
    args = []
    for operand, indices in zip(operands, subs):
        integers = operand.astype(np.int64, copy=False)
        if not np.array_equal(operand, integers):
            raise ValueError('colour probes must have integral entries')
        args += [integers.astype(dtype, copy=False), indices]
    return np.einsum(*args, outsub, optimize='greedy')


def kernel_bp(ops, bindings, final_tags, jets, cache):
    """Exact kernel of a composition as a bracket polynomial (Fraction coefficients) in the final tags.
    jets: dict O-tag -> list of components (one derivative each)."""
    # only graph operators carry r variables
    gops = [(k, op) for k, op in enumerate(ops) if op.kind != 'q0']
    offs, n = {}, 0
    for k, op in gops:
        offs[k] = n
        n += GRAPH_SPECS[op.kind]['nr']
    outs_all = {t for op in ops for t in op.outtags()}
    key = (tuple(op.kind for op in ops),
           tuple(tuple(t if (t in outs_all or t in jets) else None for t in b) for b in bindings),
           tuple(final_tags), tuple(sorted((k, tuple(v)) for k, v in jets.items())))
    if key in cache:
        return cache[key]
    one = Poly.const(1, n)
    final = set(final_tags)
    mom = {t: {t: one} for t in final_tags}          # momentum (LinMom) of each tag
    total = {((), ()): one}
    blocks = []
    # process from the last operator backwards
    for k in range(len(ops) - 1, -1, -1):
        op, bnd = ops[k], bindings[k]
        if op.kind == 'q0':
            p = lin_add((one, mom[(op.prefix, 1)]), (one, mom[(op.prefix, 2)]))
            mom[bnd[0]] = p
            continue
        lam = {v: mom[(op.prefix, v)] for v in GRAPH_SPECS[op.kind]['outputs']}
        J, c = density_and_X(op.kind, lam, n, offs[k])
        total = bp_mul(total, J)
        for j, tag in enumerate(bnd):
            pj = {}
            for v, cv in c.items():
                pj = lin_add((one, pj), (cv[j + 1], lam[v]))
            mom[tag] = pj
        blocks += blocks_for(op.kind, offs[k])
    for tag, comps in jets.items():
        for a in comps:
            total = bp_mul(total, lin_factor(mom[tag], a))
    res = bp_integrate(total, blocks)
    cache[key] = res
    return res


def bijections(types_final, ftypes):
    """all type-respecting bijections positions -> functional indices, with the Koszul sign of
    the odd letters (positions in final order) sorted by functional index."""
    N = len(types_final)
    pools = {}
    for i, t in enumerate(ftypes):
        pools.setdefault(t, []).append(i)
    pos_by = {}
    for p, t in enumerate(types_final):
        pos_by.setdefault(t, []).append(p)
    if sorted(types_final) != sorted(ftypes):
        return None, None
    groups = sorted(pos_by)
    perms = [list(itertools.permutations(pools[t])) for t in groups]
    S, sg = [], []
    for choice in itertools.product(*perms):
        s = [0] * N
        for t, ch in zip(groups, choice):
            for p, i in zip(pos_by[t], ch):
                s[p] = i
        odd = [s[p] for p in range(N) if PAR[types_final[p]]]
        inv = sum(1 for a in range(len(odd)) for b in range(a + 1, len(odd)) if odd[a] > odd[b])
        S.append(s)
        sg.append(-1 if inv & 1 else 1)
    return np.array(S, dtype=np.int64), np.array(sg, dtype=np.int64)


def eval_bp(K, S, final_tags, D3, U):
    """sum over keys of coeff * prod brackets, for each bijection row of S.
    Returns (integer numerators as object array, common denominator)."""
    pos = {t: i for i, t in enumerate(final_tags)}
    if not K:
        return np.zeros(len(S), dtype=object), 1
    den = 1
    for c in K.values():
        den = den * c.denominator // math.gcd(den, c.denominator)
    M3 = int(np.abs(D3).max()) or 1
    MU = int(np.abs(U).max()) or 1
    bound = sum(abs(int(c * den)) * M3 ** len(t) * MU ** len(l) for (t, l), c in K.items())
    dt = np.int64 if bound < 2 ** 62 else object
    acc = np.zeros(len(S), dtype=dt)
    for (tris, lins), c in K.items():
        v = np.full(len(S), int(c * den), dtype=dt)
        for t in tris:
            v = v * D3[S[:, pos[t[0]]], S[:, pos[t[1]]], S[:, pos[t[2]]]].astype(dt)
        for lab, comp in lins:
            v = v * U[S[:, pos[lab]], comp].astype(dt)
        acc = acc + v
    return acc.astype(object), den


def pair(ops, O, O_colours, jets, funcs, dim, f, cache=None, verbose=False):
    """<Phi, ops[-1] ... ops[0] O> for test functionals funcs = (types, A (N x dim), u (N x 3))."""
    if cache is None:
        cache = {}
    ftypes, A, u = funcs
    N = len(ftypes)
    D3 = np.zeros((N, N, N), dtype=np.int64)
    for i, j, k in itertools.product(range(N), repeat=3):
        x, y, z = (tuple(map(int, u[t])) for t in (i, j, k))
        D3[i, j, k] = (x[0] * (y[1]*z[2] - y[2]*z[1])
                       + x[1] * (y[2]*z[0] - y[0]*z[2])
                       + x[2] * (y[0]*z[1] - y[1]*z[0]))
    U = np.array(u, dtype=np.int64)
    total = Fr(0)
    nh = 0
    for bindings, terms in histories(ops, O):
        final_tags = sorted(t for t, _ in terms[0][1])
        # per term: reorder to final_tags with Koszul sign, collect types
        relevant = []
        for coef, mono in terms:
            types = dict(mono)
            tf = [types[t] for t in final_tags]
            if sorted(tf) != sorted(ftypes):
                continue
            order = [final_tags.index(t) for t, _ in mono]
            odd = [o for (t, ty), o in zip(mono, order) if PAR[ty]]
            inv = sum(1 for a in range(len(odd)) for b in range(a + 1, len(odd)) if odd[a] > odd[b])
            relevant.append((coef * (-1 if inv & 1 else 1), tf))
        if not relevant:
            continue
        nh += 1
        G = colour_network(ops, bindings, O_colours, final_tags, dim, f, A)
        K = kernel_bp(ops, bindings, final_tags, jets, cache)
        for coef, tf in relevant:
            S, sg = bijections(tf, ftypes)
            gv = G[tuple(S[:, p] for p in range(N))]
            kv, den = eval_bp(K, S, final_tags, D3, U)
            total += coef * Fr(int(np.sum((sg * gv).astype(object) * kv)), den)
    if verbose:
        print('      histories with matching final types: %d' % nh)
    return total
