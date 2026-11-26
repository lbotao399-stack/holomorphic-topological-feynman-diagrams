"""Exact evaluation of compositions of the operators hat Q_0, hat Q_Gamma of hCS on C^{3|2}
(S2 conventions, code of S2: s2code/salg.py, ops.py, ops2.py) paired with plane-wave test
functionals.

hat Q_Gamma (graph with m spokes):  int d^2theta  T_Gamma  delta_{A_1} ... delta_{A_m},   kernel I_Gamma,
  T_Gamma = prod_v C_{legs(v)}  prod_{outputs v} C^{D_v}(x + X_v)   (vertex order),
  legs(v) = (incoming edges in edge order, output D_v if deg v = 2, outgoing edges in edge order),
  spoke j carries A_j, tree edges carry summed indices.  For the box, (A), (B) this is S2 (1.4), (1.6), (1.7).
  Kernel I_Gamma = S2 eq-kernel-def (signed): box +det, I_A = +int J_A, I_B = -int J_B, three loops: q3lib.

Test functional Phi = (a_i, u_i), i = 1..5: a final letter c^D with jet m in a polynomial P maps to
sum_i a_i^D u_i^m eps_i (odd eps_i); <P, Phi> = coefficient of eps_1 ... eps_5.  On compositions the
letters are 'slot' letters (t, slot, 0, 0, 0) of salg: the Grassmann and theta signs are computed by the
S2 algebra itself; colours by contraction of the networks of structure constants; jets by the
momentum rule (input letter bound to delta_j of an operator with displacement map X: lambda ->
sum_v (dX_v/dz_j) lambda_v; hat Q_0: lambda -> lambda_n1 + lambda_n2).
"""
import itertools
import math
import os
import random
import sys
from fractions import Fraction as Fr

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 's2code'))
sys.path.insert(0, HERE)

import numpy as np
import sympy as sp

from salg import add, mul, scale, csuper, sderiv, theta_coeff, lderiv, letters_of, gen, const, structure_constants
import q3lib
from q3lib import Reduction, LAM, Z, R, degree

SIGNS = {'box': 1, 'A': -1, 'B': 1, 'G1': 1, 'G2': -1, 'G3': 1, 'G4': 1, 'G5': -1}   # check_kernels.py
AUT = {'G1': 2, 'G2': 2, 'G3': 1, 'G4': 8, 'G5': 2, 'A': 2, 'B': 2, 'box': 2}


def perm_sign(seq):
    seq = list(seq)
    s = 1
    for i in range(len(seq)):
        while seq[i] != i:
            j = seq[i]
            seq[i], seq[j] = seq[j], seq[i]
            s = -s
    return s


# ------------------------------------------------------------------------------ bracket algebra
def vec_add(a, b, cb=1):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0) + cb * v
    return {k: v for k, v in out.items() if v != 0}


def det_expand(V1, V2, V3):
    """det(V1, V2, V3) for V = {label: coefficient}:  {sorted label triple: coefficient}."""
    out = {}
    for l1, c1 in V1.items():
        for l2, c2 in V2.items():
            if l2 == l1:
                continue
            for l3, c3 in V3.items():
                if l3 == l1 or l3 == l2:
                    continue
                trip = (l1, l2, l3)
                key = tuple(sorted(trip))
                s = perm_sign([key.index(x) for x in trip])
                out[key] = out.get(key, 0) + s * c1 * c2 * c3
    return {k: sp.expand(v) for k, v in out.items() if sp.expand(v) != 0}


def bp_mul(P, Q):
    out = {}
    for k1, v1 in P.items():
        for k2, v2 in Q.items():
            key = tuple(sorted(k1 + k2))
            out[key] = out.get(key, 0) + v1 * v2
    return out


def bp_from_det(D):
    return {(k,): v for k, v in D.items()}


def integrate_dom(expr, dom):
    """dom: list of (kind, symbols) with kinds cube, ord2, ord3."""
    syms = [s for _, ss in dom for s in ss]
    P = sp.Poly(sp.expand(expr), *syms) if syms else None
    if P is None:
        return sp.Rational(expr)
    idx = {s: i for i, s in enumerate(syms)}
    tot = sp.Integer(0)
    for mon, c in P.terms():
        w = sp.Integer(1)
        for kind, ss in dom:
            ks = [mon[idx[s]] for s in ss]
            if kind == 'cube':
                w *= sp.Rational(1, ks[0] + 1)
            elif kind == 'ord2':
                a, b = ks
                w *= sp.Rational(1, (b + 1) * (a + b + 2))
            else:
                a, b, c3 = ks
                w *= sp.Rational(1, (c3 + 1) * (b + c3 + 2) * (a + b + c3 + 3))
        tot += c * w
    return tot


# ------------------------------------------------------------------------------ operator specs
class GraphOp:
    """Operator data of a graph: legs, outputs, density factors (Vecs over output vertices),
    domain, displacement coefficients c[v][j] = d X_v / d z_j."""

    def __init__(self, name):
        red = Reduction(name)
        self.name = name
        self.red = red
        root = red.root
        E = red.edges
        self.m = len(red.spokes)
        self.outs = [v for v in range(root) if degree(E, v) == 2]
        spoke_of = {e: j for j, (e, s) in enumerate(red.spokes)}
        self.legs = []
        for v in range(root):
            ins = [e for e, (a, b) in enumerate(E) if b == v]
            outs_ = [e for e, (a, b) in enumerate(E) if a == v]
            L = []
            for e in ins:
                L.append(('A', spoke_of[e]) if e in spoke_of else ('R', e))
            if v in self.outs:
                L.append(('D', v))
            for e in outs_:
                L.append(('A', spoke_of[e]) if e in spoke_of else ('R', e))
            assert len(L) == 3
            self.legs.append(L)
        zero = {LAM[v][a]: 0 for v in range(root) if v not in self.outs for a in range(3)}

        def lin(vexpr):
            out = {}
            for u in self.outs:
                c = sp.expand(sp.diff(vexpr[0].subs(zero), LAM[u][0]))
                if c != 0:
                    out[u] = c
            return out
        self.factors = []          # (prefactor, [Vec, Vec, Vec]) with Vec over output vertices
        for st in red.steps:
            if st['kind'] == 'box':
                self.factors.append((sp.Integer(1), [lin(l) for l in st['lam']]))
            else:
                r1 = st['r'][0]
                pre = r1 ** 2 if st['kind'] == 'chain' else r1 * (1 - r1)
                self.factors.append((pre, [lin(st['Y']), lin(st['lv']), lin(st['lw'])]))
        self.sign = SIGNS[name]
        self.nr = 3 * len(red.steps)
        self.dom = red.domain()
        X, _ = red.displacement()
        self.c = {v: [sp.expand(sp.diff(X[v][0], Z[j][0])) for j in range(self.m)] for v in self.outs}

    def instance(self, k):
        """Ratio symbols renamed for the k-th application in a composition."""
        sub = {R[i]: sp.Symbol('r%d_%d' % (i + 1, k)) for i in range(self.nr)}
        facs = [(pre.subs(sub), [{l: c.subs(sub) for l, c in V.items()} for V in Vs]) for pre, Vs in self.factors]
        dom = [(kind, tuple(sub[s] for s in ss)) for kind, ss in self.dom]
        cc = {v: [x.subs(sub) for x in row] for v, row in self.c.items()}
        return facs, dom, cc


OPS = {}


def op(name):
    if name not in OPS:
        OPS[name] = GraphOp(name)
    return OPS[name]


# ------------------------------------------------------------------------------ skeleton algebra
def sk_letter(t, slot):
    return (t, slot, 0, 0, 0)


def apply_graph_skel(P, bind, outslots):
    """hat Q_Gamma on the skeleton polynomial P: delta_j bound to slot bind[j] (j = 1..m; delta_m acts
    first), outputs csuper(outslots[k]) in vertex order, then the theta_1 theta_2 coefficient."""
    d = P
    for j in reversed(range(len(bind))):
        s = bind[j]
        nd = {}
        for l in letters_of(d):
            if l[0] != 0 and l[1] == s:
                nd = add(nd, sderiv(d, l))
        d = nd
        if not d:
            return {}
    T = const(1)
    for s in outslots:
        T = mul(T, csuper(s))
    return theta_coeff(mul(T, d), '12')


def apply_q0_skel(P, s, n1, n2):
    """hat Q_0 acting on the letter in slot s: component of csuper(n1) csuper(n2) (colour (1/2) f(s,n1,n2))."""
    Rr = mul(csuper(n1), csuper(n2))
    comp = {1: theta_coeff(Rr, '0'), 3: scale(theta_coeff(Rr, '1'), -1),
            4: scale(theta_coeff(Rr, '2'), -1), 2: theta_coeff(Rr, '12')}
    out = {}
    for l in letters_of(P):
        if l[0] != 0 and l[1] == s:
            out = add(out, mul(comp[l[0]], lderiv(P, l)))
    return out


# ------------------------------------------------------------------------------ histories
class History:
    """A composition applied to the skeleton operator O.  Records the skeleton polynomial, the colour
    network, and the momentum structure."""

    def __init__(self, P, slots_O, counter):
        self.P = P
        self.tensors = []          # list of (coefficient, legs)  (f-tensors)
        self.pref = Fr(1)
        self.counter = counter
        self.apps = []             # list of (kind, name, bind, outslots, instance index)
        self.producer = {}         # slot -> (app index, output vertex) ; O slots absent
        self.consumer = {}         # slot -> (app index, j)
        self.Oslots = list(slots_O)

    def new(self):
        self.counter += 1
        return self.counter

    def copy(self):
        h = History(dict(self.P), self.Oslots, self.counter)
        h.tensors = list(self.tensors)
        h.pref = self.pref
        h.apps = list(self.apps)
        h.producer = dict(self.producer)
        h.consumer = dict(self.consumer)
        return h

    def live_slots(self):
        return sorted({l[1] for l in letters_of(self.P) if l[0] != 0})


def step_q0(hs):
    out = []
    for h in hs:
        for s in h.live_slots():
            g = h.copy()
            n1, n2 = g.new(), g.new()
            g.P = apply_q0_skel(h.P, s, n1, n2)
            if not g.P:
                continue
            k = len(g.apps)
            g.apps.append(('Q0', None, [s], [n1, n2], k))
            g.tensors.append((s, n1, n2))
            g.pref *= Fr(1, 2)
            g.consumer[s] = (k, 0)
            g.producer[n1] = (k, 0)
            g.producer[n2] = (k, 1)
            out.append(g)
    return out


def step_graph(hs, name):
    G = op(name)
    out = []
    for h in hs:
        live = h.live_slots()
        for bind in itertools.permutations(live, G.m):
            g = h.copy()
            outslots = [g.new() for _ in G.outs]
            P = apply_graph_skel(h.P, list(bind), outslots)
            if not P:
                continue
            g.P = P
            k = len(g.apps)
            g.apps.append(('G', name, list(bind), outslots, k))
            omap = dict(zip(G.outs, outslots))
            internal = {}
            for v, L in enumerate(G.legs):
                ids = []
                for typ, x in L:
                    if typ == 'A':
                        ids.append(bind[x])
                    elif typ == 'D':
                        ids.append(omap[x])
                    else:
                        if x not in internal:
                            internal[x] = g.new()
                        ids.append(internal[x])
                g.tensors.append(tuple(ids))
            for j, s in enumerate(bind):
                g.consumer[s] = (k, j)
            for v, s in omap.items():
                g.producer[s] = (k, v)
            out.append(g)
    return out


# ------------------------------------------------------------------------------ kernels of histories
_KCACHE = {}


def structure_key(h):
    """Momentum structure: for each application, which outputs are consumed by which (app, j)."""
    key = []
    for (kind, name, bind, outslots, k) in h.apps:
        cons = tuple(h.consumer.get(s) for s in outslots)
        key.append((kind, name, cons))
    return tuple(key)


def history_kernel(h):
    """Composite kernel as a bracket polynomial in the final labels ('o', app, output index)."""
    key = structure_key(h)
    if key in _KCACHE:
        return _KCACHE[key]
    # momentum of every output slot, as Vec over final labels, computed from the last application back
    mom = {}
    inst = {}
    for (kind, name, bind, outslots, k) in h.apps:
        if kind == 'G':
            inst[k] = op(name).instance(k)

    def m_of(k, idx):
        """Vec of output idx of application k."""
        kind, name, bind, outslots, _ = h.apps[k]
        s = outslots[idx]
        if s not in h.consumer:
            return {('o', k, idx): sp.Integer(1)}
        k2, j = h.consumer[s]
        kind2, name2, bind2, outslots2, _ = h.apps[k2]
        if kind2 == 'Q0':
            return vec_add(m_of(k2, 0), m_of(k2, 1))
        G2 = op(name2)
        cc = inst[k2][2]
        V = {}
        for i, v in enumerate(G2.outs):
            V = vec_add(V, {lab: cc[v][j] * c for lab, c in m_of(k2, i).items()})
        return V

    BP = {(): sp.Integer(1)}
    doms = []
    for (kind, name, bind, outslots, k) in h.apps:
        if kind != 'G':
            continue
        G = op(name)
        facs, dom, cc = inst[k]
        doms += dom
        moms = {v: m_of(k, i) for i, v in enumerate(G.outs)}
        for pre, Vs in facs:
            Ws = []
            for V in Vs:
                W = {}
                for v, c in V.items():
                    W = vec_add(W, {lab: c * cm for lab, cm in moms[v].items()})
                Ws.append(W)
            D = det_expand(*Ws)
            BP = bp_mul(BP, {kk: pre * vv for kk, vv in bp_from_det(D).items()})
        BP = {kk: vv * G.sign for kk, vv in BP.items()}
    res = {}
    for kk, vv in BP.items():
        val = integrate_dom(vv, doms)
        if val != 0:
            res[kk] = Fr(int(val.p), int(val.q))
    _KCACHE[key] = res
    return res


# ------------------------------------------------------------------------------ evaluation
class TestData:
    def __init__(self, dim, seed, amax=2, umax=3, n=5):
        rng = random.Random(seed)
        self.n = n
        self.a = np.array([[rng.randint(-amax, amax) for _ in range(dim)] for _ in range(n)], dtype=np.int64)
        self.u = [[rng.randint(-umax, umax) for _ in range(3)] for _ in range(n)]
        self.br = {}
        for i, j, k in itertools.permutations(range(n), 3):
            M = sp.Matrix([self.u[i], self.u[j], self.u[k]]).T
            self.br[(i, j, k)] = int(M.det())


def eval_history(h, f, Ovecs, td):
    """<h.P with colour network and kernel, Phi>."""
    P = h.P
    if not P:
        return Fr(0)
    assert len(P) == 1, 'final skeleton polynomial has %d monomials' % len(P)
    (odd, even), coef = next(iter(P.items()))
    assert not even and all(l[0] == 1 for l in odd), 'final letters are not all c'
    finals = [l[1] for l in odd]                         # slot order = letter order
    # colour: contract the network, final legs open (in letter order), then with the a_i
    letters = {}
    def L(x):
        if x not in letters:
            letters[x] = len(letters)
        return letters[x]
    ops_, subs = [], []
    for legs in h.tensors:
        ops_.append(f)
        subs.append([L(x) for x in legs])
    for s, vvec in Ovecs.items():
        ops_.append(vvec)
        subs.append([L(s)])
    args = []
    for o, s in zip(ops_, subs):
        args += [o, s]
    T = np.einsum(*args, [L(s) for s in finals], optimize='greedy')
    argsf = [x.astype(float) if isinstance(x, np.ndarray) else x for x in args]
    Tf = np.einsum(*argsf, [L(s) for s in finals], optimize='greedy')
    for _ in finals:
        T = np.tensordot(T, td.a.T, axes=([0], [0]))     # legs cycle to the end, in order
        Tf = np.tensordot(Tf, td.a.T.astype(float), axes=([0], [0]))
    assert np.max(np.abs(Tf)) < 2.0 ** 52 and np.allclose(T, Tf), 'integer overflow in the colour network'
    # kernel labels of the finals
    labs = [('o',) + h.producer[s] for s in finals]
    labs = [(x[0], x[1], op(h.apps[x[1]][1]).outs.index(x[2])) if h.apps[x[1]][0] == 'G' else x for x in labs]
    K = history_kernel(h)
    tot = Fr(0)
    assert len(finals) == td.n
    for sig in itertools.permutations(range(td.n)):
        col = int(T[sig])
        if col == 0:
            continue
        w = dict(zip(labs, sig))
        kv = Fr(0)
        for key, c in K.items():
            p = 1
            for trip in key:
                p *= td.br.get((w[trip[0]], w[trip[1]], w[trip[2]]), 0)
                if p == 0:
                    break
            if p:
                kv += c * p
        if kv:
            tot += perm_sign(sig) * col * kv
    return tot * coef * h.pref


def evaluate(seq, O_types, f, Ovecs_list, tds):
    """sum over histories of the composition seq (applied right to left: seq[-1] first) on
    O = prod_k letter(O_types[k], slot k); returns [value for each test data]."""
    slots = list(range(1, len(O_types) + 1))
    P = const(1)
    for t, s in zip(O_types, slots):
        P = mul(P, gen(sk_letter(t, s)))
    hs = [History(P, slots, 100)]
    for name in reversed(seq):
        hs = step_q0(hs) if name == 'Q0' else step_graph(hs, name)
    out = []
    for td, Ovecs in zip(tds, Ovecs_list):
        tot = Fr(0)
        for h in hs:
            tot += eval_history(h, f, dict(zip(slots, Ovecs)), td)
        out.append(tot)
    return out, len(hs)


def fdense(name):
    n, fd = structure_constants(name)
    f = np.zeros((n, n, n), dtype=np.int64)
    for (A, B, C), v in fd.items():
        assert v.denominator == 1
        f[A, B, C] = int(v)
    return n, f, fd


def open_tensor(h, f, Ovecs):
    """Colour network of h contracted with the O vectors; final legs open, in letter order."""
    (odd, even), coef = next(iter(h.P.items()))
    finals = [l[1] for l in odd]
    letters = {}
    def L(x):
        if x not in letters:
            letters[x] = len(letters)
        return letters[x]
    args = []
    for legs in h.tensors:
        args += [f, [L(x) for x in legs]]
    for s, vvec in Ovecs.items():
        args += [vvec, [L(s)]]
    T = np.einsum(*args, [L(s) for s in finals], optimize='greedy')
    Tf = np.einsum(*[x.astype(float) if isinstance(x, np.ndarray) else x for x in args],
                   [L(s) for s in finals], optimize='greedy')
    assert np.max(np.abs(Tf)) < 2.0 ** 52 and np.allclose(T, Tf), 'integer overflow in the colour network'
    return T, finals, coef


def eval_with_tensor(h, T, finals, coef, td):
    assert len(finals) == td.n
    Tt, Tf = T, T.astype(float)
    for _ in finals:
        Tt = np.tensordot(Tt, td.a.T, axes=([0], [0]))
        Tf = np.tensordot(Tf, td.a.T.astype(float), axes=([0], [0]))
    assert np.max(np.abs(Tf)) < 2.0 ** 52 and np.allclose(Tt, Tf), 'integer overflow'
    labs = [('o',) + h.producer[s] for s in finals]
    labs = [(x[0], x[1], op(h.apps[x[1]][1]).outs.index(x[2])) if h.apps[x[1]][0] == 'G' else x for x in labs]
    K = history_kernel(h)
    tot = Fr(0)
    for sig in itertools.permutations(range(td.n)):
        col = int(Tt[sig])
        if col == 0:
            continue
        w = dict(zip(labs, sig))
        kv = Fr(0)
        for key, c in K.items():
            p = 1
            for trip in key:
                p *= td.br.get((w[trip[0]], w[trip[1]], w[trip[2]]), 0)
                if p == 0:
                    break
            if p:
                kv += c * p
        if kv:
            tot += perm_sign(sig) * col * kv
    return tot * coef * h.pref


def histories(seq, O_types):
    slots = list(range(1, len(O_types) + 1))
    P = const(1)
    for t, s in zip(O_types, slots):
        P = mul(P, gen(sk_letter(t, s)))
    hs = [History(P, slots, 100)]
    for name in reversed(seq):
        hs = step_q0(hs) if name == 'Q0' else step_graph(hs, name)
    hs = [h for h in hs if h.P]
    for h in hs:
        assert len(h.P) == 1
        (odd, even), c = next(iter(h.P.items()))
        assert not even and all(l[0] == 1 for l in odd), 'final letters are not all c'
    return hs, slots


def evaluate2(seq, O_types, f, Ovecs_list, tds):
    """values[iO][itd] = <seq O_{iO}, Phi_{itd}>."""
    hs, slots = histories(seq, O_types)
    vals = [[Fr(0)] * len(tds) for _ in Ovecs_list]
    for h in hs:
        for iO, Ov in enumerate(Ovecs_list):
            T, finals, coef = open_tensor(h, f, dict(zip(slots, Ov)))
            if not T.any():
                continue
            for it, td in enumerate(tds):
                vals[iO][it] += eval_with_tensor(h, T, finals, coef, td)
    return vals, len(hs)


# ------------------------------------------------------------------------------ vectorized evaluation
class TestDataV(TestData):
    def __init__(self, dim, seed, amax=2, umax=3, n=6):
        TestData.__init__(self, dim, seed, amax, umax, n)
        self.S = np.array(list(itertools.permutations(range(n))), dtype=np.int64)
        self.sgn = np.array([perm_sign(s) for s in self.S], dtype=np.int64)
        self.BR = np.zeros((n, n, n), dtype=np.int64)
        for (i, j, k), v in self.br.items():
            self.BR[i, j, k] = v


def _lcm(a, b):
    return a * b // math.gcd(a, b)


def eval_fast(h, T, finals, coef, td):
    """Same value as eval_with_tensor; sum over the n! assignments vectorized."""
    n = td.n
    assert len(finals) == n
    Tt, Tf = T, T.astype(float)
    for _ in finals:
        Tt = np.tensordot(Tt, td.a.T, axes=([0], [0]))
        Tf = np.tensordot(Tf, td.a.T.astype(float), axes=([0], [0]))
    assert np.max(np.abs(Tf)) < 2.0 ** 52 and np.allclose(Tt, Tf), 'integer overflow'
    col = Tt[tuple(td.S.T)]
    if not col.any():
        return Fr(0)
    labs = [('o',) + h.producer[s] for s in finals]
    labs = [(x[0], x[1], op(h.apps[x[1]][1]).outs.index(x[2])) if h.apps[x[1]][0] == 'G' else x for x in labs]
    pos = {l: i for i, l in enumerate(labs)}
    K = history_kernel(h)
    if not K:
        return Fr(0)
    den = 1
    for c in K.values():
        den = _lcm(den, c.denominator)
    kv = np.zeros(len(td.S), dtype=object)
    for key, c in K.items():
        num = int(c * den)
        prod = np.ones(len(td.S), dtype=np.int64)
        for trip in key:
            prod = prod * td.BR[td.S[:, pos[trip[0]]], td.S[:, pos[trip[1]]], td.S[:, pos[trip[2]]]]
        kv = kv + prod.astype(object) * num
    tot = int(np.sum(td.sgn.astype(object) * col.astype(object) * kv))
    return Fr(tot, den) * coef * h.pref


def evaluate3(seq, O_types, f, Ov, tds):
    """values[itd] = <seq O, Phi_itd> for one colour assignment Ov of O."""
    hs, slots = histories(seq, O_types)
    vals = [Fr(0)] * len(tds)
    for h in hs:
        T, finals, coef = open_tensor(h, f, dict(zip(slots, Ov)))
        if not T.any():
            continue
        for it, td in enumerate(tds):
            vals[it] += eval_fast(h, T, finals, coef, td)
    return vals, len(hs)


# ------------------------------------------------------------------------------ input derivatives on O
# O letters with jets: a letter d^k C bound to delta_j of an operator inserts prod_a (mom_a)^{k_a},
# mom = sum_v (dX_v/dz_j) lambda_v (graph) or lambda_n1 + lambda_n2 (hat Q_0).  Kernel keys are tuples
# of items ('b', l1, l2, l3) (bracket) and ('c', l, a) (component a of lambda_l).
_KCACHE2 = {}


def _slot_mom(h, s, inst):
    """Vec of the momentum carried away from slot s by its consumer (final: own label)."""
    if s not in h.consumer:
        return {('o',) + h.producer[s]: sp.Integer(1)}
    k2, j = h.consumer[s]
    kind2, name2, bind2, outslots2, _ = h.apps[k2]
    if kind2 == 'Q0':
        return vec_add(_slot_mom(h, outslots2[0], inst), _slot_mom(h, outslots2[1], inst))
    G2 = op(name2)
    cc = inst[k2][2]
    V = {}
    for i, v in enumerate(G2.outs):
        V = vec_add(V, {lab: cc[v][j] * c for lab, c in _slot_mom(h, outslots2[i], inst).items()})
    return V


def _relabel(h, V):
    out = {}
    for lab, c in V.items():
        k, x = lab[1], lab[2]
        if h.apps[k][0] == 'G':
            x = op(h.apps[k][1]).outs.index(x)
        out[('o', k, x)] = out.get(('o', k, x), 0) + c
    return out


def history_kernel2(h, Ojets):
    key = (structure_key(h), tuple(sorted((s, h.consumer.get(s), k) for s, k in Ojets.items())))
    if key in _KCACHE2:
        return _KCACHE2[key]
    inst = {}
    for (kind, name, bind, outslots, k) in h.apps:
        if kind == 'G':
            inst[k] = op(name).instance(k)
    BP = {(): sp.Integer(1)}
    doms = []
    for (kind, name, bind, outslots, k) in h.apps:
        if kind != 'G':
            continue
        G = op(name)
        facs, dom, cc = inst[k]
        doms += dom
        moms = {v: _relabel(h, _slot_mom(h, outslots[i], inst)) for i, v in enumerate(G.outs)}
        for pre, Vs in facs:
            Ws = []
            for V in Vs:
                W = {}
                for v, c in V.items():
                    W = vec_add(W, {lab: c * cm for lab, cm in moms[v].items()})
                Ws.append(W)
            D = det_expand(*Ws)
            BP = bp_mul(BP, {(('b',) + kk,): pre * vv for kk, vv in D.items()})
        BP = {kk: vv * G.sign for kk, vv in BP.items()}
    for s, jet in Ojets.items():
        V = _relabel(h, _slot_mom(h, s, inst))
        for a, ka in enumerate(jet):
            for _ in range(ka):
                BP = bp_mul(BP, {(('c', lab, a),): c for lab, c in V.items()})
    res = {}
    for kk, vv in BP.items():
        val = integrate_dom(vv, doms)
        if val != 0:
            res[kk] = Fr(int(val.p), int(val.q))
    _KCACHE2[key] = res
    return res


def eval_fast2(h, T, finals, coef, td, Ojets):
    n = td.n
    Tt, Tf = T, T.astype(float)
    for _ in finals:
        Tt = np.tensordot(Tt, td.a.T, axes=([0], [0]))
        Tf = np.tensordot(Tf, td.a.T.astype(float), axes=([0], [0]))
    assert np.max(np.abs(Tf)) < 2.0 ** 52 and np.allclose(Tt, Tf), 'integer overflow'
    col = Tt[tuple(td.S.T)]
    if not col.any():
        return Fr(0)
    labs = [('o',) + h.producer[s] for s in finals]
    labs = [(x[0], x[1], op(h.apps[x[1]][1]).outs.index(x[2])) if h.apps[x[1]][0] == 'G' else x for x in labs]
    pos = {l: i for i, l in enumerate(labs)}
    K = history_kernel2(h, Ojets)
    if not K:
        return Fr(0)
    U = np.array(td.u, dtype=np.int64)
    den = 1
    for c in K.values():
        den = _lcm(den, c.denominator)
    kv = np.zeros(len(td.S), dtype=object)
    for key, c in K.items():
        num = int(c * den)
        prod = np.ones(len(td.S), dtype=np.int64)
        for it in key:
            if it[0] == 'b':
                prod = prod * td.BR[td.S[:, pos[it[1]]], td.S[:, pos[it[2]]], td.S[:, pos[it[3]]]]
            else:
                prod = prod * U[td.S[:, pos[it[1]]], it[2]]
        kv = kv + prod.astype(object) * num
    tot = int(np.sum(td.sgn.astype(object) * col.astype(object) * kv))
    return Fr(tot, den) * coef * h.pref


def evaluate4(seq, O_types, Ojets_list, f, Ov, tds):
    """As evaluate3 with jets on the letters of O: Ojets_list[k] = jet of letter k (a 3-tuple)."""
    hs, slots = histories(seq, O_types)
    Ojets = {s: tuple(j) for s, j in zip(slots, Ojets_list) if any(j)}
    vals = [Fr(0)] * len(tds)
    for h in hs:
        T, finals, coef = open_tensor(h, f, dict(zip(slots, Ov)))
        if not T.any():
            continue
        for it, td in enumerate(tds):
            vals[it] += eval_fast2(h, T, finals, coef, td, Ojets)
    return vals, len(hs)


# ------------------------------------------------------------------------------ fast path: one graph application
_LOCAL = {}


def local_kernel(name, j=None, a=None):
    """int_D I-density(r; lambda_v) [ (sum_v dX_v/dz_j lambda_v)_a ]  in local labels v (output vertices):
    {tuple of items ('b', v1, v2, v3) / ('c', v, a): Fraction}."""
    key = (name, j, a)
    if key in _LOCAL:
        return _LOCAL[key]
    G = op(name)
    BP = {(): sp.Integer(1)}
    for pre, Vs in G.factors:
        D = det_expand(*Vs)
        BP = bp_mul(BP, {(('b',) + kk,): pre * vv for kk, vv in D.items()})
    BP = {kk: vv * G.sign for kk, vv in BP.items()}
    if j is not None:
        BP = bp_mul(BP, {(('c', v, a),): G.c[v][j] for v in G.outs if G.c[v][j] != 0})
    res = {}
    for kk, vv in BP.items():
        val = integrate_dom(vv, G.dom)
        if val != 0:
            res[kk] = Fr(int(val.p), int(val.q))
    _LOCAL[key] = res
    return res


def substitute(K, moms):
    """Replace local labels by Vecs of final labels (numeric coefficients) in a kernel dict."""
    out = {}
    for key, c in K.items():
        acc = {(): Fr(1) * c}
        for it in key:
            if it[0] == 'b':
                D = det_expand(moms[it[1]], moms[it[2]], moms[it[3]])
                part = {(('b',) + kk,): vv for kk, vv in D.items()}
            else:
                part = {(('c', lab, it[2]),): cc for lab, cc in moms[it[1]].items()}
            new = {}
            for k1, v1 in acc.items():
                for k2, v2 in part.items():
                    kk = tuple(sorted(k1 + k2))
                    new[kk] = new.get(kk, 0) + v1 * Fr(int(sp.Rational(v2).p), int(sp.Rational(v2).q))
            acc = new
        for kk, vv in acc.items():
            out[kk] = out.get(kk, 0) + vv
    return {k: v for k, v in out.items() if v != 0}


def kernel_single_graph(h, Ojets):
    """Kernel of a history with one graph application and hat Q_0's (final outputs or Q_0 outputs)."""
    gk = [k for k, ap in enumerate(h.apps) if ap[0] == 'G']
    assert len(gk) == 1
    k = gk[0]
    kind, name, bind, outslots, _ = h.apps[k]
    G = op(name)
    inst = {}            # no ratio symbols needed: moms of final labels have numeric coefficients

    def mom_out(s):
        if s not in h.consumer:
            return _relabel(h, {('o',) + h.producer[s]: sp.Integer(1)})
        k2, j = h.consumer[s]
        assert h.apps[k2][0] == 'Q0'
        o2 = h.apps[k2][3]
        return vec_add(mom_out(o2[0]), mom_out(o2[1]))
    moms = {v: mom_out(s) for v, s in zip(G.outs, outslots)}
    # jet factors: list over jet letters of lists of (kind, data): ('final', Vec, a) or ('moment', j, a)
    jet_terms = []
    for s, jet in Ojets.items():
        for a, ka in enumerate(jet):
            for _ in range(ka):
                terms = []
                def flow(slot):
                    if slot not in h.consumer:
                        terms.append(('final', _relabel(h, {('o',) + h.producer[slot]: sp.Integer(1)}), a))
                        return
                    k2, j2 = h.consumer[slot]
                    if h.apps[k2][0] == 'Q0':
                        for o in h.apps[k2][3]:
                            flow(o)
                    else:
                        assert k2 == k
                        terms.append(('moment', j2, a))
                flow(s)
                jet_terms.append(terms)
    assert len(jet_terms) <= 1, 'fast path: at most one input derivative'
    if not jet_terms:
        return substitute(local_kernel(name), moms)
    res = {}
    for t in jet_terms[0]:
        if t[0] == 'moment':
            part = substitute(local_kernel(name, t[1], t[2]), moms)
        else:
            base = substitute(local_kernel(name), moms)
            part = {}
            for kk, vv in base.items():
                for lab, cc in t[1].items():
                    key = tuple(sorted(kk + (('c', lab, t[2]),)))
                    part[key] = part.get(key, 0) + vv * Fr(int(cc))
        for kk, vv in part.items():
            res[kk] = res.get(kk, 0) + vv
    return {k_: v for k_, v in res.items() if v != 0}


def evaluate5(seq, O_types, Ojets_list, f, Ov, tds):
    """evaluate4 with the fast kernel path when the composition has one graph application."""
    hs, slots = histories(seq, O_types)
    Ojets = {s: tuple(j) for s, j in zip(slots, Ojets_list) if any(j)}
    single = sum(1 for x in seq if x != 'Q0') == 1
    vals = [Fr(0)] * len(tds)
    for h in hs:
        T, finals, coef = open_tensor(h, f, dict(zip(slots, Ov)))
        if not T.any():
            continue
        if single:
            K = kernel_single_graph(h, Ojets)
            _KCACHE2[(structure_key(h), tuple(sorted((s, h.consumer.get(s), k) for s, k in Ojets.items())))] = K
        for it, td in enumerate(tds):
            vals[it] += eval_fast2(h, T, finals, coef, td, Ojets)
    return vals, len(hs)
