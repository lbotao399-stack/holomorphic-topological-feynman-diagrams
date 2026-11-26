"""V3: BV Wick signs (verifier's own code).

Rule under test (explorer b5/c1): sigma = sgn(source -> target) [* (-1)^{#spokes O->v} for Q_n graphs],
  source = vertex fields, local order per vertex, vertices in vertex order;
  target = [(e,0),(e,1) per internal edge in edge order] + [externals/outputs in order] + [spoke fields in delta order].

Part A  obstruction thetas in engine labelling, local order L' (anomaly-structure colour_weight): sigma; and
        sigma * C_{L'} == W^Gamma of anomaly-structure claim c3 (index strings), tensors for so(3), so(4).
Part B  relabelling invariance of sigma*K*C (own relabelling code), mod sum p = 0; control with a wrong rule.
Part C  Q_n graphs built from the run-1 / S2 source strings (S2 eq-res-Q1, eq-res-TA, eq-res-TB; six-d-Q3 C1, C6;
        msym-odd theta_rootN + arm tensor of opsk.py docstring): sigma vs the nilpotency-fixed signs.
Part D  rule-free cross-check of the relative sign: Q_2 graph data with the nilpotency-fixed coefficients (msym-odd
        rootN/arm data (1,1); S2 data kappa_A = kappa_B) acting on the vertex operator O_N = (1/6) f C C C give the
        colour-kernel structure of the obstruction graphs; compare with K^eng W^Gamma (the convention behind [O_2]).
"""
import sys
import random
import itertools
import numpy as np
import sympy as sp
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/normalization/verify')
import htk
from htk import Graph
import functools
print = functools.partial(print, flush=True)
from vext import perm_sign

rng = random.Random(4242)


# ------------------------------------------------------------------ Lie algebras with kappa = delta
def so_n(N):
    basis = []
    for i in range(N):
        for j in range(i + 1, N):
            M = np.zeros((N, N), dtype=np.int64)
            M[i, j], M[j, i] = 1, -1
            basis.append(M)
    d = len(basis)
    f = np.zeros((d, d, d), dtype=np.int64)
    for a in range(d):
        for b in range(d):
            Cm = basis[a] @ basis[b] - basis[b] @ basis[a]
            for c in range(d):
                v = -np.trace(Cm @ basis[c])
                assert v % 2 == 0
                f[a, b, c] = v // 2
    # check total antisymmetry and kappa = delta
    for a in range(d):
        for b in range(d):
            assert -np.trace(basis[a] @ basis[b]) == 2 * (a == b)
    assert (f == -np.transpose(f, (1, 0, 2))).all() and (f == -np.transpose(f, (0, 2, 1))).all()
    return f


LIE = {'so3': so_n(3), 'so4': so_n(4)}
LET = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'


def sigma_rule(vertices, local, pairs, outs, spokes=(), spokes_OtoV=0):
    """vertices: list; local[v]: list of half-edge labels; pairs: list of (h0, h1) in edge order;
    outs: list of external/output half-edge labels in order; spokes: spoke-field labels in delta order."""
    source = [h for v in vertices for h in local[v]]
    target = [h for pr in pairs for h in pr] + list(outs) + list(spokes)
    assert sorted(map(repr, source)) == sorted(map(repr, target))
    pos = {repr(h): i for i, h in enumerate(target)}
    return perm_sign([pos[repr(h)] for h in source]) * (-1) ** spokes_OtoV


def contract(f, local_lists, pairs, free):
    """einsum of f_{local order} per vertex, delta on pairs, free indices in the order `free`."""
    hs = sorted({repr(h) for L in local_lists for h in L})
    let = {h: LET[i] for i, h in enumerate(hs)}
    terms = [''.join(let[repr(h)] for h in L) for L in local_lists]
    ops = [f] * len(local_lists)
    for (a, b) in pairs:
        terms.append(let[repr(a)] + let[repr(b)])
        ops.append(np.eye(f.shape[0], dtype=np.int64))
    return np.einsum(','.join(terms) + '->' + ''.join(let[repr(h)] for h in free), *ops, optimize='greedy')


def index_string_tensor(f, string, free, sign=1):
    """string like 'xyz|a1 y2 u|...' : per-vertex index names; names ending in "'" are contracted with
    the unprimed name (kappa^-1 = delta).  free: names of free indices in order."""
    verts = [v.split() for v in string.split('|')]
    names = sorted({n.rstrip("'") for v in verts for n in v})
    let = {n: LET[i] for i, n in enumerate(names)}
    let2 = {n: LET[len(names) + i] for i, n in enumerate(names)}
    terms, ops = [], []
    primed = set()
    for v in verts:
        s = ''
        for n in v:
            if n.endswith("'"):
                s += let2[n[:-1]]
                primed.add(n[:-1])
            else:
                s += let[n]
        terms.append(s)
        ops.append(f)
    for n in sorted(primed):
        terms.append(let[n] + let2[n])
        ops.append(np.eye(f.shape[0], dtype=np.int64))
    return sign * np.einsum(','.join(terms) + '->' + ''.join(let[n] for n in free), *ops, optimize='greedy')


print("V3  BV Wick signs")
# ------------------------------------------------------------------ Part A
print("Part A: obstruction thetas, engine labelling, local order L' (externals first, then incident half-edges by (edge, end))")
W_strings = {   # anomaly-structure claim c3 (digest), vertex order 0..5 of htk.theta, overall sign -1
    (1, 3, 3): "x y z|a1 y' u|a2 u' w|a3 z' s|a4 s' r|x' w' r'",
    (2, 2, 3): "x y z|a1 x' u|a2 y' w|a3 z' s|a4 s' r|u' w' r'",
}


def engine_data(pqr):
    G, inter = htk.theta(*pqr)
    ext = [v for P in inter for v in P]
    local = {}
    for v in range(G.nv):
        hs = [('x', v)] if v in ext else []
        hs += sorted([('e', e, s) for e, ab in enumerate(G.edges) for s in (0, 1) if ab[s] == v], key=lambda h: (h[1], h[2]))
        local[v] = hs
    pairs = [(('e', e, 0), ('e', e, 1)) for e in range(G.ne)]
    outs = [('x', v) for v in ext]
    return G, ext, local, pairs, outs


SIG_ENG = {}
for pqr in [(1, 3, 3), (2, 2, 3)]:
    G, ext, local, pairs, outs = engine_data(pqr)
    sig = sigma_rule(list(range(G.nv)), local, pairs, outs)
    SIG_ENG[pqr] = sig
    res = []
    for nm, f in LIE.items():
        C = contract(f, [local[v] for v in range(G.nv)], pairs, outs)
        W = index_string_tensor(f, W_strings[pqr], ['a1', 'a2', 'a3', 'a4'], sign=-1)
        res.append((nm, bool((sig * C == W).all()), int(np.abs(W).sum()) > 0))
    print("  theta%s: edges %s ext %s ; sigma = %+d ; sigma*C_{L'} == W^Gamma (c3 strings), nonzero: %s"
          % (pqr, G.edges, ext, sig, res))
print("  sigma_133 * sigma_223 = %+d" % (SIG_ENG[(1, 3, 3)] * SIG_ENG[(2, 2, 3)]))

# |Aut| with unlabelled external legs (leaf vertices), networkx
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher
for pqr in [(1, 3, 3), (2, 2, 3)]:
    G, ext, *_ = engine_data(pqr)
    M = nx.MultiGraph()
    M.add_nodes_from(range(G.nv))
    for (a, b) in G.edges:
        M.add_edge(a, b)
    for k, v in enumerate(ext):
        M.add_edge(v, 'leaf%d' % k)
    Gs = nx.Graph(M)
    naut = sum(1 for _ in GraphMatcher(Gs, Gs).isomorphisms_iter())
    print("  |Aut theta%s| (external legs unlabelled) = %d" % (pqr, naut))

# ------------------------------------------------------------------ Part B
print("Part B: relabelling invariance of sigma*K*C, (H,T) = (2,1), so(4), mod sum p = 0 (own relabelling code)")
H = 2
Psym = [[sp.Symbol('p%d_%d' % (i, c)) for c in range(H)] for i in range(4)]


def modcons(e):
    sub = {Psym[3][c]: -(Psym[0][c] + Psym[1][c] + Psym[2][c]) for c in range(H)}
    return sp.expand(sp.sympify(e).subs(sub, simultaneous=True))


def physical_quantity(G0, ext0, local0, perm, order, flips, lperm, f, wrong=False):
    """Relabel: vertex u -> perm[u]; new edge k = old edge order[k]; old edge e reversed if flips[e];
    local order at u permuted by lperm[u].  Returns sigma, K(sympy, momenta of physical externals), C tensor."""
    newid = {e: k for k, e in enumerate(order)}
    edges = []
    for k, e in enumerate(order):
        a, b = G0.edges[e]
        edges.append((perm[b], perm[a]) if flips[e] else (perm[a], perm[b]))
    G = Graph(G0.nv, edges, 'relab')

    def nh(h):
        if h[0] == 'x':
            return ('x', perm[h[1]])
        e, s = h[1], h[2]
        return ('e', newid[e], 1 - s if flips[e] else s)
    inv = {perm[u]: u for u in range(G0.nv)}
    local = {}
    for vn in range(G0.nv):
        u = inv[vn]
        L = [local0[u][i] for i in lperm[u]]
        local[vn] = [nh(h) for h in L]
    wpairs = [(('e', k, 0), ('e', k, 1)) if edges[k][0] < edges[k][1] else (('e', k, 1), ('e', k, 0)) for k in range(G.ne)]
    pairs = [(('e', k, 0), ('e', k, 1)) for k in range(G.ne)]
    outs = [('x', perm[u]) for u in ext0]
    sig = sigma_rule(list(range(G.nv)), local, pairs, outs)
    sig_wrong = sigma_rule(list(range(G.nv)), local, wpairs, outs)
    lam = {}
    for vn in range(G.nv - 1):
        u = inv[vn]
        lam[vn] = list(Psym[ext0.index(u)]) if u in ext0 else [0] * H
    K = htk.kernel_exact(G, lam, T=1, check=1)
    C = contract(f, [local[v] for v in range(G.nv)], [(('e', k, 0), ('e', k, 1)) for k in range(G.ne)], outs)
    return sig, sp.sympify(K), C, sig_wrong


def same(q1, q2):
    s1, K1, C1 = q1
    s2, K2, C2 = q2
    nz = np.argwhere((C1 != 0) | (C2 != 0))
    for idx in nz[:40]:
        idx = tuple(idx)
        if modcons(s1 * K1 * int(C1[idx]) - s2 * K2 * int(C2[idx])) != 0:
            return False
    return True


f4 = LIE['so4']
for pqr in [(1, 3, 3), (2, 2, 3)]:
    G0, ext0, local0, pairs0, outs0 = engine_data(pqr)
    ident = (list(range(G0.nv)), list(range(G0.ne)), [False] * G0.ne, {u: [0, 1, 2] for u in range(G0.nv)})
    ref = physical_quantity(G0, ext0, local0, *ident, f4)
    inv_ok, ctrl, ctrl2 = [], [], []
    for trial in range(5):
        perm = list(range(G0.nv)); rng.shuffle(perm)
        order = list(range(G0.ne)); rng.shuffle(order)
        flips = [rng.random() < 0.5 for _ in range(G0.ne)]
        lperm = {u: rng.sample([0, 1, 2], 3) for u in range(G0.nv)}
        q = physical_quantity(G0, ext0, local0, perm, order, flips, lperm, f4)
        inv_ok.append(same(ref[:3], q[:3]))
        ctrl.append(same((ref[3],) + ref[1:3], (q[3],) + q[1:3]))
        ctrl2.append(same((1,) + ref[1:3], (1,) + q[1:3]))
    print("  theta%s: K_ref = %s ; invariant under 5 random relabellings: %s ; controls: orientation-blind pairs (smaller label first): %s ; sigma dropped: %s"
          % (pqr, sp.factor(modcons(ref[1])), inv_ok, ctrl, ctrl2))

# ------------------------------------------------------------------ Part C
print("Part C: Q_n graphs (local orders parsed from the source T-strings)")


def qgraph(edges_str, outputs, T_string, vorder):
    """edges_str: list of (a,b) with 'O'; outputs: vertices with an output; T_string: per-vertex index names
    (R<a><b> = edge a->b, A<j> = j-th spoke in edge order, D<v> = output at v), vertices in vorder."""
    spokes = [k for k, (a, b) in enumerate(edges_str) if 'O' in (a, b)]
    name = {}
    for j, k in enumerate(spokes):
        a, b = edges_str[k]
        name['A%d' % (j + 1)] = ('e', k, 1) if a == 'O' else ('e', k, 0)
    verts = T_string.split('|')
    local = {}
    for v, s in zip(vorder, verts):
        L = []
        for tok in s.split():
            if tok.startswith('A'):
                h = name[tok]
                assert edges_str[h[1]][h[2]] == v, (tok, v)
            elif tok.startswith('D'):
                assert tok[1:] == str(v)
                h = ('x', v)
            else:   # R<a><b> or R-alias given as (a,b) dict below
                a, b = tok[1], tok[2]
                cand = [k for k, (x, y) in enumerate(edges_str) if str(x) == a and str(y) == b]
                assert len(cand) == 1, tok
                k = cand[0]
                h = ('e', k, 0) if str(v) == a else ('e', k, 1)
                assert str(edges_str[k][h[2]]) == str(v)
            L.append(h)
        local[v] = L
    pairs = [(('e', k, 0), ('e', k, 1)) for k, (a, b) in enumerate(edges_str) if 'O' not in (a, b)]
    outs = [('x', v) for v in vorder if v in outputs]
    spk = [name['A%d' % (j + 1)] for j in range(len(spokes))]
    nOv = sum(1 for k in spokes if edges_str[k][0] == 'O')
    return sigma_rule(vorder, local, pairs, outs, spk, nOv), local


def E(s):
    out = []
    for tok in s.strip('()').split(')('):
        out.append(tuple('O' if c == 'O' else int(c) for c in tok))
    return out


# S2 box (eq-res-Q1), A (eq-res-TA), B (eq-res-TB); R indices renamed to edge names
S2c = {
    'box': (E('(O1)(12)(23)(3O)'), [1, 2, 3], "A1 D1 R12|R12 D2 R23|R23 D3 A2", [1, 2, 3]),
    'A': (E('(O1)(12)(23)(O3)(34)(45)(5O)'), [1, 2, 4, 5], "A1 D1 R12|R12 D2 R23|R23 A2 R34|R34 D4 R45|R45 D5 A3", [1, 2, 3, 4, 5]),
}
sig_S2 = {}
for nm, (ed, outs, Ts, vo) in S2c.items():
    sig_S2[nm], _ = qgraph(ed, outs, Ts, vo)
# B with letters u,w,q,a,b -> 1..5: edges (O u)(u q)(O w)(w q)(q a)(a b)(b O); T_B = A1 Du Ru|A2 Dw Rw|Ru Rw Rq|Rq Da R'|R' Db A3
edB = [('O', 1), (1, 3), ('O', 2), (2, 3), (3, 4), (4, 5), (5, 'O')]
TB = "A1 D1 R13|A2 D2 R23|R13 R23 R34|R34 D4 R45|R45 D5 A3"
sig_S2['B'], _ = qgraph(edB, [1, 2, 4, 5], TB, [1, 2, 3, 4, 5])
print("  S2 (3,0):  sigma(box, A, B) = %s ;  S2: box sign derived (eq-box-amp) +, kappa_A = kappa_B = +1/2 (eq-Q2, eq-kAkB)"
      % [sig_S2[k] for k in ('box', 'A', 'B')])
Q3 = {
    'G1': ('(O1)(12)(23)(O3)(34)(45)(O5)(56)(67)(7O)', [1, 2, 4, 6, 7], 'A1 D1 R12|R12 D2 R23|R23 A2 R34|R34 D4 R45|R45 A3 R56|R56 D6 R67|R67 D7 A4', 2),
    'G2': ('(O1)(12)(23)(O3)(34)(45)(56)(6O)(57)(7O)', [1, 2, 4, 6, 7], 'A1 D1 R12|R12 D2 R23|R23 A2 R34|R34 D4 R45|R45 R56 R57|R56 D6 A3|R57 D7 A4', 2),
    'G3': ('(O1)(12)(23)(O3)(34)(45)(5O)(46)(67)(7O)', [1, 2, 5, 6, 7], 'A1 D1 R12|R12 D2 R23|R23 A2 R34|R34 R45 R46|R45 D5 A3|R46 D6 R67|R67 D7 A4', 1),
    'G4': ('(O1)(13)(O2)(23)(34)(45)(56)(6O)(57)(7O)', [1, 2, 4, 6, 7], 'A1 D1 R13|A2 D2 R23|R13 R23 R34|R34 D4 R45|R45 R56 R57|R56 D6 A3|R57 D7 A4', 8),
    'G5': ('(O1)(13)(O2)(23)(34)(45)(5O)(46)(67)(7O)', [1, 2, 5, 6, 7], 'A1 D1 R13|A2 D2 R23|R13 R23 R34|R34 R45 R46|R45 D5 A3|R46 D6 R67|R67 D7 A4', 2),
}
kappa_Q3 = [sp.Rational(-1, 2), sp.Rational(1, 2), -1, sp.Rational(1, 8), sp.Rational(-1, 2)]
pred = []
for nm, (es, outs, Ts, aut) in Q3.items():
    s, _ = qgraph(E(es), outs, Ts, list(range(1, 8)))
    pred.append(sp.Rational(s, aut))
print("  six-d-Q3 (n=3): sigma/|Aut| = %s ; nilpotency (C7) kappa = %s ; equal: %s" % (pred, kappa_Q3, pred == kappa_Q3))


# msym-odd theta_rootN + arm tensor (opsk.py docstring)
def rootN(p, q, r):
    nv = p + q + r - 1
    N, Np = nv - 1, nv - 2
    edges, paths, nxt = [], [], 0
    for n_ in (p, q, r):
        prev, pv = N, []
        for j in range(n_):
            cur = Np if j == n_ - 1 else nxt
            if j < n_ - 1:
                pv.append(nxt)
                nxt += 1
            edges.append((prev, cur))
            prev = cur
        paths.append(pv)
    return nv, edges, paths


def rootN_local(nv, edges, paths):
    """arm convention: interior vertex: (incoming half-edge, output, outgoing half-edge); N': heads of last edges A,B,C."""
    N, Np = nv - 1, nv - 2
    local = {}
    for v in range(nv - 1):
        ins = [('e', k, 1) for k, (a, b) in enumerate(edges) if b == v]
        outs_ = [('e', k, 0) for k, (a, b) in enumerate(edges) if a == v]
        if v == Np:
            local[v] = ins            # edge order = path order A, B, C
        else:
            assert len(ins) == 1 and len(outs_) == 1
            local[v] = [ins[0], ('x', v), outs_[0]]
    return local


sig_rootN = {}
for pqr in [(1, 3, 3), (2, 2, 3)]:
    nv, edges, paths = rootN(*pqr)
    N = nv - 1
    local = rootN_local(nv, edges, paths)
    spokes = [k for k, (a, b) in enumerate(edges) if a == N]
    pairs = [(('e', k, 0), ('e', k, 1)) for k, (a, b) in enumerate(edges) if a != N and b != N]
    outs = [('x', v) for P in paths for v in P]
    s = sigma_rule(list(range(nv - 1)), local, pairs, outs, [('e', k, 1) for k in spokes], len(spokes))
    sig_rootN[pqr] = s
print("  msym-odd rootN + arm tensor: sigma(133), sigma(223) = %s ; nilpotency null space (k_A,k_B) = (1,1) (run-1 verified)"
      % [sig_rootN[(1, 3, 3)], sig_rootN[(2, 2, 3)]])

# ------------------------------------------------------------------ Part D
print("Part D: rule-free relative sign: nilpotency-fixed Q_2 data acting on O_N = (1/6) f C C C vs K^eng W^Gamma")
# prod_j delta/delta C^{A_j} (A_1 A_2 A_3 left to right) of (1/6) f_abc C^a C^b C^c = f_{A3 A2 A1} = -f_{A1 A2 A3}
# check this Grassmann identity by brute force on so(3)
f3 = LIE['so3']


def check_delta_identity(f):
    d = f.shape[0]
    # represent cubic polynomial coefficients c[a,b,c] of C^a C^b C^c (antisymmetrized)
    for A1, A2, A3 in itertools.product(range(d), repeat=3):
        # left derivative delta_{A3}: sum over positions with sign (-1)^{pos}
        val = 0
        for (a, b, c) in itertools.product(range(d), repeat=3):
            coef = sp.Rational(int(f[a, b, c]), 6)
            word = [a, b, c]
            # apply delta_A3, then delta_A2, then delta_A1 (rightmost acts first)
            terms = [(coef, word)]
            for D in (A3, A2, A1):
                new = []
                for cf, w in terms:
                    for p, x in enumerate(w):
                        if x == D:
                            new.append((cf * (-1) ** p, w[:p] + w[p + 1:]))
                terms = new
            val += sum(cf for cf, w in terms)
        if val != -f[A1, A2, A3]:
            return False
    return True


print("  prod_j delta/delta C^{A_j} [(1/6) f C C C] = -f_{A1 A2 A3} (so(3), brute-force Grassmann): %s" % check_delta_identity(f3))


def arm_tensor_contracted(f, pqr):
    """X^col_{D...} = -f_{A_A A_B A_C} T^arm_{A_A A_B A_C; D...}, outputs in path order."""
    nv, edges, paths = rootN(*pqr)
    local = rootN_local(nv, edges, paths)
    N = nv - 1
    spokes = [k for k, (a, b) in enumerate(edges) if a == N]
    # T^arm as contraction: vertices 0..nv-2 with local orders; internal pairs; free = outputs, spoke fields
    pairs = [(('e', k, 0), ('e', k, 1)) for k, (a, b) in enumerate(edges) if a != N and b != N]
    outs = [('x', v) for P in paths for v in P]
    spk = [('e', k, 1) for k in spokes]
    T = contract(f, [local[v] for v in range(nv - 1)], pairs, outs + spk)
    return -np.einsum('ijklabc,abc->ijkl', T, f)


def S2_tensor_contracted(f, which):
    """S2 T_A / T_B with spokes contracted into O_N; outputs reordered to engine path order (odd letters)."""
    if which == 'A':
        ed = E('(O1)(12)(23)(O3)(34)(45)(5O)'); Ts = S2c['A'][2]; outsv = [1, 2, 4, 5]
        # engine identification: O->N=0, S2 3 -> N'=5, path O-1-2-3 -> B (eng 1,2), path O-5-4-3 -> C (eng 3,4)
        phys = {1: 1, 2: 2, 5: 3, 4: 4}
    else:
        ed = edB; Ts = TB; outsv = [1, 2, 4, 5]
        # B: u=1, w=2, q=3, a=4, b=5. theta(2,2,3): paths O-u-q (A, eng 1), O-w-q (B, eng 2), O-b-a-q (C, eng 3,4)
        phys = {1: 1, 2: 2, 5: 3, 4: 4}
    _, local = qgraph(ed, outsv, Ts, [1, 2, 3, 4, 5])
    spokes = [k for k, (a, b) in enumerate(ed) if 'O' in (a, b)]
    spk = [('e', k, 1) if ed[k][0] == 'O' else ('e', k, 0) for k in spokes]
    pairs = [(('e', k, 0), ('e', k, 1)) for k, (a, b) in enumerate(ed) if 'O' not in (a, b)]
    outs = [('x', v) for v in outsv]
    T = contract(f, [local[v] for v in [1, 2, 3, 4, 5]], pairs, outs + spk)
    X = -np.einsum('ijklabc,abc->ijkl', T, f)
    # reorder odd output letters from S2 vertex order to engine ext order
    order_eng = sorted(outsv, key=lambda v: phys[v])          # S2 vertices in engine order
    perm = [outsv.index(v) for v in order_eng]
    sgn = perm_sign(perm)
    return sgn * np.transpose(X, perm), phys


@functools.lru_cache(maxsize=None)
def K_engine(pqr):
    G, inter = htk.theta(*pqr)
    ext = [v for P in inter for v in P]
    lam = {v: (list(Psym[ext.index(v)]) if v in ext else [0, 0]) for v in range(G.nv - 1)}
    return sp.sympify(htk.kernel_exact(G, lam, T=1))


@functools.lru_cache(maxsize=None)
def K_rootN(pqr):
    nv, edges, paths = rootN(*pqr)
    G = Graph(nv, edges, 'rootN')
    ext = [v for P in paths for v in P]
    lam = {v: (list(Psym[ext.index(v)]) if v in ext else [0, 0]) for v in range(nv - 1)}
    return sp.sympify(htk.kernel_exact(G, lam, T=1))


@functools.lru_cache(maxsize=None)
def K_S2lab(which):
    G = htk.graph_A() if which == 'A' else htk.graph_B()
    # htk vertex k = S2 vertex k+1 (A: 1..5; B: u,w,q,a,b = 1..5); phys map S2 -> engine ext position
    phys = {1: 1, 2: 2, 5: 3, 4: 4}
    lam = {}
    for k in range(5):
        s2v = k + 1
        lam[k] = list(Psym[phys[s2v] - 1]) if s2v in phys else [0, 0]
    return sp.sympify(htk.kernel_exact(G, lam, T=1))


def ratio(Kx, Xc, Ky, Yc):
    """common scalar r with Kx*Xc = r*Ky*Yc mod sum p (entrywise); None if inconsistent."""
    kx, ky = modcons(Kx), modcons(Ky)
    rk = sp.simplify(kx / ky)
    if rk.free_symbols:
        return None
    nz = np.argwhere(Yc != 0)
    rs = set()
    for idx in nz:
        idx = tuple(idx)
        rs.add(sp.Rational(int(Xc[idx]), int(Yc[idx])))
    if (Xc[Yc == 0] != 0).any() or len(rs) != 1:
        return None
    return rk * rs.pop()


for nm, f in LIE.items():
    rN, rS = {}, {}
    for pqr, which in [((1, 3, 3), 'A'), ((2, 2, 3), 'B')]:
        Y = index_string_tensor(f, W_strings[pqr], ['a1', 'a2', 'a3', 'a4'], sign=-1)
        Ke = K_engine(pqr)
        rN[pqr] = ratio(K_rootN(pqr), arm_tensor_contracted(f, pqr), Ke, Y)
        Xs, _ = S2_tensor_contracted(f, which)
        rS[pqr] = ratio(K_S2lab(which), Xs, Ke, Y)
    print("  %s: [K^rootN X^arm] / [K^eng W] = %s ; [K^S2lab X^S2] / [K^eng W] = %s"
          % (nm, [rN[(1, 3, 3)], rN[(2, 2, 3)]], [rS[(1, 3, 3)], rS[(2, 2, 3)]]))
    print("       equal ratios for the two thetas (=> relative sign +1 of K W): rootN data %s ; S2 data %s"
          % (rN[(1, 3, 3)] == rN[(2, 2, 3)], rS[(1, 3, 3)] == rS[(2, 2, 3)]))
print("  K^eng: 133 = %s ; 223 = %s" % (sp.factor(modcons(K_engine((1, 3, 3)))), sp.factor(modcons(K_engine((2, 2, 3))))))
