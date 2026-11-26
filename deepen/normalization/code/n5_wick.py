"""(b),(c) BV Wick signs: obstruction thetas of CS on R x C^2, and the Q_n graphs of S2 / six-d-Q3 / msym-odd.

Part A  obstruction thetas, engine labelling, explorer local order L': sigma (wick.py) vs sgn(L -> L').
Part B  labelling invariance of sigma * K * C for the obstruction thetas ((H,T) = (2,1), so(3) colour tensor,
        momenta at the physical external vertices, compared modulo sum p = 0): random vertex relabellings
        (root moved), edge reversals, edge-order permutations, random local leg orders.
Part C  Q_n graphs with their decorations (S2 leg pattern: in-edges in edge order, output, out-edges; spoke j -> A_j
        in edge order; outputs in vertex order): sigma for S2 box, A, B; six-d-Q3 G1..G5; msym-odd S2 data and
        theta_rootN data (arm tensor = S2 pattern in that labelling, checked below).
Part D  labelling invariance of sigma * K * T for S2 A, B at (H,T) = (2,1) and (3,0) (O fixed as root).
"""
import sys
import random
import itertools
import numpy as np
import sympy as sp
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/normalization/code')
import htk
from htk import Graph
from wick import wick_sign, pattern_local, explorer_LLp_sign, explorer_local, perm_sign

rng = random.Random(2026)
eps3 = np.zeros((3, 3, 3), dtype=np.int64)
for (a, b, c), s in [((0, 1, 2), 1), ((1, 2, 0), 1), ((2, 0, 1), 1), ((0, 2, 1), -1), ((2, 1, 0), -1), ((1, 0, 2), -1)]:
    eps3[a, b, c] = s
kinv = np.eye(3, dtype=np.int64)
LET = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'


def colour_tensor(verts, edges, local, free, O='O'):
    """einsum of f_{local order} at each vertex in verts and kinv on internal edges; free: list of half-edges
    (outputs, spoke fields) giving the order of the free indices."""
    hs = sorted({str(h) for v in verts for h in local[v]})
    let = {h: LET[i] for i, h in enumerate(hs)}
    terms, ops = [], []
    for v in verts:
        terms.append(''.join(let[str(h)] for h in local[v]))
        ops.append(eps3)
    for e, (a, b) in enumerate(edges):
        if a != O and b != O:
            terms.append(let[str(('e', e, 0))] + let[str(('e', e, 1))])
            ops.append(kinv)
    out = ''.join(let[str(h)] for h in free)
    return np.einsum(','.join(terms) + '->' + out, *ops)


# ------------------------------------------------------------------ Part A
print("Part A: obstruction thetas, engine labelling (N = 0, N' = root), local order L' of anomaly-structure")
P = [[sp.Symbol('p%d_%d' % (i, c)) for c in range(2)] for i in range(4)]
for pqr in [(1, 3, 3), (2, 2, 3)]:
    G, inter = htk.theta(*pqr)
    ext = [v for Pp in inter for v in Pp]
    loc = explorer_local(G.nv, G.edges, ext)
    sig = wick_sign(list(range(G.nv)), G.edges, loc, [('x', v) for v in ext])
    sLL, _ = explorer_LLp_sign(G.nv, G.edges, ext)
    print("  theta%s: sigma (Wick) = %+d ; sgn(L -> L') = %+d ; equal: %s" % (pqr, sig, sLL, sig == sLL))


# ------------------------------------------------------------------ Part B
def conserve(expr):
    sub = {P[3][c]: -(P[0][c] + P[1][c] + P[2][c]) for c in range(2)}
    return sp.expand(sp.sympify(expr).subs(sub, simultaneous=True))


def relabelled_quantity(G0, ext0, H, T, perm, order, flips, locperm, O=None, outputs_phys=None, spokes_phys=None):
    """returns (sigma, K, Ctensor) for the relabelled graph; physical quantities are referred to G0 labels."""
    nv = G0.nv
    edges = []
    newid = {}
    for k, e in enumerate(order):
        a, b = G0.edges[e]
        a2, b2 = perm[a], perm[b]
        edges.append((b2, a2) if flips[e] else (a2, b2))
        newid[e] = k
    G = Graph(nv, edges, "relab")
    inv = {perm[u]: u for u in range(nv)}

    # physical half-edge (e, end) -> new half-edge
    def newh(e, end):
        return ('e', newid[e], 1 - end if flips[e] else end)

    if O is None:
        verts_new = list(range(nv))
    else:
        verts_new = [v for v in range(nv) if v != perm[O]]
    # local orders: physical local order (explorer) permuted by locperm, translated
    loc = {}
    for vn in verts_new:
        u = inv[vn]
        hs = [('x', u)] if u in (ext0 if O is None else outputs_phys) else []
        hs += [('e', e, s) for e, ab in enumerate(G0.edges) for s in (0, 1) if ab[s] == u]
        hs = [hs[i] for i in locperm[u]]
        loc[vn] = [('x', vn) if h[0] == 'x' else newh(h[1], h[2]) for h in hs]
    if O is None:
        outs = [('x', perm[u]) for u in ext0]
        Oname = 'O'
        edges_w = edges
        sig = wick_sign(verts_new, edges_w, loc, outs)
        lam = {}
        for vn in range(nv - 1):
            u = inv[vn]
            lam[vn] = list(P[ext0.index(u)]) if u in ext0 else [0] * H
        K = htk.kernel_exact(G, lam, T=T)
        free = outs
        Cten = colour_tensor(verts_new, edges_w, loc, free)
    else:
        Onew = perm[O]
        assert Onew == nv - 1
        edges_w = [tuple('O' if x == Onew else x for x in ab) for ab in edges]
        outs = [('x', perm[u]) for u in outputs_phys]
        spokes_new = [newid[e] for e in spokes_phys]
        sig = wick_sign(verts_new, edges_w, loc, outs, spokes_new)
        lam = {}
        for vn in range(nv - 1):
            u = inv[vn]
            lam[vn] = list(P[outputs_phys.index(u)]) if u in outputs_phys else [0] * H
        K = htk.kernel_exact(G, lam, T=T)
        free = outs + [('e', k, 1) if edges_w[k][0] == 'O' else ('e', k, 0) for k in spokes_new]
        Cten = colour_tensor(verts_new, edges_w, loc, free)
    return sig, K, Cten


def same(q1, q2, cons=True):
    s1, K1, C1 = q1
    s2, K2, C2 = q2
    f = conserve if cons else sp.expand
    # compare s K C entrywise: K polynomial, C integer tensor
    nz = np.argwhere((C1 != 0) | (C2 != 0))
    if len(nz) == 0:
        return f(s1 * K1 - s2 * K2) == 0 or True
    for idx in nz:
        idx = tuple(idx)
        if f(s1 * K1 * int(C1[idx]) - s2 * K2 * int(C2[idx])) != 0:
            return False
    return True


print("Part B: labelling invariance of sigma*K*C, obstruction thetas, (H,T) = (2,1), so(3), modulo sum p = 0")
for pqr in [(1, 3, 3), (2, 2, 3)]:
    G0, inter = htk.theta(*pqr)
    ext0 = [v for Pp in inter for v in Pp]
    ident = (list(range(G0.nv)), list(range(G0.ne)), [False] * G0.ne,
             {u: list(range(3)) for u in range(G0.nv)})
    ref = relabelled_quantity(G0, ext0, 2, 1, *ident)
    res = []
    flipped_only = []
    for trial in range(6):
        perm = list(range(G0.nv)); rng.shuffle(perm)
        order = list(range(G0.ne)); rng.shuffle(order)
        flips = [rng.random() < 0.5 for _ in range(G0.ne)]
        locperm = {u: rng.sample(range(3), 3) for u in range(G0.nv)}
        q = relabelled_quantity(G0, ext0, 2, 1, perm, order, flips, locperm)
        res.append(same(ref, q))
        # control: same relabelling with the Wick sign replaced by +1 must fail for some trials
        flipped_only.append(same(ref, (1, q[1], q[2])) if q[0] == -1 else None)
    print("  theta%s: 6 random relabellings: invariant %s ; controls with sigma -> +1 where sigma = -1: %s"
          % (pqr, res, flipped_only))

# ------------------------------------------------------------------ Part C
print("Part C: Q_n graphs with their decorations")


def qn_sigma(name, vertices, edges, outputs_at):
    loc = pattern_local(vertices, edges, outputs_at)
    spokes = [e for e, (a, b) in enumerate(edges) if 'O' in (a, b)]
    outs = [('x', v) for v in vertices if v in outputs_at]
    s = wick_sign(vertices, edges, loc, outs, spokes)
    return s, loc, spokes


S2 = {
    'box': ([1, 2, 3], [('O', 1), (1, 2), (2, 3), (3, 'O')], [1, 2, 3]),
    'A': ([1, 2, 3, 4, 5], [('O', 1), (1, 2), (2, 3), ('O', 3), (3, 4), (4, 5), (5, 'O')], [1, 2, 4, 5]),
    'B': (['u', 'w', 'q', 'a', 'b'], [('O', 'u'), ('u', 'q'), ('O', 'w'), ('w', 'q'), ('q', 'a'), ('a', 'b'), ('b', 'O')],
          ['u', 'w', 'a', 'b']),
}


def E(s):
    out = []
    for tok in s.split(')('):
        tok = tok.strip('()')
        a, b = tok[0], tok[1]
        out.append(('O' if a == 'O' else int(a), 'O' if b == 'O' else int(b)))
    return out


Q3 = {
    'G1': (E('(O1)(12)(23)(O3)(34)(45)(O5)(56)(67)(7O)'), [1, 2, 4, 6, 7]),
    'G2': (E('(O1)(12)(23)(O3)(34)(45)(56)(6O)(57)(7O)'), [1, 2, 4, 6, 7]),
    'G3': (E('(O1)(12)(23)(O3)(34)(45)(5O)(46)(67)(7O)'), [1, 2, 5, 6, 7]),
    'G4': (E('(O1)(13)(O2)(23)(34)(45)(56)(6O)(57)(7O)'), [1, 2, 4, 6, 7]),
    'G5': (E('(O1)(13)(O2)(23)(34)(45)(5O)(46)(67)(7O)'), [1, 2, 5, 6, 7]),
}
T_strings = {   # six-d-Q3 claim C6 (colour index strings per vertex, vertex order 1..7)
    'G1': 'A1 D1 R12|R12 D2 R23|R23 A2 R34|R34 D4 R45|R45 A3 R56|R56 D6 R67|R67 D7 A4',
    'G2': 'A1 D1 R12|R12 D2 R23|R23 A2 R34|R34 D4 R45|R45 R56 R57|R56 D6 A3|R57 D7 A4',
    'G3': 'A1 D1 R12|R12 D2 R23|R23 A2 R34|R34 R45 R46|R45 D5 A3|R46 D6 R67|R67 D7 A4',
    'G4': 'A1 D1 R13|A2 D2 R23|R13 R23 R34|R34 D4 R45|R45 R56 R57|R56 D6 A3|R57 D7 A4',
    'G5': 'A1 D1 R13|A2 D2 R23|R13 R23 R34|R34 R45 R46|R45 D5 A3|R46 D6 R67|R67 D7 A4',
}


def index_names(vertices, edges, loc, spokes):
    """S2 index names of the half-edges: spoke j -> A_j, output at v -> D_v, internal edge (a,b) -> R_ab."""
    name = {}
    for j, e in enumerate(spokes):
        a, b = edges[e]
        name[str(('e', e, 1) if a == 'O' else ('e', e, 0))] = 'A%d' % (j + 1)
    for e, (a, b) in enumerate(edges):
        if 'O' not in (a, b):
            name[str(('e', e, 0))] = name[str(('e', e, 1))] = 'R%s%s' % (a, b)
    for v in vertices:
        name[str(('x', v))] = 'D%s' % v
    return '|'.join(' '.join(name[str(h)] for h in loc[v]) for v in vertices)


for nm, (vs, ed, outs) in S2.items():
    s, loc, sp_ = qn_sigma(nm, vs, ed, outs)
    print("  S2 %-3s: sigma = %+d   (indices: %s)" % (nm, s, index_names(vs, ed, loc, sp_)))
sig_q3 = []
for nm, (ed, outs) in Q3.items():
    vs = list(range(1, 8))
    s, loc, sp_ = qn_sigma(nm, vs, ed, outs)
    sig_q3.append(s)
    names = index_names(vs, ed, loc, sp_)
    print("  six-d-Q3 %s: sigma = %+d ; leg pattern equals claim C6 T-string: %s" % (nm, s, names == T_strings[nm]))
print("  six-d-Q3 signs (G1..G5):", sig_q3)


def rootN(p, q, r):
    nv = p + q + r - 1
    N_, Np = nv - 1, nv - 2
    edges, nxt, paths = [], 0, []
    for n_ in (p, q, r):
        prev, ints = N_, []
        for j in range(n_):
            cur = Np if j == n_ - 1 else nxt
            if j < n_ - 1:
                ints.append(nxt)
                nxt += 1
            edges.append((prev, cur))
            prev = cur
        paths.append(ints)
    return nv, edges, paths


msym = {}
for pqr in [(1, 3, 3), (2, 2, 3)]:
    nv, ed, paths = rootN(*pqr)
    O = nv - 1
    ed2 = [tuple('O' if x == O else x for x in ab) for ab in ed]
    outs = [v for Pp in paths for v in Pp]
    vs = list(range(nv - 1))
    s, loc, sp_ = qn_sigma('rootN', vs, ed2, outs)
    # arm tensor of msym-odd vops.py: arm_P = f_{A_P D_1 R_1} f_{R_1 D_2 R_2} ..., T = f_{R_A R_B R_C} arm_A arm_B arm_C
    arm_ok = True
    for Pi, Pp in enumerate(paths):
        for k, v in enumerate(Pp):
            hs = loc[v]
            arm_ok &= (len(hs) == 3 and hs[1] == ('x', v) and hs[0][2] == 1 and hs[2][2] == 0)
    Npr = nv - 2
    arm_ok &= [h[1] for h in loc[Npr]] == sorted(h[1] for h in loc[Npr])   # in-edges of N' in path (= edge) order
    msym['rootN%s' % (pqr,)] = s
    print("  msym-odd theta_rootN%s: sigma = %+d ; S2 pattern = arm tensor of vops.py: %s ; spokes O->v: %d"
          % (pqr, s, arm_ok, sum(1 for e in sp_ if ed2[e][0] == 'O')))

# ------------------------------------------------------------------ Part D
print("Part D: labelling invariance of sigma*K*T for S2 A, B (O = root fixed)")
for nm, G0 in [('A', htk.graph_A()), ('B', htk.graph_B())]:
    O = G0.nv - 1
    outputs_phys = [0, 1, 3, 4]
    spokes_phys = [e for e, ab in enumerate(G0.edges) if O in ab]
    # physical local order = S2 pattern in the reference labelling
    vs = list(range(O))
    ed2 = [tuple('O' if x == O else x for x in ab) for ab in G0.edges]
    loc_ref = pattern_local(vs, ed2, outputs_phys)
    for H, T in [(2, 1), (3, 0)]:
        Pp = [[sp.Symbol('p%d_%d' % (i, c)) for c in range(H)] for i in range(4)]
        globals()['P'] = Pp
        res = []
        base = None
        for trial in range(5):
            if trial == 0:
                perm = list(range(G0.nv)); order = list(range(G0.ne)); flips = [False] * G0.ne
            else:
                perm = list(range(O)); rng.shuffle(perm); perm.append(O)
                order = list(range(G0.ne)); rng.shuffle(order)
                flips = [rng.random() < 0.5 for _ in range(G0.ne)]
            # locperm relative to the explorer-style physical list used in relabelled_quantity:
            locperm = {}
            for u in range(O):
                hs = [('x', u)] if u in outputs_phys else []
                hs += [('e', e, s) for e, ab in enumerate(G0.edges) for s in (0, 1) if ab[s] == u]
                target = [('x', u) if h[0] == 'x' else h for h in loc_ref[u]]
                if trial == 0:
                    locperm[u] = [hs.index(h) for h in target]
                else:
                    locperm[u] = rng.sample(range(3), 3)
            q = relabelled_quantity(G0, None, H, T, perm, order, flips, locperm, O=O,
                                    outputs_phys=outputs_phys, spokes_phys=spokes_phys)
            if trial == 0:
                base = q
                print("  %s (H,T)=(%d,%d): reference sigma = %+d" % (nm, H, T, q[0]))
            else:
                res.append(same(base, q, cons=False))
        print("    4 random relabellings (tree vertices, edge reversals incl. spokes, edge order, local orders): invariant %s" % res)
