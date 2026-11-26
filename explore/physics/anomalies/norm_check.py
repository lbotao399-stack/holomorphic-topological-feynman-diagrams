"""Absolute normalization of O_n from the six-d-Q3 ordered Wick rule for Q_n (explore/six-d-Q3/README 'Ordered Wick weights';
wick_weights.wick_sign imported unchanged) and O_n = Q_n I/(2n+2).
For every bulk graph (box; both thetas; nine three-loop graphs) and every root r (all 2n+2 vertices):
  F_Q(r)  = s_Q(Gamma, r) * W^f[source words, root word = spokes in edge order + root output] * I^{S2}_Gamma[root r labelling]
            output order: bivalent non-root vertices in vertex order, then the root output;
  F_rule  = rule sign * W^f[orientation rule, ext = legs + root] * K_Gamma[standard labelling];
  rho(Gamma, r) = F_Q(r) * sgn(output order -> ext order) / F_rule  (formal gl(t) colour, exact kernels at random momenta).
Claim: rho depends on n only.  Then N_n = (-1)^{n+1+n(n+1)/2} rho_n (8 pi^3)^{-n}  (derivation in REPORT)."""
import json, os, sys, random, io, contextlib, itertools
from fractions import Fraction as Fr
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/six-d-Q3')
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/circuit-general/assembly')
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
os.environ.setdefault('NUMBA_CACHE_DIR', '/private/tmp/claude-501/-Users-libotao-Desktop/07c118b5-2480-4b3c-899c-a49c36ac4d00/scratchpad/numba_anom')
from wick_weights import wick_sign
import networkx as nx
import htk, assembly
import colour as C
rng = random.Random(53)

def det3(a, b, c):
    return (a[0]*(b[1]*c[2]-b[2]*c[1]) - a[1]*(b[0]*c[2]-b[2]*c[0]) + a[2]*(b[0]*c[1]-b[1]*c[0]))

def engine_kernel(nv, edges, src):
    """src: dict vertex -> 3-vector over all vertices (sum zero); root = nv-1"""
    g = nx.Graph(); g.add_edges_from(edges)
    full = [src.get(v, [0, 0, 0]) for v in range(nv)]
    if htk.catalog.is_sp(g):
        lam = {v: full[v] for v in range(nv - 1)}
        return {'1': Fr(str(htk.kernel_exact(htk.Graph(nv, edges), lam)))}
    with contextlib.redirect_stdout(io.StringIO()):
        res = assembly.kernel_moments(nv, edges, full, order=0)
    return {k: Fr(v) for k, v in res['bare'].items()}

def ribbon_words(words, partner, slot):
    out = {}
    nv = len(words)
    for orient in itertools.product((0, 1), repeat=nv):
        nxt, sign = {}, 1
        for v in range(nv):
            w0, w1, w2 = words[v]
            cyc = (w0, w1, w2) if orient[v] == 0 else (w0, w2, w1)
            if orient[v]:
                sign = -sign
            for k in range(3):
                nxt[cyc[k]] = cyc[(k + 1) % 3]
        seen, faces, empty = set(), [], 0
        for start in [(h, s) for w in words for h in w for s in ('L', 'R')]:
            if start in seen:
                continue
            word, node = [], start
            while node not in seen:
                seen.add(node)
                h, side = node
                if side == 'R':
                    node = (nxt[h], 'L')
                elif h[0] == 'x':
                    word.append(slot[h]); node = (h, 'R')
                else:
                    node = (partner[h], 'R')
            if word:
                faces.append(C.canon_cyc(word))
            else:
                empty += 1
        key = (empty, tuple(sorted(faces)))
        out[key] = out.get(key, 0) + sign
    return {k: v for k, v in out.items() if v}

def perm_sign(p):
    s = 1
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]:
                s = -s
    return s

def scale(d, c):
    return {k: v * c for k, v in d.items() if v * c != 0}

def F_rule(nv, edges, ext, mom):
    W = C.ribbon(nv, edges, ext)
    K = engine_kernel(nv, edges, {v: mom[v] for v in ext})
    return {lab: scale(W, kv) for lab, kv in K.items()}, ext

def F_Q(nv, edges, r, mom, ext_rule):
    # relabel: root r -> nv-1, others keep order
    others = [v for v in range(nv) if v != r]
    lab = {v: i for i, v in enumerate(others)}; lab[r] = nv - 1
    inv = {i: v for v, i in lab.items()}
    E = [(lab[a], lab[b]) for a, b in edges]
    root = nv - 1
    s = wick_sign(E, root)
    deg = [0] * nv
    for a, b in E:
        deg[a] += 1; deg[b] += 1
    words, outs = [], []
    for v in range(root):
        ins = [('e', i, 1) for i, (a, b) in enumerate(E) if b == v]
        oo = [('e', i, 0) for i, (a, b) in enumerate(E) if a == v]
        mid = [('x', v)] if deg[v] == 2 else []
        if mid:
            outs.append(v)
        words.append(ins + mid + oo)
    rw = [('e', i, 1 if b == root else 0) for i, (a, b) in enumerate(E) if root in (a, b)]
    if deg[root] == 2:
        rw.append(('x', root)); outs.append(root)
    words.append(rw)
    partner = {}
    for i in range(len(E)):
        partner[('e', i, 0)] = ('e', i, 1); partner[('e', i, 1)] = ('e', i, 0)
    # slots: ext_rule order (original vertex labels)
    slot = {('x', v): ext_rule.index(inv[v]) for v in outs}
    W = ribbon_words(words, partner, slot)
    pi = [ext_rule.index(inv[v]) for v in outs]
    K = engine_kernel(nv, E, {lab[v]: mom[v] for v in ext_rule})
    return {labk: scale(W, kv * s * perm_sign(pi)) for labk, kv in K.items()}

def ratio(a, b):
    """a = rho * b with rho = +-1 ?"""
    if not any(a.values()) and not any(b.values()):
        return 0   # both vanish identically: comparison vacuous
    for rho in (1, -1):
        if all(set(a.get(k, {})) == set(b.get(k, {})) and all(a.get(k, {})[kk] == rho * b[k][kk] for kk in b.get(k, {})) for k in set(a) | set(b)):
            return rho
    return None

cases = [(1, 'box', 4, [(0, 1), (1, 2), (2, 3), (3, 0)], [0, 1, 2, 3]),
         (2, 'theta133', 6, [(0, 5), (0, 1), (1, 2), (2, 5), (0, 3), (3, 4), (4, 5)], [1, 2, 3, 4]),
         (2, 'theta223', 6, [(0, 1), (1, 5), (0, 2), (2, 5), (0, 3), (3, 4), (4, 5)], [1, 2, 3, 4])]
for rec in json.load(open('kernels_n3.json')):
    cases.append((3, 'n3#%d' % rec['idx'], rec['nv'], [tuple(e) for e in rec['edges']], rec['legs'] + [rec['nv'] - 1]))
for n, name, nv, edges, ext in cases:
    ps = [[rng.randint(-5, 5) for _ in range(3)] for _ in range(3)]
    ps.append([-sum(p[a] for p in ps) for a in range(3)])
    mom = {v: ps[i] for i, v in enumerate(ext)}
    Fr_, _ = F_rule(nv, edges, ext, mom)
    rhos = []
    for r in range(nv):
        rhos.append(ratio(F_Q(nv, edges, r, mom, ext), Fr_))
    print('n=%d %-9s rho over all %d roots: %s' % (n, name, nv, rhos), flush=True)
