"""Labelling independence of c_Gamma * W_Gamma(x,x,x,x) (gl(t), formal trace basis) at n = 3 (and the box, n = 1):
random vertex relabellings (root = any degree-2 vertex, random leg order), random edge reversals and edge orders;
kernel recomputed exactly (assembly.kernel_moments, sources e_1,e_2,e_3 on the new legs), colour by the orientation rule."""
import json, os, sys, random
from fractions import Fraction as Fr
A = (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/circuit-general/assembly'
sys.path.insert(0, A)
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
os.environ.setdefault('NUMBA_CACHE_DIR', '/private/tmp/claude-501/-Users-libotao-Desktop/07c118b5-2480-4b3c-899c-a49c36ac4d00/scratchpad/numba_anom')
import io, contextlib
import assembly
import htk
import colour as C
rng = random.Random(17)

def kernel(nv, edges, legs):
    full = [[0] * 3 for _ in range(nv)]
    for a, v in enumerate(legs):
        full[v][a] = 1
    full[nv - 1] = [-sum(r[i] for r in full[:-1]) for i in range(3)]
    import networkx as nx
    g = nx.Graph(); g.add_edges_from(edges)
    if htk.catalog.is_sp(g):
        lam = {v: full[v] for v in range(nv - 1)}
        k = Fr(str(htk.kernel_exact(htk.Graph(nv, edges), lam)))
        return {'1': k} if k != 0 else {}
    with contextlib.redirect_stdout(io.StringIO()):
        res = assembly.kernel_moments(nv, edges, full, order=0)
    return {k: Fr(v) for k, v in res['bare'].items()}

def cw(nv, edges, legs):
    K = kernel(nv, edges, legs)
    W = C.ribbon_single(C.ribbon(nv, edges, legs + [nv - 1]))
    out = {}
    for lab, kv in K.items():
        for key, w in W.items():
            out[(lab,) + key] = out.get((lab,) + key, 0) + kv * w
    return {k: v for k, v in out.items() if v != 0}

def relabel(nv, edges):
    deg = [0] * nv
    for a, b in edges:
        deg[a] += 1; deg[b] += 1
    biv = [v for v in range(nv) if deg[v] == 2]
    root = rng.choice(biv)
    others = [v for v in range(nv) if v != root]
    rng.shuffle(others)
    new = {v: i for i, v in enumerate(others)}
    new[root] = nv - 1
    E = []
    for a, b in edges:
        a2, b2 = new[a], new[b]
        if rng.random() < 0.5:
            a2, b2 = b2, a2
        E.append((a2, b2))
    rng.shuffle(E)
    legs = [new[v] for v in biv if v != root]
    rng.shuffle(legs)
    return E, legs

G3 = json.load(open('kernels_n3.json'))
cases = [('box', 4, [(0, 1), (1, 2), (2, 3), (3, 0)], [0, 1, 2])]
for r in G3:
    cases.append(('n3#%d' % r['idx'], r['nv'], [tuple(e) for e in r['edges']], r['legs']))
for name, nv, edges, legs in cases:
    ref = cw(nv, edges, legs)
    res = []
    for _ in range(4):
        E, L = relabel(nv, edges)
        res.append(cw(nv, E, L) == ref)
    print('%-6s c*W(x^4) invariant under 4 random relabellings/reversals/edge orders: %s ; reference zero: %s' % (name, res, not ref), flush=True)
