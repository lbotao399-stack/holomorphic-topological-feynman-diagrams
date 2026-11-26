"""Three-loop obstruction cochain of hCS on C^3 from the nine exact bulk kernels.
O_n = N_n * sum_Gamma |Aut Gamma|^{-1} alt(K_Gamma W_Gamma),  K_Gamma = c_Gamma Delta^n  =>
O_n(y_1..y_4) = 24 N_n Delta^n T_n(x_1..x_4),  T_n = sum_Gamma c_Gamma |Aut|^{-1} Sym W_Gamma  (n odd).
T_n given on x_1 = .. = x_4 = x in the formal gl(t) basis  t^k p_lambda  (p_j = tr x^j), plus the box (n = 1).
Output: assemble_n3.json"""
import json, sys
from fractions import Fraction as Fr
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher
import colour as C
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
import htk

def aut(nv, edges):
    g = nx.Graph(); g.add_nodes_from(range(nv)); g.add_edges_from(edges)
    return sum(1 for _ in GraphMatcher(g, g).isomorphisms_iter())

def assemble(recs):
    tot = {}
    rows = []
    for r in recs:
        nv, edges, legs = r['nv'], [tuple(e) for e in r['edges']], r['legs']
        W = C.ribbon_single(C.ribbon(nv, edges, legs + [nv - 1]))
        A = aut(nv, edges)
        c = {k: Fr(v) for k, v in r['bare'].items()}
        rows.append((r['idx'], A, c, W))
        for lab, cv in c.items():
            for key, w in W.items():
                kk = (lab,) + key
                tot[kk] = tot.get(kk, 0) + cv * w / A
    return rows, {k: v for k, v in tot.items() if v != 0}

def show(tot):
    labs = sorted(set(k[0] for k in tot))
    for lab in labs:
        terms = sorted((k[1], k[2], v) for k, v in tot.items() if k[0] == lab)
        print('   [%s] ' % lab + ' + '.join('(%s) t^%d p%s' % (v, tp, ''.join(map(str, part))) for tp, part, v in terms))

# box
box = dict(idx='box', nv=4, edges=[(0, 1), (1, 2), (2, 3), (3, 0)], legs=[0, 1, 2])
kb = Fr(str(htk.kernel_exact(htk.Graph(4, box['edges']), {0: [1, 0, 0], 1: [0, 1, 0], 2: [0, 0, 1]})))
box['bare'] = {'1': str(kb)}
rows1, T1 = assemble([box])
print('n=1: c_box =', kb, ' |Aut| =', rows1[0][1])
print(' W_box(x^4) =', rows1[0][3])
print(' T_1 ='); show(T1)
recs = json.load(open('kernels_n3.json'))
rows3, T3 = assemble(recs)
print('n=3 graphs: idx, |Aut|, c_Gamma, W_Gamma(x^4) in (t^k, partition):')
for idx, A, c, W in rows3:
    print('  #%d |Aut|=%d c=%s W=%s' % (idx, A, {k: str(v) for k, v in c.items()}, W))
print(' T_3 ='); show(T3)
# specializations: gl(N|N) = gl(t=0); sl_N: p1 = 0; general N
def special(tot, t, kill_p1):
    out = {}
    for (lab, tp, part), v in tot.items():
        if kill_p1 and 1 in part:
            continue
        key = (lab, part)
        out[key] = out.get(key, 0) + v * Fr(t) ** tp
    return {k: v for k, v in out.items() if v != 0}
print(' T_3 at t = 0 (gl(N|N), str):', {('%s,p%s' % (k[0], ''.join(map(str, k[1])))): str(v) for k, v in special(T3, 0, False).items()})
print(' T_1 at t = 0 (gl(N|N), str):', {('%s,p%s' % (k[0], ''.join(map(str, k[1])))): str(v) for k, v in special(T1, 0, False).items()})
json.dump({'T1': {repr(k): str(v) for k, v in T1.items()}, 'T3': {repr(k): str(v) for k, v in T3.items()},
           'rows3': [[i, A, {k: str(v) for k, v in c.items()}, {repr(k): v for k, v in W.items()}] for i, A, c, W in rows3]},
          open('assemble_n3.json', 'w'), indent=1)
