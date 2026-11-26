"""VERIFY (b1), (b3) |Aut|: counting, minors mod conservation, S_4 sign of Delta, automorphisms.
Own code (sympy, networkx); test graphs from htk.catalog only as a list."""
import sys, itertools
import sympy as sp
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher

# 1. internal degree >= 2 from 2|I_S| <= 3|V_S| - 4 applied to Gamma minus a vertex of internal degree d
n, d = sp.symbols('n d', integer=True, positive=True)
I, V = 3 * n + 1, 2 * n + 2
ineq = sp.simplify((2 * (I - d)) - (3 * (V - 1) - 4))
print("2(|I|-d) - (3(|V|-1)-4) =", ineq, "(must be <= 0, so d >= 3/2, d >= 2)")
E = sp.simplify(3 * V - 2 * I)
print("external legs with cubic vertices: 3|V| - 2|I| =", E)

# 2. minors mod conservation
p = [sp.Matrix(sp.symbols('p%d_1:4' % i)) for i in range(1, 5)]
sub = {p[3][c]: -(p[0][c] + p[1][c] + p[2][c]) for c in range(3)}
br = lambda i, j, k: sp.Matrix.hstack(p[i - 1], p[j - 1], p[k - 1]).det()
Delta = br(1, 2, 3)
for (i, j, k), s in (((1, 2, 4), -1), ((1, 3, 4), 1), ((2, 3, 4), -1)):
    print("[%d%d%d] = %+d Delta mod conservation: %s" % (i, j, k, s, sp.expand((br(i, j, k) - s * Delta).subs(sub, simultaneous=True)) == 0))
ok = True
for sgm in itertools.permutations(range(4)):
    sg = sp.combinatorics.Permutation(list(sgm)).signature()
    Ds = sp.Matrix.hstack(p[sgm[0]], p[sgm[1]], p[sgm[2]]).det()
    ok &= sp.expand((Ds - sg * Delta).subs(sub, simultaneous=True)) == 0
print("Delta(p_s1, p_s2, p_s3) = sgn(s) Delta mod conservation for all 24 s:", ok)

# 3. alternation of P(p) W with P = Delta^n: alt(Delta^n K) = Delta^n sum_s sgn(s)^{n+1} K^s (symbolic check of the
#    sign bookkeeping on a generic 4-tensor with 2 colour values)
for nn in (1, 2, 3):
    Ksym = {idx: sp.Symbol('K_%d%d%d%d' % idx) for idx in itertools.product(range(2), repeat=4)}
    okn = True
    for idx in itertools.product(range(2), repeat=4):
        lhs = 0
        for sgm in itertools.permutations(range(4)):
            sg = sp.combinatorics.Permutation(list(sgm)).signature()
            Ds = sp.Matrix.hstack(p[sgm[0]], p[sgm[1]], p[sgm[2]]).det()
            lhs += sg * Ds ** nn * Ksym[tuple(idx[sgm[i]] for i in range(4))]
        rhs = Delta ** nn * sum(sp.combinatorics.Permutation(list(sgm)).signature() ** (nn + 1) * Ksym[tuple(idx[sgm[i]] for i in range(4))]
                                for sgm in itertools.permutations(range(4)))
        okn &= sp.expand((lhs - rhs).subs(sub, simultaneous=True)) == 0
    print("n = %d: alt(Delta^n K) = Delta^n sum_s sgn(s)^(n+1) K^s mod conservation (all 16 index values): %s" % (nn, okn))

# 4. automorphisms with marked external vertices
def aut(edges, ext):
    G = nx.Graph()
    G.add_edges_from(edges)
    for v in G.nodes:
        G.nodes[v]['ext'] = v in ext
    GM = GraphMatcher(G, G, node_match=lambda a, b: a['ext'] == b['ext'])
    return sum(1 for _ in GM.isomorphisms_iter())


print("|Aut| box (marked 0..3):", aut([(0, 1), (1, 2), (2, 3), (3, 0)], [0, 1, 2, 3]))
print("|Aut| theta(1,3,3) (marked 1..4):", aut([(0, 5), (0, 1), (1, 2), (2, 5), (0, 3), (3, 4), (4, 5)], [1, 2, 3, 4]))
print("|Aut| theta(2,2,3) (marked 1..4):", aut([(0, 1), (1, 5), (0, 2), (2, 5), (0, 3), (3, 4), (4, 5)], [1, 2, 3, 4]))
print("|Aut| S2 graph A with O = 5 marked as operator:", aut([(5, 0), (0, 1), (1, 2), (5, 2), (2, 3), (3, 4), (4, 5)], [5]))

# 5. tight D = 3 graphs with cubic completion (internal degrees 2 or 3, four of degree 2)
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
from htk import catalog
for nl in (1, 2, 3):
    gs = catalog.tight_graphs(3, nl)
    cub = []
    for G in gs:
        degs = [sum((a == v) + (b == v) for a, b in G.edges) for v in range(G.nv)]
        cub.append((min(degs), max(degs), sorted(degs).count(2)))
    print("D=3 n=%d: %d tight graphs; (min deg, max deg, #deg-2) = %s; with all degrees in {2,3}: %d" % (
        nl, len(gs), cub, sum(1 for c in cub if c[0] >= 2 and c[1] <= 3)))
