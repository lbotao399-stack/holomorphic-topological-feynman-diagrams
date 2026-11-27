"""(H,T) = (4,1): the four tight theta graphs, exact kernels K_htk(lambda;0) by the htk reduction layer
against the closed form s c_1 pi D_theta, independence of lambda_N, Laplace form of D_theta in 4-minors,
|Aut|, and the (5,0) companions (c_0 Q(Lbar))."""
import random
from fractions import Fraction as Fr
from itertools import combinations, permutations
from common import Log, parity
import htk
import twoloop

log = Log(__file__)
SHAPES = [(1, 5, 5), (2, 4, 5), (3, 3, 5), (3, 4, 4)]
C1 = {(1, 5, 5): Fr(1, 5184), (2, 4, 5): Fr(1, 9072), (3, 3, 5): Fr(1, 10080), (3, 4, 4): Fr(1, 22680)}
C0 = {(1, 5, 5): Fr(1, 17280), (2, 4, 5): Fr(1, 17280), (3, 3, 5): Fr(1, 17280), (3, 4, 4): Fr(1, 25920)}


def s_sign(H, q):
    import math
    return (-1) ** (math.ceil(H / 2) + (H + 1) * (q + 1))


def det(M):
    from common import det_int
    return det_int(M)


rng = random.Random(41)
H = 4
log('## kernels at (H,T) = (4,1): htk.kernel_exact vs s c_1 pi D_theta')
for shape in SHAPES:
    G, interior = htk.theta(*shape)
    th = twoloop.Theta(*twoloop.theta_graph(*shape))
    assert [list(x) for x in interior] == [th.interior[P] for P in 'ABC']
    for trial in range(2):
        lam = [[rng.randint(-5, 5) for _ in range(H)] for _ in range(G.nv - 1)]
        Kx = htk.kernel_exact(G, lam, T=1, check=1)
        lamfull = [list(map(Fr, v)) for v in lam] + [[-sum(Fr(v[i]) for v in lam) for i in range(H)]]
        Dth = twoloop.D_theta(th, lamfull)
        closed = twoloop.kernel_htk_closed(th, lamfull, H, 1)          # coefficient of pi
        s = s_sign(H, shape[1])
        ok1 = Kx.c == ({1: closed} if closed else {})
        ok2 = closed == s * C1[shape] * Dth
        # lambda_N independence: shift lambda_N, adjust nothing else (root absorbs)
        lam2 = [list(v) for v in lam]
        lam2[0] = [x + rng.randint(-3, 3) for x in lam2[0]]
        Kx2 = htk.kernel_exact(G, lam2, T=1, check=1)
        log(shape, 'trial', trial, 'K_exact/pi =', Kx.coeff(1), ' closed =', closed, ' s c1 D =', s * C1[shape] * Dth,
            ' equal:', ok1 and ok2, ' lambda_N-independent:', Kx2.c == Kx.c, ' s =', s)
        assert ok1 and ok2 and Kx2.c == Kx.c

log('## Laplace form: D_theta = sum eps [4-minor][4-minor] over the 8 interior sources (vertex labels 1..8)')


def laplace_terms(shape):
    th = twoloop.Theta(*twoloop.theta_graph(*shape))
    A, B, C = th.interior['A'], th.interior['B'], th.interior['C']
    # columns: A: (q;0), B: (q;q), C: (0;q); rows 0..3 top, 4..7 bottom
    cols = [('A', u) for u in A] + [('B', u) for u in B] + [('C', u) for u in C]
    terms = {}
    # Laplace along the top 4 rows: choose 4 columns S for the top block (must contain all A, no C)
    for S in combinations(range(8), 4):
        if any(cols[c][0] == 'C' for c in S) or any(cols[c][0] == 'A' and c not in S for c in range(8)):
            continue
        Sc = [c for c in range(8) if c not in S]
        sign = (-1) ** (sum(range(4)) + sum(S))
        key = (tuple(cols[c][1] for c in S), tuple(cols[c][1] for c in Sc))
        terms[key] = terms.get(key, 0) + sign
    return th, terms


for shape in SHAPES:
    th, terms = laplace_terms(shape)
    # verify against D_theta at random integer sources
    for _ in range(3):
        lamfull = [[Fr(rng.randint(-6, 6)) for _ in range(H)] for _ in range(th.nv)]
        val = sum(s * det([[lamfull[u][i] for u in I] for i in range(H)]) * det([[lamfull[u][i] for u in J] for i in range(H)])
                  for (I, J), s in terms.items())
        assert val == twoloop.D_theta(th, lamfull)
    log(shape, 'A,B,C interior =', th.interior['A'], th.interior['B'], th.interior['C'],
        ' D_theta =', ' + '.join('(%+d)[%s][%s]' % (s, ''.join(map(str, I)), ''.join(map(str, J))) for (I, J), s in sorted(terms.items())))

log('## |Aut| of the theta graphs (graph automorphisms, external legs unlabelled)')
import networkx as nx
for shape in SHAPES:
    G, _ = htk.theta(*shape)
    g = nx.MultiGraph()
    g.add_nodes_from(range(G.nv))
    g.add_edges_from(G.edges)
    n_aut = sum(1 for _ in nx.algorithms.isomorphism.MultiGraphMatcher(g, g).isomorphisms_iter())
    from math import factorial, prod
    formula = 2 * prod(factorial(shape.count(n)) for n in set(shape))
    log(shape, '|Aut| =', n_aut, ' (2 prod m_l! =', formula, ')')
    assert n_aut == formula

log('## (H,T) = (5,0) companions: htk.kernel_exact vs (sign) c_0 Q(Lbar)')
for shape in SHAPES:
    G, interior = htk.theta(*shape)
    th = twoloop.Theta(*twoloop.theta_graph(*shape))
    lam = [[rng.randint(-4, 4) for _ in range(5)] for _ in range(G.nv - 1)]
    Kx = htk.kernel_exact(G, lam, T=0, check=1)
    lamfull = [list(map(Fr, v)) for v in lam] + [[-sum(Fr(v[i]) for v in lam) for i in range(5)]]
    closed = twoloop.kernel_htk_closed(th, lamfull, 5, 0)
    Q = twoloop.Q_form(th, lamfull, twoloop.Lbar(th, lamfull))
    log(shape, 'K_exact =', Kx, ' closed =', closed, ' K/(c0 Q(Lbar)) =', closed / (C0[shape] * Q) if Q else None, ' equal:', Kx == closed)
    assert Kx == closed
log.close()
