"""Direct numerics (htk.oracle lattice, no reduction) for the C^2 three-loop non-SP graph."""
import sys, os, time
from fractions import Fraction as Fr
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'code'))
from htk.graph import Graph
from htk import oracle
N3 = [(0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 4)]
G = Graph(5, N3)
SA, SB = Fr(1, 3), Fr(2, 5)
cfgs = {'[3;0,2,4] SA=1/3 SB=2/5': [[0, 1], [0, 0], [SA, SB], [1, 0], [-1 - SA, -1 - SB]],
        '[3;0,2]': [[0, 1], [0, 0], [-1, -1], [1, 0], [0, 0]]}
for name, lam in cfgs.items():
    t0 = time.time()
    v, err = oracle.kernel(G, [[float(x) for x in r] for r in lam[:4]], 0)
    print(name, 'K =', v, '+-', err, ' time %.0fs' % (time.time() - t0))
K = -(1 + 3 * SA - 3 * SA ** 2 + SB + 12 * SA * SB - 5 * SB ** 2) / 576
print('exact [3;0,2,4] at SA=1/3, SB=2/5:', K, float(K), '   exact [3;0,2]:', -1 / 576)
