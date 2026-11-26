"""Three-loop (n = 3, odd) tight graphs at (H,T) = (2,1): brute-force T = 1 density at random rational points with
non-degenerate sources (entries in [-40,40]); control: the two-loop thetas with the same generator."""
import sys, random
from fractions import Fraction as Fr
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
import htk
random.seed(21)
def rl(nv):
    return [[Fr(random.randint(-40, 40), random.randint(1, 7)), Fr(random.randint(-40, 40), random.randint(1, 7))] for _ in range(nv - 1)]
def rt(ne):
    return [Fr(random.randint(1, 30), random.randint(1, 30)) for _ in range(ne)]
for pqr in [(1, 3, 3), (2, 2, 3)]:
    G, _ = htk.theta(*pqr)
    print("control theta%s: nonzero at 3 points:" % (pqr,), [htk.brute.density(G, rt(G.ne), rl(G.nv), T=1).coeff != 0 for _ in range(3)])
gs = htk.catalog.tight_graphs(3, 3)
zs = []
for i, g in enumerate(gs):
    zs.append(all(htk.brute.density(g, rt(g.ne), rl(g.nv), T=1).coeff == 0 for _ in range(2)))
print("three-loop tight graphs (%d): density zero at 2 random points each: %d/%d" % (len(gs), sum(zs), len(gs)))
