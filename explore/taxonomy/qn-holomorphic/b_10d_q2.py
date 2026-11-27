"""Fork B: 10d hCS on C^5 (H=5) Q_2 graphs: rooted trees and Gamma\\O-unicyclic graphs, |Aut|, reduction plan, kernel
I = (-1)^{n(n+1)/2} K at 3 random integer sources (htk.kernel_exact, exact rationals), first moment with one input jet.
Log: b_10d_q2.log"""
import random, time
import kern_h
import rgraphs
import htk
from kern_h import GraphKernel

rng = random.Random(5)
H, n = 5, 2
gs = [(g, 'tree') for g in rgraphs.rooted_trees(H, n)] + [(g, 'unicyclic') for g in rgraphs.rooted_unicyclic(H, n)]
for k, (g, kind) in enumerate(gs, 1):
    gk = GraphKernel(H, g['nb'], g['edges'])
    vals = []
    t0 = time.time()
    for _ in range(3):
        lam = {v: [rng.randint(-3, 3) for _ in range(H)] for v in gk.outputs}
        vals.append(gk.sign_n * htk.kernel_exact(gk.G, lam, check=0))
    e0 = gk.spokes[0]
    mom = gk.sign_n * htk.kernel_exact(gk.G, lam, powers={(e0, 0): 1}, check=0)
    nout = len(gk.outputs)
    print('#%d %s %s |Aut|=%d outputs %d plan %s: I = %s; I[(p_0)_1] = %s [%.1fs]' % (
        k, kind, g['edges'], g['aut'], nout, 'tight' if gk.red.tight_plan else 'non-tight',
        ', '.join(str(v) for v in vals), mom, time.time() - t0), flush=True)
print('trees %d, unicyclic %d, sum 1/|Aut| trees = %s' % (
    sum(1 for _, kd in gs if kd == 'tree'), sum(1 for _, kd in gs if kd == 'unicyclic'),
    sum(__import__('fractions').Fraction(1, g['aut']) for g, kd in gs if kd == 'tree')), flush=True)
