"""Fork B: timing of symbolic densities (kern_h.GraphKernel.density_factors) for H=4 n=2 trees.  Log: b_time_dens.log"""
import sys, time
import rgraphs
from kern_h import GraphKernel
H, n = int(sys.argv[1]), int(sys.argv[2])
for k, g in enumerate(rgraphs.rooted_trees(H, n), 1):
    gk = GraphKernel(H, g['nb'], g['edges'])
    moves = [(len(mv.chains), [len(c) for c in mv.chains]) for mv in gk.red.moves]
    t0 = time.time()
    if len(sys.argv) > 3 and int(sys.argv[3]) != k:
        print('#%d %s tight %s moves %s' % (k, g['edges'], gk.red.tight_plan, moves), flush=True)
        continue
    d = gk.density_factors()
    nterms = [len(f.keys) for f in d] if isinstance(d, list) else 'blocktri'
    print('#%d %s tight %s moves %s: factor sizes %s [%.1fs]' % (k, g['edges'], gk.red.tight_plan, moves, nterms, time.time() - t0), flush=True)
