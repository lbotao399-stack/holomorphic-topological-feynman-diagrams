"""Fork B: symbolic densities with the b_fastdet patch; timing and agreement with spred numeric.  Log: b_time_dens2.log"""
import sys, time, random
import kern_h
import b_fastdet
import rgraphs, htk
from kern_h import GraphKernel, modq, P
import numpy as np
H, n = int(sys.argv[1]), int(sys.argv[2])
which = [int(x) for x in sys.argv[3:]]
rng = random.Random(1)
for k, g in enumerate(rgraphs.rooted_trees(H, n), 1):
    if which and k not in which:
        continue
    gk = GraphKernel(H, g['nb'], g['edges'])
    if not gk.red.tight_plan:
        print('#%d non-tight, skipped' % k, flush=True)
        continue
    t0 = time.time()
    d = gk.density_factors()
    t1 = time.time()
    lam = {v: tuple(rng.randint(-3, 3) for _ in range(H)) for v in gk.outputs}
    arr = {v: np.array([[lam[v][a] % P] for a in range(H)], dtype=np.int64) for v in gk.outputs}
    val = gk.kernel_numeric(arr, 1)
    ref = gk.sign_n * htk.kernel_exact(gk.G, {v: list(lam[v]) for v in lam}, check=0)
    print('#%d moves %s: factor sizes %s [%.1fs]; XP == exact %s' % (
        k, [[len(c) for c in mv.chains] for mv in gk.red.moves], [len(f.keys) for f in d], t1 - t0,
        int(val[0]) == modq(ref)), flush=True)
