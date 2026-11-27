"""Fork B: 8d order-3 kernels as used by b_run_nilp_f8.py (tight plans: XP symbolic densities with b_fastdet;
>= 3-chain moves: block-triangular numeric density) against htk.kernel_exact at 2 random sources each.
Log: b_check_order3_kernels.log"""
import random, time
import numpy as np
import kern_h
import b_fastdet
import htk
import rgraphs, b_theta333
from kern_h import GraphKernel, modq, P
rng = random.Random(77)
H = 4
for k, g in enumerate(rgraphs.rooted_trees(H, 3), 1):
    t0 = time.time()
    gk = GraphKernel(H, g['nb'], g['edges'])
    mode = 'tight-symbolic'
    if any(len(mv.chains) >= 3 for mv in gk.red.moves):
        if gk.red.tight_plan:
            gk, _ = b_theta333.make_gk(H, g)
        mode = 'block-triangular'
    oks = []
    for _ in range(2):
        lam = {v: tuple(rng.randint(-3, 3) for _ in range(H)) for v in gk.outputs}
        arr = {v: np.array([[lam[v][a] % P] for a in range(H)], dtype=np.int64) for v in gk.outputs}
        val = gk.kernel_numeric(arr, 1)
        ref = gk.sign_n * htk.kernel_exact(gk.G, {v: list(lam[v]) for v in lam}, check=0)
        oks.append(int(val[0]) == modq(ref))
    print('#%d |Aut| %d %s: XP == exact %s [%.1fs]' % (k, g['aut'], mode, oks, time.time() - t0), flush=True)
