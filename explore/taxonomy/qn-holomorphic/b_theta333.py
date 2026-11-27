"""Fork B: theta(3,3,3) rooted at a trivalent vertex (8d Q_2 graph G2.3, |Aut| = 6): the tight plan is one 3-chain move
whose symbolic 8x8 local determinant in 24 output-momentum parameters does not finish (> 5 min); use the generic plan
with the block-triangular numeric-lambda density (blocktri_h) and check it against spred at numeric sources.
Log: b_theta333.log"""
import random, time
import numpy as np
import rgraphs
import kern_h
from kern_h import GraphKernel, modq, P
from htk.spred import plan_generic
import htk


def make_gk(H, g):
    gk0 = GraphKernel(H, g['nb'], g['edges'])
    if any(len(mv.chains) >= 3 for mv in gk0.red.moves) and gk0.red.tight_plan:
        gk = GraphKernel(H, g['nb'], g['edges'], plan=plan_generic(gk0.G))
        gk.red.tight_plan = False
        return gk, True
    return gk0, False


if __name__ == '__main__':
    rng = random.Random(9)
    for g in rgraphs.rooted_trees(4, 2):
        gk, swapped = make_gk(4, g)
        if not swapped:
            continue
        print('graph', g['edges'], 'aut', g['aut'], 'generic plan moves', [[len(c) for c in mv.chains] for mv in gk.red.moves], flush=True)
        for trial in range(3):
            lam = {v: tuple(rng.randint(-3, 3) for _ in range(4)) for v in gk.outputs}
            t0 = time.time()
            L = 1
            arr = {v: np.array([[lam[v][a] % P] for a in range(4)], dtype=np.int64) for v in gk.outputs}
            val = gk.kernel_numeric(arr, L)
            t1 = time.time()
            ref = gk.sign_n * htk.kernel_exact(gk.G, {v: list(lam[v]) for v in lam}, check=0)
            print('  trial %d: blocktri I = %d (mod P), htk.kernel_exact I = %s -> equal %s  [%.1fs, %.1fs]' % (
                trial, int(val[0]), ref, int(val[0]) == modq(ref), t1 - t0, time.time() - t1), flush=True)
