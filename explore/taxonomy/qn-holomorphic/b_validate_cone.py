"""Fork B: 8d Q_3 cone kernels (H=4, n=3) vs circuit-general/atlas/e1_cone_kernels.log, compared as residues mod P
(validate_kernels.part_c aborted: rational reconstruction of large rationals; symbolic 32-parameter density too slow).  Then validate_kernels parts d, e.
Log: b_validate_cone.log"""
import os, sys, random, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'circuit-general', 'atlas'))
import kern_h
from kern_h import GraphKernel, modq, P
import rgraphs
import validate_kernels as vk
from htk.graph import from_networkx
from atlas import set_cone

logv = [Fr(34835017, 41472000), Fr(8447, 2073600), Fr(-18845633, 3888000), Fr(862993709, 18662400),
        Fr(2201, 6912), Fr(-143507, 248832), Fr(7143371, 37324800), Fr(47569949, 7464960),
        Fr(-2937451, 4976640), Fr(-146710723, 7776000), Fr(-11557, 259200), Fr(7914137, 11197440)]
rng = random.Random(20261123)
items = set_cone(4, 3)
print('(c) set_cone(4,3): %d graphs; rgraphs.rooted_trees(4,3): %d' % (len(items), len(rgraphs.rooted_trees(4, 3))), flush=True)
neq = 0
for idx, (G, O, S0) in enumerate(items):
    Gh, lab = from_networkx(G, root=O)
    lam = {lab[v]: [rng.randint(-3, 3) for _ in range(4)] for v in S0 if v != O}
    t0 = time.time()
    r = Gh.nv - 1
    edges = [(a, b) if (b == r or a != r) else (b, a) for a, b in Gh.edges]   # spokes -> root; (-1)^{Hk}=+1 at H=4
    gk = GraphKernel(4, r, edges)
    assert all(v in gk.outputs for v in lam)
    lam_o = {v: tuple(lam.get(v, [0] * 4)) for v in gk.outputs}
    Kf = vk.exact_I(gk, lam_o) * gk.sign_n          # spred exact on the GraphKernel plan (numeric sources)
    K = modq(Kf)
    ok = K == modq(logv[idx])
    neq += ok
    print('  #%d plan %s K = %s, log %s: equal %s [%.1fs]' % (idx, 'tight' if gk.red.tight_plan else 'non-tight', Kf, logv[idx], Kf == logv[idx], time.time() - t0), flush=True)
print('  equal %d/12' % neq, flush=True)
vk.part_d()
vk.part_e()
