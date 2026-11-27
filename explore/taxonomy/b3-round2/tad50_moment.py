"""Exact order-1 moment of tadpole(5,0), (H,T) = (1,1): M_e = int D y_e (htk conventions) with source
l0 = 1 at vertex 0, -1 at O = 5, insertion y_e on root edge e = (v,5), by circuit-general/cores.kernel
(exact periods; same engine as atlas e3_tadpole50.log).  Edge order: sorted (as e3); converted to the
(-v,w) order of q4_deriv.py by sgn(perm)^{D+1}.  Usage: python3 tad50_moment.py v   Log: tad50_moment_v<v>.log"""
import sys, os, time
HERE = os.path.dirname(os.path.abspath(__file__))
CG = os.path.join(HERE, "..", "..", "circuit-general")
sys.path.insert(0, os.path.join(CG, "atlas"))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "..", "code")))
sys.path.insert(0, os.path.join(CG, "cores"))
import cores
sys.path.insert(0, os.path.join(HERE, "..", "qn-ht"))
from hteng import perm_sign
v = int(sys.argv[1])
log = open(os.path.join(HERE, "tad50_moment_v%d.log" % v), "w")
edges = sorted([(0, 1), (1, 2), (2, 3), (3, 4), (0, 4), (1, 5), (2, 5), (3, 5), (4, 5)])
mt = sorted(edges, key=lambda e: (-e[0], e[1]))
sg = perm_sign(edges, mt) ** 3
lam = [[1], [0], [0], [0], [0], [-1]]
t0 = time.time()
ins = () if v == 0 else ((edges.index((v, 5)), 0),)
val = cores.kernel(6, edges, lam, ins, verbose=False, T=1)
s = "tadpole(5,0) insertion %s (edge %s): K = %s ; sign to (-v,w) order %d  (%.0fs)" % (ins, (v, 5) if v else None, val, sg, time.time() - t0)
print(s); log.write(s + "\n")
