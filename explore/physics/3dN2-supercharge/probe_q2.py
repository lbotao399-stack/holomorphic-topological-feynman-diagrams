import time
from fractions import Fraction as Fr
from ht3d import Engine, mono, pstr, psum, pscale, load_kernels
from theories import *
K = load_kernels(4)
th = su2_adj(1)
E = Engine(th, K)
A = E.A
for lets in [[('B',0,0),('B',1,0)], [('B',0,0),('B',1,0),('B',2,0)], [('B',0,0),('s0',1,0)], [('s0',0,0),('s0',1,0),('s0',2,0)],
             [('B',0,0),('B',0,0),('B',1,0)], [('B',0,0),('s0',1,0),('p0',2,0)]]:
    O = mono(A, *lets)
    t = time.time()
    R2 = E.Q2(O)
    Rf = E.Q2(O, 'fan'); Rt = E.Q2(O, 'tad')
    C = psum(E.Q0(R2), E.Q2(E.Q0(O)))
    print(pstr(O), "(%.1fs)" % (time.time()-t))
    print("   Q2 fan :", pstr(Rf))
    print("   Q2 tad :", pstr(Rt))
    print("   {Q0,Q2}:", pstr(C))
