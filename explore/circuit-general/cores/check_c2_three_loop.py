import sys, time
sys.path.insert(0,'.')
import cores
from fractions import Fraction as F
edges = [(0,1),(0,2),(0,3),(1,4),(4,2),(1,3),(2,3)]
lam = [(1,2),(2,-1),(-1,3),(-3,-4),(1,0)]
for ins, ref in [((), {'1':F(-71,12)}), (((0,0),), {'1':F(222317,25200),'z3':F(-3021,700)}), (((0,0),(0,0)), {'1':F(-4003781,302400),'z3':F(3772,525)})]:
    t=time.time(); v = cores.kernel(5, edges, lam, ins, verbose=True)
    print(ins, v, 'match' if v==ref else ('MISMATCH ref %s'%ref), round(time.time()-t,1),'s', flush=True)
