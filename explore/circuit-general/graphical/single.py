"""Single period I_q(G; alpha) with plan details.
usage: python3 single.py 'EDGES' 'ALPHAS' Q      e.g.  "[(0,1),(0,2),(1,2),(0,3),(1,3),(2,3)]" "[1,1,1,2,2,2]" 3"""
import ast
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import seq
import mzv
import sv

if __name__ == '__main__':
    E = ast.literal_eval(sys.argv[1])
    A = ast.literal_eval(sys.argv[2])
    q = int(sys.argv[3])
    t = time.time()
    cands = seq.plans_for(E, A, q)
    print('edges %s alpha %s q=%d: %d plans; SV_WMAX=%d' % (E, A, q, len(cands), sv.WMAX))
    val, pl = seq.period(E, A, q, verbose=True)
    print('route %s hubs %s path %s' % (pl['route'], pl['hubs'], pl['path']))
    print('I_q = %s   (%.1fs)' % (mzv.tostr(val), time.time() - t))
    print('numeric %s' % mzv.numeric(val, 25))
