"""Exact holomorphic kernels and core periods.

Kernel (T = 0 or 1, zero topological source, htk conventions):
  python3 evaluate.py kernel NV 'EDGES' 'SOURCES' ['INSERTIONS'] [T]
     EDGES       [[u,v],...] oriented, root = last vertex
     SOURCES     list over all vertices (sum 0) or over non-root vertices, H-vectors
     INSERTIONS  [[e,i],...]: multiply by component i of the current on edge e (= -d/dz_e^i)
  e.g. python3 evaluate.py kernel 5 '[[0,1],[0,2],[0,3],[1,4],[4,2],[1,3],[2,3]]' \
            '[[1,2],[2,-1],[-1,3],[-3,-4],[1,0]]' '[[0,0]]'
Period:
  python3 evaluate.py period 'EDGES' 'POWERS' Q
     int prod t^(a-1) / U^Q over RP^{N-1}_{>0}, any graph whose series-parallel reduction
     ends in W_n (n>=3, K4 = W3), the prism or K33; Q integer or half-integer ('9/2').
Output labels: '1', 'z3', 'z5', ... (T = 0);  monomials in sqrtpi, z3, ... (half-integer Q).
"""
import json
import sys
import time
from fractions import Fraction as F

import cores


def fmt(v):
    return {cores.enc(k): str(x) for k, x in sorted(v.items(), key=lambda kv: str(kv[0]))}


def main(argv):
    if not argv:
        print(__doc__)
        return
    t = time.time()
    if argv[0] == 'kernel':
        nv = int(argv[1])
        edges = [tuple(e) for e in json.loads(argv[2])]
        lam = json.loads(argv[3])
        ins = [tuple(x) for x in json.loads(argv[4])] if len(argv) > 4 else []
        T = int(argv[5]) if len(argv) > 5 else 0
        val = cores.kernel(nv, edges, lam, ins, verbose=True, T=T)
    elif argv[0] == 'period':
        edges = [tuple(e) for e in json.loads(argv[1])]
        powers = json.loads(argv[2])
        q = F(argv[3])
        if q.denominator == 1:
            q = int(q)
        val = cores.period(edges, powers, q)
    else:
        raise SystemExit(__doc__)
    print(json.dumps(fmt(val)), ' (%.1fs)' % (time.time() - t))


if __name__ == '__main__':
    main(sys.argv[1:])
