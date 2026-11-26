"""Build prism / K33 period tables at level q into assembly/tables:  python3 build_tables.py prism 8"""
import resource
import sys
import time

import periods

name, q = sys.argv[1], int(sys.argv[2])
t = time.time()
tab = periods.table(name, q)
print('%s q=%d: %d orbits, %.1fs, maxrss %.2f GB' % (name, q, len(tab), time.time() - t,
      resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9), flush=True)
