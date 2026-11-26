"""Build and cache prism / K33 period tables: python3 build_tables.py prism 7"""
import sys
import time

import cores

name, q = sys.argv[1], int(sys.argv[2])
t = time.time()
tab = cores.table(name, q)
print('%s q=%d: %d orbits, %.1fs' % (name, q, len(tab), time.time() - t), flush=True)
